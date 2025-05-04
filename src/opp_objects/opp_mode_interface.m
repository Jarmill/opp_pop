classdef (Abstract) opp_mode_interface
    %OPP_MODE_INTERFACE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        mode;          %the mode m (defines the id)
        opts;          %relevant options for the mode
        levels;        %locations for each inverter voltage level
        L;             %levels of the inverter        
        transition;    %guard measures for the partition staying within the level (no switching)
        vars;          %basic variable type
        
        Z_load;
        Symmetry; 
        f0;
        
    end
    
    methods
        %constructor
        function obj = opp_mode_interface(m, lsupp_ref, objective_mode, opts)
            %OPP_MODE_INTERFACE Construct an instance of this class
            %   Detailed explanation goes here
                        obj.mode = m;
            obj.L = opts.L;
            obj.f0 = opts.f0;
            obj.Symmetry = opts.Symmetry;
            obj.Z_load = opts.Z_load;
             

           lsupp_base = obj.refine_support(opts, m, lsupp_ref);
           obj.vars = lsupp_base.vars;

            %TODO: define grid-side filter dynamics
            %define the dynamics within the mode                        

            f = obj.all_dynamics(obj.vars, opts);
            
            %create locations for each level     
            [prefix] = obj.get_prefix();

            obj.levels = obj.make_level_locs(m, opts, lsupp_base, objective_mode, f, prefix);                  
            
            obj.transition = obj.make_transitions(m, opts,  lsupp_base, prefix);
        end

        function levels =  make_level_locs(obj, m, opts, lsupp_base, objective_mode, f, prefix)
            %create locations for each level
            vars = lsupp_base.vars;       
            N = size(opts.L, 2);
            P = opts.partition;
            levels = cell(N, P);
            X_partition = support_partition(opts.partition, vars.x, opts.Symmetry);
            for n = 1:N
                curr_f = f(:, n);
                curr_objective = objective_mode(n, :);                
                for p = 1:P
                    % curr_id = loc_id(n, j);
                    curr_id = strcat(prefix, sprintf('m%d_n%d_p%d', m, n, p));

                    curr_lsupp = lsupp_base;

                    if p>1 || obj.cannot_start(opts, n)
                        curr_lsupp.X_init = [];
                    end
                    if p<P || obj.cannot_end(opts, n)
                        curr_lsupp.X_term = [];
                    end

                    if opts.partition > 1
                        curr_lsupp.X = [curr_lsupp.X; X_partition(p)>=0];
                    end

                    if ~isempty(opts.allowed_levels)                                                
                        if obj.cannot_reach(opts, m, n)
                            curr_lsupp.X = [];
                            curr_lsupp.X_term = [];
                            curr_lsupp.X_init = [];
                        end
                    end

                    cell_info = struct('mode', m, 'partition', p, 'level', n, 'L', opts.L(n), 'id', curr_id);
                    levels{n, p} = opp_location(curr_lsupp, curr_f, curr_objective, cell_info);

                end
            end
        end

        function transition = make_transitions(obj, m, opts, lsupp_base, prefix)
            %create the transitions internal to the mode (no switching)                      
            N = size(opts.L, 2);
            P = opts.partition;
            transition = cell(N, P-1);
            vars = lsupp_base.vars;  

            for n=1:N
                for p=1:P-1
                    curr_trans_id = strcat(prefix, sprintf('trans_m%d_n%d_p%d', m, n, p));
                     if ~isempty(opts.allowed_levels) && ~opts.allowed_levels(m+1, n)
                         curr_supp = [];
                     else
                        curr_supp = lsupp_base.X(2:end-1);
                        RotAngle = 2*pi/double(P*2^opts.Symmetry);
                        dp = double(p);
                        new_con = (vars.x(1:2)==[cos(dp*RotAngle); sin(dp*RotAngle)]);
                        curr_supp= [new_con; curr_supp];
                    end
                    transition{n, p} = guard(curr_trans_id, lsupp_base.vars, ...
                        obj.levels{n, p}, obj.levels{n, p+1}, curr_supp, vars.x);
                end
            end
        end

        %dynamics
        function f_load = load_dynamics(obj, x_load, opts)
            %dynamics in the mode
            %(trig spinning around a circle, clock increasing, load
            %charging/modifying the current)
            %vars: variables (t, x)
            %opts: options from opp_options
            %n: level of the inverter

            %scaled inverter value
            Lscale = 2*opts.L/max(max(abs(opts.L)));
            % u_curr = opts.L(n)/max(abs(opts.L)); 

            N = size(opts.L, 2);

            %dynamics of the load
            if (isempty(x_load) || (imag(opts.Z_load) == 0))                     
                %purely resistive
                f_load = [];
            elseif (imag(opts.Z_load) >= 0)
                %inductive load
                %i' = -(R/L)i + (1/L) v
                %per-unit system, ignore the L value
                inductance = imag(opts.Z_load)/(2*pi*opts.f0);
                resistance= real(opts.Z_load);                
                f_load = -((resistance)/(inductance))*x_load*ones(1, N) + Lscale;
            else
                 %vc' = (v-vc)/(R*C)
                 %per-unit, ignore (R*C) factor
                 %TODO: v is from the voltage source. Modify when it is 
                 %filtered by a grid-side filter
                 capacitance= -imag(opts.Z_load)*(2*pi*opts.f0);
                 resistance= real(opts.Z_load);  
                 f_load = Lscale - x_load*ones(1, N)*(1/(resistance*capacitance));
            end

        end       
       
        %% describing constraints
        function liou = flow(obj, d)
            %return the continuity equation within the mode
            [N, P] = size(obj.levels);
            liou = cell(N, P);
            %start with continuity within the location
            %positive: incoming, negative: outgoing
            for n=1:N
                for p = 1:P
                    liou{n, p} = obj.levels{n, p}.liou_con(d);
                end
            end

            %now handle transitions within the mode
            for n=1:N
                for p = 1:P-1
                    trans_loss = obj.transition{n, p}.reset_push(d);
                    liou{n, p} = liou{n, p} - trans_loss;
                    liou{n, p+1} = liou{n, p+1} + trans_loss;
                end
            end

            %return liou
            %the manager will sew together jumps between multiple modes.

        end

        %% fetching moments
        function [mass_init_mode, mass_sum]= initial_mass(obj)
            %return the mass of the initial measure in this mode
            [N, P] = size(obj.levels);
            mass_init_mode = zeros(N, P)*mom(obj.vars.x(1));
            mass_sum = 0;

            for n=1:N 
                p=1; %initial measure will only be found at the first partition index
                    mass_curr = obj.levels{n, p}.mass_init();
                    if ~isnumeric(mass_curr)
                        mass_init_mode(n, p) = mass_curr;
                    end
                    mass_sum = mass_sum + mass_curr;
                    % mass_init_mode = mass_init_mode + obj.levels{n, p}.mass_init();
                % end
            end
        end

        function mass_term_mode = terminal_mass(obj)
            %return the mass of the initial measure in this mode
            mass_term_mode = [];
            [N, P] = size(obj.levels);
            mass_term_mode = zeros(N, P)*mom(obj.vars.x(1));

            for n=1:N
                for p = 1:P
                    mass_term_mode(n, p) = obj.levels{n, p}.mass_term();
                    % mass_init_mode = mass_init_mode + obj.levels{n, p}.mass_init();
                end
            end
        end

        function imon = init_monom(obj, d, NTRIG)
            %moments of the initial measure
            %NTRIG: ignore trigonometric variables
            if nargin < 3
                NTRIG = false;
            end
            [N, P] = size(obj.levels);
            imon = cell(N, P);
            for n=1:N
                for p = 1:P
                    if ~isempty(obj.levels{n, p}.init)
                        if NTRIG
                            [~, imon{n, p}] = obj.levels{n, p}.non_trig_monom_init(d);                            
                        else
                            imon{n, p} = obj.levels{n, p}.init.mom_monom(d);
                        end
                    else
                        imon{n, p} = 0;
                    end
                end
            end
        end

        function imon = sel_init_monom(obj, d, ind)
            %moments of the initial measure

            n = length(obj.vars.x);
            if nargin < 2
                ind = 1:n;
            end
            
            [N, P] = size(obj.levels);
            imon = cell(N, P);
            for n=1:N
                for p = 1:P
                    if ~isempty(obj.levels{n, p}.init)                        
                        [~, imon{n, p}] = obj.levels{n, p}.select_monom_init(d, ind);                            
                    else
                        imon{n, p} = 0;
                    end
                end
            end
        end

        function imon = sel_term_monom(obj, d, ind)
            %moments of the terminalmeasure

            n = length(obj.vars.x);
            if nargin < 2
                ind = 1:n;
            end
            
            [N, P] = size(obj.levels);
            imon = cell(N, P);
            for n=1:N
                for p = 1:P
                    if ~isempty(obj.levels{n, p}.term)                        
                        [~, imon{n, p}] = obj.levels{n, p}.select_monom_term(d, ind);                            
                    else
                        imon{n, p} = 0;
                    end
                end
            end
        end

        function imon = sel_occ_monom(obj, d, ind)
            %moments of the occupation measure

            n = length(obj.vars.x);
            if nargin < 2
                ind = 1:n;
            end
            
            [N, P] = size(obj.levels);
            imon = cell(N, P);
            for n=1:N
                for p = 1:P
                    if ~isempty(obj.levels{n, p}.supp.X)                        
                        [~, imon{n, p}] = obj.levels{n, p}.select_monom_occ(d, ind);                            
                    else
                        imon{n, p} = 0;
                    end
                end
            end
        end

        function [trmon, trmon_sum] = trig_occ_monom(obj, d, level_mult)
            %get moments of the occupation measure
            %for the (c, s) marginal
            [N, P] = size(obj.levels);
            trmon = cell(N, P);    
            trmon_sum = 0;
            
            for n=1:N
                for p = 1:P                    
                    %TODO: the trig monom can be reduced by the algebraic
                    %dependence (c^2+s^2=1)
                    [~, tr_base] = obj.levels{n, p}.trig_monom(d);
                    switch obj.Symmetry                                                    
                        case 0
                            tr_curr = tr_base;
                        case 1
                            %half-wave
                            [~, tr_alt] = obj.levels{N-n+1, p}.trig_monom(d, [-1, -1]);
                            tr_curr = (tr_base+ tr_alt);

                        case 2
                            %quarter wave
                            %this is probably bugged. check it.
                            [~, tr_refl] = obj.levels{n, p}.trig_monom(d, [-1, 1]);
                            [~, tr_alt] = obj.levels{N-n+1, p}.trig_monom(d, [-1, -1]);
                            [~, tr_alt_refl] = obj.levels{N-n+1, p}.trig_monom(d, [1, -1]);
                            tr_curr = (tr_base+ tr_alt+ tr_refl+ tr_alt_refl);
                    end
                    % if level_mult
                    %     tr_curr = obj.L(n);
                    % end
                    trmon{n, p} = tr_curr;
                    trmon_sum = trmon_sum + tr_curr;
                end
            end
        end

        function tmon = term_monom(obj, d, NTRIG, flip_load)
            %moments of the terminal measure
                        %NTRIG: ignore trigonometric variables
            if nargin < 3
                NTRIG = false;
            end
            if nargin < 4
                flip_load = false;
            end
            [N, P] = size(obj.levels);
            tmon = cell(N, P);
            for n=1:N
                for p = 1:P                    
                      if isempty(obj.levels{n, p}.term)
                          tmon{n, p} = 0;
                      else
                          if NTRIG
                              [~, tmon{n, p}] = obj.levels{n, p}.non_trig_monom_term(d, flip_load);
                          else    
                            [~, tmon{n, p}] = obj.levels{n, p}.term.mom_monom(d);
                          end
                      end                                                          
                end
            end
        end

        %TODO: three phase balance constraint
        function tmon = mom_sub(obj, vars, vref, level_mult)
            %three-phase balanced symmetry in the current
            
            if nargin < 4
                level_mult = false;
            end
            [N, P] = size(obj.levels);
            tmon = cell(N, P);
            for n=1:N
                for p = 1:P                    
                      if isempty(obj.levels{n, p}.sys{1}.supp)
                          tmon{n, p} = 0;
                      else
                          tmon{n, p} = obj.levels{n, p}.mom_occ_sub(vars, vref);
                      end                                    

                      if level_mult
                          tmon{n, p} = tmon{n, p}*obj.L(n);
                      end

                end
            end


        end

        function obj_min = objective(obj)
            %fetch the objective (THD) from all of the mode locations
            [N, P] = size(obj.levels);
            obj_min = 0;
            % sym_scale = 2^(-2*double(obj.Symmetry));
            sym_scale = 1;
            for n=1:N
                for p = 1:P               
                    [om_curr, ~, ~] = obj.levels{n, p}.objective_con();
                      obj_min = obj_min + sym_scale*om_curr;
                end
            end
        end

        function supp_con_out = supp_con(obj)
            %fetch all support constraints
            [N, P] = size(obj.levels);
            supp_con_out = [];
            for n=1:N
                for p = 1:P                    
                    curr_loc = obj.levels{n, p}.supp_con();
                    if p<P
                        curr_trans = obj.transition{n, p}.supp;
                    else
                        curr_trans = [];
                    end
                    supp_con_out = [supp_con_out; curr_loc; curr_trans];
                end
            end
        end

        %% recovery
        function [m_out, tr_out] = mmat(obj)
            [N, P] = size(obj.levels);
            % obj_min = 0;
            m_out = cell(N, P);
            tr_out = cell(N, P-1);
            for n=1:N
                for p = 1:P               
                    m_out{n, p} = obj.levels{n, p}.mmat();
                    if p < P
                        tr_out{n, p} = obj.transition{n, p}.mmat();
                    end
                end
            end
        end

        function [m_out, tr_out]= mmat_corner(obj)
            [N, P] = size(obj.levels);
            % obj_min = 0;
            m_out = cell(N, P);
            tr_out = cell(N, P-1);
            for n=1:N
                for p = 1:P               
                    m_out{n, p} = obj.levels{n, p}.mmat_corner();
                    if p < P
                        tr_out{n, p} = obj.transition{n, p}.mmat_corner();
                    end
                end
            end
        end


    end

    %% abstract methods (to be overloaded)
    methods (Abstract)
        %constructor
        get_prefix(obj);
        all_dynamics(obj);                      
        refine_support(obj, opts, m,  lsupp_ref)


        %reachable state logic
        cannot_start(opts, n)
        cannot_end(opts, n)
        cannot_reach(opts, m, n)
    end
end

