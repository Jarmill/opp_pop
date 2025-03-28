classdef opp_mode
    %OPP_MODE Measures describing the trajectory at mode m 
    %(before switch m, or at the end of the sequence)
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
        function obj = opp_mode(m, lsupp_ref, objective_mode, opts)
            %OPP_MODE Construct an instance of this class
            %   create locations
            obj.mode = m;
            obj.L = opts.L;
            obj.f0 = opts.f0;
            obj.Symmetry = opts.Symmetry;
            obj.Z_load = opts.Z_load;
            
            N = length(opts.L);
            P = opts.partition;
            obj.levels = cell(N, P);
            obj.transition = cell(N, P-1);

            %define the generic support set
            lsupp_base = lsupp_ref;

            vars = lsupp_base.vars;           
            obj.vars = vars;

            %define the terminal set
            %terminate
            %ignore the trig constraint (beginning) and support arc
            %constraint (end)
            start_pt = [vars.x(1)==1; vars.x(2)==0];
            switch obj.Symmetry
                case 0
                    %full-wave symmetry: end at 2pi
                    stop_pt = start_pt;
                case 1
                    %half-wave symmetry: end at pi
                    stop_pt = [vars.x(1)==-1; vars.x(2)==0];
                case 2
                    %quarter-wave symmetry: end at pi/2
                    Delta_scale = opts.f0*opts.Ts*2^(-double(opts.Symmetry));

                    stop_pt= [vars.x(1)==0; vars.x(2)==1; ...
                        vars.x(3)>=Delta_scale/2];
                    if imag(opts.Z_load)>0 && real(opts.Z_load)==0
                        stop_pt = [stop_pt; vars.x(4)==0];
                    end

                    
            end
            Xstop = [stop_pt; lsupp_base.X(2:end-1)];
            Xstart = [start_pt; lsupp_base.X(2:end-1)];

            mode_end = opts.k/(2^opts.Symmetry);
            if m==0                
               lsupp_base.X_init = Xstart;
            elseif m==mode_end || (opts.early_stop && (mod(m, 2)==0))
                lsupp_base.X_term = Xstop;
            end

            %TODO: define grid-side filter dynamics
            %define the dynamics within the mode

            X_partition = support_partition(opts.partition, vars.x, opts.Symmetry);
            

            f = obj.all_dynamics(vars, opts);

            % objective_mode = obj.all_objective_mode(vars, opts);
            

            % loc_id = N*opts.partition*m + reshape(1:(N*opts.partition), N, []);
            %create locations for each level
            for n = 1:N
                curr_f = f(:, n);
                curr_objective = objective_mode(n, :);                
                for p = 1:P
                    % curr_id = loc_id(n, j);
                    curr_id = sprintf('m%d_n%d_p%d', m, n, p);

                    curr_lsupp = lsupp_base;

                    if p>1 || ((opts.start_level~=0) && (n~= opts.start_level))
                        curr_lsupp.X_init = [];
                    end
                    if p<P || (opts.start_level~=0) && ...
                            ((opts.Symmetry==0) && (n~= opts.start_level) ||...
                            (opts.Symmetry==1) && ((N-n+1)~= opts.start_level))
                        curr_lsupp.X_term = [];
                    end

                    if opts.partition > 1
                        curr_lsupp.X = [curr_lsupp.X; X_partition(p)>=0];
                    end

                    if ~isempty(opts.allowed_levels) 
                        
                        flag = false;
                        if ~opts.allowed_levels(m+1, n)
                            flag = true;
                        elseif m>0
                            can_down = (n<N) && opts.allowed_levels(m, n+1);
                            can_up = (n>1) && opts.allowed_levels(m, n-1);

                            flag = ~(can_up || can_down);
                        end


                        % 
                        % 
                        % && ~opts.allowed_levels(m+1, n) ...
                        % || ((m>0) && (((n>1)&& ~opts.allowed_levels(m, n-1)) || ((n<N)&& ~opts.allowed_levels(m, n+1))))

                        if flag
                            curr_lsupp.X = [];
                            curr_lsupp.X_term = [];
                            curr_lsupp.X_init = [];
                        end
                    end


                    cell_info = struct('mode', m, 'partition', p, 'level', n, 'L', opts.L(n), 'id', curr_id);
                    obj.levels{n, p} = opp_location(curr_lsupp, curr_f, curr_objective, cell_info);

                end
            end

            %create transitions between the partition 
            %when advancing the angle theta
            
            % gtop =    guard(1, vars, loc1, loc2, Xgtop, x);
            for n=1:N
                for p=1:P-1
                    curr_trans_id = sprintf('trans_m%d_n%d_p%d', m, n, p);
                    
                    curr_supp = Xstop;
                    RotAngle = 2*pi/double(P*2^opts.Symmetry);
                    dp = double(p);
                    new_con = (vars.x(1:2)==[cos(dp*RotAngle); sin(dp*RotAngle)]);
                    curr_supp(1:2) = new_con;
                    if ~isempty(opts.allowed_levels) && ~opts.allowed_levels(m+1, n)
                        curr_supp = [];
                    end
                    obj.transition{n, p} = guard(curr_trans_id, lsupp_base.vars, ...
                        obj.levels{n, p}, obj.levels{n, p+1}, curr_supp, vars.x);
                end
            end
        end

        %% functions used in the constructor                
        function f = all_dynamics(obj, vars, opts)
            %create the dynamics as a matrix f
            %row: each level
            %column: each state
            f_trig = 2*pi*[-vars.x(2); vars.x(1)] / (2^double(obj.Symmetry));
            % f_phi = vars.x(3);
            f_clock = 1;

            %TODO: check for symmetry scaling in the load
            f_load = load_dynamics(obj, vars, opts) / (2^double(obj.Symmetry));
            
            % f = [f_trig; f_phi; f_load];
            N = length(opts.L);
            f_basic = [f_trig; f_clock] * ones(1, N);

            f = [f_basic; f_load];
        end


        function f_load = load_dynamics(obj, vars, opts)
            %dynamics in the mode
            %(trig spinning around a circle, clock increasing, load
            %charging/modifying the current)
            %vars: variables (t, x)
            %opts: options from opp_options
            %n: level of the inverter

            %scaled inverter value
            Lscale = 2*opts.L/max(abs(opts.L));
            % u_curr = opts.L(n)/max(abs(opts.L));           

            %dynamics of the load
            if (length(vars.x)==3) || (imag(opts.Z_load) == 0)                      
                %purely resistive
                f_load = [];
            elseif (imag(opts.Z_load) >= 0)
                %inductive load
                %i' = -(R/L)i + (1/L) v
                %per-unit system, ignore the L value
                inductance = imag(opts.Z_load)/(2*pi*opts.f0);
                resistance= real(opts.Z_load);                
                f_load = -(resistance)/(inductance)*vars.x(4) + Lscale;
            else
                 %vc' = (v-vc)/(R*C)
                 %per-unit, ignore (R*C) factor
                 %TODO: v is from the voltage source. Modify when it is 
                 %filtered by a grid-side filter
                 capacitance= -imag(opts.Z_load)*(2*pi*opts.f0);
                 resistance= real(opts.Z_load);  
                 f_load = Lscale - vars.x(4)/(resistance*capacitance);
            end

        end
    
        %% functions used to describe constraints
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


        %fetching moments

        function harm = voltage_harmonics_mom(obj, vars, harm_mon)
            %voltage harmonics constraints
            harm= mom(vars.x(1))*zeros;
            
            [N, P] = size(obj.levels);

            for n=1:N
                for p = 1:P            
                    [~, harm_base] = obj.levels{n, p}.voltage_harmonics_mom(vars, harm_mon);                           
                    
                    switch obj.Symmetry                            
                        case 0
                            harm_mom = harm_base;
                        case 1
                            %half-wave
                            [~, tr_alt] = obj.levels{n, p}.voltage_harmonics_mom(vars, harm_mon, [-1, -1]);
                            harm_mom = (harm_base- tr_alt)*(1/4);
                        case 2
                            %quarter wave
                            [~, tr_refl] = obj.levels{n, p}.voltage_harmonics_mom(vars, harm_mon, [-1, 1]);
                            [~, tr_alt] = obj.levels{n, p}.voltage_harmonics_mom(vars, harm_mon, [-1, -1]);
                            [~, tr_alt_refl] = obj.levels{n, p}.voltage_harmonics_mom(vars, harm_mon, [1, -1]);
                            harm_mom = (harm_base + tr_refl- tr_alt- tr_alt_refl)*(1/16);
                    end
                    
                    harm = harm+harm_mom;                    
                end
            end
        end


        function harm = load_harmonics_mom(obj, vars, harm_mon, harm_in)
            %voltage harmonics constraints
            % harm= mom(p)*0;
            Lmax = max(obj.L);

            % Z_type = 0;
            
            if (length(vars.x)==3) || (imag(obj.Z_load) == 0)                      
                %purely resistive
                % harm = obj.voltage_harmonics_mom(vars, harm_in);
                Z_type = 0;
                Z_scale = 1;
            else                
                if (imag(obj.Z_load) >= 0)
                    
                    %inductive load
                    %i' = -(R/L)i + (1/L) v
                    inductance = imag(obj.Z_load)/(2*pi*obj.f0);
                    Z_type = 1;
                    Z_scale = (Lmax/inductance);
                    % harm = (harm_eval.*vars.x(4)) *  .*(obj.opts.L);
                else
                     %vc' = (v-vc)/(R*C)
                     %per-unit, ignore (R*C) factor
                     %TODO: v is from the voltage source. Modify when it is 
                     %filtered by a grid-side filter                    
                     capacitance= -imag(obj.Z_load)*(2*pi*obj.f0);
                     resistance= real(obj.Z_load);                       
                     RC = resistance*capacitance;
                     Z_type = 2;
                     Z_scale = (Lmax/RC);                     
                end
            end

            for n=1:N
                for p = 1:P            
                    [~, harm_mom] = obj.levels{n, p}.load_harmonics_mom(obj, vars, harm_mon, harm_in, Z_type, Z_scale);       
                    harm = harm+harm_mom;                    
                end
            end
        end

        %TODO: current harmonics

        function [mass_init_mode, mass_sum]= initial_mass(obj)
            %return the mass of the initial measure in this mode
            [N, P] = size(obj.levels);
            mass_init_mode = zeros(N, P)*mom(obj.vars.x(1));
            mass_sum = 0;

            for n=1:N 
                p=1; %initial measure will only be found at the first partition index
                % for p = 1:P
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

        function [trmon, trmon_sum] = trig_occ_monom(obj, d, sym)
            %get moments of the occupation measure
            %for the (c, s) marginal
            if nargin < 3
                sym = 0;
            end
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
                            [~, tr_refl] = obj.levels{n, p}.trig_monom(d, [-1, 1]);
                            [~, tr_alt] = obj.levels{N-n+1, p}.trig_monom(d, [-1, -1]);
                            [~, tr_alt_refl] = obj.levels{N-n+1, p}.trig_monom(d, [1, -1]);
                            tr_curr = (tr_base+ tr_alt+ tr_refl+ tr_alt_refl);
                    end
                    trmon{n, p} = tr_curr;
                    trmon_sum = trmon_sum + tr_curr;
                end
            end
        end

        function tmon = term_monom(obj, d, NTRIG)
            %moments of the terminal measure
                        %NTRIG: ignore trigonometric variables
            if nargin < 3
                NTRIG = false;
            end
            [N, P] = size(obj.levels);
            tmon = cell(N, P);
            for n=1:N
                for p = 1:P                    
                      if isempty(obj.levels{n, p}.term)
                          tmon{n, p} = 0;
                      else
                          if NTRIG
                              [~, tmon{n, p}] = obj.levels{n, p}.non_trig_monom_term(d);
                          else    
                            [~, tmon{n, p}] = obj.levels{n, p}.term.mom_monom(d);
                          end
                      end                                    
                end
            end
        end

        %TODO: three phase balance constraint
        function tmon = mom_sub(obj, vars, vref)
            %three-phase balanced symmetry in the current
            
            [N, P] = size(obj.levels);
            tmon = cell(N, P);
            for n=1:N
                for p = 1:P                    
                      if isempty(obj.levels{n, p}.sys)
                          tmon{n, p} = 0;
                      else
                          tmon{n, p} = obj.levels{n, p}.mom_occ_sub(vars, vref);
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
end

