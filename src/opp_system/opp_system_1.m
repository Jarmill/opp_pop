classdef opp_system_1 < opp_system_interface
    %OPP_SYSTEM_1 Storage for all single-phase terms
    %   includes the N-switching constraints
    
    % properties
    %     Property1
    % end
    
    methods
        %% constructor
        function obj = opp_system_1(opts)
            %OPP_SYSTEM_1 Construct an instance of this class
            %   Detailed explanation goes here
            obj@opp_system_interface(opts)
        end
        
        function [vars] = create_vars(obj, opts)
              %declare the variables
            load_state = imag(opts.Z_load)~=0;
            mpol('c', 1, 1);
            mpol('s', 1, 1)
            mpol('phi', 1, 1)
            mpol('I', 1, 1)
            x = [c; s; phi];
            if load_state
                x = [x; I];
            end

            if opts.TIME_INDEP
               t = [];
            else
                mpol('t', 1, 1)
                % t = t;
            end
            vars = struct('t', t, 'x', x);
        end
       
        function [vars, jumps, modes] = create_system(obj, opts)
            %used in the constructor

            k = opts.k/(2^opts.Symmetry);
            jumps = cell(k, 1);
            modes = cell(k+1, 1);

          
            %x = [c; s; phi; l] -> [cos(theta), sin(theta), clock, load
            %state (current of load inductor/voltage of load capacitor)

            
            
            %create the basic location structure
            if ~opts.TIME_INDEP
                mpol('t', 1, 1)
            else
                t = [];
            end

            %create the basic support set
            lsupp_base = loc_support();
            lsupp_base.vars.x = x;
            lsupp_base.vars.t = t;
            vars = lsupp_base.vars;

            lsupp_base.TIME_INDEP = opts.TIME_INDEP;
            lsupp_base.FREE_TERM = 0;
            lsupp_base.Tmax = double(2^opts.Symmetry);

            %create the support set

            Theta = opts.f0*opts.Ts;

            %BUG HERE BUG HERE BUG HERE
            %X_trig = 1-x(1)^2 + x(2)^2;
            X_trig = 1-x(1)^2 - x(2)^2;

            %clock and rescaled load
            if length(x)<4
                X_load = [];
            else
                X_load = 1-x(4)^2;                
            end

            Theta_scale = Theta*2^(double(obj.opts.Symmetry));

            X_clock_mode = x(3)*(1-2*Theta_scale- x(3));    
            X_clock_jump = (x(3)-Theta_scale)*(1-2*Theta_scale- x(3));    
            X = [X_trig==0; X_clock_mode>=0; X_load>=0];             
            X_jump = [X_trig==0; X_clock_jump>=0; X_load>=0];
            
  
            lsupp_base.X = X;
            %define the reset law for the jump
            

            %declare the jumps and modes
            objective_mode = obj.objective_level(vars, opts);
            
            for m=0:k
                lsupp_curr = lsupp_base;
                arc_curr = support_arc(m, x, Theta, opts.Symmetry);
                lsupp_curr.X = [lsupp_curr.X; arc_curr>=0];
                modes{m+1} = opp_mode(m, lsupp_curr, objective_mode, opts);

                if m~=0
                    jumps{m} = opp_jump(m, opts, vars, X_jump);
                end
            end
        end

        %% constraints
        function [mom_con, supp_con_one] = cons(obj, d)

            %generate the constraints
            supp_con_one = obj.supp_con();
            
              

            %dynamics constraints
            %mass of initial measure = 1
            con_prob = obj.con_prob_dist();
             
            %initial = sum of terminal measure
            con_preserve = obj.con_return(d);

            %flow +jump continuity constraints
            con_liou = obj.con_flow(d);
            
            %trig is uniformly distributed over circle
            con_leb = obj.con_uni_circ(d);
            

            %harmonics constraints
            [con_harm, ~] = obj.con_harmonics();                           

            %without harmonics
            % mom_con = [con_prob; con_preserve; con_leb; con_threephase; con_harm];
            
            %only dynamics
            % mom_con = [con_prob; con_preserve; con_liou];

            %without dynamics
            % mom_con = [con_prob; con_leb;
            %     con_floating;
            %     con_harm;
            %     con_liou; con_preserve
            %     ];

            %with harmonics and dynamics
            mom_con = [con_prob; con_leb; %fixed probability/arc measures                
                con_harm; %harmonics constraints                
                con_preserve;  con_liou %flow
                ];                  

        end
  

        function flow_con = con_flow(obj, d)
            %the flow conservation constraint (the big one)

            %start the storage structure 
            Nmodes = length(obj.modes);
            liou_cell = cell(Nmodes, 1);
            jump_src = cell(Nmodes-1, 1);
            jump_dst = cell(Nmodes-1, 1);

            %compute all terms
            for m=1:(Nmodes)
                liou_cell{m} = obj.modes{m}.flow(d);
            end

            %add the jump to the cell terms
            for m=1:(Nmodes-1)                
                [jump_src{m}, jump_dst{m}] = obj.jumps{m}.liou_reset(d);
            end

            [N, P] = size(liou_cell{1});
            
            %iterate over all cells
            flow_con = [];
            flow_con_cell = cell(Nmodes, N, P);
            for m = 1:Nmodes
                liou_curr = liou_cell{m};
                for n=1:N
                    if isempty(obj.opts.allowed_levels) || obj.opts.allowed_levels(m, n)
                        for p = 1:P
                            if m==1
                                dst_curr = 0;
                            else
                                dst_curr = jump_dst{m-1}{n, p};
                            end
                            if m==Nmodes
                                src_curr = 0;
                            else
                                src_curr = jump_src{m}{n, p};
                            end
                            
                            flow_con_cell{m, n, p} = liou_curr{n, p} + src_curr + dst_curr==0;
                            
                            %stack them into a giant vector: flow_con
                            flow_con = [flow_con; flow_con_cell{m, n, p}];
                        end
                    end
                end
            end

            
        end

        function mass_con_eq = con_prob_dist(obj)
            %initial measure is a probability distribution (mass 1)
            
            [~, mass_init_sum] = obj.modes{1}.initial_mass();
        
            mass_con_eq = (mass_init_sum==1);
           
        end

        function return_con = con_return(obj, d)
            %conservation of position between the initial and final measure
            
            % mass_con = obj.modes{1}.mass_init_mode();
            init_monom = obj.modes{1}.init_monom(d, true);

            N = length(obj.opts.L);
            %TODO: 
            %This is for full-wave. generalize for other symmetry
            %structures

            if obj.opts.Symmetry == 2
                return_con = [];
            else
                return_mom = {init_monom{:, 1}};

                %index the terminal destination levels based on the
                %applied symmetry
                %unconstrained for quarter-wave symmetry (here at least)
                flip_load = obj.opts.Symmetry==1;                
                if flip_load
                    stop_order = N:-1:1;   
                    stop_range = 2:length(obj.modes);
                else
                    stop_order = 1:N;
                    stop_range = 3:2:length(obj.modes);
                end

                

                if obj.opts.early_stop                
                    for m = stop_range
                        % mass_con = mass_con - obj.modes{m}.mass_term_mode();
                        stop_monom = obj.modes{m}.term_monom(d, true, flip_load);
                        return_mom = madd_cell_mom(return_mom, {stop_monom{stop_order, end}}, -1);
                    end
                else
                    % mass_con = mass_con - obj.modes{end}.mass_term_mode();
                    stop_monom = obj.modes{end}.term_monom(d, true, flip_load);
                    return_mom = madd_cell_mom(return_mom, {stop_monom{stop_order, end}}, -1);
                end
    
                [N, P] = size(return_mom);
                return_con = [];
                for n = 1:N
                    for p = 1:P
                        if ~isnumeric(return_mom{n, p})
                            return_con = [return_con; return_mom{n, p}==0];
                        end
                    end
                end
            end

            % return_con = (return_mom==0);
        
        end

        function supp_con_all = supp_con(obj)
            %fetch support constraints from the model
            supp_con_all = [];

            for i = 1:length(obj.jumps)
                supp_con_all = [supp_con_all; obj.jumps{i}.supp_con()];
            end%cell of opp_switch
            for i = 1:length(obj.modes)
                supp_con_all = [supp_con_all; obj.modes{i}.supp_con()];
            end

            supp_con_all = [supp_con_all; obj.sys3.supp_con()];
          %cell of opp_mode(), contains initial/terminal/occupation measures           
         end

        %% harmonics
        function [harm_con, harm_source] = con_harmonics(obj)
            %collect harmonics constraints on the voltage source and the load
            harm = obj.harm_eval(obj.vars, obj.opts.harmonics);
            harm_load_data = obj.harm_eval(obj.vars, obj.opts.harmonics_load);
            
            % harmonics on voltage source
            if ~isempty(harm)
                harm_source = 0;
                for m = 1:length(obj.modes)   
                    harm_mom = obj.modes{m}.voltage_harmonics_mom(obj.vars, harm);
                    harm_source = harm_source + harm_mom;
                end
                harm_source_con = harmonics_process(obj.opts.harmonics, harm_source);
            else
                harm_source_con = [];                
            end


            %TODO: grid side filters

            %TODO: symmetries with harmonics on load side
            %harmonics on the load side
            % if ~isempty(harm_load_data)          
            %     harm_load = 0;
            %     for m = 0:length(obj.modes)     
            %         harm_mom_load = obj.modes{m}.load_harmonics_mom(obj.vars, harm_load, obj.opts.harmonics_load);
            %         harm_load = harm_load + harm_mom_load;
            %     end
            %     harm_load_con = harmonics_process(obj.opts.harmonics_load, harm_load);
            % else
                harm_load_con = [];
            % end
            
            harm_con = [harm_load_con; harm_source_con];

            %TODO: finish this. Sum up the harmonics over all components
        end

        function harm_monom = harm_eval(obj, vars, harm_in)            
            % [vars.x(1).^harm_in.index_cos; vars.x(2).^harm_in.index_sin]/pi;

            %compute the chebyshev moments. then divide by pi
            %cos(n theta) = T_n(cos(theta))
            %sin(n theta) = sin(theta) U_{n-1}(cos(theta))

            %TODO: modify this for symmetry
            
            c = vars.x(1);
            s = vars.x(2);

            %the *1, *2, *4 multiplications are already taken into 
            %account by the mass of the occupation measure

            if isempty(harm_in)
                harm_monom = [];
            else
                if ~isempty(harm_in.index_cos)
                    %chebyshev of the first kind
                    c_ind_max = max(max(harm_in.index_cos), 1);
    
                    T = zeros(c_ind_max+1, 1)*c;
                    T(1) = 1+0*c;
                    T(2) = c;
                    for p = 2:c_ind_max
                        T(p+1) = 2*c*T(p) - T(p-1);
                    end

                    switch obj.opts.Symmetry
                        case 0
                            cos_scale = 1;                        
                        case 1
                            % cos_scale = 2*mod(harm_in.index_cos, 2);
                            cos_scale = mod(harm_in.index_cos, 2);
                        case 2
                            cos_scale = 0;
                    end
    
                    harm_cos = cos_scale.*T(harm_in.index_cos+1);
                else
                    harm_cos = [];
                end
                
                if ~isempty(harm_in.index_sin)
                    %chebyshev of the second kind
                    s_ind_max = max(max(harm_in.index_sin), 1);
    
                    U = zeros(s_ind_max, 1)*c;
                    U(1) = 1+0*c;
                    U(2) = 2*c;
                    for p = 2:s_ind_max
                        U(p+1) = 2*c*U(p) - U(p-1);
                    end
    
                    switch obj.opts.Symmetry
                        case 0
                            sin_scale = 1;                        
                        case 1
                            % sin_scale = 2*mod(harm_in.index_sin, 2);
                            sin_scale = mod(harm_in.index_sin, 2);
                        case 2
                            % sin_scale= 4*mod(harm_in.index_sin, 2);
                            sin_scale = mod(harm_in.index_sin, 2);
                    end

                    harm_sin = sin_scale .* (s*U(harm_in.index_sin));
                else
                    harm_sin = [];
                end
    
                %package up the harmonics
                % harm_monom = [harm_cos; harm_sin]/pi;    

                %divide by pi to get the harmonic
                %then multiply by 2pi because time is scaled to [0, 1]?
                %figure this out
                harm_monom = [harm_cos; harm_sin]*2;
                % harm_monom = [harm_cos; harm_sin]/pi;

                %process the symmetry
                % harm_monom = obj.symmetry_eval(harm_monom, [c; s]);
            end
        end

        %% objective
        function objective = objective_level(obj, vars, opts)
            %return the mode-objective at each level
            %
            %TODO: 
            %the scaling factors may be wrong for inductance and
            %capacitance. check this later.
            %
            N = length(opts.L);
            Lmax = max(abs(opts.L));
            % sym_factor = double(2^opts.Symmetry);
            sym_factor = 1;
            %Another TODO: quarter-wave symmetry may break the
             %characterization of the current in the inductor/capacitor
            if (length(vars.x)==3) || (imag(opts.Z_load) == 0)                      
                %purely resistive
                %sym_factor: replicate the square according to the symmetry
                %2*pi: because time is normalized to [0, 1]
                objective = (2*pi)*(opts.L.^2);
            elseif (imag(opts.Z_load) >= 0)
                
                if real(opts.Z_load)==0
                    objective = pi^2 * (2*pi)^2*vars.x(4)^2*ones(size(opts.L));
                else
                    %inductive load
                    %i' = -(R/L)i + (1/L) v
                    %per-unit system, ignore the L value
                    inductance = imag(opts.Z_load)/(2*pi*opts.f0);
                    % resistance= real(opts.Z_load);                
                    % f_load = -(resistance)/(inductance)*vars.x(4) + Lscale;
                    objective = (2*pi)*vars.x(4)^2*(Lmax/inductance)^2*(opts.L).^2;
                end
            else
                 %vc' = (v-vc)/(R*C)
                 %per-unit, ignore (R*C) factor
                 %TODO: v is from the voltage source. Modify when it is 
                 %filtered by a grid-side filter
                 %
                 %                 
                 %
                 capacitance= -imag(opts.Z_load)*(2*pi*opts.f0);
                 resistance= real(opts.Z_load);  
                 RC = resistance*capacitance;
                 % f_load = Lscale - vars.x(4)/(resistance*capacitance);
                 objective = (opts.L.^2) + 2*(opts.L)*vars.x(4)*(Lmax/RC) + (Lmax/RC)^2*vars.x(4)^2;
            end
            objective = objective'*sym_factor/(2*pi);

        end
   

        function objective_out = objective(obj)
            %objective in the single-phase setting
            objective_out = 0;
            for i = 1:length(obj.sys1.modes)
                objective_out = objective_mode + obj.modes{i}.objective();
            end
        end      

        %% other helper functions
        function mom_I = current_mom(obj, d)
            %get the marginals of (c, s, I) from the occupation measure
            
            vars_inv = obj.vars.x([1, 2, 4]);
            p_in = mmon(vars_inv, d);
            p_in_sym = obj.symmetry_eval_current(p_in, vars_inv);

            bmom = 0;
                %dispatch into the measures
                for m = 1:length(obj.modes)
                    curr_mom = obj.modes{m}.mom_sub(obj.get_vars(), p_in_sym);
                    bmom = madd_cell_mom(bmom, curr_mom, 1);
                end

                %collect all differential moments together
                % bmom_all = 0*mom(mon_3);
                mom_I = bmom(:, 1);
                for i =1:size(bmom, 1)
                    for j = 2:size(bmom, 2)
                        mom_I{i, j} = mom_I{i, j} +  bmom{i, j};
                    end
                end
        end

    function bmom_all = three_phase_current_mom(obj, d)
            %get moments of the differntial-mode current

            % K = sqrt(2/3)*[1 -0.5 -0.5;
            % 0 sqrt(3)/2 -sqrt(3)/2] /sqrt(3);


            % Kt = [K; ones(1, 3)/3];

            vars_inv = obj.vars.x([1, 2, 4]);
            % p_in = mmon(vars_inv(1:2), d-1)*vars_inv(3);
            p_in = mmon(vars_inv, d);

            %TODO: Debug this
            p_in_sym = obj.symmetry_eval_current(p_in, vars_inv);

            mon_3 = obj.three_phase_rotate(p_in_sym, vars_inv);


            %TODO: testing only (do the whole current)
            bmom = 0;
            %dispatch into the measures
            for m = 1:length(obj.modes)
                curr_mom = obj.modes{m}.mom_sub(obj.get_vars(), mon_3);
                bmom = madd_cell_mom(bmom, curr_mom, 1);
            end

            %collect all differential moments together
            bmom_all = 0*mom(mon_3);
            for i =1:size(bmom, 1)
                for j = 1:size(bmom, 2)
                    bmom_all = bmom_all +  bmom{i, j};
                end
            end
        end        

        %% recovery
    end
end

