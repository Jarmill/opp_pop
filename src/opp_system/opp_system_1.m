classdef opp_system_1 < opp_system_interface
    %OPP_SYSTEM_1 Storage for all single-phase terms
    %   includes the N-switching constraints
    
    % properties
    %     name = '';
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
            mpol(['c', opts.name], 1, 1);
            mpol(['s', opts.name], 1, 1);

            c = eval(['c', opts.name]);
            s = eval(['s', opts.name]);

            if opts.clock
                mpol(['phi', opts.name], 1, 1);
                phi = eval(['phi', opts.name]);
            else
                phi = [];
            end

            mpol(['I', opts.name], 1, 1)
            I = eval(['I', opts.name]);
            
            x = [c; s; phi];
            if load_state
                x = [x; I];
            end

            if opts.TIME_INDEP
               t = [];
            else
                mpol(['t', opts.name], 1, 1)
                t = eval(['t', opts.name]);
                % t = t;
            end
            vars = struct('t', t, 'x', x);
        end
       
        function [modes, jumps, opts] = create_system(obj, opts)
            %used in the constructor

            k = opts.k/(2^opts.Symmetry);
            jumps = cell(k, 1);
            modes = cell(k+1, 1);

          
            
            %x = [c; s; phi; l] -> [cos(theta), sin(theta), clock, load
            %state (current of load inductor/voltage of load capacitor)

            
            
            %create the basic location structure

            x = obj.vars.x;
            t = obj.vars.t;

            %create the basic support set
            
            lsupp_base = loc_support(obj.vars);
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
            % X_trig = 1-x(1)^2 - x(2)^2;

            %clock and rescaled load
            % X = (X_trig==0); 

            %GROEBNER REDUCTION !!!!!!!
            X = (x(2)^2 == 1 - x(1)^2);
            

            I_ind = 3 + opts.clock;
            if imag(opts.Z_load)==0
                X_load = [];
            else
                %LOAD LIMIT
                X_load = 1-x(I_ind)^2;                
            end

            X = [X; X_load >= 0];
            Theta_scale = Theta*2^(double(opts.Symmetry));

            if opts.clock
                X_clock_mode = x(3)*(1-2*Theta_scale- x(3)) >= 0;    
                X_clock_jump = (x(3)-Theta_scale)*(1-2*Theta_scale- x(3)) >= 0; 
                X = [X; X_clock_mode];  
                X_jump = [X; X_clock_jump];
            else
                X_jump = X;
            end
                       
            
            
  
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

            con_dwell = obj.con_dwell_time();
            

            %harmonics constraints
            [con_harm, ~] = obj.con_harmonics();                           

            
            %quarter-matching constraints
            con_match = obj.con_quarter_match(d);


            con_power_dissipation = obj.con_power_loss();

            con_power_match = obj.power_match_con(d);
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
                con_match; %quarter-matching
                con_dwell; %soft dwell-time constraint
                con_power_dissipation; % power losses
                con_power_match; %power matching in measures
                ];                  

        end

        function [mom_con, supp_con_one] = cons_limited(obj, d)

            %generate the constraints
            supp_con_one = obj.supp_con();                          

            %initial = sum of terminal measure
            con_preserve = obj.con_return(d);

            %flow +jump continuity constraints
            con_liou = obj.con_flow(d);
            

            mom_con = [con_preserve; con_liou];              

        end


        
  
  

        function flow_con = con_flow(obj, d)
            %the flow conservation constraint (the big one)

            %start the storage structure 
            Nmodes = length(obj.mode);
            liou_cell = cell(Nmodes, 1);
            jump_src = cell(Nmodes-1, 1);
            jump_dst = cell(Nmodes-1, 1);

            %compute all terms
            for m=1:(Nmodes)
                liou_cell{m} = obj.mode{m}.flow(d);
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
                            if ~islogical(flow_con_cell{m, n, p})                               
                                flow_con = [flow_con; flow_con_cell{m, n, p}];
                            end
                        end
                    end
                end
            end

            
        end

        function mass_con_eq = con_prob_dist(obj)
            %initial measure is a probability distribution (mass 1)
            
            [~, mass_init_sum] = obj.mode{1}.initial_mass();
        
            mass_con_eq = (mass_init_sum==1);
           
        end

        function return_con = con_return(obj, d)
            %conservation of position between the initial and final measure
            
            % mass_con = obj.mode{1}.mass_init_mode();
            init_monom = obj.mode{1}.init_monom(d, true);

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
                    stop_range = 2:length(obj.mode);
                else
                    stop_order = 1:N;
                    stop_range = 3:2:length(obj.mode);
                end

                

                if obj.opts.early_stop                
                    for m = stop_range
                        % mass_con = mass_con - obj.mode{m}.mass_term_mode();
                        stop_monom = obj.mode{m}.term_monom(d, true, flip_load);
                        return_mom = madd_cell_mom(return_mom, {stop_monom{stop_order, end}}, -1);
                    end
                else
                    % mass_con = mass_con - obj.mode{end}.mass_term_mode();
                    stop_monom = obj.mode{end}.term_monom(d, true, flip_load);
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

        function con_match = con_quarter_match(obj, d)
                        
            %only use when there are an even number of jumps
            con_match = [];
            if (obj.opts.Symmetry==1) && (obj.opts.quarter_match)
            % if false
                %only match quarters 
                
                kh = length(obj.jumps);
                N = length(obj.opts.L);
                % jump_up_before = cell(N, kh);
                % jump_down_before = cell(N, kh);
                % 
                % jump_up_after = cell(N, kh);
                % jump_down_after = cell(N, kh);
                
                %moments for the occupation measures
                % for i = 1:((kh/2))
                %     oub = obj.mode{i}.trig_occ_sign_monom(d, [1; 1]);
                %     oua = obj.mode{kh - i + 1}.trig_occ_sign_monom(d, [-1; 1]);
                %     con_curr = [];
                %     for n = 1:(N)
                %         if ~isnumeric(oub{n})
                %             con_curr = [con_curr; oub{n}-oua{n}==0];
                %         end                                               
                %     end
                %     con_match = [con_match; con_curr];
                % end

                %get moments for the jumps before pi
                for i = 1:(kh/2)
                    [jub, jdb] = obj.jumps{i}.sel_monom_jump_summarize(d, [1; 2], [1; 1]);
                    % jump_up_before(:, i) = jub;                        
                    % jump_down_before(:, i) = jdb;

                    [jua, jda] = obj.jumps{kh-i+1}.sel_monom_jump_summarize(d, [1; 2], [-1; 1]);
                    % jump_up_after(:, kh - i + 1) = jua;                        
                    % jump_down_after(:, kh - i + 1) = jda;

                    con_curr = [];
                    for n = 1:(N-1)
                        if ~isnumeric(jub{n})
                            con_curr = [con_curr; jub{n}==jda{n}];
                        end
                        if ~isnumeric(jdb{n})
                            con_curr = [con_curr;  jdb{n} == jua{n}];
                        end                                                
                    end
                    
                    con_match = [con_match; con_curr];
                end

                %get flipped moments for the jumps before pi

                

                % jump_up = 

            end
        end

        

        function supp_con_all = supp_con(obj)
            %fetch support constraints from the model
            supp_con_all = [];

            for i = 1:length(obj.jumps)
                supp_con_all = [supp_con_all; obj.jumps{i}.supp_con()];
            end%cell of opp_switch
            for i = 1:length(obj.mode)
                supp_con_all = [supp_con_all; obj.mode{i}.supp_con()];
            end           
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
                for m = 1:length(obj.mode)   
                    harm_mom = obj.mode{m}.voltage_harmonics_mom(obj.vars, harm);
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
            %     for m = 0:length(obj.mode)     
            %         harm_mom_load = obj.mode{m}.load_harmonics_mom(obj.vars, harm_load, obj.opts.harmonics_load);
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


        function con_dwell = con_dwell_time(obj)
            con_dwell = [];
            if ~obj.opts.clock && (obj.opts.Ts > 0)
                for i = 1:length(obj.mode)
                    Theta_scale = obj.opts.f0*obj.opts.Ts*2^(double(obj.opts.Symmetry));        
                    con_dwell = [con_dwell; obj.mode{i}.occ_mass() >= Theta_scale];            
                end
            end
        end 



        function con_power = con_power_loss(obj)
            %power dissipation constraint
            %bounded power budget per fundamental period
            
            % [~, mass_init_sum] = obj.mode{1}.initial_mass();
            
            if obj.opts.power_budget == Inf
                con_power = [];
            else
                power_use = obj.power_dissipated(obj.opts.dispatch);
                
                con_power = [power_use <= obj.opts.power_budget];
            end

           
        end



        function con_power_match = power_match_con(obj, d)
            %how much power is dissipated over the fundamental period?
            %accumulate the conduction losses and the switching losses

            %conduction losses
            con_power_match = [];
             if obj.opts.power_budget < Inf
               
                for m = 1:length(obj.mode)   
                    con_power_match= [con_power_match; obj.mode{m}.power_match_con(d)];                                
                end
                
                for i = 1:length(obj.jumps)
                    con_power_match= [con_power_match; obj.jumps{i}.power_match_con(d)];                  
                end  
             end
        end

        function power_use = power_dissipated(obj, dispatch)
            %how much power is dissipated over the fundamental period?
            %accumulate the conduction losses and the switching losses
            if nargin == 1
                dispatch = obj.opts.dispatch;
            end

            %conduction losses
            power_cond = 0;
            for m = 1:length(obj.mode)   
                power_curr= obj.mode{m}.power_dissipated(dispatch);                
                power_cond = power_cond + power_curr;
            end

            %scale to [0, 2pi] time range
            power_use = power_cond * 2*pi * (2^(-double(obj.opts.Symmetry)));

            %switching losses
            M = length(obj.jumps);
            for i = 1:M
                power_curr = obj.jumps{i}.switching_losses(dispatch);
                power_use = power_use + power_curr;                
            end  

            %symmetry: increase power losses
            %QW: already doing quarter matching if RL ~= 0
            %and also QW would require extra trackign anyways
            if obj.opts.Symmetry > 0
                power_use = power_use * 2;
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
            I_ind = 3+ obj.opts.clock;
            %Another TODO: quarter-wave symmetry may break the
             %characterization of the current in the inductor/capacitor
            if (length(vars.x)==3) || (imag(opts.Z_load) == 0)                      
                %purely resistive
                %sym_factor: replicate the square according to the symmetry
                %2*pi: because time is normalized to [0, 1]
                objective = (2*pi)*(opts.L.^2);
            elseif (imag(opts.Z_load) >= 0)
                
                % if real(opts.Z_load)==0
                objective_pulse = pi^2 *vars.x(I_ind)^2*ones(size(opts.L));
                
                if opts.uext ~= 0
                    inductance = imag(opts.Z_load)/(2*pi*opts.f0);
                    resistance= real(opts.Z_load); 
                    tau = resistance/inductance;               
                    trig_wave = [real(opts.uext), imag(opts.uext)]*vars.x(1:2)/sqrt(tau^2 + 1);
                    objective_mix = trig_wave *  (pi * vars.x(I_ind)*ones(size(opts.L)));                    
                    objective_wave = pi*abs(opts.uext)/(tau^2+1) / (2*pi)^2;
                else
                    objective_mix =0;
                    objective_wave = 0;
                end


                    objective =  (2*pi)^2*(objective_pulse - 2*objective_mix + objective_wave);
                % else
                    %inductive load
                    %i' = -(R/L)i + (1/L) v
                    %per-unit system, ignore the L value
                    % inductance = imag(opts.Z_load)/(2*pi*opts.f0);
                    % resistance= real(opts.Z_load);                
                    % f_load = -(resistance)/(inductance)*vars.x(4) + Lscale;
                    % objective = (2*pi)*vars.x(4)^2*(Lmax/inductance)^2*(opts.L).^2;
                % end
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
                 objective = (opts.L.^2) + 2*(opts.L)*vars.x(I_ind)*(Lmax/RC) + (Lmax/RC)^2*vars.x(I_ind)^2;
            end
            objective = objective'*sym_factor/(2*pi);

        end
       
        %% other helper functions

        function [sel_up, sel_down] = sel_jump_monom(obj, d, ind)
            %used for alignment

            %get moments for all jump measures
            N = length(obj.opts.L);
            sel_up   = cell(N-1, 1);
            sel_down = cell(N-1, 1);           

            M = length(obj.jumps);

            for i = 1:M
                [curr_up, curr_down] = obj.jumps{i}.sel_monom_jump(d, ind);
                sel_up = madd_cell_mom(sel_up,curr_up, 1);
                sel_down = madd_cell_mom(sel_down,curr_down, 1);
            end                        
        end


        function mom_I = current_mom(obj, d)
            %get the marginals of (c, s, I) from the occupation measure
            
            vars_inv = obj.vars.x([1, 2, 4]);
            p_in = mmon(vars_inv, d);
            p_in_sym = obj.symmetry_eval_current(p_in, vars_inv);

            bmom = 0;
                %dispatch into the measures
                for m = 1:length(obj.mode)
                    curr_mom = obj.mode{m}.mom_sub(obj.get_vars(), p_in_sym);
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
            for m = 1:length(obj.mode)
                curr_mom = obj.mode{m}.mom_sub(obj.get_vars(), mon_3);
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

        function [m_out] = mmat(obj)
            %get the moment matrix of all measure variables
            m_out = struct;
            K = length(obj.mode);
            % [N, P] = size(obj.mode{1}.levels);
            % m_out.levels = cell(K, N, P);
            m_out.modes = cell(K, 1);
            m_out.transition = cell(K, 1);

            for i = 1:K
                [m_out.modes{i}, m_out.transition{i}] = obj.mode{i}.mmat();
            end

            m_out.jump = cell(K-1, 1);
            for i=1:(K-1)
                m_out.jump{i} = obj.jumps{i}.mmat();
            end            
        end
   
        function [m_out] = mmat_corner(obj)
            %get the moment matrix of all measure variables
            m_out = struct;
            K = length(obj.mode);
            % [N, P] = size(obj.mode{1}.levels);
            % m_out.levels = cell(K, N, P);
            m_out.modes = cell(K, 1);
            m_out.transition = cell(K, 1);

            for i = 1:K
               [m_out.modes{i}, m_out.transition{i}] = obj.mode{i}.mmat_corner();
            end

            m_out.jump = cell(K-1, 1);
            for i=1:(K-1)
                m_out.jump{i} = obj.jumps{i}.mmat_corner();
            end
        end
    
        function cs = current_square_summary(obj)
            % cs = [];
            K = length(obj.mode);
            [N, P] = size(obj.mode{1}.levels);

            cs = zeros(K, N);
            for m=1:K
                for n=1:N
                    for p = 1:P
                        if ~isempty(obj.mode{m}.levels{n, p}.sys{1}.meas_occ.supp)
                            currsq = obj.mode{m}.levels{n, p}.sys{1}.meas_occ.vars.x(end)^2;
                            cs(m, n) = cs(m, n) +  double(mom(currsq));
                        end
                    end
                end
            end
            cscale = pi^3 * 2^double(obj.opts.Symmetry);
            cs = cs * cscale;
        end

        function ms = mass_summary(obj)
            %collect the masses of the occupation measure into a neat array
            ms = struct;
            K = length(obj.mode);
            [N, P] = size(obj.mode{1}.levels);

            ms.mode = zeros(K, N, P);
            ms.trans = zeros(K, N, P-1);
            for m=1:K
                for n=1:N
                    for p = 1:P
                        ms.mode(m, n, p) = double(obj.mode{m}.levels{n, p}.sys{1}.meas_occ.mass());
                        if p < P
                            ms.trans(m, n, p) = double(obj.mode{m}.transition{n, p}.mass());
                        end
                    end
                end
            end

            ms.jump_up = zeros(K-1, N-1, P);
            ms.jump_down = zeros(K-1, N-1, P);
            for m=1:K-1
                for n=1:N-1
                    for p = 1:P
                        ms.jump_up(m, n, p) = double(obj.jumps{m}.jump_up{n, p}.mass());
                        ms.jump_down(m, n, p) = double(obj.jumps{m}.jump_down{n, p}.mass());
                    end
                end
            end

            [c, mom_harm] = obj.con_harmonics();
            ms.harm = double(mom_harm);
        end
        function [load, load_candidate] = recover_load(obj)
            %
            %load:              initial current in the load
            %load_candidate:    load current upon entering the state
            %leaving the state
                Mc = obj.mmat_corner();
                ms = obj.mass_summary();
                [N, P] = size(obj.mode{1}.levels);
                Nmodes = length(obj.mode);
                load_candidate = zeros(Nmodes+1, N)*NaN;
 
                %get the initial current
                I_ind = 4+obj.opts.clock;
                for n =1:N
                    Mcurr= obj.mode{1}.levels{n, 1}.mmat_corner();
                    init_curr = Mcurr.init;
                    if ~isempty(init_curr) && (init_curr(1, 1) > 0.99)
                        load_candidate(1, n) = init_curr(I_ind, 1)/init_curr(1, 1);
                    end
                end
                %track along the jumps
                for m = 1:(Nmodes-1)
                    for n = 1:N-1
                       
                        for p = 1:P
                            if ms.jump_up(m, n, p) > 0.99
                                jump_curr = Mc.jump{m}.up{n, p};
                                load_candidate(m+1, n+1) = jump_curr(I_ind, 1)/jump_curr(1, 1);
                            elseif ms.jump_down(m, n, p) > 0.99
                                jump_curr = Mc.jump{m}.down{n, p};
                                load_candidate(m+1, n) = jump_curr(I_ind, 1)/jump_curr(1, 1);
                            end
                        end
                    end
                end
    
                %track the exit
                %get the initial current
                    for n =1:N
                        Mcurr= obj.mode{end}.levels{n, end}.mmat_corner();
                        term_curr = Mcurr.term;
                        if ~isempty(term_curr) && (term_curr(1, 1) > 0.99)
                            load_candidate(end, n) = term_curr(I_ind, 1)/term_curr(1, 1);
                        end
                    end
                load =load_candidate(1, ~isnan(load_candidate(1, :)));
    
    
            end
    
    end
end

