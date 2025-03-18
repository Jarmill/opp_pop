classdef opp_manager 
    %OPP_MANAGER Manager for the MATLAB optimal pulse pattern synthesis
    %task
    %
    %based on meas_class/manager_interface, but not totally the same
    %
    %   Detailed explanation goes here
    
    properties
        opts;   %
        jumps;  %cell of opp_switch
        modes;  %cell of opp_mode(), contains initial/terminal/occupation measures

        sdp_settings;
        vars;
    end
    
    methods
        function obj = opp_manager(opts) 
            %OPP_MANAGER Construct an instance of this class
            %   Detailed explanation goes here
            %
            %opts: opp_options structure
            obj.opts = opts;
            obj.sdp_settings = sdpsettings('solver', opts.solver, 'mosek.MSK_DPAR_BASIS_TOL_S', 1e-8, ...
                'mosek.MSK_DPAR_BASIS_TOL_X', 1e-8, 'mosek.MSK_DPAR_INTPNT_CO_TOL_MU_RED', 1e-9, ...
                'mosek.MSK_DPAR_INTPNT_TOL_PATH', 1e-6);            

            [obj.vars, obj.jumps, obj.modes] = obj.create_system(opts);
        end

        %% construct everything
        function [vars, jumps, modes] = create_system(obj, opts)
            %used in the constructor

            k = opts.k/(2^opts.Symmetry);
            jumps = cell(k, 1);
            modes = cell(k+1, 1);

            %declare the variables
            load_state = imag(opts.Z_load)~=0;
            mpol('x', 3 + load_state);
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

            Delta = opts.f0*opts.Ts;

            %BUG HERE BUG HERE BUG HERE
            %X_trig = 1-x(1)^2 + x(2)^2;
            X_trig = 1-x(1)^2 - x(2)^2;

            %clock and rescaled load
            if length(x)<4
                X_load = [];
            else
                X_load = 1-x(4)^2;                
                        end

            X_clock_mode = x(3)*(1-2*Delta - x(3));    
            X_clock_jump = (x(3)-Delta)*(1-2*Delta - x(3));    
            X = [X_trig==0; X_clock_mode>=0; X_load>=0];             
            X_jump = [X_trig==0; X_clock_jump>=0; X_load>=0];
            

            %BUG CHECKING
            % X_clock_mode = x(3);
            % X_clock_jump = (x(3))*(1 - x(3));    
            % X = [X_trig==0; X_clock_mode==0; X_load>=0];
            % 
            % X_jump = [X_trig==0; X_clock_jump>=0; X_load>=0];
            % 

            

                   
            lsupp_base.X = X;
            %define the reset law for the jump
            

            %declare the jumps and modes

            %get the per-mode objectives
            %TODO: generalize this towards three-phase objectives
            objective_mode = obj.objective_level(vars, opts);
            
            for m=0:k
                lsupp_curr = lsupp_base;
                arc_curr = support_arc(m, x, Delta, opts.Symmetry);
                lsupp_curr.X = [lsupp_curr.X; arc_curr>=0];
                modes{m+1} = opp_mode(m, lsupp_curr, objective_mode, opts);

                if m~=0
                    jumps{m} = opp_jump(m, opts, vars, X_jump);
                end
            end
        end

        %% the main routine
        function [sol, obj] = run(obj, order)
            %RUN the main call, the full peak program at the target order
            

            d = 2*order;
            
            %formulate constraints
            [mom_con, supp_con, len_dual] = obj.cons(d);
                        
            [objective, obj_con] = opp_objective(obj);
            %solve the program
            sol = obj.solve(objective, mom_con,supp_con);
            
            %process dual variables            
            % obj = obj.dual_process(d, sol.dual_rec, len_dual);
        end
        
        %% pose the problem
        function sol = solve(obj, objective, mom_con,supp_con)
            %SOLVE formulate and solve peak estimation program from
            %constraints and objective    
            
            %TODO: incorporate minquery into maximin (minimax) formulation

            mset('yalmip',true, 'verbose', false);
            %make sure the solution is precise
            mset(obj.sdp_settings);
            % https://docs.mosek.com/9.2/pythonfusion/parameters.html
            
            
            P = msdp(min(objective), mom_con, supp_con);

            sol = struct;
            tic;
            [sol.status,sol.obj_rec, ~,sol.dual_rec]= msol(P);     
            sol.solver_time = toc;

            sol.mass = obj.mass_summary();

            [~, mom_harm] = obj.con_harmonics();
            sol.harmonics = double(mom_harm);
            
        end  

        %% form the constraints
        function [mom_con, supp_con, len_dual] = cons(obj, d)

            %generate the constraints
            supp_con = obj.supp_con();
            
              

            % %mass of initial measure = 1
            con_prob = obj.con_prob_dist();
            % 
            % %initial = sum of terminal measure
            con_preserve = obj.con_return(d);

            %flow +jump continuity constraints
            % 
            con_liou = obj.con_flow(d);
            % 
            % %harmonics constraints
            [con_harm, ~] = obj.con_harmonics();
            % 
            con_leb = obj.con_uni_circ(d);
            % 
            con_threephase = obj.con_balance(d);

            % mom_con = [con_prob; con_preserve; con_leb];
            % mom_con = [con_prob; con_preserve; con_leb; con_harm];
            % mom_con = [con_prob; con_preserve; con_leb(1:end)];
            % mom_con = [con_prob; con_leb];

            % mom_con = [con_prob; con_preserve; con_harm; con_leb; con_threephase];
            
            %without harmonics
            % mom_con = [con_prob; con_liou; con_leb; con_preserve; con_threephase];


            %with harmonics
            mom_con = [con_prob; con_preserve; con_liou; con_harm; con_leb; con_threephase];


            %TODO: objective constraints as well
            
            %TODO: index the dual variables
            len_dual = 0;

        end

        function uni_con = con_uni_circ(obj, d)
            %the signal encircles the unit disk entirely in one period
            %this means that the (c, s) marginal of the occupation measures
            %is a uniform distribution (because time has been scaled back)
            
            %get the lebesgue distribution
            if obj.opts.uniform_arc               
                pw = genPowGlopti(2, d);
                leb_circ = leb_sphere(pw,1);
    
                %time is scaled, should be uniform moments
                uni_circ = leb_circ/(2*pi);
    
                %get moments of the (c, s)-marginal
                trmon_sum = 0;
                %TODO: change the lebesgue constraint for symmetry
                for m = 1:length(obj.modes)
                    tr_curr = obj.modes{m}.trig_occ_monom(d);
                    trmon_sum = trmon_sum + cell_sum(tr_curr);
                end
    
                uni_con = (trmon_sum == uni_circ);
            else
                uni_con = [];
            end
            
        end

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

            %harmonics on the load side
            if ~isempty(harm_load_data)          
                harm_load = 0;
                for m = 0:length(obj.modes)     
                    harm_mom_load = obj.modes{m}.load_harmonics_mom(obj.vars, harm_load, obj.opts.harmonics_load);
                    harm_load = harm_load + harm_mom_load;
                end
                harm_load_con = harmonics_process(obj.opts.harmonics_load, harm_load);
            else
                harm_load_con = [];
            end
            
            harm_con = [harm_load_con; harm_source_con];

            %TODO: finish this. Sum up the harmonics over all components
        end

        function bcon = con_balance(obj, d)
            %ensure that the inductive current is three-phase balanced
            %
            %figure out if this is possible to do without inductive current
            %
            if length(obj.vars.x)>3 && imag(obj.opts.Z_load)>0 ...
                    && obj.opts.three_phase == opp_three_phase.Balanced
                
                vars_inv = obj.vars.x([1, 2, 4]);
                p_in = mmon(vars_inv, d-1)*vars_inv(3);

                mon_3 = obj.three_phase_rotate(p_in, vars_inv);

                width_3 = size(mon_3, 2);
                mon_3_sum = mon_3*ones(width_3, 1);

                %process the symmetry (if valid)
                mon_3_sum = obj.symmetry_eval(mon_3_sum, obj.vars.x([1; 2]));

                bmom = 0;
                for m = 1:length(obj.modes)
                    curr_mom = obj.modes{m}.mom_sub(obj.get_vars(), mon_3_sum);
                    bmom = madd_cell_mom(bmom, curr_mom, 1);
                end

                
                bcon = [];
                for i =1:size(bmom, 1)
                    for j = 1:size(bmom, 2)
                        bcon = [bcon; bmom{i, j}==0];
                    end
                end
                
            else
                bcon = [];
            end
        end


        function var_stack = get_vars(obj)
            %return variables
            var_stack = [obj.vars.t; obj.vars.x];
        end

        function w_sym = symmetry_eval(obj, w_in, vars_trig)
            %compensate for the symmetry structure in the problem
            %
            %original: w(c, s) over the occupation measure
            %
            %Full-Wave w(c, s)
            %Half-Wave w(c, s) - w(-c, -s)
            %Quarter-Wave w(c, s) + w(-c, s) - w(-c, -s) - w(c, -s)

            if obj.opts.Symmetry==0
                w_sym = w_in;
            else
                w_refl = subs(w_in, vars_trig, -vars_trig);
                if obj.opts.Symmetry==1
                    w_sym = w_in - w_refl;
                else
                    Rp = [-1, 0; 0, 1];
                    w_q_pos = subs(w_in, vars_trig, Rp*vars_trig);
                    w_q_neg = subs(w_in, vars_trig, -Rp*vars_trig);
                    w_sym = w_in + w_q_pos - w_q_neg - w_refl;
                end
            end
        end


        function mon_3 = three_phase_rotate(obj, p_in, vars_inv)            
            %return a vector of polynomials in vars_inv
            %rotated as [u(theta), u(theta-2pi/3), u(theta-4pi/3)]
            %variable 1 and 2 are trigonometrically related (cos and sin)
            %the others are along for the ride

            %TODO: use this in constructing three-phase symmetry

            R3 = [cos(2*pi/3), -sin(2*pi/3); sin(2*pi/3), cos(2*pi/3)];

            vars_inv_trig = vars_inv(1:2);
            %monomials times the current
            va = p_in;
            vb = subs(va, vars_inv_trig, R3*vars_inv_trig);
            vc = subs(va, vars_inv_trig, (R3*R3)*vars_inv_trig);

            mon_3 = [va, vb, vc];

            % if obj.opts.Symmetry==1
            %     R2 = [-1, 0; 0, 1];
            % 
            % 
            %     mon_3_flip = subs(mon_3, vars_inv_trig, R2*vars_inv_trig);
            % 
            %     mon_3 = [mon_3, -mon_3_flip];
            % end

        end

    
        function harm_monom = harm_eval(obj, vars, harm_in)            
            % [vars.x(1).^harm_in.index_cos; vars.x(2).^harm_in.index_sin]/pi;

            %TODO: this repeated computation (in opp_locations) is 
            %inefficient, fix later. Not the most pressing issue.

            %compute the chebyshev moments. then divide by pi
            %cos(n theta) = T_n(cos(theta))
            %sin(n theta) = sin(theta) U_{n-1}(cos(theta))

            %TODO: modify this for symmetry
            
            c = vars.x(1);
            s = vars.x(2);
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
    
                    harm_cos = T(harm_in.index_cos+1);
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
    
                    harm_sin = s*U(harm_in.index_sin);
                else
                    harm_sin = [];
                end
    
                %package up the harmonics
                % harm_monom = [harm_cos; harm_sin]/pi;    

                %divide by pi to get the harmonic
                %then multiply by 2pi because time is scaled to [0, 1]?
                %figure this out
                harm_monom = [harm_cos; harm_sin]*2;

                %process the symmetry
                harm_monom = obj.symmetry_eval(harm_monom, [c; s]);
            end
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
            for m = 1:Nmodes
                liou_curr = liou_cell{m};
                for n=1:N
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
                        
                        flow_curr = liou_curr{n, p} + src_curr + dst_curr==0;
                        
                        %stack them into a giant vector: flow_con
                        flow_con = [flow_con; flow_curr];
                    end
                end
            end

            
        end

        function mass_con_eq = con_prob_dist(obj)
            %initial measure is a probability distribution (mass 1)
            
            [~, mass_init_sum] = obj.modes{1}.initial_mass();

            % mass_init = 

            % mass_init_summary = sum(sum(mass_init_all));

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
                if obj.opts.Symmetry == 0
                    stop_order = 1:N;
                else
                    stop_order = N:-1:1;                
                end

                if obj.opts.early_stop                
                    for m = 3:2:length(obj.modes)
                        % mass_con = mass_con - obj.modes{m}.mass_term_mode();
                        stop_monom = obj.modes{m}.term_monom(d, true);
                        return_mom = madd_cell_mom(return_mom, {stop_monom{stop_order, end}}, -1);
                    end
                else
                    % mass_con = mass_con - obj.modes{end}.mass_term_mode();
                    stop_monom = obj.modes{end}.term_monom(d, true);
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
          %cell of opp_mode(), contains initial/terminal/occupation measures           
        end

        %% process the objective

        function objective = objective_level(obj, vars, opts)
            %return the mode-objective at each level
            %
            %TODO: 
            %the scaling factors may be wrong for inductance and
            %capacitance. check this later.
            %
            N = length(opts.L);
            Lmax = max(abs(opts.L));
            sym_factor = double(2^opts.Symmetry);
            %Another TODO: quarter-wave symmetry may break the
             %characterization of the current in the inductor/capacitor
            if (length(vars.x)==3) || (imag(opts.Z_load) == 0)                      
                %purely resistive
                %sym_factor: replicate the square according to the symmetry
                %2*pi: because time is normalized to [0, 1]
                objective = sym_factor * (2*pi)*(opts.L.^2);
            elseif (imag(opts.Z_load) >= 0)
                
                %inductive load
                %i' = -(R/L)i + (1/L) v
                %per-unit system, ignore the L value
                inductance = imag(opts.Z_load)/(2*pi*opts.f0);
                % resistance= real(opts.Z_load);                
                % f_load = -(resistance)/(inductance)*vars.x(4) + Lscale;
                objective = vars.x(4)*(Lmax/inductance)^2*(opts.L).^2;
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
            objective = objective'*sym_factor;


            if opts.three_phase == "Floating"
                objective = obj.three_phase_rotate(objective, opts.vars.x);
            end
        end

        function [objective, obj_con] = opp_objective(obj)
            %OPP_OBJECTIVE Form the objective of the OPP problem
            %
            %Output: 
            % objective: type mom, objective to minimize
            % obj_con:   constraints used to define the objective
            
            %will need to modify this for the three-phase quadratic program
            obj_con = [];
            objective = 0;
            if obj.opts.null_objective
                %for testing only
                % m0 = obj.modes{1}.initial_mass();
                % [~, mass_init_sum] = ones(1, length(obj.opts.L))*m0(:, 1);
                [~, mass_init_sum] = obj.modes{1}.initial_mass();
                objective = mass_init_sum;
            else
                for i = 1:length(obj.jumps)
                    objective = objective + obj.jumps{i}.objective();
                end
                
                for i = 1:length(obj.modes)
                    objective = objective + obj.modes{i}.objective();
                end
            end

        end

        function opp_out = recover(obj)
            %process and recover the solution

            opp_out = [];
        end

        function [m_out] = mmat(obj)
            %get the moment matrix of all measure variables
            m_out = struct;
            K = length(obj.modes);
            % [N, P] = size(obj.modes{1}.levels);
            % m_out.levels = cell(K, N, P);
            m_out.modes = cell(K, 1);

            for i = 1:K
                m_out.modes{i} = obj.modes{i}.mmat();
            end

            m_out.jump = cell(K-1, 1);
            for i=1:(K-1)
                m_out.jump{i} = obj.jumps{i}.mmat();
            end
        end

        function ms = mass_summary(obj)
            %collect the masses of the occupation measure into a neat array
            ms = struct;
            K = length(obj.modes);
            [N, P] = size(obj.modes{1}.levels);

            ms.mode = zeros(K, N, P);
            ms.trans = zeros(K, N, P-1);
            for m=1:K
                for n=1:N
                    for p = 1:P
                        ms.mode(m, n, p) = double(obj.modes{m}.levels{n, p}.sys{1}.meas_occ.mass());
                        if p < P
                            ms.trans(m, n, p) = double(obj.modes{m}.transition{n, p}.mass());
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

        function [pattern] = recover_pattern(obj)
            %RECOVER_PATTERN recover a switching sequence from a solution 
            %(mass of occupation measures)
            %
            %Output: struct pattern
            % occ:      percentage of time spent in each location
            % alpha:    switching angles
            % u:        voltage levels
            % energy:   energy in the switching sequence
            ms = obj.mass_summary;
            % L = obj.L;
            
            mocc = sum(ms.mode, 3);
            

            pattern = struct;
            pattern.occ = mocc;

            %TODO: generalize for other symmetry structures

            %TODO: recovery for early stopping
            [~, ind] = max(mocc, [], 2);
            ang = sum(mocc, 2);
            pattern.alpha = 2*pi*cumsum(ang(1:end-1))';
            pattern.u = obj.opts.L(ind);
            pattern.energy = (obj.opts.L(ind).^2)*(2*pi*ang);

            

        end

    end
end

