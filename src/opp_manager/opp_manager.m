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
            lsupp_base.Tmax = 1;

            %create the support set

            Delta = opts.f0*opts.Ts;
            X_trig = 1-x(1)^2 + x(2)^2;

            %clock and rescaled load
            X_clock_mode = x(3)*(1-2*Delta - x(3));    
            X_clock_jump = (x(3)-Delta)*(1-2*Delta - x(3));    
            
            if length(x)<4
                X_load = [];
            else
                X_load = 1-x(4)^2;                
            end

            

            X = [X_trig==0; X_clock_mode>=0; X_load>=0];
            X_jump = [X_trig==0; X_clock_jump>=0; X_load>=0];
                       
            lsupp_base.X = X;
            %define the reset law for the jump
            

            %declare the jumps and modes
            for m=0:k
                lsupp_curr = lsupp_base;
                arc_curr = support_arc(m, x, Delta);
                lsupp_curr.X = [lsupp_curr.X; arc_curr>=0];
                modes{m+1} = opp_mode(m, lsupp_curr, opts);

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

            mset('yalmip',true);
            %make sure the solution is precise
            mset(obj.sdp_settings);
            % https://docs.mosek.com/9.2/pythonfusion/parameters.html
            
            
            P = msdp(min(objective), mom_con, supp_con);

            sol = struct;
            tic; 
            [sol.status,sol.obj_rec, ~,sol.dual_rec]= msol(P);     
            sol.solver_time = toc;
        end  

        %% form the constraints
        function [mom_con, supp_con, len_dual] = cons(obj, d)

            %generate the constraints
            supp_con = obj.supp_con();
            
            %mass of initial measure = 1
            con_prob = obj.con_prob_dist();
            
            %initial = sum of terminal measure
            con_preserve = obj.con_return();

            %flow +jump continuity constraints
            
            con_liou = obj.con_flow(d);

            %harmonics constraints
            con_harm = obj.con_harmonics();

            con_leb = obj.con_lebesgue_circ(d);

            mom_con = [con_prob; con_preserve; con_liou; con_harm; con_leb];


            %TODO: objective constraints as well
            
            %TODO: index the dual variables
            len_dual = 0;

        end

        function leb_con = con_lebesgue_circ(obj, d)
            %the signal encircles the unit disk entirely in one period
            %this means that the (c, s) marginal of the occupation measures
            %is a Lebesgue distribution
            
            %get the lebesgue distribution
            pw = genPowGlopti(2, d);
            leb_circ = LebesgueMomSphere(pw,1);

            %get moments of the (c, s)-marginal
            trmon_sum = 0;
            for m = 1:length(obj.modes)
                tr_curr = obj.modes{m}.trig_occ_monom(d);
                trmon_sum = trmon_sum + cell_sum(tr_curr);
            end

            leb_con = (trmon_sum == leb_circ);
            
        end

        function harm_con = con_harmonics(obj)
            %collect harmonics constraints on the voltage source and the load
            harm = obj.opts.harmonics;
            harm_load_data = obj.opts.harmonics_load;

            % harmonics on voltage source
            if ~isempty(harm)
                harm_source = 0;
                for m = 1:length(obj.modes)   
                    harm_mom = obj.modes{m}.voltage_harmonics_mom(obj.vars, harm);
                    harm_source = harm_source + harm_mom;
                end
                harm_source_con = harmonics_process(harm, harm_source);
            else
                harm_source_con = [];                
            end


            %TODO: grid side filters

            %harmonics on the load side
            if ~isempty(harm_load_data)          
                harm_load = 0;
                for m = 0:length(obj.modes)     
                    harm_mom_load = obj.modes{m}.load_harmonics_mom(obj.vars, harm_load_data);
                    harm_load = harm_load + harm_mom_load;
                end
                harm_load_con = harmonics_process(harm, harm_load);
            else
                harm_load_con = [];
            end
            
            harm_con = [harm_load_con; harm_source_con];

            %TODO: finish this. Sum up the harmonics over all components
        end

        function flow_con = con_flow(obj)
            %the flow conservation constraint (the big one)

            %start the storage structure 
            Nmodes = length(obj.modes);
            liou_cell = cell(Nmodes, 1);
            jump_src = cell(Nmodes-1, 1);
            jump_dst = cell(Nmodes, 1);

            %compute all terms
            for m=1:(Nmodes)
                liou_cell{m} = obj.modes.flow(d);
            end

            %add the jump to the cell terms
            for m=1:(Nmodes)                
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

            %TODO: 
            %This is for full-wave. generalize for other symmetry
            %structures

            return_mom = {init_monom{:, 1}};
            if obj.opts.early_stop                
                for m = 3:2:length(obj.modes)
                    % mass_con = mass_con - obj.modes{m}.mass_term_mode();
                    stop_monom = obj.modes{m}.term_monom(d, true);
                    return_mom = madd_cell_mom(return_mom, {stop_monom{:, end}}, -1);
                end
            else
                % mass_con = mass_con - obj.modes{end}.mass_term_mode();
                stop_monom = obj.modes{end}.term_monom(d, true);
                return_mom = madd_cell_mom(return_mom, {stop_monom{:, end}}, -1);
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

        function [objective, obj_con] = opp_objective(obj)
            %OPP_OBJECTIVE Form the objective of the OPP problem
            %
            %Output: 
            % objective: type mom, objective to minimize
            % obj_con:   constraints used to define the objective
            
            %will need to modify this for the three-phase quadratic program
            obj_con = [];
            objective = 0;
            for i = 1:length(obj.jumps)
                objective = objective + obj.jumps{i}.objective();
            end
            
            for i = 1:length(obj.modes)
                objective = objective + obj.modes{i}.objective();
            end

        end

        function opp_out = recover(obj)
            %process and recover the solution

            opp_out = [];
        end

    end
end

