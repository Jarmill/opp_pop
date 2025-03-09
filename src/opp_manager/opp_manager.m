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

            [jumps, modes] = obj.create_system(opts);
        end

        %% construct everything
        function [jumps, modes] = create_system(opts)
            %used in the constructor

            k = opts.k;
            jumps = cell(k, 1);
            modes = cell(k+1, 1);

            %declare the variables
            load_state = imag(opts.Z_load)~=0;
            mpol('x', 3 + load_state);
            %x = [c; s; phi; l] -> [cos(theta), sin(theta), clock, load
            %state (current of load inductor/voltage of load capacitor)

            
            %create the basic location structure
            if opts.TIME_INDEP
                mpol('t', 1, 1)
            else
                t = [];
            end

            %create the basic support set
            lsupp_base = loc_support();
            lsupp_base.vars.x = x;
            lsupp_base.vars.t = t;

            lsupp_base.TIME_INDEP = opts.TIME_INDEP;
            lsupp_base.FREE_TERM = 0;
            lsupp_base.Tmax = 1;

            %create the support set

            Delta = opts.f0*opts.Ts;
            X_trig = x(1)^2 + x(2)^2 ==1;

            %clock and rescaled load
            X_clock_mode = x(3)*(1-2*Delta - x(3));    
            X_clock_jump = (x(3)-Delta)*(1-2*Delta - x(3));    
            
            if length(vars.x)<4
                X_load = [];
            else
                X_load = 1-x(4)^2;                
            end

            lsupp_base.X = X;

            X = [X_trig==0; X_clock_mode>=0; X_load>=0];
            X_jump = [X_trig==0; X_clock_jump>=0; X_load>=0];
                       
            %define the reset law for the jump
            

            %declare the jumps and modes
            for m=0:k
                lsupp_curr = lsupp_base;
                arc_curr = support_arc(m, x, Delta);
                lsupp_curr.X = [lsupp_curr.X; arc_curr>=0];
                modes{m+1} = opp_mode(m, lsupp_curr, opts);
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
            
            %mass of initial measure = 1

            %initial = sum of terminal measure

            %flow +jump continuity constraints


        end

        
        function [objective, obj_con] = opp_objective(obj)
            %OPP_OBJECTIVE Form the objective of the OPP problem
            %
            %Output: 
            % objective: type mom, objective to minimize
            % obj_con:   constraints used to define the objective
            
            obj_con = [];
            objective = 0;
            for i = 1:length(obj.jumps)
                objective = objective + obj.jumps{i}.objective();
            end
            
            for i = 1:length(obj.modes)
                objective = objective + obj.jumps{i}.modes();
            end

        end

        function opp_out = recover(obj)
            %process and recover the solution

            opp_out = [];
        end

    end
end

