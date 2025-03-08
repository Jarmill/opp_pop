classdef distance_manager_hy < manager_hy_interface
    %PEAK_MANAGER_HY Hybrid peak estimation manager
    %   Detailed explanation goes here
    
    %TODO: rework this with subsystems and samplers
        
    methods
        function obj = distance_manager_hy(locations_in,guards_in)
            %PEAK_MANAGER_HY Construct an instance of this class
            %   Detailed explanation goes here
%             obj.Property1 = inputArg1 + inputArg2;

            obj@manager_hy_interface(locations_in, guards_in);
        end

        function sol = solve(obj, objective, mom_con,supp_con)
            %SOLVE formulate and solve peak estimation program from
            %constraints and objective    
            
            %TODO: incorporate minquery into maximin (minimax) formulation

            mset('yalmip',true);
            %make sure the solution is precise
            mset(obj.sdp_settings);
            % https://docs.mosek.com/9.2/pythonfusion/parameters.html
            
            
            tic
            P = msdp(min(objective), mom_con, supp_con);
            timerec = toc;
            sol = struct;
            [sol.status,sol.obj_rec, ~,sol.dual_rec]= msol(P);
            sol.solver_time = timerec;
        end  
        
        %% Formulating and solving program
        function [objective, mom_con, supp_con, len_dual] = cons(obj,d, Tmax)
            %add in the wasserstein marginal constraints

            if nargin <3
                Tmax = 1;
            end
        

            [objective, mom_con, supp_con, len_dual] = cons@manager_hy_interface(obj, d, Tmax);

            len_dual.w = zeros(length(obj.loc));
            for i = 1:length(obj.loc)
                %include the csp clique overlap constraints
                [marg_cons_curr, len_dual.w(i)] = obj.loc{i}.marg_cons(d);
                mom_con= [mom_con; marg_cons_curr];
               
            end
        end

        function obj = dual_process(obj, d, dual_rec, len_dual)           
            %DUAL_PROCESS dispatch the dual variables from solution to
            %locations and measures, turn the variables into nonnegative
            %functions along trajectories
            
            %v: coefficients of liouville (per loc)
            %beta: coefficients of cost (per loc)
            %alpha: dual of zeno gaps (per guard)
            %gamma: initial mass
            %gamma_occ: occ. mass (if time t not a variable) 
            
            %unfortunate copy-paste

            rec_eq = dual_rec{1};
            rec_ineq = dual_rec{2};
  
            gamma = rec_eq(1);

            liou_offset = 1;
            cost_con_offset = 0;
            marg_con_offset = sum(len_dual.v)+1;

            %locations
            for i = 1:length(obj.loc)                               
                
                %liouville                
                liou_len_curr = len_dual.v(i);
                
                v_coeff = rec_eq((1:liou_len_curr) + liou_offset);
                                
                liou_offset = liou_offset + liou_len_curr;
                
                %maximin cost duals

                if len_dual.beta(i)
                    beta_curr = rec_ineq((1:len_dual.beta(i)) + cost_con_offset);
                    cost_con_offset = cost_con_offset + len_dual.beta(i);
                else
                    beta_curr = 1;
                end
                
                marg_len_curr = len_dual.w(i);
                w_coeff = rec_eq((1:marg_len_curr) + marg_con_offset);
                marg_con_offset = marg_con_offset + marg_len_curr;
                
                obj.loc{i}.len_dual = struct('v', liou_len_curr, 'zeta', [], 'beta', len_dual.beta(i), 'w', marg_len_curr);
                obj.loc{i} = obj.loc{i}.dual_process(d, [v_coeff;w_coeff], beta_curr, gamma);
                
            end            
%           
            %TODO: time-independent formulation
            if len_dual.gamma_occ
                cost_con_offset = 1;
            else
%                 mass_occ_con
            end

            %guards
            for i = 1:length(obj.guards)                
                obj.guards{i} = obj.guards{i}.dual_process(rec_ineq(cost_con_offset + i));
            end                        
        end
        
    end
    
end

