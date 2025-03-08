classdef manager_hy_interface < manager_interface
    %PEAK_MANAGER_HY Hybrid peak estimation manager
    %   Detailed explanation goes here
    
    %TODO: rework this with subsystems and samplers
    
    properties        
        guards;                        
    end
    
    methods
        function obj = manager_hy_interface(locations_in,guards_in)
            %PEAK_MANAGER_HY Construct an instance of this class
            %   Detailed explanation goes here
%             obj.Property1 = inputArg1 + inputArg2;

            if ~iscell(locations_in)
                locations_in = {locations_in};
            end
            
            obj@manager_interface(locations_in);
            
            if ~iscell(guards_in)
                obj.guards = {guards_in};                        
            else
                obj.guards = guards_in;                        
            end           
            
            obj.solver = 'mosek';
        end
        
        %% Formulating and solving program
        
        function [objective, mom_con, supp_con, len_dual] = cons(obj,d, Tmax)
            %CONS formulate support and measure constraints for peak
            %program at degree d
            %Input:
            %   d:      Monomials involved in relaxation (2*order)
            %   Tmax:   Maximum time (only when time-independent)
            %
            %Output:
            %   objective:  target to maximize  (@mom)
            %   mom_con:    moment constraints  (@momcon)
            %   supp_con:   support constraints (@supcon)
            %   len_liou:   number of liouville constraints (uint32)

            supp_con = [];       %support constraint     
            mass_init_sum = 0;   %mass of initial measure should be 1
            objective = 0;                                   
            
            liou_con = cell(length(obj.loc), 1);
            obj_con = cell(length(obj.loc), 1);
            
            mass_occ_sum = 0;

            len_dual = struct('v', zeros(length(obj.loc), 1), 'beta', zeros(length(obj.loc), 1), 'alpha', ones(length(obj.guards), 1), 'gamma', 1, 'gamma_occ', 0);
            
            %process the location measures
            for i = 1:length(obj.loc)
                loc_curr = obj.loc{i};
                
                %accumulate support constraints
                supp_con = [supp_con; loc_curr.supp_con()];

                %find moment constraints of current location
                [obj_curr_max, obj_con_curr_ineq, obj_con_curr_eq] = loc_curr.objective_con();                                
                
                %deal with possible chance (mean/quantile) SDE objective
                %modifications over here
                [obj_curr, cons_eq_obj] = obj.objective_process(obj_curr_max);
                
                if ~isempty(obj_curr)
                    objective = objective + obj_curr;
                end
                
                %hack the SDE content over here
                
                obj_con{i} = [obj_con_curr_ineq; obj_con_curr_eq; cons_eq_obj];
                
                liou_con{i} = loc_curr.liou_con(d);  
                
                %track the lengths of dual variables
                len_dual.v(i) = length(liou_con{i});
                len_dual.beta(i) = length(obj_con_curr_ineq);

                %initial measure has mass 1
                if ~isempty(loc_curr.init)
                    mass_init_sum = mass_init_sum + loc_curr.mass_init();
                end
                
                if isempty(loc_curr.vars.t)
                    mass_occ_sum = mass_occ_sum + loc_curr.mass_occ();
                end
            end
            
            %process the guards
            zeno_con = []; %zeno mass constraints
            
            for i = 1:length(obj.guards)
                g_curr = obj.guards{i};
                
                %support constraints
                supp_con = [supp_con; g_curr.supp];
                
                %zeno
                zeno_con = [zeno_con; g_curr.zeno_con()];
                
                
                %liouville constraints
                [mom_src, mom_dest] = g_curr.liou_reset(d);
                
                liou_con{g_curr.src.id}  = liou_con{g_curr.src.id}  + mom_src;
                liou_con{g_curr.dest.id} = liou_con{g_curr.dest.id} + mom_dest;
            end
            
            %finalize moment constraints
            
            %mass of initial measure sums to one
            if isnumeric(mass_init_sum)
                mass_init_con = [];
            else
                % mass_init_sum == 0 eliminates the constraint when there
                % is one location. keep it in.
                mass_init_con = (mass_init_sum - 1 == 0);
            end
            
            %TIME-INDEPENDENT mass of occupation measures are less than
            %Tmax. Only when t is not included as a variable;
            if isnumeric(mass_occ_sum) %&& (nargin > 2)
                
                %TODO: access Tmax
                mass_occ_con = [];
            else
                % mass_init_sum == 0 eliminates the constraint when there
                % is one location. keep it in.
                %TODO: confirm this behavior, check how this impacts
                %indexing and dual recovery
                mass_occ_con = (mass_occ_sum  <= Tmax);
                len_dual.gamma_occ = 1;
            end
            
            
            loc_con = [];
            liou_con_all = [];
            obj_con_all = [];           
            for i = 1:length(liou_con)
                %TODO: figure out the sign convention for dual variables
                %should the negative sign be there on liou_con?
%                 loc_con = [loc_con; -liou_con{i} == 0; obj_con{i}];
                liou_con_all = [liou_con_all; -liou_con{i} == 0];
                obj_con_all  = [obj_con_all; obj_con{i}];
            end
                
            len_liou = length(liou_con_all);
            
%             len_dual = [];

            
                        
            mom_con = [mass_init_con; liou_con_all; mass_occ_con; obj_con_all; zeno_con];

        end                    
    
        
        function [objective, cons_eq_obj] = objective_process(obj, obj_max)
            %OBJECTIVE_PROCESS implement the possible chance-constraints in
            %case of stochastic systems.
            
            %by default, consider only the mean-maximization task.
            
            %other code could involve chance-bounds (concentration bounds
            %such as cantelli and vp)
            objective = obj_max(1);
            cons_eq_obj = [];
        end
 
        function s_out = mmat_corner(obj)
            %get the top corner of the moment matrix for all measures
%             s_out = struct('locations',  cell(length(obj.loc), 1), ...
%                            'guards', cell(length(obj.loc), 1));
                       
            s_out = {};
            s_out.locations = cell(length(obj.loc), 1);
            s_out.guards = cell(length(obj.guards), 1);
            
            for i = 1:length(obj.loc)
                s_out.locations{i} = obj.loc{i}.mmat_corner();
            end
            
            for i = 1:length(obj.guards)
                s_out.guards{i} = obj.guards{i}.mmat_corner();
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
            
            rec_eq = dual_rec{1};
            rec_ineq = dual_rec{2};
  
            gamma = rec_eq(1);

            liou_offset = 1;
            cost_con_offset = 0;

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
                
                obj.loc{i}.len_dual = struct('v', liou_len_curr, 'zeta', [], 'beta', len_dual.beta(i));
                obj.loc{i} = obj.loc{i}.dual_process(d, v_coeff, beta_curr, gamma);
                
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
        

        
        %% Recovery
        
        function [optimal, mom_out, corner] = recover(obj, tol)
            %RECOVER if top corner of the moment matrix is rank-1, then
            %return approximate optimizer
            
            if nargin < 2
                tol = 5e-4;
            end
            
            optimal = zeros(length(obj.loc), 1);
            mom_out = cell(length(obj.loc), 1);
            corner = cell(length(obj.loc), 1);
            for i = 1:length(obj.loc)
                [optimal(i), mom_out{i}, corner{i}] = obj.loc{i}.recover(tol);                 
            end            
        end
        
        
    end
    
    
end

