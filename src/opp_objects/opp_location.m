classdef opp_location < location_interface
    %OPP_LOCATION Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        mode;
        partition; 
    end
    
    methods
        function obj = opp_location(loc_supp, f, objective, id)
            %OPP_LOCATION Construct an instance of this class
            %   Detailed explanation goes here
            if nargin < 3             
                %by default, no objective
                objective = [];                
            end
            
            if nargin < 4
                id = [];            
            end
            obj@location_interface(loc_supp, f, objective, id);
            
            %TODO: make id the last argument                                       
            obj.sys  = subsystem_base(obj.supp, obj.f, 1, id);
        end
        
        function vars_out = get_vars_end(obj)
            %GET_VARS_END variables at endpoint measures
            %   initial and terminal, without time-dependent
            vars_out = [obj.vars.t; obj.vars.x];
        end
        
        function vars_out = get_vars(obj)
            %GET_VARS add more variables as necessary
            vars_out = [obj.vars.t; obj.vars.x];
        end

        function [objective, cons_eq, cons_ineq, len_dual] = all_cons(obj, d)
            %ALL_CONS all constraints involving solely location
            %does not include sum of mass of initial measures
            %Output:
            %   cons_eq: equality constraints
            %   cons_ineq: inequality constraints (objective)
            
            %gather all constraints together
            liou = obj.liou_con(d);
            len_liou = length(liou);
               
            [objective, cons_ineq] = obj.objective_con();
            
            %package up the output
            len_dual = struct;
            len_dual.v = len_liou;
            len_dual.beta = length(cons_ineq);
            
            %ensure this iss the correct sign
            cons_eq = (-liou==0);                        
        end      

        function [obj_min, obj_con_ineq, obj_con_eq] = objective_con(obj, objective)
            %OBJECTIVE_CON deal with the objective, which may be maximin
                                    
            %TODO: This should maybe go in the manager
            %The current implementation is only for peak estimation
            
            %TODO: include support for putting objectives on initial and
            %occupation measures as well as the terminal measure
            if nargin == 1
                objective = obj.get_objective();
            end
                                    
            obj_con_eq = [];
            obj_con_ineq = [];
            
            var_end = obj.var_index(obj.vars, {'t', 'x'});
            if isempty(objective)
                obj_min = 0;
            elseif length(objective) == 1    
                obj_subs = obj.sys{1}.meas_occ.var_sub_mom(var_end, objective);
                obj_min = (obj_subs);                            
            % else
            %     obj_subs = obj.term.var_sub_mom(var_end, objective);
            %     q_name = ['q_', num2str(obj.id)];
            %     mpol(q_name, 1, 1);
            %     q = eval(q_name);
            %     muq = meas(q);
            %     obj.cost_q = q;
            % 
            %     obj_min = mom(q);
            %     obj_con_eq = [mass(q) == 1];
            %     obj_con_ineq= (mom(q) <= obj_subs);
            end            
        end

    end
end

