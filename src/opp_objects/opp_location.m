classdef opp_location < location_interface
    %OPP_LOCATION Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        mode;
        partition; 
        level;
    end
    
    methods
        function obj = opp_location(loc_supp, f, objective, info)
            %OPP_LOCATION Construct an instance of this class
            %   Detailed explanation goes here
          
            id = info.id;

            obj@location_interface(loc_supp, f, objective, id);
            
            obj.mode = info.mode;
            obj.partition = info.partition;
            obj.level= info.level;
            % obj.id = info.id;
            %TODO: make id the last argument                                       
            obj.sys  = subsystem_base(obj.supp, obj.f, [], id);
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

        
        function [v_trig, mon_trig] = trig_monom(obj, d)
            %moments of [c, s] (trigonometric lift, used for Lebesgue
            %constraint)
            x_curr = obj.sys{1}.meas_occ.vars.x;
            x_trig = x_curr(1:2);
            v_curr = mmon(x_curr_trig, 0, d);

            mon_trig = mom(v_curr);

        end

        function [harm_poly, harm_mom] = voltage_harmonics_mom(obj, vars, harm_in)
            %voltage harmonics evaluation 
            %equivalent to a resistive load
            harm_eval = [vars.x(1).^harm_in.index_cos; vars.x(1).^harm_in.index_sin]/pi;
            sub_eval = obj.sys{1}.meas_occ.var_sub(vars, harm_eval);

            harm_poly = obj.level*sub_eval;
            harm_mom = mom(harm_poly);            
        end
        function [harm_poly, harm_mom] = load_harmonics_mom(obj, vars, harm_in, Z_type, Z_scale)
            %voltage harmonics evaluation 
            %equivalent to a resistive load

            if Z_type == 0
                [harm_poly, harm_mom]  = obj.voltage_harmonics_mom(vars, harm_in);
            elseif Z_type ==1
                [harm_poly, harm_mom]  = obj.capacitance_harmonics_mom(vars, harm_in, Z_scale);
            else
                [harm_poly, harm_mom]  = obj.inductance_harmonics_mom(vars, harm_in, Z_scale);
            end           
        end

        function [harm_poly, harm_mom] = capacitance_harmonics_mom(obj, vars, Z_scale)
            %current evaluation for a capacative load
            harm_eval = [vars.x(1).^harm_in.index_cos; vars.x(1).^harm_in.index_sin]/pi;
            harm_poly = harm_eval.*(obj.level-vars.x(4));

            harm_mom = mom(harm_poly);
        end

        function [harm_poly, harm_mom] = inductance_harmonics_mom(obj, vars, Z_scale)
            %current evaluation for a capacative load
            harm_eval = [vars.x(1).^harm_in.index_cos; vars.x(1).^harm_in.index_sin]/pi;
            harm_poly = harm_eval.*(Z_scale*vars.x(4));

            harm_mom = mom(harm_poly);
        end

        %holdovers from abstract class
        function dual_out = dual_process(obj)
            dual_out = []
        end

        function leq= len_eq_cons(obj)
            leq= []
        end
    end
end

