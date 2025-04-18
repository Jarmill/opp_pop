classdef opp_location < location_interface
    %OPP_LOCATION Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        mode;
        partition; 
        level;
        L;
    end
    
    methods
        function obj = opp_location(loc_supp, f, objective, info)
            %OPP_LOCATION Construct an instance of this class
            %   Detailed explanation goes here
          
            id = info.id;

            obj@location_interface(loc_supp, f, [], id);
            
            obj.mode = info.mode;
            obj.partition = info.partition;
            obj.level= info.level;
            obj.L = info.L;
            obj.f = f;
            obj.objective = objective;
            % obj.id = info.id;
            %TODO: make id the last argument                                       
            obj.sys  = {subsystem_base(obj.supp, obj.f, [], id)};
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
            
            %ensure this is the correct sign
            cons_eq = (-liou==0);                        
        end      

        %
        %TODO: need to deal with the objective
        %


        function [obj_min, obj_con_ineq, obj_con_eq] = objective_con(obj, objective)
            %OBJECTIVE_CON deal with the objective, which may be maximin

            %TODO: This should maybe go in the manager
            %The current implementation is only for peak estimation

            %TODO: include support for putting objectives on initial and
            %occupation measures as well as the terminal measure
            if nargin == 1
                objective = obj.get_objective();
            end
            % 
            % obj_con_eq = [];
            % obj_con_ineq = [];
            % 
            var_end = obj.var_index(obj.vars, {'t', 'x'});
            % if isempty(objective)
            %     obj_min = 0;
            % elseif isscalar(objective)

            %this allows for three-phase considerations
            if isnumeric(objective)
                obj_subs = objective*obj.sys{1}.meas_occ.mass();
            else
                if isempty(obj.supp.X)
                    obj_subs = 0;
                else
                    obj_subs = mom(obj.sys{1}.meas_occ.var_sub([var_end.t; var_end.x], objective));
                end
            end

                obj_min = (obj_subs);   
                obj_con_ineq = [];
                obj_con_eq = [];
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
        

        
            function [v_trig, mom_trig] = trig_monom(obj, d, signs)
            %moments of [c, s] (trigonometric lift, used for Lebesgue
            %constraint)
            if nargin < 3
                signs = [1; 1];
            end
            if isempty(obj.supp.X)
                v_trig = 0;
                mom_trig = 0;
            else
                x_curr = obj.sys{1}.meas_occ.vars.x;
                x_trig = x_curr(1:2);
                             
                v_trig = mmon(x_trig, 0, d);
                
                v_trig =  subs(v_trig, x_trig, diag(signs)*x_trig);                
                               

    
                mom_trig = mom(v_trig);
            end

        end

        function [v_ntrig, mon_ntrig] = non_trig_monom_init(obj, d)
            %moments of all other variables [phi, l] 
            if isempty(obj.supp.X)
                v_ntrig = 0;
                mon_ntrig = 0;
            else
                x_curr = obj.init.meas{1}.vars.x;
                x_ntrig = x_curr(3:end);
                v_ntrig = mmon(x_ntrig, 0, d);
    
                mon_ntrig = mom(v_ntrig);
            end
        end

        function [v_ntrig, mon_ntrig] = non_trig_monom_term(obj, d)
            %moments of all other variables [phi, l] 
            x_curr = obj.term.meas{1}.vars.x;
            x_ntrig = x_curr(3:end);
            v_ntrig = mmon(x_ntrig, 0, d);

            mon_ntrig = mom(v_ntrig);

        end
      

        function mom_out = mom_occ_sub(obj, vars, vref)
            v_sub = obj.sys{1}.meas_occ.var_sub(vars, vref);
            mom_out = mom(v_sub);            
        end

        function [harm_poly, harm_mom] = voltage_harmonics_mom(obj, vars, harm_mon, signs)
            if nargin < 4
                signs = [1, 1];
            end
            
            harm_mon = subs(harm_mon, vars.x(1:2), diag(signs)*vars.x(1:2));
            %voltage harmonics evaluation 
            %equivalent to a resistive load
            if obj.L ~= 0 && ~isempty(obj.supp.X)
                % harm_monom = obj.harm_eval(vars, harm_in);
                sub_eval = obj.sys{1}.meas_occ.var_sub([vars.t; vars.x], harm_mon);

                harm_poly = obj.L*sub_eval;
                harm_mom = mom(harm_poly);            
            else
                harm_poly = 0;
                harm_mom = 0;
            end
        end
        function [harm_poly, harm_mom] = load_harmonics_mom(obj, vars, harm_mon, Z_type, Z_scale)
            %voltage harmonics evaluation 
            %equivalent to a resistive load
            if isempty(obj.supp.X)
                harm_poly = 0;
                harm_mom = 0;
            else
                if Z_type == 0
                    [harm_poly, harm_mom]  = obj.voltage_harmonics_mom(vars, harm_mon);
                elseif Z_type ==1
                    [harm_poly, harm_mom]  = obj.capacitance_harmonics_mom(vars, harm_in, Z_scale);
                else
                    [harm_poly, harm_mom]  = obj.inductance_harmonics_mom(vars, harm_mon, Z_scale);
                end           
            end
        end

        function [harm_poly, harm_mom] = capacitance_harmonics_mom(obj, vars, harm_mon, Z_scale)
            %current evaluation for a capacative load
            
            harm_poly = harm_mon.*(obj.L-vars.x(4));

            harm_mom = mom(harm_poly);
        end

        function [harm_poly, harm_mom] = inductance_harmonics_mom(obj, vars, harm_mon, Z_scale)
            %current evaluation for a capacative load
 
            harm_poly = harm_mon.*(Z_scale*vars.x(4));

            harm_mom = mom(harm_poly);
        end        

        %recover the solution
        function m_out = mmat(obj)
            m_out = struct('init', [], 'term', [], 'occ', []);
            if ~isempty(obj.init)
                m_out.init = obj.init.mmat();
            end

            if ~isempty(obj.term)
                m_out.term = obj.term.mmat();
            end
            
            m_out.occ = obj.sys{1}.meas_occ.mmat();
            
        end

        function m_out = mmat_corner(obj)
            m_out = struct('init', [], 'term', [], 'occ', []);
            if ~isempty(obj.init)
                m_out.init = obj.init.mmat_corner();
            end

            if ~isempty(obj.term)
                m_out.term = obj.term.mmat_corner();
            end
            
            m_out.occ = obj.sys{1}.meas_occ.mmat_corner();
            
        end

        function [optimal, mom_out, corner] = recover(obj, tol)
            %RECOVER if top corner of the moment matrix is rank-1, then
            %return approximate optimizer
            
            if nargin < 2
                tol = 5e-4;
            end
                        
            if isempty(obj.init)
                opt_init = 1;
                mom_init.t = []; mom_init.x = [];
                corner_init = 0;
            else
                [opt_init, mom_init, corner_init] = obj.init.recover(tol);
            end
            if isempty(obj.term)
                opt_term = 1;
                mom_term.t = []; mom_term.x = [];
                corner_term = 0;
            else
                [opt_term, mom_term, corner_term] = obj.term.recover(tol);
            end
            
            optimal = opt_init && opt_term;
            
            mom_out = struct('t0', mom_init.t, 'x0', mom_init.x, ...
                             'tp', mom_term.t, 'xp', mom_term.x);     
            corner = struct('init', corner_init, 'term', corner_term);
        end


        %holdovers from abstract class
        function dual_out = dual_process(obj)
            dual_out = [];
        end

        function leq= len_eq_cons(obj)
            leq= [];
        end
    end
end

