classdef opp_location < location_interface
    %OPP_LOCATION Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        mode;
        partition; 
        level;
        L;
        
        I_split = [];        
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
            
            %split the current into positive and negative 
            if info.I_split 
                obj.I_split = opp_current_split(loc_supp, id);
            end
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
               
            [objective, cons_obj] = obj.objective_con();


            con_power = obj.power_match_con();

            cons_ineq = [cons_obj];
            
            %package up the output
            len_dual = struct;
            len_dual.v = len_liou;
            len_dual.beta = length(cons_ineq);
            
            %ensure this is the correct sign
            cons_eq = [-liou==0; con_power];                        
        end      

        function power_use = conduction_losses(obj, dispatch)
            %power used in a fundamental cycle dissipated by conduction
            m = obj.mode;
            % cond = [0, 0];
            power_use = 0;
            ind_pos = dispatch.conduction_pos{m+1};
            ind_neg = dispatch.conduction_neg{m+1};
            % I_lin = obj.I_split.mom_lin();

            %scale by the maximum rated current
            Ipos = dispatch.IT * obj.I_split.pos.vars.x;
            Ineg = dispatch.IT * obj.I_split.neg.vars.x;
            for i = 1:length(ind_pos)
                cond_curr = dispatch.topology{ind_pos(i)}.conduction;
                power_use = power_use + mom(Ipos * (cond_curr(1) + cond_curr(2)*Ipos));
            end
            for i = 1:length(ind_neg)
                cond_curr = dispatch.topology{ind_neg(i)}.conduction;
                power_use = power_use + mom(Ineg * (cond_curr(1) + cond_curr(2)*Ineg));
            end            

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

        function con_power = power_match_con(obj, d)
    
            %splitting the current into positive and negative components
            %for bounding the power losses
            if isempty(obj.I_split)
                con_power = [];
            else
                con_power = [obj.sys{1}.mom_monom(d) - obj.I_split.mom_monom(d) == 0];
            end
        
        end

        function om = occ_mass(obj)
            %mass of the occupation measure in this level
            om = obj.sys{1}.meas_occ.mass();
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

        function [v_sel, mon_sel] = select_monom_init(obj, d, ind)
            %moments of all other variables [phi, l] 
            
            if isempty(obj.supp.X)
                v_sel = 0;
                mon_sel= 0;
            else
                x_curr = obj.init.meas{1}.vars.x;
                x_sel = x_curr(ind);                
                v_sel= mmon(x_sel, 0, d);
                
    
                mon_sel = mom(v_sel);
            end
        end

        function [v_sel, mon_sel] = select_monom_term(obj, d, ind)
            %moments of all other variables [phi, l] 
            
            if isempty(obj.supp.X)
                v_sel = 0;
                mon_sel= 0;
            else
                x_curr = obj.term.meas{1}.vars.x;
                x_sel = x_curr(ind);                
                v_sel= mmon(x_sel, 0, d);
                
    
                mon_sel = mom(v_sel);
            end
        end

        function [v_sel, mon_sel] = select_monom_occ(obj, d, ind)
            %moments of all other variables [phi, l] 
            
            if isempty(obj.supp.X)
                v_sel = 0;
                mon_sel= 0;
            else
                x_curr = obj.sys{1}.meas_occ.vars.x;
                x_sel = x_curr(ind);                
                v_sel= mmon(x_sel, 0, d);
                
    
                mon_sel = mom(v_sel);
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

        function [v_ntrig, mon_ntrig] = non_trig_monom_term(obj, d, flip_load)
            %moments of all other variables [phi, l] 
            if nargin < 3
                flip_load = 0;
            end
            x_curr = obj.term.meas{1}.vars.x;
            x_ntrig = x_curr(3:end);
            v_ntrig = mmon(x_ntrig, 0, d);

            if flip_load && (length(x_ntrig) > 1)
                %flip the load currents in case of half-wave symmetry
                if flip_load == 1
                    %single phase (keep the clock)
                    v_ntrig = subs(v_ntrig, x_ntrig, [1; -1*ones(length(x_ntrig)-1, 1)] .* x_ntrig);
                else
                    %three-phase (no clock)
                    v_ntrig = subs(v_ntrig, x_ntrig, -x_ntrig);
                end
            end
            
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

