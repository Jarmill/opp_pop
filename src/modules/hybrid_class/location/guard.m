classdef guard < meas_base
    %GUARD guard measure governing transitions
    %   Detailed explanation goes here
    
    properties
        id = [];
        src = [];
        dest = [];
        
        reset = [];
        reset_null = [];
        reset_identity = 0;
        
        zeno_cap = 4; %maximum number of transitions along guard
        dual = struct('zeno', 0, 'solved', 0);  %dual variable to zeno constraints        
    end
    
    properties(Access=private)
        TIME_INDEP= 0;
    end
    
    methods
        function obj = guard(id, vars_old, src,dest,supp_old, reset_old)
            %GUARD Construct an instance of this class
            %   Detailed explanation goes here
%             obj.Property1 = inputArg1 + inputArg2;
            
            obj@meas_base(vars_old, supp_old);
            obj.id = id;
            if isnumeric(obj.id)
                obj.id = num2str(obj.id);
            end
            obj.src = src;
            obj.dest = dest;
%             obj.reset = reset_old;
            
            obj = obj.var_def('g', vars_old, supp_old, reset_old);
                        
        end
        
        function obj = set_zeno(zeno_new)
            obj.zeno_cap = zeno_new;
        end
        
        function obj = var_def(obj, suffix, vars_old, supp_old, reset_old)
            %VAR_DEF create new variables 't[suffix]_id',
            %'x[suffix]_id
            
            
            if isempty(obj.vars.t)
                t_new = [];
                obj.TIME_INDEP = 1;
            else
                tname = ['t', suffix, '_', obj.id];                       
                mpol(tname, 1, 1);
                t_new = eval(tname);                   
            end
            
            xname = ['x', suffix, '_', obj.id];                        
            mpol(xname, length(obj.vars.x), 1);
            x_new = eval(xname);
            
            obj.vars = struct('t', t_new, 'x', x_new);
            
            obj.supp = subs_vars(supp_old, [vars_old.t; vars_old.x], ...
                                obj.get_vars());
            obj.reset = subs_vars(reset_old, [vars_old.t; vars_old.x], ...
                obj.get_vars());
            
            if isequal(obj.reset, obj.vars.x)
                obj.reset_identity = 1;
            else
                %check if some of the variables are nulled to zero
                obj.reset_null = false(length(obj.vars.x), 1);
                for i = 1:length(obj.vars.x)
                    if coef(obj.reset(i))==0 && pow(obj.reset(i))==0
                        obj.reset_null(i) = true;
                    end
                end
            end


            
            obj.meas = meas(obj.get_vars());
        end
        
        function con = zeno_con(obj)
            %zeno execution constraint
            %at most (zeno_cap) transitions will occur on the guard
            if (obj.zeno_cap > 0) && (obj.zeno_cap < Inf)               
                con = [obj.mass() <= obj.zeno_cap];
            else
                con = [];
            end
        end
        
        function mom_out = reset_push(obj, d)            
%             f_curr = obj.var_sub(vars_old, f_old);
            if obj.reset_identity
                v = obj.monom(d);
                %trivial reset map (local measures)
                mom_out = mom(v);
            else
                %reset only involves x, t stays the same
                if any(obj.reset_null)
                    %TODO: make this faster. Don't generate all monomials
                    %if some of them are zero.
                    v = obj.monom(d);
                    Rv = v;
                    for i = 1:length(v)
                        pv = get_representation(v(i), obj.vars.x);
                        if (pv.pow * obj.reset_null) > 0 
                            Rv(i) = 0*v(i);
                        end
                    end
                    % Rv = subs(v, obj.vars.x(obj.reset_null), zeros(sum(obj.reset_null), 1));
                    
                    if ~isequal(obj.vars.x(~obj.reset_null), obj.reset(~obj.reset_null))
                        Rv = subs(Rv, obj.vars.x(~obj.reset_null), [obj.reset(~obj.reset_null)]);
                    end
                else
                    v = obj.monom(d);
                    Rv = subs(v, obj.vars.x, [obj.reset]);
                end
                
                mom_out = mom(Rv);
            end
        end
        
        function [mom_src, mom_dest] = liou_reset(obj, d)
            %LIOU_RESET Liouville expressions from the transition
            %mom_src:  from the source location
            %mom_dest: to the destination location
            if isempty(obj.supp)
                mom_src = 0;
                mom_dest = 0;
            else
                mom_src = -obj.mom_monom(d);
                mom_dest = obj.reset_push(d);
            end
        end
        
        function reset_out = reset_eval(obj, x)
            %output of the reset map
            reset_out =  (eval(obj.reset, obj.vars.x, x));
        end
        
        function supp_out = supp_eval(obj, t, x, tol)
            %is (t, x) in the support of the guard?
            if nargin < 4
                tol = 1e-6;
            end
            if obj.TIME_INDEP
                supp_out =  all(eval(obj.supp, obj.get_vars(), x, tol));            
            else
                supp_out =  all(eval(obj.supp, obj.get_vars(), [t; x], tol));            
            end
        end
        
        function nn_out  = nonneg(obj, t, x)
            %nonnegative dual function at this state transition
            %x comes from the space of src
            
            %check the sign convention
            %v should decrease along the jump
            if obj.reset_identity
                Rx = x;
            else
                Rx = eval(obj.reset, obj.vars.x, x);
            end
            
            vsrc = obj.src.v_eval(t, x');
            vdest = obj.dest.v_eval(t, Rx');
            
            nn_out = vsrc - vdest - obj.dual.zeno;
            
        end
        
        function obj = dual_process(obj, zeno_dual)
            %store the dual variable
            obj.dual = struct('zeno', zeno_dual, 'solved', 1);
        end
                
        
    end
end

