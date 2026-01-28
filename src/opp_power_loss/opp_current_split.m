classdef opp_current_split
    %OPP_CURRENT_SPLIT Store the positive and negative parts of the load
    %current. used for semiconductor losses.
    
    properties
        % Property1
        % id;
        pos;
        neg;
    end
    
    methods
        function obj = opp_current_split(loc_supp, prefix)
            %CURRENT_SPLIT Construct an instance of this class
            %   Detailed explanation goes here



            %get the load current, and only the load current
            if nargin < 2
                prefix = [];
            end

            if ~isempty(loc_supp)                
                obj.pos = obj.meas_def([prefix, '_pos'], 1);
                obj.neg = obj.meas_def([prefix, '_neg'], -1);
            end
        end

        %% measure access
        function mmmon_out = mom_monom(obj, dmin, dmax)
            %MOM_MONOM moments of monomials
            if nargin < 3
                dmax = dmin;
                dmin = 0;
            end
            
            if isempty(obj.pos)
                mmmon_out = 0;
            else
                mmmon_out = obj.pos.mom_monom(dmin, dmax) + obj.neg.mom_monom(dmin, dmax);
            end
        end  

        function mmmon_out = mom_lin(obj)
            %MOM_MONOM moments of monomials           
            if isempty(obj.pos)
                mmmon_out = 0;
            else
                mmmon_out = [obj.pos.mom_monom(1,1), obj.neg.mom_monom(1,1)];
            end
        end  


        function supp_out = supp(obj)
            %return supports of the measures
            if isempty(obj.pos)
                supp_out = [];
            else
                supp_out = [obj.pos.supp; obj.neg.supp];
            end
        end

         %% Measure Creation
         function [meas_sgn] = meas_def(obj, suffix, sign)           
            %MEAS_DEF Define the measures in the collection
            %declare a variable for each measure (index ind in the union)

            vars_new = struct;
            old_stack =[];
            new_stack = [];

            new_name = ['I_', suffix];

            mpol(new_name, 1, 1);   
            I_var = eval(new_name);
            vars = struct('x', I_var);

            %define the measure
            
            meas_sgn = meas_base(vars, [sign*vars.x >= 0]);
        end  
    end
end

