classdef opp_diff_current_split
    %OPP_DIFF_CURRENT_SPLIT processes the differential-mode tdd objective    
    %
    %
    %creates 3 different measures 
    %tau1: [c, s, Ia, Ib]
    %tau2: [c, s, Ib, Ic]
    %tau3: [c, s, Ia, Ic]
    %
    %
    %aligns marginals between the measures
    %does quadratic objectives. 
    %Only 4 states, so more computationally friendly
    %

    properties
        x = {};
        tau = {};        
        testing = 0;
    end
    
    methods
        function obj = opp_diff_current_split(floating)
            %OPP_DIFF_3 Construct an instance of this class
            %   Detailed explanation goes here
            obj.tau = cell(3, 1);
            obj.x = cell(3, 1);
            if nargin >= 1 || floating == false
                for i =1:3
                    xname = sprintf('x_tau_%d', i);
                    mpol(xname, 4, 1);
                    x_tau = eval(xname); %bad bad bad bad bad

                    obj.tau{i} = meas(x_tau);         
                    obj.x{i} = x_tau;
                end
            end
        end
        
        function m_out = mass(obj)
            m_out = sum(cellfun(@mass, obj.tau));
        end


        function sc = supp_con(obj)
            %support constraint
            sc = [];
            for i = 1:3
                xcurr = obj.x{i};
                sc_curr = [sum(xcurr(1:2).^2)==1; xcurr(3:4).^2 <= 1];
                sc = [sc; sc_curr];
            end
        end     

        function [objective] = objective_diff(obj)
            %create the common-mode current
            mquad = 0;
            for i = 1:3
                xcurr = obj.x{i};

                if obj.testing
                    mquad = mquad + mom(sum(xcurr(3:4).^2)*(1/6));
                else
                    quad_diag = mom(sum(xcurr(3:4).^2)*(1/9));
                    
                    mquad = mquad + quad_diag;
                    
                    quad_off = mom(xcurr(3)*xcurr(4))*(2/9);
                    mquad = mquad - quad_off;
                
                end
                    
            end


            objective = (2*pi) * (pi)^2 * mquad;            

            % xab = obj.x{i};
            % Ia = xab(3);
            % objective = (2*pi) * (pi)^2 * mom(Ia.^2);        
        end

        function marg_con = con_diff(obj, d, three_phase_mom)
            %form the three-phase 

            % ind = [1 2; 1 3; 2 3];
            ind = [1 2; 2 3; 1 3];

            marg_con = [];
            for i = 1:3
                xcurr = obj.x{i};
                v1 = mmon(xcurr([1, 2, 3]), d);
                v2 = mmon(xcurr([1, 2, 4]), d);

                marg_v = [mom(v1), mom(v2)];

                three_phase_curr = three_phase_mom(:, ind(i, :));

                marg_con = [marg_con, (marg_v == three_phase_curr)];
            end
            % 
            % va = mmon(obj.x([1, 2, 3]), d);
            % vb = mmon(obj.x([1, 2, 4]), d);
            % vc = mmon(obj.x([1, 2, 5]), d);
            % 
            % marg_v = [mom(va), mom(vb), mom(vc)];
            % 
            % marg_con = (marg_v == three_phase_mom);
            marg_con = reshape(marg_con, [], 1);
        end
    
        %fetch the recovered entries
        function M = mmat(obj)
            if isempty(obj.tau)
                M = [];
            else
                M = cell(3, 1);
                for i = 1:3
                    M{i} = double(mmat(obj.tau{i}));
                end
                
            end
        end

        function Mc = mmat_corner(obj)
            if isempty(obj.tau)
                Mc = [];
            else                
                Mc = cell(3, 1);
                for i = 1:3
                    xcurr = obj.x{i};
                    Mc{i} = double(mom([1;xcurr]*[1; xcurr]'));
                end
                
            end
        end
    end
end

