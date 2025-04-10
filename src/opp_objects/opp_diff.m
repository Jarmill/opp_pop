classdef opp_diff
    %OPP_DIFF processes the differential-mode tdd objective    
    
    properties
        xtrig = [];
        tau = [];
        soc = {};
        testing = 0;
    end
    
    methods
        function obj = opp_diff(floating)
            %OPP_DIFF Construct an instance of this class
            %   Detailed explanation goes here
            if nargin >= 1 || floating == false
                mpol('x_tau', 2 + obj.testing, 1);
                obj.tau = meas(x_tau);         
                obj.xtrig = x_tau;
            end
        end
        
        function m_out = mass(obj)
            m_out = mass(obj.tau);
        end


        function sc = supp_con(obj)
            %support constraint
            sc = [sum(obj.xtrig.^2) == 1];
        end

        function [obj_con] = objective_diff(obj, d, mom_diff)
            %create the soc-bounding for differential-mode distortion
            %mom_diff: moments of I*(w(c, s)) across the measures

            Nsoc = size(mom_diff, 1);
            

            tau_mom = mom(mmon(obj.xtrig, 0, d-1));
            
            obj.soc = cell(Nsoc, 1);
            
            obj_con = [];
            for i = 1:Nsoc
                xname= ['x_soc_', num2str(i)];                
                mpol(xname, 2 + obj.testing, 1);
                xcurr = eval(xname);
                meas_curr = meas(xcurr);
                obj.soc{i} = meas_curr;

                %form an SOC-type constraint               
                %[<w, tau>, <w Ia, mu>, <w Ib, mu>;
                % <w Ia, mu>, 1, 0;
                % <w Ib, mu>, 0, 1] >= 0

                %<w, tau> >= <w Ia, mu>^2 + <w Ib, mu>^2
                %
                %warning! achtung! bad code! bad code!
                
                obj_con_curr = [mass(meas_curr)==tau_mom(i);
                    mom(xcurr(1)) == mom_diff(i, 1); 
                    mom(xcurr(2)) == mom_diff(i, 2);
                    mom([xcurr(1)^2; xcurr(2)^2; xcurr(1)*xcurr(2)]) == [1; 1; 0];
                ];

                if obj.testing
                    obj_con_curr = [obj_con_curr;
                        mom(xcurr(3)) == mom_diff(i, 3);
                        mom(xcurr(3)^2) == 1;
                        mom(xcurr(3)*xcurr(1)) == 0;
                        mom(xcurr(3)*xcurr(2)) == 0];
                end

                obj_con = [obj_con; obj_con_curr];
            end
        end
    end
end

