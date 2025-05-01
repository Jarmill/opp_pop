classdef opp_diff_dyn
    %OPP_DIFF_DYN Storage of three-phase terms
    %   imposes that the three-phase terms obey dynamics as well
    %   should hopefully reduce conservatism
    
    properties
        x;
        modes;
        jumps;
        G;
        opts;
        testing = 0;
        objective = 0;
    end
    
    methods
        function obj = opp_diff_dyn(opts)
            %OPP_DIFF_DYN Construct an instance of this class
            %   Detailed explanation goes here
            if opts.three_phase ~= "Ignore"
                mpol('x_tau', 5, 1);                        
                obj.x = x_tau;
            end           
            obj.opts = opts;

            [obj.modes, obj.jumps] = obj.create_modes();
            %generate the modes and jumps
        end

        %create the modes and jumps for the three-phase system
        function [modes, jumps] = create_modes(obj)
            [levels, G] = obj.mode_switch_graph();

            %partitions of the three-phase system
            %if full-wave only, break up the symmetry structure
            P = 1+(obj.opts.Symmetry==0);

            modes = [];
            jumps = [];


        end

        %determine the switching logic
        function [V, G] = mode_switch_graph(obj)

            %determine the allowable set of voltages
            %Output:
            %   V:  three-phase voltage level configurations
            %   G:  reduced switching graph between levels
            L = obj.opts.L;
            N = length(L);
            N1 = ones(1, N);
            

            %all possible voltage triples
            va = kron(kron(L, N1), N1);
            vb = kron(kron(N1, L), N1);
            vc = kron(kron(N1, N1), L);

            V = [va; vb; vc];

            %prune the voltage set

            %unimodal: no negative va in first half-cycle
            if obj.opts.unipolar && obj.opts.Symmetry > 0
                V = V(:, va >= 0);
            end

            %common-mode constraints (typically 1/3)
            Vcm = sum(V, 1)/3;
            cm_valid = abs(Vcm) <= obj.opts.common_mode;
            V = V(:, cm_valid);

            %now form the graph of allowable switches
            %
            
            N3 = size(V, 2);

            
            G = zeros(N3);
            for i = 1:N3
                for j = 1:(i-1)
                    dv = abs(V(:, i) - V(:, j));      
                            %if max(dv) <= 1
                    %reduce the set of switches: one at a time

                    %not restrictive, because the guards are not in play on
                    %the three-phase signal w.r.t. switching on multiple
                    %phases at a time

                    if nnz(dv)==1
                        G(i, j) = 1;
                    end
                end
            end   

            G = G + G';

        end


        function [objective] = objective_diff(obj)
            %create the common-mode current

            if obj.testing==0
                Q = (eye(3) - ones(3)/3);
            else
                Q = eye(3);
            end
            
            
            xi = obj.x(3:5);
            
            % Q = ones(3);

           
            quad = (xi'*Q*xi)*(1/3);

            objective = (2*pi) * (pi)^2 * mom(quad);
            % objective = 0;
        end
        
        
        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
end

