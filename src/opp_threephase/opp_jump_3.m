classdef opp_jump_3
    %OPP_JUMP_3 describes jumps between the modes of a three-phase
    %dynamical system.
    % 
    %this implementation ignores the clock
    
    properties
       mode;          %the mode m (defines the id)
       opts;          %relevant options for the mode
       jump;       %guards for transiting up a level       
       src;           %source mode 
       dst;           %destination mode
       L;             %levels of the inverter  
       L_single;      %three-phase levels of the inverter
       vars;          %basic variable type    
       orbits;        %cyclic symmetry for the guards
    end
    
    methods
        %% constructor and related methods
        function obj = opp_jump_3(lsupp_ref, opts_3, G)
            %OPP_JUMP_3 Construct an instance of this class
            %   Detailed explanation goes here
            obj.L = opts_3.L;
            obj.L_single = opts_3.L_single;
            % obj.G = G;
            [obj.src, obj.dst] = find(G);     

            


            obj.opts = opts_3;
            obj.vars = lsupp_ref.vars;
            
            L_aug = [obj.L(:, obj.src); obj.L(:, obj.dst)];
            obj.orbits = get_orbits(L_aug);

            obj.jump = obj.create_jumps(lsupp_ref);


            
        end
        
        function jumps = create_jumps(obj, lsupp_ref)
            %CREATE_JUMPS fill in the cell of possible jumps
            
            reset_law = obj.vars.x;

            P = obj.opts.partition;
            N = size(obj.L, 2);

            X_partition = support_partition(P, obj.vars.x, 0);

            % E = nnz(G);
            E = length(obj.src);

            jumps = cell(E, P);
            X_jump = lsupp_ref.X;

            for e = 1:E
                for p = 1:P
                    if P > 1
                        X_p = X_partition(p)>=0;
                    else
                        X_p = [];
                    end

                    supp_curr = [X_jump; X_p];

                    name = sprintf('tau_jump_n_%d_%d_p_%d', obj.src(e), obj.dst(e), p);

                    jumps{e, p} = guard(name, obj.vars, [], [], ...
                        supp_curr, reset_law);

                end
            end
        end

        %% liouville constraints
        function [mom_src, mom_dst] = liou_reset(obj, d)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            [E, P] = size(obj.jump);


            N3 = size(obj.L, 2);
            

            E = length(obj.src);

            %moments of the guards going up and down in the transition
            mom_src = cell(N3, P);
            mom_dst = cell(N3, P);


            for e =1:N3
                for p=1:P
                    mom_src{e, p} = 0;
                    mom_dst{e, p} = 0;
                end
            end
            
            for p = 1:P
                for e=1:E
                    [mom_src_curr, mom_dst_curr] = obj.jump{e, p}.liou_reset(d);
                                                            
                    mom_src{obj.src(e), p} = mom_src{obj.src(e), p} + mom_src_curr;                   
                    mom_dst{obj.dst(e), p} = mom_dst{obj.dst(e), p} + mom_dst_curr;
                    % end
                end
            end
        end

        function supp_con_out = supp_con(obj)
            %fetch all support constraints
            [E, P] = size(obj.jump);
            supp_con_out = [];
            for n=1:E
                for p = 1:P                    
                    curr= obj.jump{n, p}.supp;                    
                    supp_con_out = [supp_con_out; curr];
                end
            end
        end

        %% recovery
         function m_out = mmat(obj)
            %get moment matrices of the jumps
            m_out = struct;
            [Np, P] = size(obj.jump);
            % N = Np+1;
            m_out.jumps= cell(Np, P);
            
            for n=1:Np
                for p = 1:P                    
                    m_out.jumps{n, p} = obj.jump{n, p}.mmat();                    
                end
            end
        end


        function m_out = mmat_corner(obj)
            %get moment matrices of the jumps
            m_out = struct;
            [Np, P] = size(obj.jump);
            % N = Np+1;
            m_out.jumps= cell(Np, P);            
            for n=1:Np
                for p = 1:P                                        
                    m_out.jumps{n, p} = obj.jump{n, p}.mmat_corner();                    
                end
            end
        end


    end
end

