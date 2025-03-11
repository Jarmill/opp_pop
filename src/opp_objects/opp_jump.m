classdef opp_jump < handle
    %OPP_JUMP measures describing the switch from m-1 to m
    %(m=1 to k)
    
    properties
       mode;          %the mode m (defines the id)
       opts;          %relevant options for the mode
       jump_up;       %guards for transiting up a level
       jump_down;     %guards for transiting down a level
       src;           %source mode 
       dst;           %destination mode
       L;             %levels of the inverter               
       vars;          %basic variable type       
    end
    
    methods
        function obj = opp_jump(m, opts, vars, X_jump)
            %OPP_JUMP Construct an instance of this class
            %   Detailed explanation goes here
            
            obj.vars = vars;
            obj.L = opts.L;

            reset_law = vars.x;
            reset_law(3) = 0*vars.x(3);

            P = opts.partition;
            N = length(obj.L);

            obj.mode = m;
            X_partition = support_partition(opts.partition, vars.x, opts.Symmetry);
            
            obj.jump_up = cell(N-1, P);
            obj.jump_down = cell(N-1, P);

            for n=1:N-1
                for p = 1:P
                    if P > 1
                        X_p = X_partition(p)>=0;
                    else
                        X_p = [];
                    end
                    supp_curr = [X_jump; X_p];
                    % curr_name = sprintf('jump_m%d_n%d_p%d', m, n, p);
                    name_down = sprintf('down_m%d_n%d_p%d', m, n, p);
                    name_up = sprintf('up_m%d_n%d_p%d', m, n+1, p);
                    obj.jump_up{n, p} = guard(name_up, vars, [], ...
                        [], supp_curr, reset_law);
                    obj.jump_down{n, p} = guard(name_down, vars, [], ...
                        [], supp_curr, reset_law);
                end
            end
            
        end
        
        function [mom_src, mom_dst] = liou_reset(obj, d)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            [Np, P] = size(obj.jump_up);
            N = Np+1;

            %moments of the guards going up and down in the transition
            mom_src = cell(N, P);
            mom_dst = cell(N, P);


            for n =1:N
                for p=1:P
                    mom_src{n, p} = 0;
                    mom_dst{n, p} = 0;
                end
            end
            
            for p = 1:P
                for n=1:Np
                    [mom_src_up, mom_dst_up] = obj.jump_up{n, p}.liou_reset(d);
                    [mom_src_down, mom_dst_down] = obj.jump_down{n, p}.liou_reset(d);
                    
                    % if n<Np
                        mom_src{n, p} = mom_src{n, p} + mom_src_up;
                        mom_dst{n+1, p} = mom_dst{n+1, p} + mom_dst_up;
                    % end
                    % if n>1
                        mom_src{n+1, p} = mom_src{n+1, p} + mom_src_down;
                        mom_dst{n, p} = mom_dst{n, p} + mom_dst_down;
                    % end
                end
            end
        end

        % function 
        %TODO: add switching losses as a constraint

        function supp_con_out = supp_con(obj)
            %fetch all support constraints
            [Np, P] = size(obj.jump_up);
            supp_con_out = [];
            for n=1:Np
                for p = 1:P                    
                    curr_up = obj.jump_up{n, p}.supp;
                    curr_down = obj.jump_down{n, p}.supp;
                    supp_con_out = [supp_con_out; curr_up; curr_down];
                end
            end
        end

        function m_out = mmat(obj)
            %get moment matrices of the jumps
            m_out = struct;
            [Np, P] = size(obj.jump_up);
            % N = Np+1;
            m_out.up = cell(Np, P);
            m_out.down = cell(Np, P);
            for n=1:Np
                for p = 1:P                    
                    m_out.up = obj.jump_up{n, p}.mmat();
                    m_out.down = obj.jump_down{n, p}.mmat();                    
                end
            end
        end

        function obj_out = objective(obj)
            %objective for the jump
            %currently null
            obj_out = 0;
        end
    end
end

