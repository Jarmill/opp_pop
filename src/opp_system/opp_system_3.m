classdef opp_diff_dyn < opp_system_interface
    %OPP_DIFF_DYN Storage of three-phase terms
    %   imposes that the three-phase terms obey dynamics as well
    %   should hopefully reduce conservatism
    
    properties
        t; 
        x;
        mode = [];
        jumps = [];
        G;
        opts;
        testing = 0;
        objective = 0;
        correspond = [];
        DYNAMICS = true;
    end
    
    methods
        function obj = opp_diff_dyn(opts)
            %OPP_DIFF_DYN Construct an instance of this class
            %   Detailed explanation goes here

            obj@opp_system_interface(opts)

            if (opts.three_phase ~= "Ignore") || (opts.common_mode < Inf)
                obj.correspond = (1:length(obj.opts.L_single)) * (obj.opts.L_single' ==obj.opts.L(1, :));
                obj.mode.correspond = obj.correspond;               
            end

            %generate the modes and jumps
        end

        function [t, x] = create_vars(obj, opts)
            if (opts.three_phase ~= "Ignore") || (opts.common_mode < Inf)
                mpol('c_tau', 1, 1)
                mpol('s_tau', 1, 1)
                mpol('I_tau', 3, 1)
                x_tau = [c_tau; s_tau; I_tau];
                x = x_tau;

                if opts.TIME_INDEP
                   t = [];
                else
                    mpol('t_tau', 1, 1)
                    t = t_tau;
                end
            else
                t = [];
                x = [];
            end
        end

        function I_marg = get_I_marginals(obj, d)
            %maybe this isn't needed?
            I_marg = obj.mode.get_I_marginals(d);            
        end

        function I_marg = get_I_marginal(obj, d)
            %maybe this isn't needed?
            I_marg = obj.mode.get_I_marginal(d);            
        end

        %create the modes and jumps for the three-phase system
        function [mode, jumps, opts_3] = create_system(obj)
            if (obj.opts.three_phase ~= "Ignore") || (obj.opts.common_mode < Inf)
                [L3, G] = obj.mode_switch_graph();
    
                %partitions of the three-phase system
                %if full-wave only, break up the symmetry structure
                
                %create the mode
                lsupp_base = obj.get_support();                                    
    
                opts_3 = obj.opts;
                opts_3.partition = 1 + (opts_3.Symmetry==0);   
                opts_3.L_single = opts_3.L;
                opts_3.L = L3;
                opts_3.allowed_levels = ones(1, size(L3, 2));
                N = size(L3, 2);
    
                mode = opp_mode_3(lsupp_base, obj.objective*ones(N, 1), opts_3);
               
                %create the jump
                jumps = opp_jump_3(lsupp_base, opts_3, G);
            else
                mode = [];
                jumps = [];
                opts_3 = obj.opts;

            end
        end


        %% support constraints

        function sc = supp_con_base(obj)
            %support constraint
            sc = [sum(obj.x(1:2).^2) == 1; 
                obj.x(3:5).^2 <= 1];

            if obj.opts.three_phase == "Balanced"
                %produce a balanced current
                sc = [sc; obj.x(5) == (-obj.x(3) - obj.x(4))];
            end
        end
        
        %% determine the switching logic
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

                    if obj.opts.common_mode == 0
                        edge_con = max(dv) <= 1;
                    else
                        edge_con = nnz(dv)==1;
                    end
                    if edge_con
                        G(i, j) = 1;
                    end
                end
            end   

            G = G + G';

        end


        

        %% moment material
        function [objective] = create_objective(obj)
            %create the common-mode current

            if obj.testing==0
                Q = (eye(3) - ones(3)/3);
            else
                Q = eye(3);
            end
            
            
            xi = obj.x(3:5);
            
            % Q = ones(3);

           
            quad = (xi'*Q*xi)*(1/3);

            % objective = (2*pi) * (pi)^2 * mom(quad);
            objective = (2*pi) * (pi)^2 * quad;            
            % objective = 0;
        end
        

        

        %% moment constraints (internal)
        %TODO: fill this in

        function [mom_con, supp_con] = cons(obj, d)
            %get all constraints involving only this structure

            if isempty(obj.mode)
                mom_con = [];
                supp_con = [];
            else
                supp_con = obj.supp_con();
    
                
                % con_liou = [];
                con_preserve = obj.con_return(d);
                con_prob = obj.con_prob_dist();
                con_liou = obj.con_flow(d);

                con_sym = obj.con_rotate_symmetry(d);
    
                %TODO: internal marginal constraints (ensure three-phase
                %symmetry in the occupation and jump measures)
    
                mom_con = [con_liou; con_preserve; con_prob];            
            end
        end

       function mass_con_eq = con_prob_dist(obj)
            %initial measure is a probability distribution (mass 1)
            
            [~, mass_init_sum] = obj.mode.initial_mass();

            mass_con_eq = (mass_init_sum==1);
           
        end

       function flow_con = con_flow(obj, d)
            %the flow conservation constraint for dynamics

            %compute all terms
            [jump_src, jump_dst] = obj.jumps.liou_reset(d);
            liou = obj.mode.flow(d);

            %add the constraints
            [N, P] = size(liou);
            
            %iterate over all cells
            flow_con = [];
            flow_con_cell = cell(N, P);
            for n=1:N
                    
                for p = 1:P                                                       
                    flow_con_cell{n, p} = liou{n, p} + jump_src{n, p} + jump_dst{n, p}==0;
                    
                    %stack them into a giant vector: flow_con
                    flow_con = [flow_con; flow_con_cell{n, p}];
                end                
            end
        end
           
        function con_sym = con_rotate_symmetry(obj, d)
            %ensure that the three-phase occupation measures satisfy
            %symmetries under rotations

            %TODO: figure out the math of this. then implement it.
            %should hopefully reduce conservatism of the approach
            
            con_sym =[];
        end

       function return_con = con_return(obj, d)
            %conservation of position between the initial and final measure
            
            % mass_con = obj.modes{1}.mass_init_mode();
            init_monom = obj.mode.init_monom(d, true);
            
            %TODO:  fix the quarter-wave structure
                                      
                return_mom = init_monom;

                %index the terminal destination levels based on the
                %applied symmetry
                %unconstrained for quarter-wave symmetry (here at least)
                flip_load = 2*(obj.opts.Symmetry>0);  
                              

                stop_monom = obj.mode.term_monom(d, true, flip_load);
                return_mom = madd_cell_mom(return_mom, stop_monom, -1);
                
                return_mom_1 = return_mom(:, 1);
                [N, P] = size(return_mom);
                return_con = [];

                for n = 1:N
                    for p = 2:P
                        return_mom_1{n, 1} = return_mom_1{n, 1} + return_mom{n, p};
                    end
                end

                for n = 1:N
                    if ~isnumeric(return_mom{n})
                        return_con = [return_con; return_mom_1{n}==0];
                    end
                end

            % return_con = (return_mom==0);
        
        end
       
        %% moment constraints (external)
        %constraints to ensure alignment with the k-and-clock-constrained
        %switching sequence

        function sel_out = sel_init_monom(obj, d, ind)
            %used for alignment
            sel_orig = obj.mode.sel_init_monom(d, ind);
            N = length(obj.opts.L_single);
            sel_out = cell(N, 1);            
            for i = 1:N
                sel_out{i} = 0;
            end
            NN = size(obj.opts.L, 2);
            for v = 1:NN
                corr_i = obj.correspond(v);
                sel_out{corr_i} = sel_out{corr_i} + sel_orig{v, 1};
            end
        end


        function sel_out = sel_term_monom(obj, d, ind)
            %used for alignment
            sel_orig = obj.mode.sel_term_monom(d, ind);
            N = length(obj.opts.L_single);
            sel_out = cell(N, 1);            
            for i = 1:N
                sel_out{i} = 0;
            end
            NN = size(obj.opts.L, 2);
            for v = 1:NN
                corr_i = obj.correspond(v);
                sel_out{corr_i} = sel_out{corr_i} + sel_orig{v, end};
            end
        end

        %% recovery
        function [m_out, tr_out, j_out] = mmat_corner(obj)

            if ~isempty(obj.mode)
                [m_out, tr_out]= obj.mode.mmat_corner();
    
                mj = obj.jumps.mmat_corner();
                j_out = mj.jumps();
            else
                m_out = [];
                tr_out = [];
                j_out = [];
            end
        end

       function [m_out, tr_out, j_out] = mmat(obj)

            if ~isempty(obj.mode)
                [m_out, tr_out]= obj.mode.mmat();
    
                mj = obj.jumps.mmat();
                j_out = mj.jumps();
            else
                m_out = [];
                tr_out = [];
                j_out = [];
            end
       end

       function ms = mass_summary(obj)
           ms = struct;
           [N, P] = size(obj.mode.levels);
           ms.mode = zeros(N, P);
           ms.trans = zeros(N, P-1);
            for n=1:N
                for p = 1:P
                    ms.mode(n, p) = double(obj.mode.levels{n, p}.sys{1}.meas_occ.mass());
                    if p < P
                        ms.trans(n, p) = double(obj.mode.transition{n, p}.mass());
                    end
                end
            end       

            E = length(obj.jumps.src)
            ms.jump = zeros(E, 2);
            
            
            for e = 1:E
                for p = 1:2
                    ms.jump(n, p) = double(obj.jumps.jump{n, p}.mass());                    
                end
            end
        
% 
            % [c, mom_harm] = obj.con_harmonics();
            % ms.harm = double(mom_harm);

       end

    end
end

