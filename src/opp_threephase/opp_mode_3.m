classdef opp_mode_3 < opp_mode_interface
    %OPP_MODE_3 Summary of this class goes here
    %   Detailed explanation goes here       
    properties
        L_single; 
    end

    methods
        function obj = opp_mode_3(lsupp_ref, objective_mode, opts_3)
            %OPP_MODE_3 Construct an instance of this class
            %   Detailed explanation goes here                      

            obj@opp_mode_interface(0, lsupp_ref, objective_mode, opts_3)
            obj.L = opts_3.L;
            obj.L_single = opts_3.L_single;

        end

        function [prefix] = get_prefix(obj)
            prefix = 'tau_';                
        end


        function f = all_dynamics(obj, vars, opts)
            %create the dynamics as a matrix f
            %row: each level
            %column: each state
            f_trig = 2*pi*[-vars.x(2); vars.x(1)];
            
            %TODO: check for symmetry scaling in the load            
            x_load = vars.x(3:5);

            f_load = load_dynamics(obj, x_load, opts);
            %check the dimensions here
                        
            N = size(opts.L, 2);
            f_basic = [f_trig] * ones(1, N);

            f = [f_basic; f_load];
        end
        

        %pruning levels for reachability
        %flag=1: can't reach
        function flag = cannot_start(obj, opts, n)
            %can the sequence start at level n?
            va = opts.L(1, n);
            
            flag = ((opts.start_level~=0) && (va~= opts.L_single(1, opts.start_level)));
        end


        function flag = cannot_end(obj, opts, n)
            %can the sequence end at level n?
            N = size(opts.L, 2);
            va = opts.L(1, n);            
            
            flag = (opts.start_level~=0) && ...
                            ((opts.Symmetry==0) && (va~=opts.L_single(1, opts.start_level)) ||...
                            (opts.Symmetry==1) && (va~= opts.L_single(1, N-opts.start_level+1)));
        end

        function flag = cannot_reach(obj, opts, m, n)
            %can the sequence reach at level n at mode m?
            flag = false;            
        end

        function lsupp_base = refine_support(obj, opts, m,  lsupp_ref)
             %define the generic support set
            lsupp_base = lsupp_ref;

            vars = lsupp_base.vars;           
         
            %define the terminal set
            %terminate
            %ignore the trig constraint (beginning) and support arc
            %constraint (end)
            start_pt = [vars.x(1)==1; vars.x(2)==0];
            switch obj.Symmetry
                case 0
                    %full-wave symmetry: end at 2pi
                    stop_pt = start_pt;
                case 1
                    %half-wave symmetry: end at pi
                    stop_pt = [vars.x(1)==-1; vars.x(2)==0];
                case 2
                    %quarter-wave symmetry: end at pi/2
                    stop_pt = [vars.x(1)==-1; vars.x(2)==0];                  
            end
            Xstop = [stop_pt; lsupp_base.X(2:end-1)];
            Xstart = [start_pt; lsupp_base.X(2:end-1)];

            
            %three-phase: go all the way around
            lsupp_base.X_init = Xstart;           
            lsupp_base.X_term = Xstop;
        
        end
    end
end

