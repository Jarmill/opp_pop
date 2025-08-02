classdef opp_mode < opp_mode_interface
    %OPP_MODE Measures describing the trajectory at mode m 
    %(before switch m, or at the end of the sequence)
    %   Detailed explanation goes here        
    
    methods
        function obj = opp_mode(m, lsupp_ref, objective_mode, opts)
            %OPP_MODE Construct an instance of this class
            %   create locations

            obj@opp_mode_interface(m, lsupp_ref, objective_mode, opts);
        end

        

        %% functions used in the constructor   

        function [prefix] = get_prefix(obj)
                prefix = '';                
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
                    Theta_scale = opts.f0*opts.Ts*2^(double(opts.Symmetry));                    
                    stop_pt= [vars.x(1)==0; vars.x(2)==1];
                    if opts.clock
                        start_pt = [start_pt; vars.x(3)>=Theta_scale/2];
                        stop_pt = [stop_pt; vars.x(3)>=Theta_scale/2];
                    end
                    if imag(opts.Z_load)>0 && real(opts.Z_load)==0
                        I_ind = 3 + opts.clock;
                        stop_pt = [stop_pt; vars.x(I_ind)==0];
                    end                    
            end
            Xstop = [stop_pt; lsupp_base.X(2:end-1)];
            Xstart = [start_pt; lsupp_base.X(2:end-1)];

            mode_end = opts.k/(2^opts.Symmetry);
            if m==0                
               lsupp_base.X_init = Xstart;
            elseif m==mode_end || (opts.early_stop && ((opts.Symmetry==2) || (mod(m, 2)==0)))
                lsupp_base.X_term = Xstop;
            end
        end
 

        function flag = cannot_start(obj, opts, n)
            %can the sequence start at level n?
            flag = ((opts.start_level~=0) && (n~= opts.start_level));
        end


        function flag = cannot_end(obj, opts, n)
            %can the sequence end at level n?
            N = size(opts.L, 2);
            flag = (opts.start_level~=0) && ...
                            ((opts.Symmetry==0) && (n~= opts.start_level) ||...
                            (opts.Symmetry==1) && ((N-n+1)~= opts.start_level));
        end

        function flag = cannot_reach(obj, opts, m, n)            
            %can the sequence reach at level n at mode m?
            N = size(opts.L, 2);
            flag = false;
            if ~opts.allowed_levels(m+1, n)
                flag = true;
            elseif m>0
                can_down = (n<N) && opts.allowed_levels(m, n+1);
                can_up = (n>1) && opts.allowed_levels(m, n-1);

                flag = ~(can_up || can_down);
            end
        end

        

        function f = all_dynamics(obj, vars, opts)
            %create the dynamics as a matrix f
            %row: each level
            %column: each state
            f_trig = 2*pi*[-vars.x(2); vars.x(1)] / (2^double(obj.Symmetry));
            % f_phi = vars.x(3);
            if opts.clock
                f_clock = 1;
            else
                f_clock = [];
            end

            %TODO: check for symmetry scaling in the load
            if imag(opts.Z_load) == 0
                x_load = [];
            else
                x_load = vars.x((3 + opts.clock):end);
            end
            f_load = obj.load_dynamics(x_load, opts) / (2^double(obj.Symmetry));
            
            % f = [f_trig; f_phi; f_load];
            N = size(opts.L, 2);
            f_basic = [f_trig; f_clock] * ones(1, N);

            f = [f_basic; f_load];
        end


 
    
        %% harmonics constraints
        function harm = voltage_harmonics_mom(obj, vars, harm_mon)
            %voltage harmonics constraints
            harm= mom(vars.x(1))*zeros;
            
            [N, P] = size(obj.levels);

            for n=1:N
                for p = 1:P            
                    [~, harm_base] = obj.levels{n, p}.voltage_harmonics_mom(vars, harm_mon);                                                                   
                    
                    harm = harm + harm_base;                    
                end
            end
        end


        function harm = load_harmonics_mom(obj, vars, harm_mon, harm_in)
            %voltage harmonics constraints
            % harm= mom(p)*0;
            Lmax = max(obj.L);

            % Z_type = 0;
            
            if (length(vars.x)==3) || (imag(obj.Z_load) == 0)                      
                %purely resistive
                % harm = obj.voltage_harmonics_mom(vars, harm_in);
                Z_type = 0;
                Z_scale = 1;
            else                
                if (imag(obj.Z_load) >= 0)
                    
                    %inductive load
                    %i' = -(R/L)i + (1/L) v
                    inductance = imag(obj.Z_load)/(2*pi*obj.f0);
                    Z_type = 1;
                    Z_scale = (Lmax/inductance);
                    % harm = (harm_eval.*vars.x(4)) *  .*(obj.opts.L);
                else
                     %vc' = (v-vc)/(R*C)
                     %per-unit, ignore (R*C) factor
                     %TODO: v is from the voltage source. Modify when it is 
                     %filtered by a grid-side filter                    
                     capacitance= -imag(obj.Z_load)*(2*pi*obj.f0);
                     resistance= real(obj.Z_load);                       
                     RC = resistance*capacitance;
                     Z_type = 2;
                     Z_scale = (Lmax/RC);                     
                end
            end

            for n=1:N
                for p = 1:P            
                    [~, harm_mom] = obj.levels{n, p}.load_harmonics_mom(obj, vars, harm_mon, harm_in, Z_type, Z_scale);       
                    harm = harm+harm_mom;                    
                end
            end
        end

       function om = occ_mass(obj)
            %mass of the occupation measure in this mode
            om = 0;
            [N, P] = size(obj.levels);
            for n=1:N
                for p = 1:P 
                    om = om + obj.levels{n, p}.occ_mass();
                end
            end
        end

        
    end


end

