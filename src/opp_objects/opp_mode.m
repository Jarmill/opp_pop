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
            I_ind = 3 + opts.clock;
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
                        
                        stop_pt = [stop_pt; vars.x(I_ind)==0];
                    end                    
            end


            Xstop = [stop_pt; lsupp_base.X(2:end-1)];
            Xstart = [start_pt; lsupp_base.X(2:end-1)];

            if opts.hard_stage_costs
                x = vars.x;
                stage_ind = (I_ind+1):length(x);

                %constrain the harmonics at the final time
                %must fulfill the stage cost constraints (hard)
                [nharm, bounds, types] = harm_screen(opts);
                X_harm = [];
                sym_scale = 2^double(opts.Symmetry);
                for i = 1:nharm
                    bndcurr = bounds(i, :) / sym_scale;
                    if bndcurr(1) == bndcurr(2);
                        X_harm_new = [x(stage_ind(i)) == bndcurr(1)];
                    else
                        xcurr = x(stage_ind(i));
                        X_harm_new = [(xcurr - bndcurr(1)) * (bndcurr(2) - xcurr) >= 0];
                    end
                    %TODO: URGENT: for debugging
                    % X_harm_new = [x(stage_ind)==0];
                    X_harm = [X_harm; X_harm_new];
                end
                

                Xstart = [Xstart; x(stage_ind)==0];
                Xstop = [Xstop; X_harm];
            end


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
            I_ind = (3 + opts.clock);

            if opts.hard_stage_costs
                % [n, bounds, types] = harm_screen(obj.opts);
                stage_ind = (I_ind+1):length(vars.x);
                nstage = length(stage_ind);

                
                f_stage = obj.stage_dynamics(opts, vars.x([1, 2, I_ind]));


            else
                nstage = 0;
                f_stage = [];
            end
            
            

            if imag(opts.Z_load) == 0
                x_load = [];
            else
                x_load = vars.x(I_ind:(end-nstage));
            end
            f_load = obj.load_dynamics(x_load, opts) / (2^double(obj.Symmetry));
            
            %apply the external voltage
            if ~isempty(f_load) && (opts.uext ~= 0)
                uext = opts.uext/max(max(abs(opts.L)));
                f_ext = 2*pi*[real(uext), imag(uext)]*vars.x(1:2)/ (2^double(obj.Symmetry));
                f_load = f_load + f_ext;
            end

            
            N = size(opts.L, 2);
            f_basic = [f_trig; f_clock] * ones(1, N);

            f = [f_basic; f_load; f_stage];


            %stage costs
            
        end


 
    
        %% harmonics constraints

        function f_stage = stage_dynamics(obj, opts, x)
            %running (stage) costs. the harmonics constraints and power
            %loss constraints, if hard constraints are used.
            
            %a hard version of opp_system_1.con_harmonics()

            f_stage = [];
            %x: [1, 2, 3]: [cos, sin, I]
            if opts.hard_stage_costs
            
                [n, bounds, types] = harm_screen(opts);
            
                %terms from opp_system_1.harm_eval(obj, vars, harm_in)            

                
                nmax = max(n);
                c = x(1);
                s = x(2);

                %cosine harmonics
                T = zeros(nmax+1, 1)*c;
                T(1) = 1+0*c;
                T(2) = c;
                for p = 2:nmax
                    T(p+1) = 2*c*T(p) - T(p-1);
                end

                %sin harmonics
                U = zeros(nmax, 1)*c;
                U(1) = 1+0*c;
                U(2) = 2*c;
                for p = 2:nmax
                    U(p+1) = 2*c*U(p) - U(p-1);
                end

                % Lrescale = (max(max(abs(opts.L))));
                % Lscale = 2*opts.L/Lrescale;

                Lz = zeros(size(opts.L));
                % fscale = 1;
                fscale = 2*pi * (2^double(-obj.Symmetry));


                %having nonzero dynamics here causes infeasibility
                %why?

                % f_stage = ones(length(n), length(opts.L));

                % f_stage = zeros(length(n), length(opts.L));

                for i = 1:length(n)
                    if types(i) == 1
                        %sine harmonic
                        % f_stage = [f_stage; fscale * s*U(n(i)) * opts.L];
                        f_stage = [f_stage; fscale * s*U(n(i)) * (opts.L)];

                    else
                        %cos harmonic
                        f_stage = [f_stage; fscale * T(n(i)+1) * opts.L];
                    end

                    % f_stage = [f_stage; Lz];
                end

            end


        end

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

