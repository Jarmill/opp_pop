classdef (Abstract) opp_system_interface
    %OPP_SYSTEM_INTERFACE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        t; 
        x;
        mode; 
        jumps;
        opts;
        objective;
    end
    
    methods
        function obj = opp_system_interface(opts)
            %OPP_SYSTEM_INTERFACE Construct an instance of this class
            %   Detailed explanation goes here
            obj.opts = opts;
            [obj.t, obj.x] = obj.create_vars(opts);
            obj.objective = obj.create_objective();
            [obj.mode, obj.jumps, obj.opts] = obj.create_system();    
        end         

    %% constraints
        function sc = supp_con(obj)
            %support of all measures in the three-phase assembly

            if isempty(obj.mode)
                sc = [];
            else                
                sc_mode = obj.mode.supp_con();
                sc_jump = obj.jumps.supp_con();
                sc = [sc_mode; sc_jump];
            end
        end

        function lsupp_base = get_support(obj)
            %get the generic support of the mode measures
            lsupp_base = loc_support();
            lsupp_base.vars.x = obj.x;
            lsupp_base.vars.t = obj.t;

            lsupp_base.TIME_INDEP = obj.opts.TIME_INDEP;
            lsupp_base.FREE_TERM = 0;
            lsupp_base.Tmax = 1;

            lsupp_base.X = obj.supp_con_base(); 
            
        end

    end

    %% abstract methods (to be overloaded)
    methods (Abstract)
        %constructor
        create_system(obj)
        create_objective(obj)

        %constraints
        

        %recovery
        mmat(obj)
        mmat_corner(obj)
        mass_summary(obj)
    end
end

