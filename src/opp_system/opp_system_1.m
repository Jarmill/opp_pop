classdef opp_system_1 < opp_system_interface
    %OPP_SYSTEM_1 Storage for all single-phase terms
    %   includes the N-switching constraints
    
    properties
        Property1
    end
    
    methods
        function obj = opp_system_1(opts)
            %OPP_SYSTEM_1 Construct an instance of this class
            %   Detailed explanation goes here
            obj@opp_system_interface(opts)
        end
        
        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
end

