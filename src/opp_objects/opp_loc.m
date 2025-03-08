classdef opp_loc
    %OPP_LOC A location at a specific 
    %   Detailed explanation goes here
    
    properties
        Property1
    end
    
    methods
        function obj = opp_loc(inputArg1,inputArg2)
            %OPP_LOC Construct an instance of this class
            %   Detailed explanation goes here
            obj.Property1 = inputArg1 + inputArg2;
        end
        
        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
end

