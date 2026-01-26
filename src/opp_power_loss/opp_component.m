classdef opp_component
    %OPP_COMPONENT Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        id; 
        type
        conduction;
        loss;
    end
    
    methods
        function obj = opp_component(id, type, conduction, loss)
            %OPP_COMPONENT Construct an instance of this class
            %   Detailed explanation goes here
            obj.id = id;
            obj.type = type;
            obj.conduction = conduction; 
            obj.loss = loss;
        end
        
        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
end

