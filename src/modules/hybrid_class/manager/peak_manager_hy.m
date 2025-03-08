classdef peak_manager_hy < manager_hy_interface
    %PEAK_MANAGER_HY Hybrid peak estimation manager
    %   Detailed explanation goes here
    
    %TODO: rework this with subsystems and samplers
        
    methods
        function obj = peak_manager_hy(locations_in,guards_in)
            %PEAK_MANAGER_HY Construct an instance of this class
            %   Detailed explanation goes here
%             obj.Property1 = inputArg1 + inputArg2;

            obj@manager_hy_interface(locations_in, guards_in);
        end
        
        %% Formulating and solving program
        
    end
    
end

