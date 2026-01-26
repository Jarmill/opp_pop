classdef opp_component
    %OPP_COMPONENT Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        id; 
        GCT; %true: GCT, false: freewheeling diode
        conduction; 
        switch_on; 
        switch_off; 
    end
    
    methods
        function obj = opp_component(id,  GCT, conduction,  switch_on, switch_off)
            %OPP_COMPONENT Construct an instance of this class
            %   Detailed explanation goes here
            obj.id = id;
            opts.GCT = GCT;
            obj.conduction = conduction; 
            obj.switch_on = switch_on;
            obj.switch_off = switch_off;
        end                

    end
end

