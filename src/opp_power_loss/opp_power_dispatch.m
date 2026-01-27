classdef opp_power_dispatch
    %OPP_POWER_DISPATCH Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        topology;

        %default to the 3-level topology
        conduction_pos = {[7, 8], [2, 9], [1, 2]};
        conduction_neg = {[3, 4], [3, 10], [5, 6]};
        
        
        jump_up_pos_on = {[2],[1]};
        jump_up_pos_off = {[8],[9]};
        
        jump_up_neg_on = {[], []};
        jump_up_neg_off = {[4], [3]};
        
        jump_down_pos_on = {[],[]};
        jump_down_pos_off = {[2],[1]};
        
        jump_down_neg_on = {[4], [3]};
        jump_down_neg_off = {[10], [5]};
        

        Vdc = 5000;
        vT = 2400; %max rated voltage (V)
        IT = 4500; %max rated current (A)
        I_rated = 2200; %nominal rated current (A)
    end
    
    methods
        function obj = opp_power_dispatch(topology_in)
            %OPP_POWER_DISPATCH Construct an instance of this class
            %   Detailed explanation goes here
            
            if nargin == 0
                cond_gct = [0.97, 0.245e-3];
                cond_diode = [1.19, 0.395e-3];
                
              
                on_gct = [0, 0.095e-3];
                off_gct = [0, 2.6e-3];
                on_diode = 0;
                off_diode = [15.2, 0, 0, 0]; %reverse recovery, make this fancier later
                
                
                
                obj.topology = cell(10, 1);
                %figure 4
                for i = 1:4
                    obj.topology{i} = opp_component(i, true, cond_gct, on_gct, off_gct);
                end
                for i = 5:10
                    obj.topology{i} = opp_component(i, false, cond_diode, on_diode, off_diode);
                end
            else
                obj.topology = topology_in;
            end

        end

        function [coeff_pos, coeff_neg] = get_jump_coeff(m)
            %m is the mode/index for the jump

            coeff_pos = [0,0];

            coeff_neg = [0,0];
        end
    end
end

