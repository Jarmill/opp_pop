classdef opp_harmonics
    %OPP_HARMONICS Summary of this class goes here
    %   Detailed explanation goes here
    
%example (cos(0t)=0, cos(1t)=0, sin(1t)=1, 0<=sin(2t)<=0.1)
%index_cos = [0; 1]
%bound_cos = [0, 0; 0, 0]
%index_sin = [1; 2]
%bound_sin = [1, 1; 0, 0.1]

    properties
        index_cos = [0; 1];
        bound_cos = [0, 0; 0, 0];
        index_sin = [1];       
        bound_sin = [1, 1];
    end  
end

