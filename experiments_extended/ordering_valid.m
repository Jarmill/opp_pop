function [valid] = ordering_valid(pc, Theta)
%ORDERING_VALID Summary of this function goes here
%   Detailed explanation goes here
% outputArg1 = inputArg1;
valid = false;
if ~isempty(pc.warm)
    af = pc.warm.alpha_q;
    acirc = [af; 2*pi + af(1)];
    dacirc = diff(acirc);
    con_order = dacirc >= Theta;

    valid = all(con_order);
end
end

