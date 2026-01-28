function [power_diss] = pulse_power_draw(out, dispatch)
%PULSE_POWER_DRAW compute the power consumption of the pulse pattern
% outputArg1 = inputArg1;
% outputArg2 = inputArg2;

top = dispatch.topology;

%conduction losses
power_cond = zeros(length(top), 1);



%switching losses
power_switch = zeros(length(top), 1);



power_diss = power_cond + power_switch;

% a = out.alpha;
% u = out.u;

% I_val = out.I_
end

