function [out] = pulse_current_voltage(u, alpha, sym, I0)
%Find the energy of the pulse function.
%
%energy = integrate x(th)^2 dth for th in [0, 2pi]

%replicate the voltages and angles
if sym == 0
    %full-wave symmetry
    uf = u;
    af = alpha;
elseif sym==1
    %half-wave symmetry
    uf = [u; -u(2:end)];
    af = [alpha; alpha+pi];
else
    %quarter-wave symmetry
    urev = u(end-1:-1:1);
    uf = [u; urev; -u(2:end); -urev];
    af= [alpha; pi - (alpha(end:-1:1)); pi + alpha; 2*pi - alpha(end:-1:1)];
end

ah = [0; af; 2*pi];
da = diff(ah);

out = struct;
out.alpha = ah;
out.voltage = uf;

%compute the current
% I0 = 0;
I_step = uf.*da;
I0_val = cumsum([0; I_step]);
if nargin < 4
    mean_I = mean(I0_val);
    I_val = I0_val - mean_I;
else
    I_val = I_val + I0;
end

out.current = I_val;

nalpha = length(alpha);

%compute the energy
% energy = 0;
% for i = 1:nalpha
%     da = (ah(i+1) - ah(i));
%     dx = u(i)^2;
% 
%     energy = energy + dx*da;        
% end

end

