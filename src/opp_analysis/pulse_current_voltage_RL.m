function [out] = pulse_current_voltage_RL(u, alpha, sym, N, tau, I0)
%Find the energy of the current when applied to an RL circuit.
%
%tau = R/L ratio (default to 1)
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
% I_step = uf.*da;
% I0_val = cumsum([0; I_step]);

I_s = I0*ones(size(ah));

alpha_val = linspace(0, 2*pi, N);
I_val = zeros(1, N);
I_val(1) = I0;

Na = (length(I_s)-1);

%compute the current
for i = 1:Na
    ucurr = uf(i);
    Iprev = I_s(i);
    dt = da(i);
    I_s(i+1) = ucurr*(1-exp(-tau*dt))/tau + Iprev*exp(-tau*dt);

    if N > 0
        a_range = (alpha_val >= ah(i)) & (alpha_val <= ah(i+1));
        dt_range = alpha_val(a_range) - ah(i);
        I_val(a_range) = ucurr*(1-exp(-tau*dt_range))/tau + Iprev*exp(-tau*dt_range);
    end
end

%compute the energy
E_s = zeros(Na, 1);

for i = 1:Na
    ucurr = uf(i);
    Iprev = I_s(i);
    dt = da(i);

    E0 = (ucurr-tau*Iprev)*(3*ucurr + tau*Iprev);
    E_denom = 2*tau^3;
    E_num_1 = 2*ucurr^2*tau*dt;
    E_num_2 = exp(-2*tau*dt)*(ucurr - tau*Iprev) * ...
        (tau*Iprev + ucurr*(4*exp(tau*dt)-1));
    
    E_s(i) = (E_num_1 + E_num_2 - E0)/E_denom;
end

% energy = da'*E_s;
energy = sum(E_s);
out.energy_break = E_s;
out.energy = energy;
% out.I_avg = I_avg;
out.I = I_s;
if N
out.I_val = I_val;
out.alpha_val = alpha_val;
end
out.u = uf;

%compute the energy
% energy = 0;
% for i = 1:nalpha
%     da = (ah(i+1) - ah(i));
%     dx = u(i)^2;
% 
%     energy = energy + dx*da;        
% end

end

