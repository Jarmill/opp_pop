function [energy] = pulse_energy(u, alpha)
%Find the energy of the pulse function
%
%energy = integrate x(th)^2 dth for th in [0, 2pi]

a0 = [0, alpha, 2*pi];
nalpha = length(alpha);

energy = 0;
for i = 1:nalpha
    da = (a0(i+1) - a0(i));
    dx = u(i)^2;

    energy = energy + dx*da;        
end

end

