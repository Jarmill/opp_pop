function [na, nb] = pulse_harmonics(nmax, u, alpha)
%Find the nth harmonic of the pulse function
%
% na = int cos(n th) x(th) d th
% nb = int sin(n th) x(th) d th
%energy = integrate x(th)^2 dth for th in [0, 2pi]

na = zeros(nmax+1, 1);
nb = zeros(nmax+1, 1);


for n = 0:nmax

    [na(n+1), nb(n+1)] = pulse_harmonic(n, u, alpha);
end

end