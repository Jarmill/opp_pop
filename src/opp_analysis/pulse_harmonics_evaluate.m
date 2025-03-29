function [harm_pattern, valid] = pulse_harmonics_evaluate(pattern, harm, tol)
%Find the nth harmonic of the pulse function
%
% na = int cos(n th) x(th) d th
% nb = int sin(n th) x(th) d th
%energy = integrate x(th)^2 dth for th in [0, 2pi]
if nargin < 3
    tol = 1e-4;
end

nmax = max(max(harm.index_cos), max(harm.index_sin));
            [na, nb] = pulse_harmonics(nmax, pattern.u, pattern.alpha);

harm_pattern = [na(harm.index_cos+1); 
    nb(harm.index_sin+1)];

valid_check = harmonics_process(harm, harm_pattern, tol);
valid = all(valid_check);

end