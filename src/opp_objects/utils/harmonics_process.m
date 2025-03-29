function [mom_con] = harmonics_process(harm, harm_moments, tol)
%HARMONICS_PROCESS create moment constraints from the harmonics structure

if nargin < 3
    tol = 0;
end

bnd = [harm.bound_cos; harm.bound_sin];

ind_eq = bnd(:, 1) == bnd(:, 2);
ind_ineq = bnd(:, 1) ~= bnd(:, 2);

% if ~isempty(ind_eq)
if tol==0
    con_eq = harm_moments(ind_eq)==bnd(ind_eq, 1);
else
    con_eq = [(harm_moments(ind_eq)<=bnd(ind_eq, 1) + tol);  (harm_moments(ind_eq)>=bnd(ind_eq, 1) - tol)];
end
% else
    % con_eq = [];
% end

if ~isempty(ind_ineq)
    
    ind_geq = ind_ineq & (bnd(:, 2)< Inf);
    ind_leq = ind_ineq & (bnd(:, 1)> -Inf);
    bnd_geq = bnd(ind_geq, 2);
    bnd_leq = bnd(ind_leq, 1);
    con_ineq = [harm_moments(ind_geq) <= bnd_geq; harm_moments(ind_leq)>= bnd_leq];

end

mom_con = [con_eq; con_ineq];

end

