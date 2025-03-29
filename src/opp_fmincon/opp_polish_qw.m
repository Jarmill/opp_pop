function [out] = opp_polish_qw(osc)
%OPP_POLISH_QW take a pattern solution from the SDP and try to generate a 
%feasible  solution using fmincon (under a fixed switching sequence)
%
%Input: 
%   osc:    the output structure of opp_manager.recover()
%Output:
%   out:    A structure with fields 'warm' and 'cold'
%           'warm' does a warm-start from the reference pattern, 'cold' does not
%           fields:
%               alpha:      The switching angles
%               u:          The levels
%               I:          The current in over the switching sequence
%               b:          Fourier coefficients (sin) for the harmonics
%               objective:  The signal energy in the current
%               tdd:        current tdd sqrt(objective/pi - b1^2)
L = osc.opts.L;

[dd, N] = size(osc.pattern.levels);
d = dd-1;
u = osc.pattern.u(1:dd)';

%declare the variables
a = sdpvar(d, 1);
af = [0; a; pi/2];
da = diff(af);
I0 = sdpvar(1, 1);

%prepare the warm start (but don't assign just yet)
I0_start = osc.pattern.I(1);
alpha_start = osc.pattern.alpha(1:d)';
% alpha_start = [0.01; 0.02; 1];
af0 = [0; alpha_start; pi/2];

modulation = osc.opts.harmonics.bound_sin(1, 1);
harm_tol = 1e-7;

%% form the energy objective
I = cumsum([0; da.*u])+I0;
E = zeros(length(I)-1, 1, 'like', sdpvar);

for k = 1:length(I)-1
    if u(k)==0
        E(k) = I(k)^2 * da(k);
    else
        E(k) = (I(k+1)^3 - I(k)^3)/(3*u(k));
    end
end
E_all = 4*sum(E);
objective = E_all;

%% form the harmonics constraints

%TODO: generalize these for other types of symmetries
sinind = osc.opts.harmonics.index_sin;
b = zeros(length(sinind), 1, 'like', sdpvar);
b0 = zeros(length(sinind), 1); %for sanity checks
for l  = 1:length(sinind)
    lc = sinind(l);
    for m = 1:d        
        % costerm =  (u(m)/lc)*(cos(lc*af(m+1)) - cos(lc*af(m)))
        % b(l) = b(l) + costerm;
        chopblock = (u/lc) .* (-cos(lc*af(2:end)) + cos(lc*af(1:end-1)));
        b(l) = 4*sum(chopblock)/pi;

        chopblock0 = (u/lc) .* (-cos(lc*af0(2:end)) + cos(lc*af0(1:end-1)));
        b0(l) = 4*sum(chopblock0)/pi;
    end
end

con_harm = harmonics_process(osc.opts.harmonics, [zeros(size(osc.opts.harmonics.index_cos)); b], harm_tol);

%% form the ordering constraints
Theta = osc.opts.f0*osc.opts.Ts*2*pi;
Theta_lim = Theta*[ones(d, 1); 0.5];
adiff = (af(2:end) - af(1:end-1) + Theta_lim);
con_order = adiff>=0;

E_start = replace(E_all, [I0; a], [I0_start; alpha_start]);


%% assemble and solve
% cons = [con_harm; con_order];

% cons = [con_harm; con_order; I<=0; diff(I) <= 0];
cons = [con_harm; con_order; I<=0];

sdpopts = sdpsettings('solver', 'fmincon', 'verbose', 0, 'usex0',1);
sdpopts_cold = sdpsettings('solver', 'fmincon', 'verbose', 0);

% sdpopts = sdpsettings('solver', 'fmincon','usex0',1);
% sdpopts_cold = sdpsettings('solver', 'fmincon');

out = struct;
sol_cold = optimize(cons, objective, sdpopts_cold);
if sol_cold.problem ==0
    a_cold = value(a);
    I_cold = value(I);    
    E_cold = value(objective);
    b_cold = value(b);
    tdd_cold = sqrt(E_cold/pi - modulation^2);
    out.cold = struct('alpha', a_cold, 'u', u,'I', I_cold, 'b', b_cold, 'objective', E_cold,  'tdd', tdd_cold);
else
    out.cold = [];
end

%start with an initial point (from the recovered pattern)
assign(a, alpha_start);
assign(I0, I0_start);

sol = optimize(cons, objective, sdpopts);
if sol.problem == 0
    a_warm = value(a);
    I_warm = value(I);    
    E_warm = value(objective);
    b_warm = value(b);
    tdd_warm = sqrt(E_warm/pi - modulation^2);
    out.warm = struct('alpha', a_warm, 'u', u, 'I', I_warm, 'b', b_warm, 'objective', E_warm,  'tdd', tdd_warm);
else
    out.warm= [];
end


end