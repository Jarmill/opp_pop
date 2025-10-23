function [out] = opp_polish_RL(osc)
%OPP_POLISH_RL take a pattern solution from the SDP and try to generate a 
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

Sym = double(osc.opts.Symmetry);


kappa = real(osc.opts.Z_load)/imag(osc.opts.Z_load) * (2 * pi * osc.opts.f0); 

[dd, N] = size(osc.pattern.occ);
% d = dd-1;

if Sym == 1 && osc.opts.quarter_match
    Sym = 2;
    dd = 1 + (dd-1)/2;
end

u = osc.pattern.u(1:dd)';

%declare the variables
a = sdpvar(dd-1, 1);
% af = [0; a; pi*2^(-Sym)];
% da = diff(af);
I0 = sdpvar(1, 1);

%prepare the warm start (but don't assign just yet)
I0_start = osc.pattern.I(1);
% alpha_start = osc.pattern.alpha(1:d)';
% % alpha_start = [0.01; 0.02; 1];
% af0 = [0; alpha_start; pi/2];

modulation = osc.opts.harmonics.bound_sin(1, 1);
harm_tol = 1e-7;

%% form the energy objective


    % function energy = curr_objective(alpha)
    %     out_r =pulse_current_voltage_RL(u, alpha, Sym, 0, kappa, I0);
    %     energy = out_r.energy;
    % end

if Sym == 0
    %full-wave symmetry
    uf = u;
    af = a;
elseif Sym==1
    %half-wave symmetry
    uf = [u; -u(2:end)];
    af = [a; a+pi];
else
    %quarter-wave symmetry
    urev = u(end-1:-1:1);
    uf = [u; urev; -u(2:end); -urev];
    af= [a; pi - a(end:-1:1); pi + a; 2*pi - a(end:-1:1)];
end

ah = [0; af; 2*pi];
da = diff(ah);

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
    I_s(i+1) = ucurr*(1-exp(-kappa*dt))/kappa + Iprev*exp(-kappa*dt);
end

%compute the energy
% E_s = zeros(Na, 1)*a(1);
energy = 0;
for i = 1:Na
    ucurr = uf(i);
    Iprev = I_s(i);
    dt = da(i);

    E0 = (ucurr-kappa*Iprev)*(3*ucurr + kappa*Iprev);
    E_denom = 2*kappa^3;
    E_num_1 = 2*ucurr^2*kappa*dt;
    E_num_2 = exp(-2*kappa*dt)*(ucurr - kappa*Iprev) * ...
        (kappa*Iprev + ucurr*(4*exp(kappa*dt)-1));
    
    energy = energy +  (E_num_1 + E_num_2 - E0)/E_denom;
end

objective = energy;

%% form the harmonics constraints
harm = osc.opts.harmonics;

if Sym > 0
    elim_sin = mod(harm.index_sin, 2)==0;
    harm.bound_sin(elim_sin, :) = [];
    harm.index_sin(elim_sin, :) = [];
    if Sym == 1
        elim_cos = mod(harm.index_cos, 2)==0;
        harm.bound_cos(elim_cos, :) = [];
        harm.index_cos(elim_cos, :) = [];
    else
        harm.bound_cos = [];
        harm.index_cos= [];
    end
end

if Sym==2
    nmax =  max(harm.index_sin);
else
    nmax = max(max(harm.index_cos), max(harm.index_sin));
end


na = zeros(nmax+1, 1, 'like', sdpvar);
nb = zeros(nmax+1, 1, 'like', sdpvar);


for n = 0:nmax
    [nacurr, nbcurr] = pulse_harmonic(n, osc.pattern.u, af');
    if ~isnumeric(nacurr)
     na(n+1) = nacurr;
    end
    if ~isnumeric(nbcurr)
     nb(n+1) = nbcurr;
    end
end

% [na, nb] = pulse_harmonics(nmax, osc.pattern.u, af');

harm_pattern = [na(harm.index_cos+1); 
    nb(harm.index_sin+1)];

harm_tol = 1e-4;
con_harm = harmonics_process(harm, harm_pattern, harm_tol);
% [harm_pattern, con_harm] = pulse_harmonics_evaluate(struct('u', osc.pattern.u, 'alpha', alpha), harm, tol);

%% form the ordering constraints
Theta = osc.opts.f0*osc.opts.Ts*2*pi;
acirc = [af; 2*pi + af(1)];
dacirc = diff(acirc);
con_order = dacirc >= Theta;
% con_order = [da >= 0 ];
% if Sym == 2
%     Theta_lim = Theta*[0.5; ones(length(da)-1, 1); 0.5];
% else
%     Theta_lim = Theta*[0; ones(length(da)-1, 1)];
% end
% Theta_lim = 0;
% adiff = (da - Theta_lim);
% con_order = adiff>=0;

% E_start = replace(E_all, [I0; a], [I0_start; alpha_start]);


%% assemble and solve
% cons = [con_harm; con_order];

% cons = [con_harm; con_order; I<=0; diff(I) <= 0];
% cons = [con_harm; con_order; I_s(1)<=0];
% cons = [con_order; con_harm];
cons = [con_harm; con_order; I0 <= 0; I0 == I_s(end)];
% cons = [con_order; I0 <= 0];
feval_max = 1000;

solver = 'fmincon';
% solver = 'ipopt';

sdpopts = sdpsettings('solver', solver, 'verbose', 1, 'usex0',1, ...
    'fmincon.maxfunevals', 1.5*feval_max, 'fmincon.MaxIter', feval_max, 'fmincon.EnableFeasibilityMode', true);
sdpopts_cold = sdpsettings('solver', solver, 'verbose', 1, ...
    'fmincon.maxfunevals', 1.5*feval_max, 'fmincon.MaxIter', feval_max, 'fmincon.EnableFeasibilityMode', true);

% sdpopts = sdpsettings('solver', 'fmincon','usex0',1);
% sdpopts_cold = sdpsettings('solver', 'fmincon');

out = struct;
sol_cold = optimize(cons, objective, sdpopts_cold);
if sol_cold.problem ==0
    a_cold = value(af);
    I_cold = value(I_s);    
    E_cold = value(objective);
    % b_cold = value(b);
    % a_cold_full = [a_cold; pi-a_cold(end:-1:1); a_cold+pi; 2*pi - a_cold(end:-1:1)];    
    % I_cold_full = [I_cold(1:end-1); -I_cold(end-1:-1:2); -I_cold(2:end-1); I_cold(end-1:-1:1)];
    tdd_cold = sqrt(E_cold/pi - modulation^2);
    out.cold = struct('alpha', a_cold,  'u', osc.pattern.u,'I', I_cold, 'objective', E_cold,  'tdd', tdd_cold,...
        'na', value(na), 'nb', value(nb), 'solvertime', sol_cold.solvertime, 'yalmiptime', sol_cold.yalmiptime);
else
    tdd_cold = Inf;
    out.cold = [];
end

%start with an initial point (from the recovered pattern)
assign(a, osc.pattern.alpha(1:dd-1)');
assign(I0, osc.pattern.I(1));

sol_warm = optimize(cons, objective, sdpopts);
if sol_warm.problem == 0
    a_warm = value(a);
    I_warm = value(I0);    
    E_warm = value(objective);
    % b_warm = value(b);
    a_warm_full = value(af);    
    I_warm_full = value(I_s);
    tdd_warm = sqrt(E_warm/pi - modulation^2);
    out.warm = struct('alpha', a_warm_full, 'alpha_q', a_warm, 'u', osc.pattern.u, ...
        'I', I_warm_full, 'I_q', I_warm, 'objective', E_warm,  'tdd', tdd_warm,...
        'solvertime', sol_warm.solvertime, 'yalmiptime', sol_warm.yalmiptime);

else
    tdd_warm = Inf;
    out.warm= [];
end

% out.tdd = min(tdd_warm, tdd_cold);


end