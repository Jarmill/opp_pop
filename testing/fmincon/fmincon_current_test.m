load('experiment_N_5_order_2.mat', 'out_std');

osc = out_std{1, 2};
osc.pattern;
L = osc.opts.L;

[dd, N] = size(osc.pattern.levels);
d = dd-1;
u = osc.pattern.u(1:dd)';
% u = [u, u(end)];

%declare the variables
a = sdpvar(d, 1);
af = [0; a; pi/2];
da = diff(af);
I0 = sdpvar(1, 1);

%start with an initial point (from the recovered pattern)
I0_start = osc.pattern.I(1);
alpha_start = osc.pattern.alpha(1:d)';
af0 = [0; alpha_start; pi/2];
assign(a, alpha_start);
assign(I0, I0_start);

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
harm_tol = 1e-5;
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

sdpopts = sdpsettings('solver', 'fmincon', 'usex0',1);
sol = optimize(cons, objective, sdpopts);

if sol.problem ==0
    a_rec = value(a);
    I_rec = value(I);    
    E_rec = value(objective);
    b_rec = value(b);
end
