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

%form the energy
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


I0_start = osc.pattern.I(1);
alpha_start = osc.pattern.alpha(1:d)';

E_start = replace(E_all, [I0; a], [I0_start; alpha_start])

cons = [];

