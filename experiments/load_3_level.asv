mset clear
yalmip('clear')

opts = opp_options;
opts.L = [-1, 0, 1];
% opts.L = [-1, -0.5, 0, 0.5, 1];
% opts.L = [-1, 1];
% opts.L = [-2, -1, 0, 1, 2];
opts.harmonics = opp_harmonics();
% opts.partition = 1;
opts.partition = 2;
% opts.partition = 3;
% opts.partition = 4;
% opts.partition = 8;
% opts.partition = 16;
% opts.TIME_INDEP = true;
opts.TIME_INDEP = true;
opts.early_stop = 0;
opts.Symmetry = 2;
opts.k = 4;
% opts.k = 8;
% opts.k = 12;
% opts.k=16
% opts.k = 24;

% modulation = 0.7;
modulation = 1;
% opts.Z_load = 0;
% opts.Z_load = 1.0j;

opts.harmonics.bound_sin = modulation*[1, 1];

%k=4 example
% opts.allowed_levels = sparse(1:5, 2+[0, 1, 0, -1, 0], ones(5, 1));

% modulation = 1;
% opts.harmonics.index_cos = [opts.harmonics.index_cos; 2; 3; 4];
% opts.harmonics.bound_cos = [opts.harmonics.bound_cos; 0, 0; 0, 0; -0.1, 0.1];
% opts.harmonics.index_sin= [1; 3];
% opts.harmonics.bound_sin = [modulation, modulation; -0.1, 0.1];

%% test a manager
MG = opp_manager(opts);
% order = 4;
% order = 3;
order = 2;
% order = 1;
d = 2*order;

sol = MG.run(order);

disp(sol)


%% diagnose the solution
if sol.status==0
    % m_out = MG.mmat();
    ms = MG.mass_summary();
    pattern_rec = MG.recover_pattern();

    Mc = MG.mmat_corner();
    M = MG.mmat();
    bound_lower = sol.obj_rec;
    if opts.Z_load==1.0j
        bound_upper = pattern_rec.energy_I;
    else
        bound_upper = pattern_rec.energy;
    end
% save('experiments/k_16_full.mat', 'sol', 'opts', 'Mc', 'M', 'pattern_rec', 'ms', 'order')
% save('experiments/k_8_full.mat', 'sol', 'opts', 'Mc', 'M', 'pattern', 'ms', 'order')

% M = MG.mmat();

%% plotting 



    %plot the signal
    N = 1000;
th = linspace(0, 2*pi, N);

%function
pu = pattern_rec.u;
pa = pattern_rec.alpha;
x = pulse_func(th, pu, pa);
I0_rec = M.modes{1}{2}.init(1,5);
%need to perform appropriate scaling
% xi = pi*(cumsum(2*x)/(N) + I0_rec);
xi = pi*(cumsum(2*x)/(N)) + I0_rec;

% [t, y] = ode45(@(t, th) pulse_func(th, pattern.u, pattern.alpha), [0, 2*pi], I0_rec*pi);


cc = linspecer(3);
figure(1)
clf
hold on
plot(th, x, 'linewidth', 3, 'color', cc(1, :))
plot(th, modulation*sin(th), 'k', 'linewidth', 3);
xlim([0, 2*pi]) 
title(sprintf('k=%d, Lower=%0.4f, Upper=%0.4f', opts.k, bound_lower, bound_upper), 'FontSize',16)
nmax = 100;
[na, nb] = pulse_harmonics(nmax, pu, pa);
% energy_L_h = pi*sum(((na(2:end).^2 + nb(2:end).^2)./(1:Nh)'.^2));
% energy_L = sum(iL.^2)/N;

figure(2)
clf
hold on
plot(th, xi, 'linewidth', 3, 'color', cc(2, :));
plot(th, -modulation*cos(th), 'k', 'linewidth', 3);

figure(3)
clf
subplot(2, 1,  1)
hold on
stem(na)
title('Cosine Harmonics')
xlabel('n')
ylabel('a_n')
subplot(2, 1, 2)
stem(nb)
title('Sine Harmonics')
xlabel('n')
ylabel('b_n')
end


































