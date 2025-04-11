mset clear
yalmip('clear')

opts = opp_options;
opts.L = [-1, 0, 1];
% opts.L = [-1, -0.5, 0, 0.5, 1];
% opts.L = [-1, 1];
% opts.L = [-2, -1, 0, 1, 2];
opts.harmonics = opp_harmonics();
% opts.partition = 1;
opts.partition = 1;
opts.TIME_INDEP = true;
% opts.start_level = 0;
% opts.start_level = 2;
% opts.start_level = 3;
opts.early_stop = 0;
% opts.start_level = 0;
% opts.start_level = 2;
% opts.null_objective = true;
opts.null_objective = false;
% opts.Symmetry = 0;
% opts.Symmetry = 1;
opts.Symmetry = 2;
% opts.three_phase = "Balanced";
opts.three_phase = "Floating";
opts.k = 4;
% opts.k = 8;
% opts.k = 12;
% opts.k = 16;
% opts.k=20;
% opts.k = 24;
% opts.k = 36;

% opts.common_mode = 1;
% opts.common_mode = 1/3;

% opts.common_mode = 1/3;

% modulation = 0.6;
modulation = 1;
% opts.Z_load = 0;
opts.Z_load = 1.0j;

opts.harmonics.bound_sin = modulation*[1, 1];
% opts.harmonics.bound_cos = [0,  0; 0.5, 0.5];

%k=4 example
% opts.allowed_levels = sparse(1:5, 2+[0, 1, 0, -1, 0], ones(5, 1));

% modulation = 1;
% opts.harmonics.index_cos = [opts.harmonics.index_cos; 2; 3; 4];
% opts.harmonics.bound_cos = [opts.harmonics.bound_cos; 0, 0; 0, 0; -0.1, 0.1];
% opts.harmonics.index_sin= [1; 2; 3; 4];
% opts.harmonics.bound_sin = [modulation, modulation; 0, 0; 0, 0; -0.1, 0.1];



%% test a manager

% k_range = 4:4:20;

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
    bn_lower = sqrt(bound_lower/pi - modulation^2);
    bn_upper = sqrt(bound_upper/pi - modulation^2);
% save('experiments/k_16_full.mat', 'sol', 'opts', 'Mc', 'M', 'pattern_rec', 'ms', 'order')
% save('experiments/k_8_full.mat', 'sol', 'opts', 'Mc', 'M', 'pattern', 'ms', 'order')

% M = MG.mmat();

%% plotting 



    %plot the signal
    N_interp = 900;
th = linspace(0, 2*pi, N_interp);

%function
pu = pattern_rec.u;
pa = pattern_rec.alpha;
x = pulse_func(th, pu, pa);
I0_rec = M.modes{1}{2}.init(1,5);
% I0_rec = M.modes{1}{3}.init(1,5);
%need to perform appropriate scaling
xi = pi*(cumsum(2*x)/(N_interp) + I0_rec);

% [t, y] = ode45(@(t, th) pulse_func(th, pattern.u, pattern.alpha), [0, 2*pi], I0_rec*pi);


cc = linspecer(4);
figure(1)
clf
tiledlayout(3, 1)
nexttile
hold on
plot(th, modulation*sin(th), 'k', 'linewidth', 3);
plot(th, x, 'linewidth', 3, 'color', cc(1, :))

ylabel('$u(\theta)$', 'Interpreter', 'latex', 'FontSize',14);

xlim([0, 2*pi]) 
title(sprintf('M=%0.1f, k=%d, Lower=%0.3f\\%%, Upper=%0.3f\\%%', modulation, opts.k, bn_lower, bn_upper), ...
    'FontSize',16, 'Interpreter', 'latex')

nexttile
hold on
plot(th, -modulation*cos(th), 'k', 'linewidth', 3);
plot(th, xi, 'linewidth', 3, 'color', cc(2, :));
ylabel('$I(\theta)$', 'Interpreter', 'latex', 'FontSize',14);
xlabel('$\theta$', 'Interpreter', 'latex', 'FontSize',14);
xlim([0, 2*pi])


% th_interp = linspace(0, 2*pi, 900);
% xi_interp = interp1(th,xi,th_interp);

xa = x;
xb = circshift(xa, N_interp/3);
xc = circshift(xa, 2*N_interp/3);

xcm = (xa + xb + xc)/3;

nexttile
hold on
plot(th, xcm, 'linewidth', 3, 'color', cc(4, :))
plot([0, 2*pi], [0, 0], ':k')
if opts.common_mode < Inf
    plot([0, 2*pi], [1, 1]*opts.common_mode, 'k')
    plot([0, 2*pi], -[1, 1]*opts.common_mode, 'k')
    ylim([-1, 1]*1.25*opts.common_mode)
end
xlim([0, 2*pi])
ylabel('$v_{cm}(\theta)$', 'Interpreter', 'latex', 'FontSize',14);
xlabel('$\theta$', 'Interpreter', 'latex', 'FontSize',14);



figure(3)
clf
nmax = 21;
[na, nb] = pulse_harmonics(nmax, pu, pa);


di = 0:nmax;
di = di(mod(di, 3) ~= 0);
di = di(2:end); %drop the first harmonic
energy_3 = sum((nb(di+1)./di').^2);
subplot(2, 1,  1)
hold on
stem(0:nmax, na)
title('Cosine Harmonics')
xlabel('n')
ylabel('a_n')
subplot(2, 1, 2)
stem(0:nmax, nb)
title('Sine Harmonics')
xlabel('n')
ylabel('b_n')
end
