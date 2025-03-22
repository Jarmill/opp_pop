mset clear

opts = opp_options;
opts.L = [-1, 0, 1];
% opts.L = [-1, 1];
% opts.L = [-2, -1, 0, 1, 2];
opts.harmonics = opp_harmonics();
opts.partition = 4;
% opts.partition = 8;
% opts.partition = 16;
% opts.TIME_INDEP = true;
opts.TIME_INDEP = true;
% opts.start_level = 0;
opts.start_level = 2;
opts.early_stop = 0;
% opts.null_objective = true;
opts.null_objective = false;
opts.Symmetry = 0;
% opts.Symmetry = 1;
% opts.three_phase = "Balanced";
opts.k = 4;

% modulation = 0.5;
modulation = 1;
% opts.Z_load = 0;
opts.Z_load = 1.0j;

opts.harmonics.bound_sin = modulation*[1, 1];

%k=4 example
opts.allowed_levels = sparse(1:5, 2+[0, 1, 0, -1, 0], ones(5, 1));

% modulation = 1;
% opts.harmonics.index_cos = [opts.harmonics.index_cos; 2; 3; 4];
% opts.harmonics.bound_cos = [opts.harmonics.bound_cos; 0, 0; 0, 0; -0.1, 0.1];
% opts.harmonics.index_sin= [1; 2; 3; 4];
% opts.harmonics.bound_sin = [modulation, modulation; 0, 0; 0, 0; -0.1, 0.1];


%% test a manager
MG = opp_manager(opts);
% order = 4;
order = 2;
% order = 1;
d = 2*order;

sol = MG.run(order);

disp(sol)


%% diagnose the solution
if sol.status==0
    % m_out = MG.mmat();
    ms = MG.mass_summary();
    pattern = MG.recover_pattern();

    Mc = MG.mmat_corner();
    M = MG.mmat();
    bound_lower = sol.obj_rec;
    if opts.Z_load==1.0j
        bound_upper = pattern.energy_I;
    else
        bound_upper = pattern.energy;
    end



% M = MG.mmat();

%% plotting 



    %plot the signal
    N = 1000;
th = linspace(0, 2*pi, N);

%function
x = pulse_func(th, pattern.u, pattern.alpha);
I0_rec = M.modes{1}{2}.init(1,5);
xi = cumsum(x)/(N) + I0_rec;



cc = linspecer(3);
figure(1)
clf
hold on
plot(th, x, 'linewidth', 3, 'color', cc(1, :))
plot(th, modulation*sin(th), 'k', 'linewidth', 3);
xlim([0, 2*pi]) 
title(sprintf('k=%d, Lower=%0.4f, Upper=%0.4f', opts.k, bound_lower, bound_upper), 'FontSize',16)
nmax = 100;
[na, nb] = pulse_harmonics(nmax, pattern.u, pattern.alpha);
% energy_L_h = pi*sum(((na(2:end).^2 + nb(2:end).^2)./(1:Nh)'.^2));
% energy_L = sum(iL.^2)/N;

figure(2)
clf
hold on
plot(th, xi, 'linewidth', 3, 'color', cc(2, :));
% plot(th, -modulation*cos(th), 'k', 'linewidth', 3);

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


































