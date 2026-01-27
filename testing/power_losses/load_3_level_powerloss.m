mset clear
yalmip('clear')

opts = opp_options;
opts.L = [-1, 0, 1];
opts.harmonics = opp_harmonics();
opts.partition = 1;
opts.TIME_INDEP = true;
opts.early_stop = 0;
opts.null_objective = false;
opts.Symmetry = 1;
opts.unipolar = 0; %need to debug this
opts.three_phase = "Ignore";


% opts.power_budget = 100;
% opts.power_budget = Inf;
opts.power_budget = 5000;

%external voltage
opts.k = 4;
% opts.k = 8;
% opts.k = 12;
% opts.k = 16;
% opts.k=20;
% opts.k = 24;
% opts.k = 36;

opts.common_mode = Inf;

modulation = 0.8;

% kappa = 0;
kappa = 0.5;
% kappa = 1.5;
% kappa = 2;

% opts.Z_load = 0;
opts.Z_load = kappa/(2*pi*opts.f0) + 1.0j;
% opts.Z_load = 4.0j;



opts.harmonics.bound_sin = modulation*[1, 1];
% opts.harmonics.bound_cos = [0,  0; 0.5, 0.5];

%% test a manager

opts.dispatch = opp_power_dispatch();

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
    [ms1, ms3] = MG.mass_summary();
    pattern_rec = MG.recover_pattern();

    [Mc1, Mc3] = MG.mmat_corner();
    if opts.clock_split
        mcc = MG.sysclock.mmat_corner();
    else
        mcc = [];
    end
    % M = MG.mmat();
    
    % Q = (eye(3) - ones(3)/3);
    % if opts.three_phase == "Floating"
    %     % Mcd = Mc.diff(4:6, 4:6);
    % 
    % else
    %     Mcd = [];
    % end


    bound_lower = sol.obj_rec;
    if imag(opts.Z_load)>0
        bound_upper = pattern_rec.energy_I;
    else
        bound_upper = pattern_rec.energy;
    end
    bn_lower = sqrt(bound_lower/pi - modulation^2/(1+kappa^2));
    bn_upper = sqrt(bound_upper/pi - modulation^2/(1+kappa^2));
% save('experiments/k_16_full.mat', 'sol', 'opts', 'Mc', 'M', 'pattern_rec', 'ms', 'order')
% save('experiments/k_8_full.mat', 'sol', 'opts', 'Mc', 'M', 'pattern', 'ms', 'order')

% M = MG.mmat();

    RESOLVE = 0;

    if RESOLVE
        opts2 = opts;
        opts2.allowed_levels = pattern_rec.levels;
        MG2 = opp_manager(opts2);
        sol2 = MG2.run(order);
        out2 = MG2.recover(sol2);

        bound_lower2 = sol2.obj_rec;
        bound_upper2 = out2.tdd_upper;
        bound_upper = out2.tdd_upper;        
        out_polish = opp_polish_RL(out2);
        bound_upper = out_polish.warm.tdd;
    end

% pattern_rec = out_polish.warm;

% if ~RESOLVE
%     pu = out.pattern.u;
%     pa = out.pattern.alpha;
%     thi = [0, pa, 2*pi];
%     xi = out.pattern.I;
%     % I0_rec = out.pattern.I(1);
% else
%     pu = out2.pattern.u;
%     pa = out2.pattern.alpha;
%     thi = [0, pa, 2*pi];
%     xi = out2.pattern.I;
%     % I0_rec = out2.pattern.I(1);
% end

%% plotting 



%plot the signal
N_interp = 900;
th = linspace(0, 2*pi, N_interp);

%function
pu = pattern_rec.u;
% pa = pattern_rec.alpha(2:end-1)';
pa = pattern_rec.alpha;
% x = pulse_func(th, pu, pa);
x = pulse_func(th, pu, pa);



I0_rec = pattern_rec.I(1);
% I0_rec = M.modes{1}{3}.init(1,5);
%need to perform appropriate scaling
xi = pi*(cumsum(2*x)/(N_interp)) + I0_rec;

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
plot(th, -modulation/sqrt(1+kappa^2)*(cos(th + atan(kappa))) - ...
    abs(opts.uext)/(kappa^2+1)*(kappa*sin(th + angle(opts.uext)) - cos(th + angle(opts.uext))), 'k', 'linewidth', 3);
if kappa > 0
    plot(pattern_rec.alpha_val, pattern_rec.I_val, 'linewidth', 3, 'color', cc(2, :));
else
    plot(th, xi, 'linewidth', 3, 'color', cc(2, :));
end
scatter([0, pattern_rec.alpha, 2*pi], pattern_rec.I, 200, cc(2, :), 'filled');
ylabel('$I(\theta)$', 'Interpreter', 'latex', 'FontSize',14);
xlabel('$\theta$', 'Interpreter', 'latex', 'FontSize',14);
xlim([0, 2*pi])


% th_interp = linspace(0, 2*pi, 900);
% xi_interp = interp1(th,xi,th_interp);

nexttile
hold on
% xi_query = interp1(thi,xi, th);
res_xi = pattern_rec.I_val + modulation/(1+kappa)*cos(pattern_rec.alpha_val+ atan(kappa));
plot(pattern_rec.alpha_val, res_xi, 'linewidth', 3, 'color', cc(3, :));
plot([0, 2*pi], [0, 0], ':k')
xlim([0, 2*pi])
ylabel('$I(\theta)-I^*(\theta)$', 'Interpreter', 'latex', 'FontSize',14);
xlabel('$\theta$', 'Interpreter', 'latex', 'FontSize',14);

% xa = x;
% xb = circshift(xa, N_interp/3);
% xc = circshift(xa, 2*N_interp/3);
% 
% xcm = (xa + xb + xc)/3;
% 
% nexttile
% hold on
% plot(th, xcm, 'linewidth', 3, 'color', cc(4, :))
% plot([0, 2*pi], [0, 0], ':k')
% if (opts.common_mode < Inf) && (opts.common_mode > 0)
%     plot([0, 2*pi], [1, 1]*opts.common_mode, 'k')
%     plot([0, 2*pi], -[1, 1]*opts.common_mode, 'k')
%     ylim([-1, 1]*1.25*opts.common_mode)
% end
% xlim([0, 2*pi])
% ylabel('$v_{cm}(\theta)$', 'Interpreter', 'latex', 'FontSize',14);
% xlabel('$\theta$', 'Interpreter', 'latex', 'FontSize',14);


% 
% figure(3)
% clf
% nmax = 21;
% [na, nb] = pulse_harmonics(nmax, pu, pa);
% 
% 
% di = 0:nmax;
% di = di(mod(di, 3) ~= 0);
% di = di(2:end); %drop the first harmonic
% energy_3 = sum((nb(di+1)./di').^2);
% subplot(2, 1,  1)
% hold on
% stem(0:nmax, na)
% title('Cosine Harmonics')
% xlabel('n')
% ylabel('a_n')
% subplot(2, 1, 2)
% stem(0:nmax, nb)
% title('Sine Harmonics')
% xlabel('n')
% ylabel('b_n')
% end


%% external voltage support
% figure(5)
% % nexttile
% clf
% hold on
% plot(th, -modulation/sqrt(1+kappa^2)*(cos(th + atan(kappa))), 'k', 'linewidth', 3);
% % plot(th, modulation/(1+kappa^2)*(kappa*sin(th + atan(kappa)) + cos(th + atan(kappa))), 'k', 'linewidth', 3);
% % plot(th, -modulation/(1+kappa)*(cos(th + atan(kappa))) - ...
% %     abs(opts.uext)/(kappa^2+1)*(kappa*sin(th + angle(opts.uext)) - cos(th + angle(opts.uext))), 'k', 'linewidth', 3);
% if kappa > 0
%     plot(pattern_rec.alpha_val, pattern_rec.I_val, 'linewidth', 3, 'color', cc(2, :));
% else
%     plot(th, xi, 'linewidth', 3, 'color', cc(2, :));
% end
% ylabel('$I(\theta)$', 'Interpreter', 'latex', 'FontSize',14);
% xlabel('$\theta$', 'Interpreter', 'latex', 'FontSize',14);
% xlim([0, 2*pi])
end