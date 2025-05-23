mset clear
yalmip('clear')

opts = opp_options;
opts.L = [-1, 0, 1];
% opts.L = [-1, -0.5, 0, 0.5, 1];
% opts.L = [-1, 1];
% opts.L = [-2, -1, 0, 1, 2];
opts.harmonics = opp_harmonics();
opts.partition = 1;
opts.TIME_INDEP = true;
opts.early_stop = 0;
opts.null_objective = false;
opts.Symmetry = 1;
% opts.Symmetry = 2;
opts.unipolar = 1; %need to debug this
% opts.three_phase = "Balanced";
% opts.three_phase = "Floating";
opts.three_phase = "Ignore";
% opts.k = 4;
% opts.k = 8;
% opts.quarter_match = true;
% opts.k = 12;
% opts.k = 16;
% opts.k=20;
% opts.k = 24;
% opts.k = 36;
opts.k = 40;

% opts.common_mode = 1;
% opts.common_mode = 1/3;
% opts.common_mode = 0;
opts.common_mode = Inf;

% opts.common_mode = 1/3;

% modulation = 0.6;
% modulation = 0.25;
modulation = 0.8;
% modulation = 1;

% kappa = 0;
kappa = 0.5;
% kappa = 2;
% kappa = 1;
% kappa = 1.5;
% kappa = 2;

Ascale = sqrt((1)/(1+kappa^2));

% opts.Z_load = 0;
opts.Z_load = kappa/(2*pi*opts.f0) + 1.0j;
% opts.Z_load = 4.0j;




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


MG = opp_manager(opts);
% order = 4;\
order = 3;
% order = 2;
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
    bn_lower = sqrt(bound_lower/pi - modulation^2*Ascale^2);
    bn_upper = sqrt(bound_upper/pi - modulation^2*Ascale^2);
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



I0_rec = pattern_rec.I(1);
% I0_rec = M.modes{1}{3}.init(1,5);
%need to perform appropriate scaling
xi = pi*(cumsum(2*x)/(N_interp)) + I0_rec;

% [t, y] = ode45(@(t, th) pulse_func(th, pattern.u, pattern.alpha), [0, 2*pi], I0_rec*pi);


cc = linspecer(4);
figure(2)
clf
tiledlayout(2, 1)
nexttile
hold on
plot(th, modulation*sin(th), 'k', 'linewidth', 3);
plot(th, x, 'linewidth', 3, 'color', cc(1, :))

ylabel('$u(\theta)$', 'Interpreter', 'latex', 'FontSize',14);

xlim([0, 2*pi]) 
title(sprintf('M=%0.1f, k=%d, Lower=%0.3f\\%%, Upper=%0.3f\\%%', modulation, opts.k, bn_lower, bn_upper), ...
    'FontSize',16, 'Interpreter', 'latex')

%find the current

%there might be an incorrect 2pi-scaling somewhere
a_jump = [];
I_jump = [];
for i = 1:length(Mc1.jump)
    for n = 1:length(opts.L)-1
        juc = Mc1.jump{i}.up{n};
        jdc = Mc1.jump{i}.down{n};
        %jump up
        if ~isempty(juc) && abs(juc(1, 1) - 1) < 1e-6
            alpha_curr = atan2(juc(3, 1), juc(2, 1));
            a_jump = [a_jump; alpha_curr];
            I_jump = [I_jump; juc(5, 1)];
        end
        %jump down
        if ~isempty(jdc) && abs(jdc(1, 1) - 1) < 1e-6
            alpha_curr = atan2(jdc(3, 1), jdc(2, 1));
            a_jump = [a_jump; alpha_curr];
            I_jump = [I_jump; jdc(5, 1)];
        end
    end
end

if opts.Symmetry == 1
    a_jump = [a_jump; a_jump + pi];
    I_jump = [I_jump; -I_jump];
end



nexttile
hold on
plot(th, -modulation*Ascale*cos(th + atan(kappa)), 'k', 'linewidth', 3);
if kappa > 0
    plot(pattern_rec.alpha_val, pattern_rec.I_val, 'linewidth', 3, 'color', cc(2, :));
    scatter([0, pattern_rec.alpha, 2*pi], pattern_rec.I, 100, 'k', 'filled')    
else
    plot(th, xi, 'linewidth', 3, 'color', cc(2, :));
end
scatter(a_jump, pi*I_jump, 100, 'g', 'filled')


%quarter-match
I0 = Mc1.modes{1}{2}.init(5, 1);
Iterm = Mc1.modes{end}{2}.term(5, 1);
scatter([0, pi, 2*pi], pi*[I0; Iterm; I0], 100, 'b', 'filled')


ylabel('$I(\theta)$', 'Interpreter', 'latex', 'FontSize',14);
xlabel('$\theta$', 'Interpreter', 'latex', 'FontSize',14);
xlim([0, 2*pi])



%% plot the harmonics
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
