mset clear
yalmip('clear')
mset('yalmip',true);
mset(sdpsettings('solver', 'mosek'));
RESOLVE = 0;

opts = opp_options;
opts.L = [-1, -0.5, 0, 0.5, 1];
% opts.L = [-1, 0, 1];

opts.harmonics = opp_harmonics();

% opts.partition = 1;
opts.partition = 2;
% opts.partition = 4;
opts.TIME_INDEP = true;
opts.early_stop = 0;
opts.null_objective = false;
% opts.Symmetry = 0;
% opts.Symmetry = 1;
opts.Symmetry = 1;
opts.unipolar = 1;
opts.quarter_match = 0;
% opts.three_phase = "Balanced";
% opts.k = 4;
% opts.k = 8;
% opts.k = 12;
% opts.k = 16;
% opts.k=20;
opts.k = 24;
% opts.k= 28;
% opts.k = 32;
% opts.k = 36;

% kappa = 0;
kappa = 0.5;
% kappa = 1;
% kappa = 1.5;
% kappa = 2;

% opts.Z_load = 0;
opts.Z_load = kappa/(2*pi*opts.f0) + 1.0j;

% modulation = 1.1;
% modulation = 0.9;
% modulation = 0.5;
modulation = 0.8;
% modulation = 0.05;
% modulation = 1.05;
% opts.Z_load = 1.0j;
opts.verbose = 0;

%k=4 example
% opts.allowed_levels = sparse(1:5, 2+[0, 1, 0, -1, 0], ones(5, 1));

% modulation = 1;
opts.harmonics.index_sin= [1;  3];
opts.harmonics.bound_sin = [modulation, modulation; -0.01, 0.01];


% opts.harmonics.index_sin= [1];
% opts.harmonics.bound_sin = [modulation, modulation];

% opts.allowed_levels = [0 0 1 0 0;
%                        0 0 0 1 0;
%                        0 0 0 0 1;
%                        0 0 0 1 0;
%                        ones(14, 5)];
% sparse(1:5, 2+[0, 1, 0, -1, 0], ones(5, 1));


%% test a manager


MG = opp_manager(opts);
% order = 1;
% order = 2;
% order = 3;
order = 4;
% order = 5;
% order = 6;
% order = 7;
d = 2*order;

sol = MG.run(order);

disp(sol)


%% diagnose the solution
if sol.status==0
    % m_out = MG.mmat();
    ms = MG.mass_summary();
    pattern_rec = MG.recover_pattern();

    % Mc = MG.mmat_corner();
    % M = MG.mmat();
    out = MG.recover(sol);

    % harm_valid = out.pattern.harm_valid;
    bound_upper = out.tdd_upper;
    %solve again
    if RESOLVE
        opts2 = opts;
        opts2.allowed_levels = out.pattern.levels;
        MG2 = opp_manager(opts2);
        sol2 = MG2.run(order);
        out2 = MG2.recover(sol2);

        bound_lower2 = sol2.tdd_lower;
        bound_upper2 = out2.tdd_upper;
        bound_upper = out2.tdd_upper;        
        % out_polish = opp_polish_qw(out2);
        out_polish = opp_polish_RL(out2);
        bound_upper = out_polish.warm.tdd;
    else
        if kappa == 0
            out_polish = opp_polish_qw(out2);
        else
            out_polish = opp_polish_RL(out);
        end
        % out_polish = out;
        % out_polish.warm = 1;
        
        
    end
    

    harm_valid = ~isempty(out_polish.warm);
    
    if harm_valid
        validstr = ' Valid';
    else
        validstr = ' Invalid';
    end


    summary_str = strcat(sprintf('M=%0.1f, k=%d, TDD Lower=%0.3e, TDD Upper=%0.3e ', ...
    modulation, opts.k, out.tdd_lower, bound_upper), ' ', validstr);



fprintf(strcat(summary_str, '\n'));

%% plotting 
    %plot the signal
    % N = 1000;
% th = linspace(0, 2*pi, N);

%function
if ~RESOLVE
    pu = out.pattern.u;
    pa = out.pattern.alpha;
    % th_switch = [0, pa, 2*pi];
    % xi_switch = out.pattern.I;
    thi = out.pattern.alpha_val;
    xi = out.pattern.I_val;
    pI = out.pattern.I;
    % I0_rec = out.pattern.I(1);
else
    pu = out2.pattern.u;
    pa = out2.pattern.alpha;
    % thi = [0, pa, 2*pi];
    % xi = out2.pattern.I;
    thi = out2.pattern.alpha_val;
    xi = out2.pattern.I_val;
    % I0_rec = out2.pattern.I(1);
end


%% plot

    pu = out_polish.warm.u';
    pa = out_polish.warm.alpha;
    pI = out_polish.warm.I;
    % thi = [0, pa, 2*pi];
    % xi = out_polish.warm.I;
    
    pat = pulse_current_voltage_RL(pu, pa, 0, 1000, kappa,-0.625416052405516);
    thi = pat.alpha_val;
    xi = pat.I_val;
puu = kron(pu, [1; 1]);
paa = [0; kron(pa, [1; 1]); 2*pi];
% x = pulse_func(th, pu, pa);
% I0_rec = M.modes{1}{2}.init(1,5);
% I0_rec = M.modes{1}{3}.init(1,5);

%need to perform appropriate scaling
% xi = 
% [t, y] = ode45(@(t, th) pulse_func(th, pattern.u, pattern.alpha), [0, 2*pi], I0_rec*pi);


cc = linspecer(3);
figure(1)
clf
tiledlayout(3, 1)
nexttile
hold on
plot(thi, modulation*sin(thi), 'k', 'linewidth', 3);
plot(paa, puu, 'linewidth', 3, 'color', cc(1, :))
plot([0, 2*pi], [0, 0], ':k')

ylabel('$u(\theta)$', 'Interpreter', 'latex', 'FontSize',14);

xlim([0, 2*pi]) 
% title(summary_str, ...
    % 'FontSize',16, 'Interpreter', 'latex')

nexttile
hold on
plot(thi, -modulation/sqrt(1+kappa^2)*cos(thi + atan(kappa)), 'k', 'linewidth', 3);
plot(thi, xi, 'linewidth', 3, 'color', cc(2, :));
plot([0, 2*pi], [0, 0], ':k')
scatter([0, pa', 2*pi], pI, 100, cc(2, :), 'filled')
ylabel('$I(\theta)$', 'Interpreter', 'latex', 'FontSize',14);
xlim([0, 2*pi]) 

nexttile
hold on
% xi_query = interp1(thi,xi, th);
res_xi = xi + modulation/sqrt(1+kappa^2)*cos(thi + atan(kappa));
plot(thi, res_xi, 'linewidth', 3, 'color', cc(3, :));
plot([0, 2*pi], [0, 0], ':k')
xlim([0, 2*pi])
ylabel('$I(\theta)-I^*(\theta)$', 'Interpreter', 'latex', 'FontSize',14);
xlabel('$\theta$', 'Interpreter', 'latex', 'FontSize',14);


nmax = 20;
[na, nb] = pulse_harmonics(nmax, pu, pa);
% energy_L_h = pi*sum(((na(2:end).^2 + nb(2:end).^2)./(1:Nh)'.^2));
% energy_L = sum(iL.^2)/N;


% figure(4)
% clf
% hold on
% % plot(th, -modulation*cos(th), 'k', 'linewidth', 3);

% 
% figure(3)
% clf
% if opts.Symmetry == 2
%     tiledlayout(1, 1)
% else
%     tiledlayout(2, 1)
%     nexttile
% % subplot(2, 1,  1)
% hold on
% stem(0:nmax, na)
% title('Cosine Harmonics')
% xlabel('n')
% ylabel('a_n')
% end
% % subplot(2, 1, 2)
% nexttile
% stem(0:nmax,nb)
% title('Sine Harmonics')
% xlabel('n')
% ylabel('b_n')
end
