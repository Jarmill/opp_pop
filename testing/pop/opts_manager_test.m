mset clear

opts = opp_options;
opts.L = [-1, 0, 1];
% opts.L = [-1, 1];
% opts.L = [-2, -1, 0, 1, 2];
opts.harmonics = opp_harmonics();
opts.partition = 2;
opp.Z_load = 1.0j;
% opts.partition = 8;
% opts.partition = 16;
opts.TIME_INDEP = true;
opts.start_level = 0;
opts.early_stop = 0;
% opts.null_objective = true;
opts.null_objective = false;
opts.Symmetry = 0;
% opts.Symmetry = 1;
% opts.three_phase = "Balanced";
opts.k = 4;
% opts.k = 8;
% opts.k = 12;
% opts.k = 16;
% opts.k=2;
% opts.Ts = (pi/8)/opts.f0;
% opts.Ts = 1e-3;
% opts.k=12;


% R = 0;
% L = 5e-3;
% X = 1.0j*(2*pi)*opts.f0*L;
% opts.Z_load = R+X;



% modulation = 1;
% opts.harmonics.index_cos = [opts.harmonics.index_cos; 2; 3; 4];
% opts.harmonics.bound_cos = [opts.harmonics.bound_cos; 0, 0; 0, 0; -0.1, 0.1];
% opts.harmonics.index_sin= [1; 2; 3; 4];
% opts.harmonics.bound_sin = [modulation, modulation; 0, 0; 0, 0; -0.1, 0.1];


%% test a manager
MG = opp_manager(opts);
order = 4;
d = 2*order;

%k=4, full-wave symmetry

%the manager can be created
%do the constraints work though?

sol = MG.run(order);

disp(sol)

%% step through the constraints

% %support constraints can be generated
% sc = MG.supp_con();
% 
% lc = MG.con_lebesgue_circ(d); %works
% pc = MG.con_prob_dist(); %works
% rc = MG.con_return(d); 
% 
% %%harmonics
% 
% hc= MG.con_harmonics();
% 
% % fc = MG.con_flow(d);
% 
% % mom_con = [fc; hc; rc; pc; lc];
% % supp_con = sc;

% vi = MG.vars.x(1:2);
% vi = MG.vars.x;
% mon_3 = MG.three_phase_rotate(vi, d);
% bc = MG.con_balance(d);
% om = MG.opp_objective();

%% diagnose the solution
if sol.status==0
    % m_out = MG.mmat();
    ms = MG.mass_summary();
    pattern = MG.recover_pattern();

    bound_lower = sol.obj_rec;
    bound_upper = pattern.energy;


    %plot the signal
    N = 1000;
th = linspace(0, 2*pi, N);

%function
x = pulse_func(th, pattern.u, pattern.alpha);


%% plotting 
figure(1)
clf
hold on
plot(th, x)
plot(th, sin(th), 'k');
xlim([0, 2*pi]) 
title(sprintf('k=%d, Lower=%0.4f, Upper=%0.4f', opts.k, bound_lower, bound_upper), 'FontSize',16)
iL = cumsum(x)/N;
nmax = 100;
[na, nb] = pulse_harmonics(nmax, pattern.u, pattern.alpha);
% energy_L_h = pi*sum(((na(2:end).^2 + nb(2:end).^2)./(1:Nh)'.^2));
% energy_L = sum(iL.^2)/N;

figure(2)
plot(th, iL);

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


































