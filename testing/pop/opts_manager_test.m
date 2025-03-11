mset clear

opts = opp_options;
opts.L = [-1, 0, 1];
% opts.L = [-2, -1, 0, 1, 2];
opts.harmonics = opp_harmonics();
opts.partition = 2;
opts.TIME_INDEP = false;
opts.start_level = 0;
opts.early_stop = 1;
opts.null_objective = true;
opts.Symmetry = 2;
% opts.Symmetry = 1;
% opts.three_phase = "Balanced";
opts.k = 4;
% opts.k=12;


% R = 0;
% L = 5e-3;
% X = 1.0j*(2*pi)*opts.f0*L;
% opts.Z_load = R+X;



% modulation = 1;
% opts.harmonics.index_cos = [opts.harmonics.index_cos; 2; 3; 4];
% opts.harmonics.bound_cos = [opts.harmonics.bound_cos; 0, 0; 0, 0; -0.1, 0.1];
% opts.harmonics.index_sin= [1; 2; 3; 4];
% opts.harmonics.bound_sin = [modulation, modulation;
%     0, 0; 0, 0; -0.1, 0.1];


%% test a manager
MG = opp_manager(opts);
order = 1;
d = 2*order;

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