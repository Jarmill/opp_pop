mset clear
yalmip('clear')

opts = opp_options;
opts.L = [-1, 0, 1];
opts.harmonics = opp_harmonics();
opts.partition = 1;
% opp.Z_load = 1.0j;
opp.Z_load = 0;
opts.TIME_INDEP = true;
opts.start_level = 0;
opts.early_stop = 0;
opts.null_objective = false;
% opts.Symmetry = 0;
opts.Symmetry = 1;
% opts.three_phase = "Balanced";
opts.k = 4;
% opts.k = 8;
% opts.k = 12;
modulation = 0.6;

opts.harmonics.bound_sin = modulation*[1, 1];

%% test a manager
MG = opp_manager(opts);
% order = 4;
order = 2;
% order = 1;
d = 2*order;

%k=4, full-wave symmetry

%the manager can be created
%do the constraints work though?

sol = MG.run(order);

disp(sol)