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
opts.k = 4;
% opts.k = 8;
% opts.k = 12;
% opts.k = 16;
% opts.k=20;
% opts.k = 24;
% opts.k = 36;
opts.quarter_match = true;

% opts.common_mode = 1;
% opts.common_mode = 1/3;
% opts.common_mode = 0;
opts.common_mode = Inf;

% opts.common_mode = 1/3;

% modulation = 0.6;
% modulation = 0.25;
modulation = 0.8;

%R/L ratio
kappa = 3;

% opts.Z_load = 0;
opts.Z_load = kappa + 1.0j*(2*pi*opts.f0);


% opts.Z_load = 4.0j;




opts.harmonics.bound_sin = modulation*[1, 1];

%% test a manager


MG = opp_manager(opts);
% order = 4;
% order = 3;
% order = 2;
order = 1;
d = 2*order;

con_match = MG.sys1.con_quarter_match(d);
con_liou = MG.sys1.con_flow(d);


% sol = MG.run(order);
% 
% disp(sol)

