
mset('clear')
opts = opp_options;
opts.L = [-1, 0, 1];
% opts.L = [-1, -0.5, 0, 0.5, 1];
% opts.L = [-1, 1];
% opts.L = [-2, -1, 0, 1, 2];
opts.harmonics = opp_harmonics();
% opts.partition = 1;
opts.partition = 2;
opts.TIME_INDEP = true;
opts.early_stop = 0;
opts.null_objective = false;
opts.Symmetry = 0;
opts.Symmetry = 1;
opts.Z_load = 1.0j;
% opts.Symmetry = 2;

%3-level inverter, unipolar
% opts.unipolar = 1;
% opts.common_mode = 0;   %5 configurations
% opts.common_mode = 1/3; %13 configurations
% opts.common_mode = 2/3; %17 configurations
% opts.common_mode = Inf; %18 configurations

%3-level inverter, multipolar
opts.unipolar = 0;
opts.common_mode = 0;     %7 configurations
% opts.common_mode = 1/3; %19 configurations
% opts.common_mode = 2/3; %25 configurations
% opts.common_mode = Inf; %27 configurations

% opts.three_phase = "Floating";
opts.three_phase = "Balanced";
opts.k = 4;

Diff = opp_diff_dyn(opts);
% [levels, G] = Diff.mode_switch_graph();