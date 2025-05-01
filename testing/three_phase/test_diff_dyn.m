
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
% opts.Symmetry = 0;
opts.Symmetry = 1;
% opts.Symmetry = 2;
opts.unipolar = 1;
opts.common_mode = 1/3;
opts.three_phase = "Floating";
opts.k = 4;

Diff = opp_diff_dyn(opts);
[levels, G] = Diff.mode_switch_graph();