load('experiment_N_5_order_2.mat', 'out_std');

osc = out_std{1, 2};
out = opp_polish_qw(osc)