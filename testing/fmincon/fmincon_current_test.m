% load('experiment_N_5_order_2.mat', 'out_std');

% osc = out_std{1, 2};

mi = 14;
ki = 8;
ordi = 3;

m = m_list(mi);

osc = out_std{mi, ki, ordi};
oscr = out_resolve{mi, ki, ordi};

alpha_orig = osc.pattern.alpha(1:k_list(ki)/4)';

tdd_lower_orig = osc.tdd_lower;
tdd_upper_orig = osc.tdd_upper;
out = opp_polish_qw(osc);

tdd_lower_r = oscr.tdd_lower;
tdd_upper_r = oscr.tdd_upper;
alpha_r = oscr.pattern.alpha(1:k_list(ki)/4)';
out_r = opp_polish_qw(oscr);

