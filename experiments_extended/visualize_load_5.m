load('harm_noR.mat')
out_L = out;
out_polish_L = out_polish;

load('harm_kappa.mat');
out_kappa = out;
out_polish_kappa = out_polish;



%suboptimality
sub_kappa = out_polish_kappa.warm.objective - out_kappa.energy_lower;
sub_L = out_polish_L.warm.objective - out_L.energy_lower;


tdd_upper = sqrt(out_polish_kappa.warm.objective/pi - 1/(1+kappa^2) * modulation^2);
tdd_lower = sqrt(out_kappa.energy_lower/pi - 1/(1+kappa^2) * modulation^2);
tdd_gap = tdd_upper - tdd_lower
%angles
a_L = out_polish_L.warm.alpha;
a_kappa = out_polish_kappa.warm.alpha;


% obj_kappa_rec = pulse_current_voltage_RL(...
%     out_polish_L.warm.u, out_polish_L.warm.alpha, 0, 1000, kappa, out_polish_kappa.warm.I_q)