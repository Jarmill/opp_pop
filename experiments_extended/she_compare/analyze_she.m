load('she_gen_08_1_valid.mat')

%one valid solution at d=7

M = root_list{end}.M;

u_base = root_list{end}.u;
alpha_base = root_list{end}.alpha;

% arev = alpha_base(end:-1:1);
% af = [alpha_base; pi-arev; pi+alpha_base; 2*pi-arev];
% uf = [u_base; u_base(end-1:-1:2); -u_base; -u_base(end-1:-1:1)];


N = 1000;
kappa = 1;
% I0 = -0.3;



I0_rec = fzero(@(I0) get_root(u_base, alpha_base, kappa, I0), -0.3);

out_she = pulse_current_voltage_RL(u_base, alpha_base,2, 1000, kappa, I0_rec);

%energy 1.005950354968238
%polish: 1.005872865526797
%order 2: 1.005715663590256
%I0_rec = -0.396822128069477
%SHE is successful, b_(1:13) = 0 and QW is engaged

[na, nb] = pulse_harmonics(19, out_she.u', out_she.alpha(2:end-1)')

function gap = get_root(u_base, alpha_base, kappa, I0)
    out = pulse_current_voltage_RL(u_base, alpha_base,2, 1000, kappa, I0);
    gap = out.I(end) - I0;  
end