%figure out periodicity and zero-average-value

%% load characteristics
% R = 1.5;
% R = 10;
% R = 0.1;
% R = 1e-1;
R = 1;
% R = 0;
% R = 1e-4;
% R = 1;
L = 1;
kappa = R/L;

%external voltage signal
Aext = -0.5;
% Aext = -1;
% Aext = 0;
% Aext = -0.25;
% phiext = 0;
% phiext = pi;
% phiext = pi/3;
phiext = pi/6;

N = 1000;

I0 = -0.473;
%% make the signal
% aq = [0.423323100643366	0.646289267844426	0.806615717344141	1.55083012024646	]';
% uq = [0 1 0  1 0]';

% aq = [0.423323100643366	0.646289267844426	0.806615717344141	1.55083012024646	]';
aq = [0.720186460835395,1.447102394507780]';
% uq = [0 1 0  1 0]';
uq = [0 1 0]';

param = struct('aq', aq, 'uq', uq, 'kappa', kappa);



% [outq] = pulse_current_voltage_RL(uq, aq, 2, N, kappa, I0);
% outq.uu = kron(outq.u, [1; 1]);    
% outq.aa = [0; kron(outq.alpha(2:end-1), [1; 1]); 2*pi];

% y0 = -0.472;
y0 =-0.3;

NA = 10000;
[y, cval, info] = fmincon(@(I) objective_Imatch(I, param), y0, [1; -1], [1; 1]);
[outq] = pulse_current_voltage_RL_ext(param.uq, param.aq, 2, NA, param.kappa, y, Aext, phiext);

e_trapz = trapz(outq.alpha_val, outq.I_val.^2)

ediff  = outq.energy - e_trapz;


outq.uu = kron(outq.u, [1; 1]);    
outq.aa = [0; kron(outq.alpha(2:end-1), [1; 1]); 2*pi];

disp(outq.I(end)-outq.I(1))

%% plot the voltage and the current
figure(1)
c = linspecer(2);
clf
hold on
% plot(reshape((0:(T-1))*2*pi + outq.aa, [], 1), reshape(repmat(outq.uu, [1, T]), [], 1),...
    plot(outq.aa, outq.uu, 'LineWidth', 3, 'Color', c(1, :));
plot(outq.alpha_val, outq.I_val, 'LineWidth', 3,  'Color', c(2, :));
scatter(outq.alpha, outq.I, 100, 'k', 'filled')
xlim([0, 2*pi]);

% figure(2)
% plot(I0_list, y_list(:, end), 'color', 'k')
% xlabel('$I(0)$', 'Interpreter', 'latex')
% ylabel('$I(2\pi)$', 'Interpreter', 'latex')

function c = objective_Imatch(I0, param)
    [outq] = pulse_current_voltage_RL(param.uq, param.aq, 2, 1000, param.kappa, I0);
    % [outq] = pulse_current_voltage_RL_ext(param.uq, param.aq, 2, NA, param.kappa, y, Aext, phiext);

    c = abs(outq.I_val(1) - outq.I_val(end));
end