%figure out periodicity and zero-average-value

%% load characteristics
% R = 1.5;
% R = 10;
% R = 0.1;
% R = 1e-1;
% R = 1.5;
% R = 1;
R = 0.4;
L = 1;
kappa = R/L;

N = 1000;

I0 = -0.473;
%% make the signal
aq = [0.423323100643366	0.646289267844426	0.806615717344141	1.55083012024646	]';
uq = [0 1 0  1 0]';

param = struct('aq', aq, 'uq', uq, 'kappa', kappa);



y0 =-0.3;

[y, cval, info] = fmincon(@(I) objective_Imatch(I, param), y0, [1; -1], [1; 1]);
[outq] = pulse_current_voltage_RL(param.uq, param.aq, 2, 1000, param.kappa, y);
[outL] = pulse_current_voltage(param.uq, param.aq, 2);

outq.uu = kron(outq.u, [1; 1]);    
outq.aa = [0; kron(outq.alpha(2:end-1), [1; 1]); 2*pi];

disp(outq.I(end)-outq.I(1))

%% plot the voltage and the current
figure(1)
clf
tiledlayout(3, 1)

c = linspecer(7);


nexttile
hold on
% plot(reshape((0:(T-1))*2*pi + outq.aa, [], 1), reshape(repmat(outq.uu, [1, T]), [], 1),...
    plot(outq.aa, outq.uu, 'LineWidth', 3, 'Color', c(2, :));
% plot(outq.alpha_val, outq.I_val, 'LineWidth', 3,  'Color', c(2, :));
% scatter(outq.alpha, outq.I, 100, 'k', 'filled')
plot([pi, pi], [-1, 1], ':k','LineWidth',2)
plot([pi/2, pi/2], [-1, 1], ':k','LineWidth',2)
plot([pi*3/2, pi*3/2], [-1, 1], ':k','LineWidth',2)
xlim([0, 2*pi]);
title('$L_{{load}}=0$ (QW-Symmetric) \\', 'interpreter', 'latex', 'fontsize', 16)
axis off


nexttile
hold on
plot([pi, pi], [-1, 1], ':k','LineWidth',2)
plot([pi/2, pi/2], [-1, 1], ':k','LineWidth',2)
plot([pi*3/2, pi*3/2], [-1, 1], ':k','LineWidth',2)
% plot([0, 2*pi], [0, 0], '--k','LineWidth',0.5)
plot(outL.alpha, outL.current, 'LineWidth',3, 'Color', c(1, :))
title('$R_{{load}}=0$ (QW-Antisymmetric) \\', 'interpreter', 'latex', 'fontsize', 16)
xlim([0, 2*pi]);
axis off


nexttile

xlim([0, 2*pi])
hold on
plot([pi, pi], [-1, 1], ':k','LineWidth',2)
plot([pi/2, pi/2], [-1, 1], ':k','LineWidth',2)
plot([pi*3/2, pi*3/2], [-1, 1], ':k','LineWidth',2)
% plot([0, 2*pi], [0, 0], '--k','LineWidth',0.5)
plot(outq.alpha_val, outq.I_val, 'LineWidth',3, 'Color', c(3, :))
title('$R_{{load}}/L_{{load}} = 0.4$ (HW-Symmetric, not QW) \\', 'interpreter', 'latex', 'fontsize', 16)
xlim([0, 2*pi]);
axis off

function c = objective_Imatch(I0, param)

    [outq] = pulse_current_voltage_RL(param.uq, param.aq, 2, 1000, param.kappa, I0);

    c = abs(outq.I_val(1) - outq.I_val(end));
end