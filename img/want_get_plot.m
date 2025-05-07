%quarter wave
aq = [0.423323100643366	0.646289267844426	0.806615717344141	1.55083012024646	]';
      uq = [0 1 0  1 0]';

aa_quad = [0; kron(aq, [1; 1]); pi/2];
u_quad = kron(uq, [1; 1]);
aa_allq = [aa_quad; pi - (aa_quad(end:-1:1)); pi + aa_quad; 2*pi - aa_quad(end:-1:1)];
u_allq = [u_quad; u_quad(end:-1:1); -u_quad; -u_quad(end:-1:1)];

modulation = 1;
th = linspace(0, 2*pi, 1200);
xdes = modulation*sin(th);
ides = -modulation*cos(th);

%% compute the current
[outq] = pulse_current_voltage(uq, aq, 2);


%% plot only the quarter-wave
figure(2)
clf
tiledlayout(2, 2)



nexttile
plot(th, xdes, 'k',  'LineWidth', 2)
box off
xlim([0, 2*pi])
ylabel('$u^*(\theta)$', 'Interpreter', 'latex', 'FontSize',14);
xlabel('$\theta$', 'Interpreter', 'latex', 'FontSize',14);

mi = max(1, max(outq.current));
c = linspecer(5);
nexttile
hold on
% plot([pi, pi], [-1, 1], ':k','LineWidth',2)
% plot([pi/2, pi/2], [-1, 1], ':k','LineWidth',2)
% plot([pi*3/2, pi*3/2], [-1, 1], ':k','LineWidth',2)
plot(aa_allq, u_allq, 'LineWidth',3, 'Color', c(1, :))

xlim([0, 2*pi])
xlabel('$\theta$', 'Interpreter', 'latex', 'FontSize',14);

% axis off
ylabel('$u(\theta)$', 'Interpreter', 'latex', 'FontSize',14);


nexttile
plot(th, ides, 'k',  'LineWidth', 2)

box off
xlim([0, 2*pi])
ylabel('$I^*(\theta)$', 'Interpreter', 'latex', 'FontSize',14);
xlabel('$\theta$', 'Interpreter', 'latex', 'FontSize',14);
ylim([-mi, mi])
nexttile
hold on 

xlim([0, 2*pi])

% plot([pi, pi], [-1, 1], ':k','LineWidth',2)
% plot([pi/2, pi/2], [-1, 1], ':k','LineWidth',2)
% plot([pi*3/2, pi*3/2], [-1, 1], ':k','LineWidth',2)
% plot([0, 2*pi], [0, 0], '--k','LineWidth',0.5)
plot(outq.alpha, outq.current, 'LineWidth',3, 'Color', c(2, :))
ylabel('$I(\theta)$', 'Interpreter', 'latex', 'FontSize',14);
xlabel('$\theta$', 'Interpreter', 'latex', 'FontSize',14);
ylim([-mi, mi])