%quarter wave
aq = [0.5764188652596557;
     1.3;          
     1.5];
      uq = [0 1 0  1]';

aa_quad = [0; kron(aq, [1; 1]); pi/2];
u_quad = kron(uq, [1; 1]);
aa_allq = [aa_quad; pi - (aa_quad(end:-1:1)); pi + aa_quad; 2*pi - aa_quad(end:-1:1)];
u_allq = [u_quad; u_quad(end:-1:1); -u_quad; -u_quad(end:-1:1)];

%half wave
ah = [0.5;     0.55;     0.85;     1.2;     1.3;     1.5];
uh = [0 1 0 -1 0 1 0]';



aa_half = [0; kron(ah, [1; 1]); pi];
u_half = kron(uh, [1; 1]);
aa_allh = [aa_half ; pi + aa_half];
u_allh = [u_half; -u_half];

%full wave
af = 2*pi*[0.1; 0.15; 0.3; 0.4; 0.6; 0.8; 0.85; 0.9];
uf = [0; 1; 0; 1; 0; -1; 0; 1; 0];
df = diff([0; af; 2*pi]);
cf = sum(df .* uf);
a_full = [0; kron(af, [1; 1]); 2*pi];
u_full = kron(uf, [1; 1]);

%% compute the current
[outf] = pulse_current_voltage(uf, af, 0);
[outh] = pulse_current_voltage(uh, ah, 1);
[outq] = pulse_current_voltage(uq, aq, 2);

%% plot the voltage
figure(1)
clf
tiledlayout(3, 2)

nexttile
c = linspecer(7);
plot([0, 2*pi], [0, 0], 'k')
plot(a_full, u_full, 'LineWidth',3, 'Color', c(6, :))
title('Voltage $u(\theta)$ \\', 'interpreter', 'latex', 'fontsize', 16)
xlim([0, 2*pi])
axis off

nexttile
plot(outf.alpha, outf.current, 'LineWidth',3, 'Color', c(6, :))
title('Current $I(\theta)$ \\', 'interpreter', 'latex', 'fontsize', 16)
xlim([0, 2*pi])
axis off

nexttile
hold on
plot(aa_allh, u_allh, 'LineWidth',3, 'Color', c(5, :))
plot([pi, pi], max(abs(u_allh))*[-1, 1], ':k','LineWidth',2)
xlim([0, 2*pi])
axis off


nexttile
hold on
plot([pi, pi], max(abs(outh.current))*[-1, 1], ':k','LineWidth',2)
plot(outh.alpha, outh.current, 'LineWidth',3, 'Color', c(5, :))
xlim([0, 2*pi])
axis off


nexttile
hold on
plot([pi, pi], [-1, 1], ':k','LineWidth',2)
plot([pi/2, pi/2], [-1, 1], ':k','LineWidth',2)
plot([pi*3/2, pi*3/2], [-1, 1], ':k','LineWidth',2)
plot(aa_allq, u_allq, 'LineWidth',3, 'Color', c(7, :))
xlim([0, 2*pi])
axis off

nexttile
hold on

xlim([0, 2*pi])

plot([pi, pi], max(abs(outq.current))*[-1, 1], ':k','LineWidth',2)
plot([pi/2, pi/2], max(abs(outq.current))*[-1, 1], ':k','LineWidth',2)
plot([pi*3/2, pi*3/2], max(abs(outq.current))*[-1, 1], ':k','LineWidth',2)

plot(outq.alpha, outq.current, 'LineWidth',3, 'Color', c(7, :))
axis off

%% plot only the quarter-wave
figure(2)
clf
tiledlayout(2, 1)

nexttile
hold on
plot([pi, pi], [-1, 1], ':k','LineWidth',2)
plot([pi/2, pi/2], [-1, 1], ':k','LineWidth',2)
plot([pi*3/2, pi*3/2], [-1, 1], ':k','LineWidth',2)
plot(aa_allq, u_allq, 'LineWidth',3, 'Color', c(2, :))
xlim([0, 2*pi])
ylabel('$u(\theta)$', 'Interpreter', 'latex', 'FontSize',14);
% axis off

nexttile
hold on

xlim([0, 2*pi])

plot([pi, pi], [-1, 1], ':k','LineWidth',2)
plot([pi/2, pi/2], [-1, 1], ':k','LineWidth',2)
plot([pi*3/2, pi*3/2], [-1, 1], ':k','LineWidth',2)
% plot([0, 2*pi], [0, 0], '--k','LineWidth',0.5)
plot(outq.alpha, outq.current, 'LineWidth',3, 'Color', c(4, :))
ylabel('$I(\theta)$', 'Interpreter', 'latex', 'FontSize',14);
xlabel('$\theta$', 'Interpreter', 'latex', 'FontSize',14);