aq = [0.423323100643366	0.646289267844426	0.806615717344141	1.55083012024646	]';
      uq = [0 1 0  1 0]';
c = linspecer(5);

% evaluate
[outq] = pulse_current_voltage(uq, aq, 2);
outq.uu = kron(outq.u, [1; 1]);    
outq.aa = [0; kron(outq.alpha(2:end-1), [1; 1]); 2*pi];


%sample
N_interp = 900;
th = linspace(0, 2*pi, N_interp);

pu = outq.u;
pa = outq.alpha;
x = pulse_func(th, pu, pa(2:end-1)');

xa = x;
xb = circshift(xa, N_interp/3);
xc = circshift(xa, 2*N_interp/3);
xcm = (xa + xb + xc)/3;

%current
I0_rec = outq.current(1);
xia = pi*(cumsum(2*x)/(N_interp)) + I0_rec;
xib = circshift(xia, N_interp/3);
xic = circshift(xia, 2*N_interp/3);
xicm = (xia + xib + xic)/3;



figure(3)
clf
tiledlayout(2, 2)
nexttile
hold on

plot([0, 2*pi], [0, 0], ':k', 'HandleVisibility','off')
plot(outq.aa, outq.uu, 'LineWidth',3, 'Color', c(1, :))
ylabel('$u(\theta)$', 'Interpreter', 'latex', 'FontSize',14);
xlabel('$\theta$', 'Interpreter', 'latex', 'FontSize',14);
% ylim([-mi, mi])
% legend({'$u(\theta)$', '$I(\theta)$'}, 'interpreter', 'latex', 'location', 'northeast', 'FontSize', 12)
xlim([0, 2*pi])
box off
ylim([-1, 1])

nexttile 
plot([0, 2*pi], [0, 0], ':k', 'HandleVisibility','off')
plot(th, xcm, 'LineWidth',3, 'Color', c(4, :))
ylabel('$u_{cm}(\theta)$', 'Interpreter', 'latex', 'FontSize',14);
xlabel('$\theta$', 'Interpreter', 'latex', 'FontSize',14);
box off
xlim([0, 2*pi])
ylim([-1, 1])

nexttile
hold on
plot([0, 2*pi], [0, 0], ':k', 'HandleVisibility','off')
plot(outq.alpha, outq.current, 'LineWidth',3, 'Color', c(2, :))
ylabel('$I(\theta)$', 'Interpreter', 'latex', 'FontSize',14);
xlabel('$\theta$', 'Interpreter', 'latex', 'FontSize',14);
box off
xlim([0, 2*pi])
ylim([-1, 1])

nexttile 
plot([0, 2*pi], [0, 0], ':k', 'HandleVisibility','off')
plot(th, xicm, 'LineWidth',3, 'Color', c(5, :))
ylabel('$I_{cm}(\theta)$', 'Interpreter', 'latex', 'FontSize',14);
xlabel('$\theta$', 'Interpreter', 'latex', 'FontSize',14);
box off
xlim([0, 2*pi])
ylim([-1, 1])