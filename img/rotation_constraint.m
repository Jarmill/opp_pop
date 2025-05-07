aq = [0.423323100643366	0.646289267844426	0.806615717344141	1.55083012024646	]';
      uq = [0 1 0  1 0]';
c = linspecer(5);
c(1, :) = [60,78,220]/255;

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
xs = [xa; xb; xc];
%current
I0_rec = outq.current(1);
xia = pi*(cumsum(2*x)/(N_interp)) + I0_rec;
xib = circshift(xia, N_interp/3);
xic = circshift(xia, 2*N_interp/3);
xicm = (xia + xib + xic)/3;

s = [1 -1 0]';
s1 = s([2, 3, 1]);
s2 = s([3, 1, 2]);
ta = double(all(xs == s, 1));
tb = double(all(xs == s1, 1));
tc = double(all(xs == s2, 1));
ta(ta==0) = NaN;
tb(tb==0) = NaN;
tc(tc==0) = NaN;

figure(3)
clf
tiledlayout(2, 1)
% nexttile
% hold on
% plot([0, 2*pi], [0, 0], ':k', 'HandleVisibility','off')
% plot(th, xa, 'LineWidth',3, 'Color', c(1, :))
% plot(th, xb, 'LineWidth',3, 'Color', c(2, :))
% plot(th, xc, 'LineWidth',3, 'Color', c(3, :))
% ylabel('$u(\theta)$', 'Interpreter', 'latex', 'FontSize',14);
% xlabel('$\theta$', 'Interpreter', 'latex', 'FontSize',14);
% box off
% xlim([0, 2*pi])
% ylim([-1, 1])
% hold on
% 
% plot([0, 2*pi], [0, 0], ':k', 'HandleVisibility','off')
% plot(outq.aa, outq.uu, 'LineWidth',3, 'Color', c(1, :))
% ylabel('$u(\theta)$', 'Interpreter', 'latex', 'FontSize',14);
% xlabel('$\theta$', 'Interpreter', 'latex', 'FontSize',14);
% % ylim([-mi, mi])
% % legend({'$u(\theta)$', '$I(\theta)$'}, 'interpreter', 'latex', 'location', 'northeast', 'FontSize', 12)
% xlim([0, 2*pi])
% box off
% ylim([-1, 1])

nexttile 
hold on
plot([0, 2*pi], [0, 0], ':k', 'HandleVisibility','off')
plot(th, xia, 'LineWidth',3, 'Color', c(1, :))
plot(th, xib, 'LineWidth',3, 'Color', c(2, :))
plot(th, xic, 'LineWidth',3, 'Color', c(3, :))
ylabel('$I(\theta)$', 'Interpreter', 'latex', 'FontSize',14);
xlabel('$\theta$', 'Interpreter', 'latex', 'FontSize',14);
box off
xlim([0, 2*pi])
ylim([-1, 1])

nexttile
c2 = linspecer(8);
% s = [1 -1 1];
hold on
plot(th, 0.5*ta, 'LineWidth',3)
plot(th, 0.5*tb, 'LineWidth',3)
plot(th, 0.5*tc, 'LineWidth',3)
xlim([0, 2*pi])
legend({'$u=[1, -1, 0]$', '$u=[-1, 0, 1]$', '$u=[0, 1, -1]$'}, 'Interpreter',...
     'latex', 'fontsize', 12)
xlabel('$\theta$', 'Interpreter', 'latex', 'FontSize',14);
ylabel('Location', 'Interpreter', 'latex', 'FontSize',14)
set(gca,'ytick',[])
set(gca,'yticklabel',[])
box off