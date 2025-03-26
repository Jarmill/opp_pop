%circle

theta = linspace(0, 2*pi, 400);
% circ = [cos(theta); sin(theta)];
c = cos(theta);
s = sin(theta);


k = 8;
kt = 5;


%parameters
Delta = pi/12;
th10 = 0;
th20 = 2*pi - (k)*Delta;

th11 = Delta;
th21 = 2*pi - (k-1)*Delta;
% th21 = 2*pi - Delta;

arc0 = angle_range_arc(c, s, th10, th20);
arc1 = angle_range_arc(c, s, th11, th21);

th15 = kt*Delta;
th25 = 2*pi - (k-kt)*Delta;
% th25 = 2*pi - mod(kt, 2)*Delta;
arc5 = angle_range_arc(c, s, th15, th25);

th1k = k*Delta;
th2k = 2*pi - (0)*Delta;
arck = angle_range_arc(c, s, th1k, th2k);

figure(1)

cc = linspecer(4);

tshift = -0.4;
clf
% tiledlayout(2, 2)
tiledlayout(1, 4)
nexttile
hold on
plot(c(arc0>=0), s(arc0>=0), 'linewidth', 3, 'color', cc(1, :))
xlim([-1, 1])
ylim([-1, 1])
scatter(1, 0, 100, 'k', 'filled')
xlabel('$\cos(\theta)$', 'Interpreter', 'latex')
ylabel('$\sin(\theta)$', 'Interpreter', 'latex')
text(tshift , 0, '$i$ = 0', 'Interpreter', 'latex', 'fontsize', 16)
axis square
axis off

nexttile
hold on
plot(c(arc1>=0), s(arc1>=0), 'linewidth', 3, 'color', cc(2, :))
xlim([-1, 1])
ylim([-1, 1])

scatter(1, 0, 100, 'k', 'filled')
xlabel('$\cos(\theta)$', 'Interpreter', 'latex')
% ylabel('$\sin(\theta)$', 'Interpreter', 'latex')
text(tshift , 0, '$i$ = 1', 'Interpreter', 'latex', 'fontsize', 16)
axis square
axis off

nexttile
hold on
plot(c(arc5>=0), s(arc5>=0), 'linewidth', 3, 'color', cc(3, :))
xlim([-1, 1])
ylim([-1, 1])
scatter(1, 0, 100, 'k', 'filled')
xlabel('$\cos(\theta)$', 'Interpreter', 'latex')
% ylabel('$\sin(\theta)$', 'Interpreter', 'latex')
text(tshift , 0, sprintf('$i$ = %d', kt), 'Interpreter', 'latex', 'fontsize', 16)
axis square
axis off

nexttile
hold on
plot(c(arck>=0), s(arck>=0), 'linewidth', 3, 'color', cc(4, :))
xlim([-1, 1])
scatter(1, 0, 100, 'k', 'filled')
ylim([-1, 1])
xlabel('$\cos(\theta)$', 'Interpreter', 'latex')
% ylabel('$\sin(\theta)$', 'Interpreter', 'latex')
text(tshift , 0, sprintf('$i$ = %d', k), 'Interpreter', 'latex', 'fontsize', 16)
axis square
axis off