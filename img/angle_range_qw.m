%circle

theta = linspace(0, 2*pi, 400);
% circ = [cos(theta); sin(theta)];
c = cos(theta);
s = sin(theta);


d = 4;
dt = 2;


%parameters
Delta = pi/12;
th10 = 0;
th20 = pi/2 - (d)*Delta;

th11 = Delta;
th21 = pi/2 - (d-1)*Delta;
% th21 = 2*pi - Delta;

arc0 = angle_range_arc(c, s, th10, th20);
arc1 = angle_range_arc(c, s, th11, th21);

th15 = dt*Delta;
th25 = pi/2 - (d-dt)*Delta;
% th25 = 2*pi - mod(kt, 2)*Delta;
arc5 = angle_range_arc(c, s, th15, th25);

th1k = d*Delta;
th2k = pi/2;
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
xlim([-0.2, 1.1])
scatter(1, 0, 100, 'k', 'filled')
scatter(0, 1, 100, 'k', 'linewidth', 2)
ylim([-0.2, 1.1])
scatter(1, 0, 100, 'k', 'filled')
xlabel('$c$', 'Interpreter', 'latex')
ylabel('$s$', 'Interpreter', 'latex')
title(sprintf('$i$ = %d', 0), 'Interpreter', 'latex', 'fontsize', 16)
axis square
% axis off

nexttile
hold on
plot(c(arc1>=0), s(arc1>=0), 'linewidth', 3, 'color', cc(2, :))
xlim([-0.2, 1.1])
scatter(1, 0, 100, 'k', 'filled')
scatter(0, 1, 100, 'k', 'linewidth', 2)
ylim([-0.2, 1.1])
scatter(1, 0, 100, 'k', 'filled')
xlabel('$c$', 'Interpreter', 'latex')
% ylabel('$s$', 'Interpreter', 'latex')
title(sprintf('$i$ = %d', 1), 'Interpreter', 'latex', 'fontsize', 16)
axis square
% axis off

nexttile
hold on
plot(c(arc5>=0), s(arc5>=0), 'linewidth', 3, 'color', cc(3, :))
xlim([-0.2, 1.1])
scatter(1, 0, 100, 'k', 'filled')
scatter(0, 1, 100, 'k', 'linewidth', 2)
ylim([-0.2, 1.1])
scatter(1, 0, 100, 'k', 'filled')
xlabel('$c$', 'Interpreter', 'latex')
% ylabel('$s$', 'Interpreter', 'latex')
title(sprintf('$i$ = %d', dt), 'Interpreter', 'latex', 'fontsize', 16)
axis square
% axis off

nexttile
hold on
plot(c(arck>=0), s(arck>=0), 'linewidth', 3, 'color', cc(4, :))
xlim([-0.2, 1.1])
scatter(1, 0, 100, 'k', 'filled')
scatter(0, 1, 100, 'k', 'linewidth', 2)
ylim([-0.2, 1.1])
xlabel('$c$', 'Interpreter', 'latex')
% ylabel('$s$', 'Interpreter', 'latex')
% text(0.2 , 0.3, sprintf('$i$ = %d', d), 'Interpreter', 'latex', 'fontsize', 14)
title(sprintf('$i$ = %d', d), 'Interpreter', 'latex', 'fontsize', 16)
axis square
% axis off