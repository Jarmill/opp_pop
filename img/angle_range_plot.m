%circle

theta = linspace(0, 2*pi, 400);
% circ = [cos(theta); sin(theta)];
c = cos(theta);
s = sin(theta);


%parameters
Delta = pi/6;
th10 = 0;
th20 = 2*pi - 2*Delta;

th11 = Delta;
th21 = 2*pi - Delta;

arc0 = angle_range_arc(c, s, th10, th20);
arc1 = angle_range_arc(c, s, th11, th21);

th15 = 5*Delta;
th25 = 2*pi;
arc5 = angle_range_arc(c, s, th15, th25);

figure(1)

cc = linspecer(3);

tshift = -0.3;
clf
tiledlayout(1, 3)
nexttile
plot(c(arc0>=0), s(arc0>=0), 'linewidth', 3, 'color', cc(1, :))
xlim([-1, 1])
ylim([-1, 1])
xlabel('$\cos(\theta)$', 'Interpreter', 'latex')
ylabel('$\sin(\theta)$', 'Interpreter', 'latex')
text(tshift , 0, '$m$ = 0', 'Interpreter', 'latex', 'fontsize', 16)
axis square
axis off

nexttile
plot(c(arc1>=0), s(arc1>=0), 'linewidth', 3, 'color', cc(2, :))
xlim([-1, 1])
ylim([-1, 1])
xlabel('$\cos(\theta)$', 'Interpreter', 'latex')
% ylabel('$\sin(\theta)$', 'Interpreter', 'latex')
text(tshift , 0, '$m$ = 1', 'Interpreter', 'latex', 'fontsize', 16)
axis square
axis off

nexttile
plot(c(arc5>=0), s(arc5>=0), 'linewidth', 3, 'color', cc(3, :))
xlim([-1, 1])
ylim([-1, 1])
xlabel('$\cos(\theta)$', 'Interpreter', 'latex')
% ylabel('$\sin(\theta)$', 'Interpreter', 'latex')
text(tshift , 0, '$m$ = 5', 'Interpreter', 'latex', 'fontsize', 16)
axis square
axis off