c = linspecer(5);
c(1, :) = [60,78,220]/255;
Nth = 1200;
NP = 2;
th = linspace(0, NP*2*pi, Nth);
xa = sin(th);
xb = sin(th - 2*pi/3);
xc = sin(th - 4*pi/3);

figure(1)
clf
hold on
plot(th, xa, 'linewidth', 3, 'color', c(1, :))
plot(th, xb, 'linewidth', 3, 'color', c(2, :))
plot(th, xc, 'linewidth', 3, 'color', c(3, :))
plot([0, 2*pi*NP], [0, 0], 'linewidth', 3, 'color', c(5, :))
axis off
