


a = [0 0.5 0.5 ]';
u = [ 0 0 1 2]';
a2 = [0 0.5 0.55]';


aa = [0; kron(a, [1; 1]); 1];
uu = kron(u, [1; 1]);
aa2 = [0; kron(a2, [1; 1]); 1];

figure(3)
clf
c = linspecer(2);
tiledlayout(1, 2)
nexttile
plot(aa, uu, 'linewidth', 3, 'color', c(1, :))
text(0.55, 1, 0, 'Prohibited', 'interpreter', 'latex', 'FontSize', 14)
axis off

nexttile
plot(aa2, uu, 'linewidth', 3, 'color', c(1, :))
text(0.45, 1.15, 0, '$\Theta$', 'interpreter', 'latex', 'FontSize', 14)
axis off