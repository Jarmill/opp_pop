load('experiment_N_7.mat')

aa = [0; kron(out_polish.warm.alpha, [1; 1]); pi*2];
uu = kron(out_polish.warm.u', [1; 1]);

figure
c = linspecer(2);
plot(aa, uu, 'LineWidth', 2, 'color', c(1, :))
axis off