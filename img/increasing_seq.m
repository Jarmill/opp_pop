load('experiment_N_5_sweep_3_std.mat')

mi = 11;
ki = [1, 3, 6, 10];

c = linspecer(2);
c = c(2, :);

figure(3)
clf
tiledlayout(2, 2)

for i = 1:4
    % kk = ki(i);
    nexttile
    hold on
    ppi = out_std_3{mi, ki(i)};
    % aa = [0; kron(ppi.pattern.alpha', [1; 1]); pi*2];
    aa = [0; ppi.pattern.alpha'; pi*2];
    II = ppi.pattern.I;
    % text(0, 0, 0, '$k=', num2str(klist(ki)), '$')
    % uu = kron(ppi.pattern.u', [1; 1]);
    tstr = strcat('$k=', sprintf('%d', k_list(ki(i))), '$');
    text(pi-0.6, 0, 0, tstr, 'interpreter', 'latex', 'fontsize', 14)

    plot(aa, II, 'LineWidth', 3, 'color', c)
    axis off

end