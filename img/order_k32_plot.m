order_list = 1:7;
tdd_bound = [NaN
    0.0114538
0.0115790
0.0115846
0.0115871
0.0115951
0.0115948];

tdd_resolve = [NaN
0.0115790
0.0115900
0.0115922
0.0115872
0.0115912
0.0115928];


solver_time = [0.2660 0.5110 2.3176 6.9334 29.2760 155.1534 607.1231];

preprocess_time = [3.7216 10.4791 25.4855 54.2122 103.8996 180.6256 305.7879];

	


polish= 1.16004e-2;

cc = linspecer(5);
figure(30)
clf
tiledlayout(2, 1)
nexttile
hold on
plot(order_list, tdd_bound, '-o', 'Color',cc(1, :), 'MarkerFaceColor','auto', 'linewidth', 2)
% plot(order_list, tdd_resolve, '-o', 'Color',cc(2, :), 'MarkerFaceColor','auto')
plot(order_list, polish*ones(size(order_list)), ':k', 'linewidth', 2)
xlabel('degree $\beta$', 'interpreter', 'latex', 'FontSize',14)
ylabel('Q', 'interpreter', 'latex', 'FontSize',14)
legend({'SDP (lower bound)', 'Recovered (upper bound)'}, 'Location','southeast', 'interpreter', 'latex', 'fontsize', 14)

nexttile
hold on
semilogy(order_list, solver_time, '-o', 'Color',cc(4, :), 'MarkerFaceColor','auto', 'linewidth', 2)
semilogy(order_list, preprocess_time, '-s', 'Color',cc(5, :), 'MarkerFaceColor','auto', 'linewidth', 2)
legend({'Preprocess', 'Solve'}, 'Location','southeast', 'interpreter', 'latex', 'fontsize', 14)
xlabel('degree $\beta$', 'interpreter', 'latex', 'FontSize',14)
ylabel('Time (seconds)', 'interpreter', 'latex', 'FontSize',14)
set(gca, 'Yscale', 'log')
