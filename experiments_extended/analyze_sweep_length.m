load('order_sweep_5_load_length.mat');

en =  cellfun(@(c) c.energy_lower, result_std.out);

modulation = 0.8;
orderlist = 1:6;
klist = 8:4:40;
kappa = 1;

en_lb = modulation.^2 * pi * 1/(1+kappa^2);

en_diff = en - en_lb;
en_diff0 = max(en_diff, zeros(size(en_diff)));
tdd_lower = sqrt(en_diff0/pi);

FS = 14;

figure(1)
clf
hold on
cmap = parula(length(klist));
for i = 1:length(klist)
    plot(orderlist, tdd_lower(i, :), '.-', 'LineWidth', 2, 'MarkerSize', 30, 'color', cmap(i, :));
end
xlabel('degree $\beta$', 'Interpreter', 'latex', 'fontsize', FS)
ylabel('TDD[$I$]', 'Interpreter', 'latex', 'fontsize', FS)
cbar = colorbar;
cbar.Ticks = linspace(0, 1, length(klist));
tlc = {};
for i = 1:length(klist)
    tlc{i} = sprintf('%d', klist(i));
end
% cbar.TickLabels = arrayfun(@(d) sprintf('%d', d), klist, 'UniformOutput', 'False');
cbar.TickLabels = tlc;
cbar.Label.String = '$k$';
cbar.Label.Interpreter= 'latex';
cbar.Label.FontSize= FS;
% text(4, 0.008, '$k=8$', 'Interpreter', 'latex', 'fontsize', FS)
% text(4, 0.055, '$k=40$', 'Interpreter', 'latex', 'fontsize', FS)


figure(2)
tl = tiledlayout(1, 2);
nexttile;
hold on
% figure(1)
% plot(orderlist, result_std.solver_time, '.-', 'LineWidth', 2, 'MarkerSize', 30)
for i = 1:length(klist)
    plot(orderlist, result_std.preprocess_time(i, :), '.-', 'LineWidth', 2, 'MarkerSize', 30, 'color', cmap(i, :));
end
ylim([0, 840])
set(gca, 'yscale', 'Log')
xlabel('degree $\beta$', 'Interpreter', 'latex', 'fontsize', FS)
ylabel('Preprocess Time (sec.)', 'Interpreter', 'latex', 'fontsize', FS)
nexttile;
hold on
for i = 1:length(klist)
    plot(orderlist, result_std.solver_time(i, :), '.-', 'LineWidth', 2, 'MarkerSize', 30, 'color', cmap(i, :));
end
ylim([0, 840])
set(gca, 'yscale', 'Log')

xlabel('degree $\beta$', 'Interpreter', 'latex', 'fontsize', FS)
ylabel('Solver Time (sec.)', 'Interpreter', 'latex', 'fontsize', FS)

cbar = colorbar;
cbar.Ticks = linspace(0, 1, length(klist));
tlc = {};
for i = 1:length(klist)
    tlc{i} = sprintf('%d', klist(i));
end
cbar.TickLabels = tlc;
cbar.Label.String = '$k$';
cbar.Label.Interpreter= 'latex';
cbar.Label.FontSize= FS;

% cbar = colorbar;
% text(4, 0.008, '$k=8$', 'Interpreter', 'latex', 'fontsize', FS)
% text(4, 0.055, '$k=40$', 'Interpreter', 'latex', 'fontsize', FS)

% figure(2)
% plot(klist, tdd_lower')