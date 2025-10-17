% load('order_sweep_5_load_mod.mat');

% modulation = 0.8;
% modulationlist = [0, 0.1, 0.5, 1, 2 5];
% klist = 1:6;
kappa = 1;
en =  cellfun(@(c) c.energy_lower, result_std.out);


% figure(5)
% clf
% plot(klist, en)

% 
en_lb = modulationlist.^2 * pi .* 1./(1+kappa.^2);
% 
en_diff = en - en_lb;
en_diff0 = max(en_diff, zeros(size(en_diff)));
tdd_lower = sqrt(en_diff0);

FS = 14;
figure(5)
clf
tiledlayout(1, 2)
colormap('copper')
cmap = copper(length(modulationlist));
%ENERGY
nexttile
hold on
for i = 1:length(modulationlist)
    plot(klist, en(:, i), '.-', 'LineWidth', 2, 'MarkerSize', 30, 'color', cmap(i, :));
end
xlim([min(klist), max(klist)])
xlabel('length $k$', 'Interpreter', 'latex', 'fontsize', FS)
ylabel('$||I||_2^2$', 'Interpreter', 'latex', 'fontsize', FS)

%TDD
nexttile
hold on

% plot(klist, tdd_lower)
for i = 1:length(modulationlist)
    plot(klist, tdd_lower(:, i), '.-', 'LineWidth', 2, 'MarkerSize', 30, 'color', cmap(i, :));
end
xlabel('length $k$', 'Interpreter', 'latex', 'fontsize', FS)
ylabel('TDD[$I$]', 'Interpreter', 'latex', 'fontsize', FS)
cbar = colorbar;
cbar.Ticks = linspace(0, 1, length(modulationlist));
tlc = {};
for i = 1:4:length(modulationlist)
    tlc{i} = sprintf('%0.2f', modulationlist(i));
end
% cbar.TickLabels = arrayfun(@(d) sprintf('%d', d), klist, 'UniformOutput', 'False');
cbar.TickLabels = tlc;
cbar.Label.String = '$M$';
cbar.Label.Interpreter= 'latex';
cbar.Label.FontSize= FS;

xlim([min(klist), max(klist)])


hold off


% figure(1)
% clf
% tl = tiledlayout(1, 2);
% nexttile;
% hold on
% colormap('autumn')
% cmap = autumn(length(modulationlist));
% for i = 1:length(modulationlist)
%     plot(klist, en(:, i), '.-', 'LineWidth', 2, 'MarkerSize', 30, 'color', cmap(i, :));
% end
% xlabel('degree $\beta$', 'Interpreter', 'latex', 'fontsize', FS)
% ylabel('$||I||_2^2$', 'Interpreter', 'latex', 'fontsize', FS)
% % cbar = colorbar;
% % % cbar.Ticks = linspace(0, 1, length(modulationlist));
% % tlc = {};
% % for i = 1:length(modulationlist)
% %     tlc{i} = sprintf('%d', modulationlist(i));
% % end
% % 
% nexttile;
% hold on
% % colormap('autumn')
% % cmap = autumn(length(modulationlist));
% for i = 1:length(modulationlist)
%     plot(klist, tdd_lower(:, i), '.-', 'LineWidth', 2, 'MarkerSize', 30, 'color', cmap(i, :));
% end
% xlabel('degree $\beta$', 'Interpreter', 'latex', 'fontsize', FS)
% ylabel('TDD[$I$]', 'Interpreter', 'latex', 'fontsize', FS)
% cbar = colorbar;
% cbar.Ticks = linspace(0, 1, length(modulationlist));
% tlc = {};
% for i = 1:length(modulationlist)
%     tlc{i} = sprintf('%0.1f', modulationlist(i));
% end
% % cbar.TickLabels = arrayfun(@(d) sprintf('%d', d), klist, 'UniformOutput', 'False');
% cbar.TickLabels = tlc;
% cbar.Label.String = '$\tau$';
% cbar.Label.Interpreter= 'latex';
% cbar.Label.FontSize= FS;
% % text(4, 0.008, '$k=8$', 'Interpreter', 'latex', 'fontsize', FS)
% % text(4, 0.055, '$k=40$', 'Interpreter', 'latex', 'fontsize', FS)
% 
% % 

%% timing
figure(2)
tl = tiledlayout(1, 2);
nexttile;
hold on
% figure(1)
% plot(klist, result_std.solver_time, '.-', 'LineWidth', 2, 'MarkerSize', 30)
for i = 1:length(modulationlist)
    plot(klist, result_std.preprocess_time(:, i), '.-', 'LineWidth', 2, 'MarkerSize', 30, 'color', cmap(i, :));
end
ylim([3, 1000])
xlim([8, 40])
set(gca, 'yscale', 'Log')
xlabel('length $k$', 'Interpreter', 'latex', 'fontsize', FS)
ylabel('Preprocess Time (sec.)', 'Interpreter', 'latex', 'fontsize', FS)
nexttile;
hold on
for i = 1:length(modulationlist)
    plot(klist, result_std.solver_time(:, i), '.-', 'LineWidth', 2, 'MarkerSize', 30, 'color', cmap(i, :));
end
% ylim([0, 840])
set(gca, 'yscale', 'Log')

xlabel('length $k$', 'Interpreter', 'latex', 'fontsize', FS)
ylabel('Solver Time (sec.)', 'Interpreter', 'latex', 'fontsize', FS)

cbar = colorbar;
cbar.Ticks = linspace(0, 1, length(modulationlist));
tlc = {};
for i = 1:4:length(modulationlist)
    tlc{i} = sprintf('%0.2f', modulationlist(i));
end
% cbar.TickLabels = arrayfun(@(d) sprintf('%d', d), klist, 'UniformOutput', 'False');
cbar.TickLabels = tlc;
cbar.Label.String = '$M$';
cbar.Label.Interpreter= 'latex';
cbar.Label.FontSize= FS;

xlim([8, 40])
ylim([3, 1000])
% 
% % cbar = colorbar;
% % text(4, 0.008, '$k=8$', 'Interpreter', 'latex', 'fontsize', FS)
% % text(4, 0.055, '$k=40$', 'Interpreter', 'latex', 'fontsize', FS)
% 
% % figure(2)
% % plot(klist, tdd_lower')