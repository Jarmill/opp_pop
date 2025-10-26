% load('order_sweep_5_load_mod.mat');
% load('order_sweep_5_load_mod_polish.mat', 'polish_out');
% modulation = 0.8;
% modulationlist = [0, 0.1, 0.5, 1, 2 5];
% klist = 1:6;
kappa = 1;
en =  cellfun(@(c) c.energy_lower, result_std.out);


% figure(5)
% clf
% plot(klist, en)

%% compute energy
en_lb = modulationlist.^2 * pi .* 1./(1+kappa.^2);
% SDP bounds
en_diff = en - en_lb;
en_diff0 = max(en_diff, zeros(size(en_diff)));
tdd_lower = sqrt(en_diff0/pi);

%polish
Nk = length(klist);
Nmod = length(modulationlist);
en_polish = NaN*ones(Nk, Nmod);
for i = 1:Nk
    for j = 1:Nmod
        if ~isempty(polish_out{i, j}.warm)
            en_polish(i, j) = polish_out{i, j}.warm.objective;
        end
    end
end

en_diff_polish = en_polish - en_lb;
% en_diff0_polish = max(en_diff_polish, zeros(size(en_diff)));
tdd_polish = sqrt(en_diff_polish/pi);

tdd_diff = tdd_polish - tdd_lower;

%SHE
load('she_output_sweep_ex2_result.mat', 'en_she', 'en_she_space', 'na_she', 'nb_she', 'drange', 'Mrange')
tdd_she = sqrt(en_she/pi - (Mrange').^2/(1+kappa^2));
tdd_she_space = sqrt(en_she_space/pi - (Mrange').^2/(1+kappa^2));


%% plots

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

%% plot
figure(7)
clf
tiledlayout(2, 1)
nexttile
colormap('parula')
% Lt = log10(tdd);
% h = imagesc(Lt);
h = imagesc(tdd_diff);
set(h, 'AlphaData', ~isnan(tdd_diff))

ylabel('$k$', 'interpreter', 'latex')
xticks(1:length(modulationlist))
xticklabels(modulationlist)

yticks(1:9)
yticklabels(4*drange)
title('$TDD[I_{loc}] - p^*_3$', 'interpreter', 'latex', 'FontSize', 18)
% cb = colorbar(); 
% cb.Ticks = linspace(0, 0.025, 6);
% yl = ylabel(cb,'$TDD[I_{SHE}] - p^*_3$','FontSize',16, 'interpreter', 'latex');

xll = diff(xlim);
yll = diff(ylim);
pbaspect([xll, yll, 1])

nexttile
h = imagesc(tdd_she_space - tdd_polish);
set(h, 'AlphaData', ~isnan(tdd_she_space))
cb = colorbar(); 
ylabel('$k$', 'interpreter', 'latex')
xlabel('$M$', 'interpreter', 'latex')
xticks(1:length(modulationlist))
xticklabels(modulationlist)
yticks(1:9)
yticklabels(4*drange)

yl = ylabel(cb,'TDD Gap','FontSize',16, 'interpreter', 'latex');

title('$TDD[I_{SHE}] - TDD[I_{loc}]$', 'interpreter', 'latex', 'FontSize', 18)
xll = diff(xlim);
yll = diff(ylim);
pbaspect([xll, yll, 1])
% cb.Ticks = linspace(-0.005, 0.025, 7);
cb.Layout.Tile = 'east';
