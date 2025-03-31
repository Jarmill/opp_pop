load('experiments\sweep_gap.mat')

figure(2)
td = tdd_diff_std;
td(td==Inf) = NaN;

Lt = log10(td);
h = imagesc(Lt');
set(h, 'AlphaData', ~isnan(Lt'))
cb = colorbar(); 
ylabel('$d$', 'interpreter', 'latex')
xlabel('$M$', 'interpreter', 'latex')
xticks(1:length(mod_list))
xticklabels(mod_list)
yl = ylabel(cb,'$log_{10}( Q^{rec} - Q_3)$','FontSize',16, 'interpreter', 'latex');

xll = diff(xlim);
yll = diff(ylim);
pbaspect([xll, yll, 1])
% yl.Position(1) = min(xlim(cb));
% yl.VerticalAlignment = 'bottom';