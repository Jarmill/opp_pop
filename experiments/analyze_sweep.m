load("experiments\experiment_N_5_sweep_3_std.mat",...
    'mod_list', 'k_list', 'objective_std_3');

energy = objective_std_3;
tdd2= (energy*(1/pi) - (mod_list').^2);
tdd2(energy==0) = NaN;
tdd = sqrt(tdd2);



%drop 1.1, it is problematic

tdd = tdd(1:end-1, :);
Nmod = length(mod_list)-1;
figure(40)
clf
% imagesc(tdd)
% linecolors = autumn(Nmod);
linecolors = parula(Nmod);
hold on 
for i=1:Nmod
    % plot(x,x.^(i/3),'color',linecolors(i,:));
    % hold on;
    semilogy(k_list/4, tdd(i, :), 'color',linecolors(i,:));
end 
set(gca, 'YScale', 'log')
xlabel('$d$', 'interpreter', 'latex','FontSize',16)
ylabel('$Q_3$ lower bound', 'interpreter', 'latex','FontSize',16)
cb = colorbar(); 
yl = ylabel(cb,'$M$','FontSize',16, 'interpreter', 'latex');
yl.Position(1) = min(xlim(cb));
yl.VerticalAlignment = 'bottom';
% cb.TickLabels=mod_list(1:end-1);
% clabel('$M$', 'interpreter', 'latex')
clim([0, 1.05]);

