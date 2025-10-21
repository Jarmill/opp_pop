load('order_sweep_5_load_mod.mat');

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
tdd_lower = sqrt(en_diff0/pi);


load('she_output_sweep_ex_2.mat')


Nd = length(drange);
NM = length(Mrange);

en_she = NaN([Nd, NM]);
NH = 2*Nd+1;
na_she = NaN([Nd, NM, NH+1]);
nb_she = NaN([Nd, NM, NH+1]);

mod_list = Mrange;

kappa = 1;
for i = 1:Nd
    for j= 1:NM
% for i = 4
    % for j = 20
        pcurr = pulse{i,j};
        if ~isempty(pcurr)
            en_min = Inf;
            for k = 1:length(pcurr)            
                %find the minimum energy
                u_base = pcurr{k}.u;
                alpha_base = pcurr{k}.alpha;

                I0_rec = fzero(@(I0) get_root(u_base, alpha_base, kappa, I0), -0.3);
                out_she = pulse_current_voltage_RL(u_base, alpha_base,2, 1000, kappa, I0_rec);
                
                [na, nb] = pulse_harmonics(NH, out_she.u', out_she.alpha(2:end-1)');
                en_she(i, j) = out_she.energy;
                na_she(i, j, :) = na;
                nb_she(i, j, :) = nb;
            end
        end
    end
end

save('she_output_sweep_ex2_result.mat', 'en_she', 'na_she', 'nb_she', 'drange', 'Mrange')


%% visualize
% load("experiments\experiment_N_5_sweep_3_std.mat",...
    % 'mod_list', 'k_list', 'objective_std_3');

tdd = sqrt(en_she/pi - (Mrange').^2/(1+kappa^2));
% tdd2= (energy*(1/pi) - (mod_list').^2);
% tdd2(energy==0) = NaN;
% tdd = sqrt(tdd2);

tdd_diff = tdd - tdd_lower;

%% plot
figure(6)
clf
colormap('parula')
% Lt = log10(tdd);
% h = imagesc(Lt);
h = imagesc(tdd_diff);
set(h, 'AlphaData', ~isnan(tdd_diff))
cb = colorbar(); 
ylabel('$k$', 'interpreter', 'latex')
xlabel('$M$', 'interpreter', 'latex')
xticks(1:length(mod_list))
xticklabels(mod_list)
yticklabels(drange)
yl = ylabel(cb,'$TDD[I_{SHE}] - p^*_3$','FontSize',16, 'interpreter', 'latex');

xll = diff(xlim);
yll = diff(ylim);
pbaspect([xll, yll, 1])

% %% stop
% %drop 1.1, it is problematic
% 
% tdd = tdd(1:end-1, :);
% Nmod = length(mod_list)-1;
% figure(40)
% clf
% % imagesc(tdd)
% % linecolors = autumn(Nmod);
% linecolors = parula(Nmod);
% hold on 
% for i=1:Nmod
%     semilogy(drange, tdd(:, i), 'color',linecolors(i,:));
% end 
% set(gca, 'YScale', 'log')
% xlabel('$d$', 'interpreter', 'latex','FontSize',16)
% ylabel('$Q_3$ lower bound', 'interpreter', 'latex','FontSize',16)
% cb = colorbar(); 
% yl = ylabel(cb,'$M$','FontSize',16, 'interpreter', 'latex');
% yl.Position(1) = min(xlim(cb));
% yl.VerticalAlignment = 'bottom';
% % cb.TickLabels=mod_list(1:end-1);
% % clabel('$M$', 'interpreter', 'latex')
% clim([0, 1.05]);




%% after
%one valid solution at d=7
% M = root_list{end}.M;
% 
% u_base = root_list{end}.u;
% alpha_base = root_list{end}.alpha;
% 
% N = 1000;
% kappa = 1;
% % I0 = -0.3;

% %energy 1.005950354968238
% %polish: 1.005872865526797
% %order 2: 1.005715663590256
% %I0_rec = -0.396822128069477
% %SHE is successful, b_(1:13) = 0 and QW is engaged
function gap = get_root(u_base, alpha_base, kappa, I0)
    out = pulse_current_voltage_RL(u_base, alpha_base,2, 1000, kappa, I0);
    gap = out.I(end) - I0;  
end