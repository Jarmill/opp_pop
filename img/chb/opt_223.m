% load('experiment_N_7.mat')

aq = [ 0.38119444187739643
 0.9753206096796204
 1.3099716822849397
 1.3659375430458054];

uq = [0 1 2 7/2 2]';

% aa = [0; kron(alpha, [1; 1]); pi*2];
% % uu = kron(u', [1; 1]);

I0 = -1;
outq = pulse_current_voltage(uq, aq, 2, I0);

aa = [0; kron(outq.alpha(2:end-1), [1; 1]); pi*2];
uu = kron(outq.u, [1; 1]);

Ipos =outq.I_val - min(outq.I_val);
Irefl = Ipos - max(Ipos)/2;


figure(4)
clf
c = linspecer(2);
tiledlayout(2, 1)

nexttile
plot(aa, uu, 'LineWidth', 2, 'color', c(1, :))
box off
ylabel('$u(\theta)$', 'FontSize', 14, 'interpreter', 'latex')
xlim([0, 2*pi])

nexttile
plot(outq.alpha, Irefl, 'LineWidth', 2, 'color', c(2, :))
xlabel('$\theta$', 'FontSize', 14, 'interpreter', 'latex')
ylabel('$I(\theta)$', 'FontSize', 14, 'interpreter', 'latex')
% axis off
box off
xlim([0, 2*pi])


%% harmonic content
NH = 25;
[a, b] = pulse_harmonics(NH, outq.u, outq.alpha(2:end-1)');

H = 1+[1, 5, 7, 11];
lb = b;
lb(1:2:NH) = NaN;
lb(H(2:end)) = 0;
figure(5)
clf

berr = b(H);
berr(1) = berr(1) - 2;

% stem(0:NH, b)
stem(0:NH, log10(abs(lb)),'filled', 'markersize', 10)
yl = ylim;
hold on

for i = 2:length(H)
    plot(H(i)*[1,1]-1, yl, '--k', 'linewidth', 2)
end

xlabel('$\ell$', 'FontSize', 14, 'interpreter', 'latex')
ylabel('$\log_{10} | b_{\ell}|$', 'FontSize', 14, 'interpreter', 'latex')
