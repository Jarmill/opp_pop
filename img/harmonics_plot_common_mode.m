aq = [0.423323100643366	0.646289267844426	0.806615717344141	1.55083012024646	]';
      uq = [0 1 0  1 0]';

I0_rec = 0;
% I0_rec = M.modes{1}{3}.init(1,5);
%need to perform appropriate scaling
% xi = pi*(cumsum(2*x)/(N_interp)) + I0_rec;
% xi = xi - max(xi)/2;

% aa_quad = [0; kron(aq, [1; 1]); pi]
% u_quad = kron(uq, [1; 1]);
% pa = [aa_quad; pi - (aa_quad(end:-1:1)); pi + aa_quad; 2*pi - aa_quad(end:-1:1)];
% pu= [u_quad; u_quad(end:-1:1); -u_quad; -u_quad(end:-1:1)];

% pa = [aq; pi - (aq(end:-1:1)); pi + aq; 2*pi - aq(end:-1:1)];

c = linspecer(5);

figure(3)
clf
% tiledlayout(2, 1)
% nexttile
% hold on
% [outq] = pulse_current_voltage(uq, aq, 2);
%  outq.uu = kron(outq.u, [1; 1]);    
%     outq.aa = [0; kron(outq.alpha(2:end-1), [1; 1]); 2*pi];
% 
% plot([0, 2*pi], [0, 0], ':k', 'HandleVisibility','off')
% plot(outq.aa, outq.uu, 'LineWidth',3, 'Color', c(1, :))
% plot(outq.alpha, outq.current, 'LineWidth',3, 'Color', c(2, :))
% ylabel('Signal', 'Interpreter', 'latex', 'FontSize',14);
% xlabel('$\theta$', 'Interpreter', 'latex', 'FontSize',14);
% % ylim([-mi, mi])
% legend({'$u(\theta)$', '$I(\theta)$'}, 'interpreter', 'latex', 'location', 'northeast', 'FontSize', 12)
% xlim([0, 2*pi])
% box off
% 
% nexttile
hold on
nmax = 41;
[na, nb] = pulse_harmonics(nmax, outq.u, outq.alpha(2:end-1)');

in = 1:(nmax+1);
in(mod(in-1, 2)==0)= [];
nb0 = nb;
nb = nb(in);
nbL = nb./((in-1).^2)';

cm3 = (mod(in, 3)==1);

% stem(in, nb, 'filled', 'Linewidth', 3, 'Color', c(1, :));
% stem(in, nbL, 'filled', 'Linewidth', 3, 'Color', c(2, :));
% scatter(2, 1, 300, 'k', 'linewidth', 3)


stem(in, log10(nbL.^2), 'filled', 'Linewidth', 3, 'Color', c(2, :));
stem(in, log10(nb.^2), 'filled', 'Linewidth', 3, 'Color', c(1, :));
stem(in(cm3), log10(nbL(cm3).^2), 'filled', 'Linewidth', 3, 'Color', c(4, :));
stem(in(cm3), log10(nb(cm3).^2), 'filled', 'Linewidth', 3, 'Color', c(5, :));
scatter(2, 0, 300, 'k', 'linewidth', 3)
% set(gca, 'Yscale', 'log')

xlabel('$\ell$', 'FontSize', 14, 'Interpreter', 'latex')
ylabel('Odd Harmonics', 'FontSize', 14, 'Interpreter', 'latex')
legend({'$\log(b_\ell^2/\ell^2)$', '$\log(b_\ell^2)$', '$\log(b_\ell^2/\ell^2)$ common-mode', '$\log(b_\ell^2)$ common-mode'}, 'interpreter', 'latex', 'location', 'southwest', 'FontSize', 12)
xlim([0, nmax])


box off
% figure(3)
