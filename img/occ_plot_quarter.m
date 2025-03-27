
aq = [0.5 1.3 1.5]';
uq = [0 1 0 1]';
aq2 = [0.3 0.6 0.9]';
uq2 = [0 -1 0 1]';

aa_quad = [0; kron(aq, [1; 1]); pi/2];
u_quad = kron(uq, [1; 1]);
aa_allq = [aa_quad; pi - (aa_quad(end:-1:1)); pi + aa_quad; 2*pi - aa_quad(end:-1:1)];
u_allq = [u_quad; u_quad(end:-1:1); -u_quad; -u_quad(end:-1:1)];

aa_quad2 = [0; kron(aq2, [1; 1]); pi/2];
u_quad2 = kron(uq2, [1; 1]);
aa_allq2 = [aa_quad2; pi - (aa_quad2(end:-1:1)); pi + aa_quad2; 2*pi - aa_quad2(end:-1:1)];
u_allq2 = [u_quad2; u_quad2(end:-1:1); -u_quad2; -u_quad2(end:-1:1)];


% a1_full = [0; kron(aa_allq, [1; 1]); 2*pi];
% a2_full = [0; kron(aa_allq2, [1; 1]); 2*pi];
% u1_full = kron(u_allq, [1; 1]);
% u2_full = kron(u_allq2, [1; 1]);
% u2_full = -u1_full;
figure(3)

da1 = diff([0; aq; pi/2]);
da2 = diff([0; aq2; pi/2]);

mass_1 = zeros(3, length(uq));
mass_2 = mass_1;
for i = 1:length(uq)
    mass_1(2-uq(i), i) = da1(i);
    mass_2(2-uq2(i), i) = da2(i);
end

c = 0.6;
mass_mixed = c*mass_1 + (1-c)*mass_2;
% mass_pure(2+uf, 1:5) = da1;

clf

cc = linspecer(5);
% tiledlayout()
subplot(2, 2, [1, 2])
hold on
plot(aa_allq, u_allq, 'linewidth', 2)
plot(aa_allq2, u_allq2, '--', 'linewidth', 2)
title('Switching Signals $U_1$, $U_2$', 'Interpreter', 'latex', 'FontSize',14)
% axis off
xl = xlabel('$\theta$', 'interpreter', 'latex');
xl.Position(2) = xl.Position(2) + 0.2;
ylabel('$u(\theta)$', 'interpreter', 'latex')

subplot(2, 2, 3)
imagesc(mass_1)
title('Pure ($U_1$) Mass', 'interpreter', 'latex')
% colorbar
clim([0, 1])
% axis off
xlabel('$n$', 'interpreter', 'latex')
ylabel('$i$', 'interpreter', 'latex')

subplot(2, 2, 4)
imagesc(mass_mixed)
title('Mixed ($0.6U_1 + 0.4U_2$) Mass', 'interpreter', 'latex')
% colorbar
clim([0, 1])
xlabel('$n$', 'interpreter', 'latex')
ylabel('$i$', 'interpreter', 'latex')
% /axis off



