a1 = [0.05 0.45 0.55 0.95]'*2*pi;
a2 = [0.3 0.5 0.7 0.9]'*2*pi;

uf = [0 1 0 -1 0]';

a1_full = [0; kron(a1, [1; 1]); 2*pi];
a2_full = [0; kron(a2, [1; 1]); 2*pi];
u1_full = kron(uf, [1; 1]);
u2_full = -u1_full;
figure(3)

da1 = diff([0; a1; 2*pi]);
da2 = diff([0; a2; 2*pi]);

mass_1 = zeros(3, 5);
mass_2 = mass_1;
for i = 1:5
    mass_1(2-uf(i), i) = da1(i);
    mass_2(2+uf(i), i) = da2(i);
end

c = 0.6;
mass_mixed = c*mass_1 + (1-c)*mass_2;
% mass_pure(2+uf, 1:5) = da1;

clf

cc = linspecer(5);
% tiledlayout()
subplot(2, 2, [1, 2])
hold on
plot(a1_full, u1_full, 'linewidth', 2)
plot(a2_full, u2_full, '--', 'linewidth', 2)
title('Switching Signals $U_1$, $U_2$', 'Interpreter', 'latex', 'FontSize',14)
% axis off
xl = xlabel('$\theta$', 'interpreter', 'latex');
xl.Position(2) = xl.Position(2) + 0.2;
ylabel('$u(\theta)$', 'interpreter', 'latex')

subplot(2, 2, 3)
imagesc(mass_1)
title('Pure ($U_1$) Mass', 'interpreter', 'latex')
% colorbar
clim([0, 2.6])
% axis off
xlabel('$n$', 'interpreter', 'latex')
ylabel('$i$', 'interpreter', 'latex')

subplot(2, 2, 4)
imagesc(mass_mixed)
title('Mixed ($0.6U_1 + 0.4U_2$) Mass', 'interpreter', 'latex')
% colorbar
clim([0, 2.6])
xlabel('$n$', 'interpreter', 'latex')
ylabel('$i$', 'interpreter', 'latex')
% /axis off



