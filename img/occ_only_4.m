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


figure(40)
clf

cc = linspecer(5);


imagesc(mass_1)
% title('Pure ($U_1$) Mass', 'interpreter', 'latex')

cb = colorbar;
clim([0, 2.6])
% axis off
ylabel('$n$', 'interpreter', 'latex', 'fontsize', 14)
xlabel('$i$', 'interpreter', 'latex', 'fontsize', 14)
ylabel(cb,'Dwell Angle','FontSize',14,'interpreter', 'latex')
% cb.Position(1) = cb.Position(1) -0.05;
xl = xlim;
% set(gca, 'XTicks', 1:5)
% set(gca, 'XTicks', 1:3)
yl = ylim;
pbaspect([diff(xl), diff(yl), 1])

