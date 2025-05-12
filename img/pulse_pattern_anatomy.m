a1 = [0.15 0.45 0.65 0.95]'*2*pi;

uf = [0 1 0 -1 0]';

a1_full = [0; kron(a1, [1; 1]); 2*pi];

da1 = diff([0; a1; 2*pi]);
ca1 = cumsum(da1);

cc = linspecer(5);
% tiledlayout()
figure(1)
clf
hold on
plot(a1_full, u1_full, 'linewidth', 2)
% axis off
xl = xlabel('$\theta$', 'interpreter', 'latex');
% xl.Position(2) = xl.Position(2) + 0.2;
ylabel('$u(\theta)$', 'interpreter', 'latex')

dx = -0.03*2*pi;
dy = 0.1;
for i = 1:5
    utext = strcat('$u^', num2str(i-1), '$');
    xc = ca1(i) - da1(i)/2;
    text(xc + dx, uf(i) + dy, 0, utext, "interpreter","latex", ...
        "FontSize", 27)
end



dx = 0.01*2*pi;
dy = 0;
for i = 1:4
    atext = strcat('$\alpha^', num2str(i), '$');
    text(a1(i) + dx, uf(i) + dy, 0, atext, "interpreter","latex", ...
        "FontSize", 27)
end

ylim([-1.15, 1.15])
xlim([0, 2*pi])
axis off