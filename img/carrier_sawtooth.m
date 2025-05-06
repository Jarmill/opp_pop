Nt = 1200;
th = linspace(0, 2*pi, Nt);

% fc = 10;
% fc = 4
% fc = 8;
% fc = 3;
% M = 0.8;
fc = 10;
M = 1;
% N = 3;
N=5;cos(
% N=7;
Nb = (N-1)/2;

xref = M*sin(th);
offset = ([-Nb:(Nb-1)]/Nb)';

% xc0 = (1/(2*Nb))*(sin(th*fc)+1);
xc0 = (1/(2*Nb))*(sawtooth(th*fc, 0.5)+1);
xc = xc0 + offset;
cc = (1/Nb)*(xc<xref) + offset;


figure(1)
clf
tiledlayout(1, 3)
nexttile
hold on
plot(th, xc, 'LineWidth', 2)
plot(th, xref, 'k', 'LineWidth', 2)
xlim([0, 2*pi]);
axis off

nexttile
plot(th, cc, 'LineWidth', 2)
xlim([0, 2*pi]);
axis off

c = linspecer(3);
nexttile
hold on;
plot(th, sum(cc), 'LineWidth',2, 'color', c(3, :))
plot(th, xref, 'k', 'LineWidth', 2)
xlim([0, 2*pi]);
axis off