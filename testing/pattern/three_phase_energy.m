%figure out energy and harmonics

%quarter wave
aq = [0.449391882818937	0.715879300295729	0.877351945248463	2.26424070834133	2.42571335329406	2.69220077077086	3.59098453640873	3.85747195388552	4.01894459883826	5.40583336193112	5.56730600688386	5.83379342436065];
uq = [0 1 0 1 0 1 0 -1 0 -1 0 -1 0];

Iq = [   -0.9599
   -0.9599
   -0.6934
   -0.6934
    0.6935
    0.6935
    0.9599
    0.9599
    0.6935
    0.6935
   -0.6934
   -0.6934
   -0.9599
   -0.9599];

N_interp = 9000;
th = linspace(0, 2*pi, N_interp);
x = pulse_func(th, uq, aq);


I0_rec =    -0.3056;
xi = pi*(cumsum(2*x)/(N_interp) + I0_rec);

xa = xi;
xb = circshift(xa, N_interp/3);
xc = circshift(xa, 2*N_interp/3);

xcm = (1/3)*(xa + xb + xc);

%% harmonics and energy
nmax = 900;
[na, nb] = pulse_harmonics(nmax, uq, aq);



di = 0:nmax;
di = di(mod(di, 3) ~= 0);
di = di(2:end); %drop the first harmonic

energy_I = sum((nb(2:end)./(1:nmax)').^2);
energy_3 = sum((nb(di+1)./di').^2);




%% plots

c = linspecer(7);

figure(3)
clf

tiledlayout(3, 1)
nexttile
hold on
plot(th, modulation*sin(th), 'k', 'linewidth', 3);
plot(th, x, 'linewidth', 3, 'color', cc(1, :))

ylabel('$u(\theta)$', 'Interpreter', 'latex', 'FontSize',14);

xlim([0, 2*pi]) 
title(sprintf('M=%0.1f, k=%d, Lower=%0.3f\\%%, Upper=%0.3f\\%%', modulation, opts.k, bn_lower, bn_upper), ...
    'FontSize',16, 'Interpreter', 'latex')

nexttile
hold on
plot(th, -modulation*cos(th), 'k', 'linewidth', 3);
plot(th, xi, 'linewidth', 3, 'color', cc(2, :));
ylabel('$I(\theta)$', 'Interpreter', 'latex', 'FontSize',14);
xlabel('$\theta$', 'Interpreter', 'latex', 'FontSize',14);
xlim([0, 2*pi])

nexttile
hold on
% plot(th, xa, 'linewidth', 3, 'color', cc(2, :));
plot(th, xcm, 'linewidth', 3, 'color', cc(3, :));
% plot(th, xb, 'linewidth', 3, 'color', cc(3, :));
% plot(th, xc, 'linewidth', 3, 'color', cc(4, :));
xlim([0, 2*pi])
ylabel('$I_{cm}(\theta)$', 'Interpreter', 'latex', 'FontSize',14);


figure(2)
clf
subplot(2, 1,  1)
hold on
stem(0:nmax, na)
title('Cosine Harmonics')
xlabel('n')
ylabel('a_n')
subplot(2, 1, 2)
stem(0:nmax, nb)
title('Sine Harmonics')
xlabel('n')
ylabel('b_n')