N = 9000;
th = linspace(0, 2*pi, N);

CM = true;

% % x = sin(th) + 0.5*cos(3*th);

% f = @(t) sin(t) - 0.5*sin(5*t) - 0.2*sin(7*t);
% f = @(t) sin(t) - 0.5*sin(5*t) - 0.2*sin(11*t);
if ~CM
    f = @(t) sin(t) - 0.5*sin(5*t) - 0.2*sin(11*t);
else
    f = @(t) sin(t) - 0.5*sin(5*t) - 0.2*sin(11*t) + 0.2*sin(9*t);
end

x = f(th);

f_cm = @(t) (f(t) + f(t+2*pi/3)+ f(t+4*pi/3))/3;

K = sqrt(2/3)*[1 -0.5 -0.5;
                0 sqrt(3)/2 -sqrt(3)/2; 
                ones(1, 3)];
        % 

% x_opp = circshift(x, N/2);
x_opp = f(th + pi);

x_cm2 = (x + x_opp)/2;

xa = x;
xb = f(th - 2*pi/3);
xc = f(th -4*pi/3);
% xb = circshift(x, N/3);
% xc = circshift(x, 2*N/3);


xclarke = K * [xa; xb; xc];

x_cm3 = (xa + xb + xc)/3;

c = linspecer(6);
figure(30)
clf
subplot(3, 1, 1)
hold on
plot(th, x, 'LineWidth', 3, 'color', c(1, :))
plot(th, x_cm2, 'LineWidth', 3, 'color', c(2, :))
plot(th, x_cm3, 'LineWidth', 3, 'color', c(6, :))
xlim([0, 2*pi])
xlabel('$\theta$', 'interpreter', 'latex')
ylabel('$x(\theta)$', 'interpreter', 'latex')

subplot(3, 1, 2)
hold on
plot(th, [xa; xb; xc], 'LineWidth', 3)
plot(th, x_cm3, 'LineWidth', 3, 'color', c(6, :))
xlabel('$\theta$', 'interpreter', 'latex')
ylabel('$x(\theta+2\pi k/3)$', 'interpreter', 'latex')
xlim([0, 2*pi])

subplot(3, 1, 3)
hold on
plot(th, xclarke(1, :), 'LineWidth', 3, 'color', c(4, :))
plot(th, xclarke(2, :), 'LineWidth', 3, 'color', c(5, :))
plot(th, xclarke(3, :), 'LineWidth', 3, 'color', c(6, :))
xlabel('$\theta$', 'interpreter', 'latex')
ylabel('$x_{\alpha \beta}(\theta)$', 'interpreter', 'latex')
xlim([0, 2*pi])