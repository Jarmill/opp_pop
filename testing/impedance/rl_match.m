%figure out periodicity and zero-average-value

%also find explicit expression for the signal

%% make the signal
aq = [0.423323100643366	0.646289267844426	0.806615717344141	1.55083012024646	]';
uq = [0 1 0  1 0]';

[outq] = pulse_current_voltage_RL(uq, aq, 2);
outq.uu = kron(outq.u, [1; 1]);    
outq.aa = [0; kron(outq.alpha(2:end-1), [1; 1]); 2*pi];

%% load characteristics

R = 4;
L = 2;
kappa = R/L;

T = 1;
N = 1000;
t = linspace(0, T*2*pi, N);

% I0 = -3;
% I0 = -0.471070901026913;
I0 = -0.473;
s = tf('s');
sys = 1/(s + kappa);

%compute the current and the energy
I_s = I0*ones(size(outq.alpha));
Na = length(outq.alpha)-1;
E_s = zeros(Na, 1);

da = diff(outq.alpha);

for i = 1:Na
    ucurr = outq.u(i);
    Iprev = I_s(i);
    dt = da(i);
    I_s(i+1) = ucurr*(1-exp(-kappa*dt))/kappa + Iprev*exp(-kappa*dt);
    

    E0 = (ucurr-kappa*Iprev)*(3*ucurr + kappa*Iprev);
    E_denom = 2*kappa^3;
    E_num_1 = 2*ucurr^2*kappa*dt;
    E_num_2 = exp(-2*kappa*dt)*(ucurr - kappa*Iprev) * ...
        (kappa*Iprev + ucurr*(4*exp(kappa*dt)-1));
    
    E_s(i) = (E_num_1 + E_num_2 - E0)/E_denom;
    % E_s(i)= (exp(-2*kappa*dt)*(ucurr*(exp(kappa*dt)-1)+(kappa*Iprev))^3 ...
        % - (kappa*Iprev)^3)/(3*kappa^3);
end

u0 = 1;
energy  = sum(E_s);
% u = ones(size(t))*u0;

% u = square(t*(2*pi))*u0;
u = pulse_func(t, outq.u, outq.alpha(2:end-1));

% y = lsim(sys, u, t, I0);
A = -R/L;
B = 1;
C = 1;
D = 0;

NI = 100;
I0_list = linspace(-5, 5, NI);
y_list = zeros(NI, N);
for i = 1:NI
    y_list(i, :) = lsim(A, B, C, D, u, t, I0_list(i))';
end

y = lsim(A, B, C, D, u, t, I0)';
% y = lsim(A, B, C, D, u, t, I0);

%% plot the voltage and the current
figure(1)
c = linspecer(2);
clf
hold on
plot(reshape((0:(T-1))*2*pi + outq.aa, [], 1), reshape(repmat(outq.uu, [1, T]), [], 1),...
    'LineWidth', 3, 'Color', c(1, :));
plot(t, y, 'LineWidth', 3,  'Color', c(2, :));
scatter(outq.alpha, I_s, 200, 'k', 'filled')

% figure(2)
% plot(I0_list, y_list(:, end), 'color', 'k')
% xlabel('$I(0)$', 'Interpreter', 'latex')
% ylabel('$I(2\pi)$', 'Interpreter', 'latex')