%figure out periodicity and zero-average-value

%also find explicit expression for the signal

%% make the signal

ah = out.pattern.alpha';
uh = out.pattern.u';
I0 = out.pattern.I(1);
[outh] = pulse_current_voltage_RL(uh, ah, 0, 1000, kappa, I0);
outh.uu = kron(outh.u, [1; 1]);    
outh.aa = [0; kron(outh.alpha(2:end-1), [1; 1]); 2*pi];

%% load characteristics

% R = 4;
% L = 2;
% kappa = R/L;

T = 1;
N = 1000;
t = linspace(0, T*2*pi, N);

% I0 = -3;
% I0 = -0.471070901026913;
% I0 = -0.473;
out
s = tf('s');
sys = 1/(s + kappa);

%compute the current and the energy
I_s = I0*ones(size(outh.alpha));
Na = length(outh.alpha)-1;
E_s = zeros(Na, 1);

da = diff(outh.alpha);

for i = 1:Na
    ucurr = outh.u(i);
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
u = pulse_func(t, outh.u, outh.alpha(2:end-1));

% y = lsim(sys, u, t, I0);
A = -kappa;
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
% plot(reshape((0:(T-1))*2*pi + outh.aa, [], 1), reshape(repmat(outh.uu, [1, T]), [], 1),...
%     'LineWidth', 3, 'Color', c(1, :));
plot(t, y, 'LineWidth', 3,  'Color', c(2, :));
scatter(outh.alpha, I_s, 200, 'k', 'filled')

% figure(2)
% plot(I0_list, y_list(:, end), 'color', 'k')
% xlabel('$I(0)$', 'Interpreter', 'latex')
% ylabel('$I(2\pi)$', 'Interpreter', 'latex')