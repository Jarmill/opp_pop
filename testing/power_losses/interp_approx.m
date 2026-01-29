N = 10000;
I = linspace(200, 3000, N)';

c = [0.72096, -0.000113, -0.0806, 0.0469];

V = c(1) + c(2) * I + c(3) * log(1 + I) + c(4) * sqrt(I);
PW = V .* I;
% basis = [ones(N, 1), I, I.^2,  sqrt(I)];
basis = [ones(N, 1), I, I.^2,  sqrt(I)];
% basis = [ones(N, 1), I, 2*I.^2 - 1, 4*I.^2 - 3,  sqrt(I)];

beta = basis \ V;

err = norm(basis * beta - V);
err_pw = norm(I .* (basis * beta - V));
V_rec = basis * beta;
PW_rec = I.* V_rec;
figure(5)
clf
subplot(1, 2, 1)
plot(I, V - V_rec)
subplot(1, 2, 2)
plot(I, I.*(V - V_rec))