function [x] = pulse_func(th, u, alpha)
%PULSE_FUNC plot a pulse train with step values u and angle changes alpha
%
%alpha in [0, 2*pi]
th = mod(th, 2*pi);
N = length(th);
x = zeros(1, N);
a0 = [0, alpha, 2*pi];
nalpha = length(alpha);

for i = 1:nalpha
    x((th<=a0(i+1)) & (th>=a0(i))) = u(i);
end
end

