function [na, nb] = pulse_harmonic(n, u, alpha)
%Find the nth harmonic of the pulse function
%
% na = int cos(n th) x(th) d th
% nb = int sin(n th) x(th) d th
%energy = integrate x(th)^2 dth for th in [0, 2pi]


a0 = [0, alpha, 2*pi];
nalpha = length(alpha);

na = 0;
nb = 0;
for i = 1:nalpha
    % da = (a0(i+1) - a0(i));
    thnext = a0(i+1);
    thprev = a0(i);

    if n==0
        na = na + u(i)*(thnext-thprev)/pi;
        nb = 0;
    else
        na=  na+ (u(i)/n)*(sin(n*thnext) - sin(n*thprev))/pi;
        nb=  nb+ (u(i)/n)*(cos(n*thnext) - cos(n*thprev))/pi;
    end


    % energy = energy + dx*da;        
end

end