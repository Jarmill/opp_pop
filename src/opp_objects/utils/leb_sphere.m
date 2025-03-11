function [y] = leb_sphere(a, r)
%LEB_UNIT_CIRC Summary of this function goes here
%   Detailed explanation goes here

if nargin < 2
    r = 1;
end

[m,n] = size(a);
R = r^(n-1);
y = zeros(m,1);
%formula from 
% https://www.johndcook.com/blog/2018/01/31/integrating-polynomials-over-a-sphere-or-ball/
for k = 1:m
    acurr = a(k, :);
    if any(mod(acurr, 2)==1)
        y(k)=0;
    else        
        bcurr = (acurr+1)/2;
        ynum = prod(gamma(bcurr));
        ydenom = gamma(sum(bcurr));

        y(k) = R*2*ynum/ydenom;
    end
end
end

