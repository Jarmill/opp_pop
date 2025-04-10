A = (3*eye(3) - ones(3))/2;
% [VV, DD] = eig(A)

K = sqrt(2)/3*[1 -0.5 -0.5;
    0 sqrt(3)/2 -sqrt(3)/2];

Kt = [K; ones(1, 3)/3]

K2 = K'*K
Kt2 = Kt'*Kt



% V =[-1/sqrt(2), -1/sqrt(2), 1/sqrt(3); 
%      0, 1/sqrt(2), 1/sqrt(3); 
%      1/sqrt(2), 0, 1/sqrt(3)];
% D = diag([3; 3; 0]);

% A2 = V'*D*V

