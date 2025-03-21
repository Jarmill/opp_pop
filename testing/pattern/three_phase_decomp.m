A = 3*eye(3) - ones(3);

V =[-1/sqrt(2), -1/sqrt(2), 1/sqrt(3); 
     0, 1/sqrt(2), 1/sqrt(3); 
     1/sqrt(2), 0, 1/sqrt(3)];
D = diag([3; 3; 0]);

A2 = V'*D*V

[VV, DD] = eig(A)

