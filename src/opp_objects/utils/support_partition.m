function [X_partition] = support_partition(partition, vars)
%SUPPORT_PARTITION segment the circle c^2 + s^2 = 1 into equal proportions
c = vars.x(1);
s = vars.x(2);

mD = 2*pi/partition;
base = cos(mD) - cos(mD)*c + sin(mD)*s;

R = [cos(mD) sin(mD); -sin(mD) cos(mD)];
b_partition = ones(partition, 1)*c;
for p = 1:(partition)
    bt = subs(base, [c; s], (R^(p-1))*[c; s]);
    b_partition(p) = bt;
end

X_partition = b_partition >= 0;



end

