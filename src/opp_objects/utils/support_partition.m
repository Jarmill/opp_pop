function [X_partition] = support_partition(partition, x)
%SUPPORT_PARTITION segment the circle c^2 + s^2 = 1 into equal proportions
%
%
%angle interval example: [0, pi/2], [pi/2, pi], [pi, 3pi/2], [3pi/2, 2pi]
%but in a trigonometrically-lifted manner
c = x(1);
s = x(2);

mD = 2*pi/double(partition);

if partition == 1
    X_partition = [];
else
    X_partition = ones(partition, 1)*c;
    for i = 1:partition
        di = double(i);
        th1 = mD*(di-1);
        th2 = mD*di;

        s1 = sin(th1);
        s2 = sin(th2);
        c1 = cos(th1);
        c2 = cos(th2);
    
        offset = (s2-s1)*c1 - (c2-c1)*s1;

        X_partition(i) = (s2-s1)*c - (c2-c1)*s - offset;

    end
end


% base = cos(mD) - cos(mD)*c + sin(mD)*s;
% 
% R = [cos(mD) sin(mD); -sin(mD) cos(mD)];
% b_partition = ones(partition, 1)*c;
% for p = 1:(partition)
%     bt = subs(base, [c; s], (R^(p-1))*[c; s]);
%     b_partition(p) = bt;
% end
% 
% X_partition = b_partition >= 0;



end

