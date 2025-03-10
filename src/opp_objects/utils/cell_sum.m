function [cs] = cell_sum(cell_in)
%CELL_SUM sum all elements of a 2d array
%useful for momcon constructions
cs = cell_in{1}*0;

[N, P] = size(cell_in);
for n = 1:N
    for p = 1:P
        cs = cs + cell_in{n, p};
    end
end
end

