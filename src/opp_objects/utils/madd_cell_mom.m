function [cell_out] = madd_cell_mom(cell1,cell2, C)
%MADD_CELL_MOM multiply and add cells of @mom data (moments)
% cell_out = cell1 + cell2*C

cell_out = cell1;
[N, P] = size(cell1);
for n = 1:N
    for p = 1:P
        cell1{n, p} = cell1{n, p} + C*cell2{n, P};
    end
end

end

