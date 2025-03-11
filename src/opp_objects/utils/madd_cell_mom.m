function [cell_out] = madd_cell_mom(cell1,cell2, C)
%MADD_CELL_MOM multiply and add cells of @mom data (moments)
% cell_out = cell1 + cell2*C

cell_out = cell2;
[N, P] = size(cell2);
for n = 1:N
    for p = 1:P
        if isnumeric(cell1)
            cell_out{n, p} = C*cell2{n, p};
        else        
            cell_out{n, p} = cell1{n, p} + C*cell2{n, p};
        end
    end
end

end

