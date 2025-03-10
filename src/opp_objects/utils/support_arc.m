function [arc] = support_arc(m, x, Delta)
%SUPPORT_ARC what angles on the unit disc are permissible for mode m
c = x(1);
s = x(2);

        if m==0        
            mD = 2*Delta;        
            %arc [0, 2pi-2alpha]
            arc = cos(mD) - cos(mD)*c + sin(mD)*s;
        elseif m==1
            mD = Delta;
            %arc [alpha, 2pi-alpha]
            arc= cos(mD) - cos(mD)*c;
        else
            mD = double(m)*Delta;        
            %arc [malpha, 2pi]
            arc = cos(mD) - cos(mD)*c - sin(mD)*s;
        end
end

