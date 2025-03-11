function [arc] = support_arc(m, x, Delta, Symmetry)
%SUPPORT_ARC what angles on the unit disc are permissible for mode m?
c = x(1);
s = x(2);

if nargin < 3
    Symmetry = 0;
end
m = double(m);

    if Symmetry==0
    %full-wave symmetry
        if m==0        
            th1 = 0;
            th2 = 2*pi-2*Delta;
            % mD = 2*Delta;        
            %arc [0, 2pi-2Delta]
            % arc = cos(mD) - cos(mD)*c + sin(mD)*s;
            
        % elseif m==1
        %     mD = Delta;
        %     %arc [Delta, 2pi-Delta]
        %     arc= cos(mD) - cos(mD)*c;
        else
            th1 = (m-1)*Delta;
            th2 = 2*pi;
            % mD = double(m-1)*Delta;        
            %arc [(m-1)Delta, 2pi]
            % arc = cos(mD) - cos(mD)*c - sin(mD)*s;
        end
    elseif Symmetry ==1
        if m==0
            %arc [0, pi-Delta]
            th1 = 0;
            th2 = pi-Delta;
        else
            %arc [mDelta, pi]
            th1 = (m-1)*Delta;
            th2 = pi;
        end
    elseif Symmetry==2
        if m==0
            %arc [Delta/2, pi/2-Delta/2]
            th1 = Delta/2;
            th2 = pi/2-Delta/2;
        else
            %arc [0, pi/2-2Delta]
            th1 = (m-1)*Delta;
            th2 = pi/2-Delta/2;
        end
    end

    arc = angle_range_arc(c, s, th1, th2);

end

