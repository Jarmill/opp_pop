function [arc] = angle_range_arc(c, s, th1, th2)
%ANGLE_RANGE_CON Summary of this function goes here


        s1 = sin(th1);
        s2 = sin(th2);
        c1 = cos(th1);
        c2 = cos(th2);
    
        offset = (s2-s1)*c1 - (c2-c1)*s1;

        arc = (s2-s1)*c - (c2-c1)*s - offset;
end

