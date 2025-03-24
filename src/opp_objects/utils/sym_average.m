function v_out= sym_average(v_in, x_trig, sym)
%SYM_AVERAGE Summary of this function goes here
%   Detailed explanation goes here

v_out = v_in;
if sym == 1
    v_out = v_out + subs(v_out, x_trig, -x_trig);                
end
if sym==2
    v_out = v_out + subs(v_out, x_trig, diag(-1, 1)*x_trig);
end

end

