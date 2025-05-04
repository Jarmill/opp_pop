mset('clear')
mpol('x', 1, 1)
mset('yalmip',true, 'verbose', true);


p = (x-2)^2;

objective = mom(p);
% mom_con = mass(x)==1;
mom_con = [mom(x^2) <= 9; mom(x^2)>= -9; mom(x) >= 3; mom(x) <= 3; mom(x)];
% supp_con = [];
% supp_con = [x^2 <= 100];
% supp_con = (x^2>=0);
% supp_con = [];

  
% P = msdp(min(objective), mom_con, supp_con);
 P = msdp(min(objective), mom_con);
sol = struct;
tic;
[sol.status,sol.obj_rec, ~,sol.dual_rec]= msol(P);     
sol.solver_time = toc;
