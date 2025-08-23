mpol('x', 2, 1);
mu = meas(x);

p = 6;

C0 = [2; 1];
f = sum((x - C0).^p);

xstar = C0/norm(C0);
fstar = sum((xstar - C0).^2);

objective = mom(f);
% supp_con = (sum(x.^2)-1==0);
supp_con = (x(2)^2 == 1 - x(1)^2);
supp_con_2 = (-x(2)^2 + 1 - x(1)^2==0);
mom_con = [mass(mu)==1];

P = msdp(min(objective), mom_con, supp_con);
P2 = msdp(min(objective), mom_con, supp_con_2);


tic;
[sol.status,sol.obj_rec, ~,sol.dual_rec]= msol(P);     
sol.solver_time = toc;

