%it is possible to impose a lebesgue constraint on (c, s)-marginals.
%so why does it fail in the OPP setting?

mset clear
mset('yalmip',true);
mpol('x', 2, 1)
mpol('y', 1, 1)
mu =meas([x; y]);

supp_con = [sum(x.^2)==1; 1-y^2>=0];

order = 4;
d = 2*order;

pw = genPowGlopti(2, d);
a = leb_sphere(pw);
vx = mmon(x, d);
mom_con = [mom(vx)==a; mom(y)==0.5];

objective = mom(y^2);


P = msdp(min(objective), mom_con, supp_con);

            sol = struct;
            tic;
            [sol.status,sol.obj_rec, ~,sol.dual_rec]= msol(P);     
            sol.solver_time = toc;

M = double(mmat(mu))