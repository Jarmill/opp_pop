opts = opp_options;

%test opp_location
mpol('x', 3, 1);
t = [];
vars = struct('t', t, 'x', x);
lsupp_base = loc_support(vars);
lsupp_base.vars.x = x;
lsupp_base.vars.t = t;

lsupp_base.TIME_INDEP = opts.TIME_INDEP;
lsupp_base.FREE_TERM = 0;
lsupp_base.Tmax = 1;
lsupp_base.X = [sum(x.^2) <= 1];

m = 2;
p = 3;
n=1;
id = '1_2_3';
cell_info = struct('mode', m, 'partition', p, 'level', n, 'id', id);

f = [-2*pi*x(2); 2*pi*x(1); 1];

objective = 1;

oloc = opp_location(lsupp_base, f, objective, cell_info);