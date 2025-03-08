opts = opp_options;
opts.solver

mset clear
mpol('t');
mpol('x', 4, 1); %(c, s, phi, l)
vars = struct('t', t, 'x', x);

lsupp1 = loc_support(vars);
lsupp1.TIME_INDEP = 1;
lsupp1.Tmax = 1;

lsupp2 = lsupp1;
lsupp2.X = [40, 50]