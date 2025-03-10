mset clear

opts = opp_options;
opts.L = [-1, 0, 1];
opts.harmonics = opp_harmonics();
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
% lsupp_base.X = [sum(x.^2) <= 1];
lsupp_base.X = [x(1)^2 + x(2)^2 == 1; x(3)-x(3)^2 >= 0];

%% test a location
% m = 0;
% p = 3;
% n=1;
% id = '1_2_3';
% cell_info = struct('mode', m, 'partition', p, 'level', n, 'id', id);
% f = [-2*pi*x(2); 2*pi*x(1); 1];
% objective = 1;
% oloc = opp_location(lsupp_base, f, objective, cell_info);

%% test a mode
% Delta = opts.f0*opts.Ts;
% arc_curr = support_arc(m, x, Delta);
% lsupp_base.X = [lsupp_base.X; (arc_curr>=0)];
% mode_test = opp_mode(m, lsupp_base, opts);
% 
% m0= mode_test.initial_mass();
% 
% % objective = 1;

%% test a manager
MG = opp_manager(opts);

%the manager can be created
%do the constraints work though?