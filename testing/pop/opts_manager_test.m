mset clear

opts = opp_options;
opts.L = [-1, 0, 1];
opts.harmonics = opp_harmonics();
opts.partition = 3;
opts.TIME_INDEP = true;
opts.start_level = 0;
opts.early_stop = 1;
opts.k = 12;


%% test a manager
MG = opp_manager(opts);
order = 2;
d = 2*order;

%the manager can be created
%do the constraints work though?

%% step through the constraints

%support constraints can be generated
% sc = MG.supp_con();

% lc = MG.con_lebesgue_circ(d); %works
% pc = MG.con_prob_dist(); %works
rc = MG.con_return(d); 