mset clear
yalmip('clear')

RESOLVE = 1;

opts = opp_options;
opts.L = [-1, -0.5, 0, 0.5, 1];
opts.harmonics = opp_harmonics();
opts.partition = 1;
opts.TIME_INDEP = true;
opts.early_stop = 0;
opts.Symmetry = 2;
opts.unipolar = 1;
opts.Z_load = 1.0j;
opts.verbose = 0;
%index

% mod_list = 0.7;
% k_list = 8;
% k_list = 4*(1:10);
% mod_list = 0.05:0.05:1.1;
% k_list = [8, 12];
% mod_list = [0.7, 1];

Nmod = length(mod_list);
Nk = length(k_list);

out_std = cell(Nmod, Nk);
out_resolve = cell(Nmod, Nk);
objective_std = zeros(Nmod, Nk);
objective_resolve = zeros(Nmod, Nk);
order = 2;
% order_list = 1:3;

%% test a manager
fname = sprintf('experiments/experiment_N_%d_order_%d.mat', length(opts.L), order);
for ki = 1:Nk
    k_curr = k_list(ki);
    for mi = 1:Nmod
        M_curr = mod_list(mi);

        opts_curr = opts;
        opts_curr.k = k_curr;
        opts_curr.harmonics.index_sin= [1;  3];
        opts_curr.harmonics.bound_sin = [M_curr, M_curr; -0.01, 0.01];

        MGc = opp_manager(opts_curr);
        sol_curr = MGc.run(order);
        
        
        if sol_curr.status==0
            out_curr = MGc.recover(sol_curr);
            out_std{mi, ki} = out_curr;
            objective_std(mi, ki) = sol_curr.obj_rec;

            if RESOLVE
                opts_2 = opts_curr;                
                opts_2.allowed_levels = out_curr.pattern.levels;
                MG2 = opp_manager(opts_2);
                sol2 = MG2.run(order);
                out_2 = MG2.recover(sol2);
                objective_resolve(mi, ki) = sol2.obj_rec;

                if sol2.status==0
                    out_resolve{mi, ki} = out_2;
                else
                    objective_resolve(mi, ki) = Inf;
                end
            end
            
        else
            objective_std(mi, ki) = Inf;
        end
        save(fname, 'out_std', 'out_resolve', 'objective_std', 'objective_resolve')
    end
end

