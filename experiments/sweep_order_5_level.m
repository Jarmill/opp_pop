mset clear
yalmip('clear')

RESOLVE = true;

opts = opp_options;
opts.L = [-1, -0.5, 0, 0.5, 1];
opts.harmonics = opp_harmonics();
opts.partition = 1;
opts.TIME_INDEP = true;
opts.early_stop = 0;
opts.null_objective = false;
opts.Symmetry = 2;
opts.unipolar = 1;

opts.k = 32;

modulation = 0.9;
opts.Z_load = 1.0j;
opts.verbose = 0;


opts.harmonics.index_sin= [1;  3];
opts.harmonics.bound_sin = [modulation, modulation; -0.01, 0.01];


%% iterate over the orders


orderlist = 1:2;
Norder = length(orderlist);

%% now run the experiment

MG = opp_manager(opts);
result_std = struct;
result_std.out = cell(1, Norder);
result_std.tdd_lower= NaN*ones(1, Norder);
result_std.solver_time  = NaN*ones(1, Norder);
result_std.preprocess_time  = NaN*ones(1, Norder);

result_resolve = struct;
result_resolve.out = cell(1, Norder);
result_resolve.tdd_lower= NaN*ones(1, Norder);
result_resolve.solver_time  = NaN*ones(1, Norder);
result_resolve.preprocess_time  = NaN*ones(1, Norder);

for i = 1:Norder
    order = orderlist(i);    
    sol = MG.run(order);    
    disp(sol)    
    
    %% diagnose the solution
    if sol.status==0
        ms = MG.mass_summary();
        pattern_rec = MG.recover_pattern();
        out = MG.recover(sol);
        result_std.out{i} = out;
    
        % harm_valid = out.pattern.harm_valid;
        result_std.tdd_lower(i) = out.tdd_lower;
        result_std.solver_time(i) = out.sol.solver_time;
        result_std.preprocess_time(i) = out.sol.preprocess_time;

        fprintf('unipolar order %d: tdd>=%0.4e', order, out.tdd_lower)
        %solve again.
        if RESOLVE
            opts2 = opts;
            opts2.allowed_levels = out.pattern.levels;
            MG2 = opp_manager(opts2);
            sol2 = MG2.run(order);
            if sol2.status == 0
                out2 = MG2.recover(sol2);
                result_resolve.out{i} = out2;
    
                % harm_valid = out.pattern.harm_valid;
                result_resolve.tdd_lower(i) = out2.tdd_lower;
                result_resolve.solver_time(i) = out2.sol.solver_time;
                result_resolve.preprocess_time(i) = out2.sol.preprocess_time;
                fprintf('restricted order %d: tdd>=%0.4e', order, out2.tdd_lower)
            end
        end
    end
        
    save('experiments/order_sweep_5.mat','result_std', 'result_resolve')
    
end
