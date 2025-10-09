mset clear
yalmip('clear')

%sweep on pulse length and order
%keep 

RESOLVE = false;

opts = opp_options;
opts.L = [-1, -0.5, 0, 0.5, 1];
opts.harmonics = opp_harmonics();
opts.partition = 2;
opts.TIME_INDEP = true;
opts.early_stop = 0;
opts.null_objective = false;
opts.Symmetry = 1;
opts.unipolar = 1;
opts.quarter_match = 0;

modulation_list = linspace(0.05, 1.25);
order = 4;

% modulation = 0.8;
kappa = 1;
opts.Z_load = 1.0j + kappa/(2*pi*opts.f0);
opts.verbose = 0;
modulation = 1;
opts.harmonics.index_sin= [1];
opts.harmonics.bound_sin = [modulation, modulation];

% opts.harmonics.index_sin= [1;  3];
% opts.harmonics.bound_sin = [modulation, modulation; -0.1, 0.1];


%% iterate over the orders
modulationlist = 0.05:0.05:1.25;
klist = 8:4:40;

% modulationlist = 0.8;
% klist = 12;

% kappalist = [0, 0.1, 0.5, 1, 2, 5, Inf];
% [mm, kk, kaka] = meshgrid(mlist, klist, kappalist);

% Nm = length(mlist);
Nk = length(klist);
Nmod = length(modulationlist);

%% now run the experiment


result_std = struct;
result_std.out = cell(Nk, Nmod);
result_std.tdd_lower= NaN*ones(Nk, Nmod);
result_std.solver_time  = NaN*ones(Nk, Nmod);
result_std.preprocess_time  = NaN*ones(Nk, Nmod);

result_resolve = struct;
result_resolve.out = cell(Nk, Nmod);
result_resolve.tdd_lower= NaN*ones(Nk, Nmod);
result_resolve.solver_time  = NaN*ones(Nk, Nmod);
result_resolve.preprocess_time  = NaN*ones(Nk, Nmod);

for jj = 1:Nk
    for i = 1:Nmod    
    
% for jj = 4:4
    % opts
    opts.k = klist(jj);
    opts.harmonics.bound_sin(1, :) = [1, 1] * modulationlist(i); 
    MG = opp_manager(opts);
    % for i = 2:2
        
        
        
        sol = MG.run(order);    
        disp(sol)    
        
        %% diagnose the solution
        if sol.status==0
            ms = MG.mass_summary();
            pattern_rec = MG.recover_pattern();
            out = MG.recover(sol);
            out_polish = opp_polish_RL(out);
            result_std.out{jj, i} = out;
            result_std.out_polish{jj, i} = out_polish;
        
            % harm_valid = out.pattern.harm_valid;
            result_std.tdd_lower(jj, i) = out.tdd_lower;
            result_std.solver_time(jj, i) = out.sol.solver_time;
            result_std.preprocess_time(jj, i) = out.sol.preprocess_time;
            result_std.energy(jj, i) = out.energy_lower;
    
            fprintf('unipolar k=%d, modulation %d: tdd>=%0.4e', klist(jj), modulationlist(i), out.tdd_lower)
            %solve again.
            if RESOLVE
                opts2 = MG.opts;
                opts2.allowed_levels = out.pattern.levels;
                MG2 = opp_manager(opts2);
                sol2 = MG2.run(order);
                if sol2.status == 0
                    out2 = MG2.recover(sol2);
                    result_resolve.out{jj, i} = out2;
        
                    % harm_valid = out.pattern.harm_valid;
                    result_resolve.tdd_lower(jj, i) = out2.tdd_lower;
                    result_resolve.solver_time(jj, i) = out2.sol.solver_time;
                    result_resolve.preprocess_time(jj, i) = out2.sol.preprocess_time;
                    fprintf('restricted order %d: energy >= %0.4e, tdd>=%0.4e \n', order, out2.energy, out2.tdd_lower)
                    save('order_sweep_5_mod.mat','result_std', 'result_resolve')
                end
            else
                save('order_sweep_5_load_mod.mat','result_std', 'modulationlist', 'klist', 'opts')
            end
        end
    end
        
    
    
end
