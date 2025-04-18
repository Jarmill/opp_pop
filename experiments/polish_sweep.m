
%% sweep over standard
load('experiment_N_5_sweep_3_std.mat');


polish_std = cell(length(mod_list), length(k_list));
for mi = 1:length(mod_list)
    for ki = 1:length(k_list)
        osc = out_std_3{mi, ki};
        if ~isempty(osc)
            polish_std{mi, ki} = opp_polish_qw(osc);
           
        end
        disp([mi, ki])
    end
end

    save('experiments/experiment_N_5_sweep_3_std.mat', 'k_list', 'mod_list', ...
        "out_std_3", "objective_std_3");

%% sweep over resolve
load('experiment_N_5_sweep_3_resolve.mat');


polish_resolve = cell(length(mod_list), length(k_list));
for mi = 1:length(mod_list)
    for ki = 1:length(k_list)
        osc = out_resolve_3{mi, ki};
        if ~isempty(osc)
            polish_resolve{mi, ki} = opp_polish_qw(osc);
        end
    end
end

    save('experiments/experiment_N_5_sweep_3_resolve.mat', 'k_list', 'mod_list', ...
        "out_resolve_3", "objective_resolve_3");


%% store the objectives
polish_resolve_tdd = Inf*ones(length(mod_list), length(k_list));
out_resolve_tdd = Inf*ones(length(mod_list), length(k_list));
polish_std_tdd = Inf*ones(length(mod_list), length(k_list));
out_std_tdd = Inf*ones(length(mod_list), length(k_list));
for mi = 1:length(mod_list)
    for ki = 1:length(k_list)
        if ~isempty(out_resolve_3{mi, ki})
            out_resolve_tdd(mi, ki) = out_resolve_3{mi, ki}.tdd_lower;
        end
        if ~isempty(polish_resolve{mi, ki}) && ~isempty(polish_resolve{mi, ki}.warm)
            polish_resolve_tdd(mi, ki) = polish_resolve{mi, ki}.warm.tdd;
        end

        if ~isempty(out_std_3{mi, ki})
            out_std_tdd(mi, ki) = out_std_3{mi, ki}.tdd_lower;
        end
        if ~isempty(polish_std{mi, ki}) && ~isempty(polish_std{mi, ki}.warm)
            polish_std_tdd(mi, ki) = polish_std{mi, ki}.warm.tdd;
        end
    end
end


%% analyze

polish_upper = min(polish_std_tdd, polish_resolve_tdd);
% tdd_diff = -out_std_tdd + polish_upper;
tdd_diff_std = polish_std_tdd - out_std_tdd;

save('experiments/sweep_gap.mat', 'tdd_diff_std', 'polish_std_tdd', 'out_std_tdd')