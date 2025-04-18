% load('experiment_N_5_sweep.mat');



% for ordi = 1:3
    % fname = sprintf('experiment_N_5_sweep_%d.mat', ordi);

    out_std_1 = out_std(1:end, 1:end, 1);
    out_resolve_1 = out_resolve(:, :, 1);
    objective_std_1 = objective_std(:, :, 1);
    objective_resolve_1 = objective_resolve(:, :, 1);
    save('experiments/experiment_N_5_sweep_1.mat', 'k_list', 'mod_list', ...
        'out_std_1', "out_resolve_1", "objective_resolve_1", "objective_std_1");



    out_std_2 = (out_std(:, :, 2));
    out_resolve_2 = (out_resolve(:, :, 2));
    objective_std_2 = objective_std(:, :, 2);
    objective_resolve_2 = objective_resolve(:, :, 2);
    save('experiments/experiment_N_5_sweep_2.mat', 'k_list', 'mod_list', ...
        'out_std_2', "out_resolve_2", "objective_resolve_2", "objective_std_2");



    out_std_3 = (out_std(:, :, 3));
    out_resolve_3 = (out_resolve(:, :, 3));
    objective_std_3 = objective_std(:, :, 3);
    objective_resolve_3 = objective_resolve(:, :, 3);
    save('experiments/experiment_N_5_sweep_3_std.mat', 'k_list', 'mod_list', ...
        'out_std_3', "objective_std_3");

    save('experiments/experiment_N_5_sweep_3_resolve.mat', 'k_list', 'mod_list', ...
        "out_resolve_3", "objective_resolve_3");
% end
% save()