load('experiment_N_5_sweep.mat');



% for ordi = 1:3
    % fname = sprintf('experiment_N_5_sweep_%d.mat', ordi);

    out_std_1 = out_std{:, :, 1};
    out_resolve_1 = out_resolve{:, :, 1};
    objective_std_1 = objective_std(:, :, 1);
    objective_resolve_1 = objective_resolve(:, :, 1);
    save('experiment_N_5_sweep_1.mat', 'k_list', 'mod_list', ...
        'out_std_1', "out_resolve_1", "objective_resolve_1", "objective_std_1");



    out_std_2 = (out_std{:, :, 2});
    out_resolve_2 = (out_resolve{:, :, 2});
    objective_std_2 = objective_std(:, :, 2);
    objective_resolve_2 = objective_resolve(:, :, 2);
    save('experiment_N_5_sweep_2.mat', 'k_list', 'mod_list', ...
        'out_std_2', "out_resolve_2", "objective_resolve_2", "objective_std_2");



    out_std_3 = (out_std{:, :, 3});
    out_resolve_3 = (out_resolve{:, :, 3});
    objective_std_3 = objective_std(:, :, 3);
    objective_resolve_3 = objective_resolve(:, :, 3);
    save('experiment_N_5_sweep_3.mat', 'k_list', 'mod_list', ...
        'out_std_3', "out_resolve_3", "objective_resolve_3", "objective_std_3");
% end
% save()