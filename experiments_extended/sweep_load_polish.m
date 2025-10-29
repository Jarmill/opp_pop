load('order_sweep_5_load_rev.mat');

modulation = 0.8;
kappalist = [0, 0.1, 0.5, 1, 2 5];
kappalist_fine = [0, logspace(-1, 1.5, 60)];
% orderlist = 1:6;

order = 4;

k = 6;

% klist = 8;
Nkappa = length(kappalist);
Nfine = length(kappalist_fine);

polish_out = cell(Nkappa, 1);
polish_out_fine = cell(Nfine, 1);

%% main loop
% for jj = 6:6
for jj = order:order
    % i = 1;
    for i = 1:Nfine   
        [mm, ii] = min(abs(kappalist_fine(i)-kappalist));
        outcurr = result_std.out{jj, ii};
        outcurr.opts.Z_load = kappalist_fine(i)/(2*pi*outcurr.opts.f0) + 1.0j;
        if kappalist_fine(i)==0
            polish_out_fine{i} = opp_polish_qw(outcurr);
        else
            polish_out_fine{i} = opp_polish_RL(outcurr);
        end

        

        fprintf('polish: kappa=%0.2f', kappalist_fine(i));

        if ~isempty(polish_out_fine{i}.warm)
            fprintf('tdd=%0.3e', polish_out_fine{i}.warm.tdd);
        end
        save('order_sweep_load_deg_fine.mat', 'polish_out', 'polish_out_fine');
    end
end


%% analyze
en_polish = cellfun(@(p) p.warm.objective, polish_out_fine);
% en_gap = en_polish - result_std.energy(1, :)

en_bound = cellfun(@(c) c.energy_lower, result_std.out);
en_gap = en_polish - en_bound';