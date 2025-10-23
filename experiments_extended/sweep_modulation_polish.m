load('order_sweep_5_load_mod.mat');

order = 3;

% modulation = 0.8;
kappa = 1;


modulationlist = 0.05:0.05:1.25;
klist = 8:4:40;
% klist = 8;
Nk = length(klist);
Nmod = length(modulationlist);

polish_out = cell(Nk, Nmod);

%% main loop
% for jj = 5:Nk
% for jj = Nk:Nk
for jj = 5:5;
    % i = 2;
    i = 1;
    % parfor i = 1:Nmod  
        outcurr = result_std.out{jj, i};
        polish_out{jj, i} = opp_polish_RL(outcurr);

        fprintf('polish: k=%d, M=%0.2f', klist(jj), modulationlist(i));
    % end
    save('order_sweep_5_load_mod_polish.mat', 'polish_out');
end


%% analyze
en_polish = cellfun(@(p) p.warm.objective, polish_out);
% en_out = cellfun(@(p) p.warm.objective, polish_out);
% en_gap = en_polish - result_std.energy(1, :)
en_gap = en_polish - result_std.energy;