load('order_sweep_5_load_mod.mat');
load('order_sweep_5_load_mod_polish.mat')
order = 3;

% modulation = 0.8;
kappa = 1;

Theta = opts.f0*opts.Ts*2*pi;


modulationlist = 0.05:0.05:1.25;
klist = 8:4:40;
% klist = 8;
Nk = length(klist);
Nmod = length(modulationlist);

% polish_out = cell(Nk, Nmod);

%% main loop
% for jj = 5:Nk
% for jj = Nk:Nk
% for jj = 5:5;
    % i = 2;
    % i = 1;
% Theta = pi/100;
for jj = 1:Nk
    for i = 1:Nmod
    % parfor i = 1:Nmod  
        af = polish_out{jj, i}.warm.alpha_q;
        acirc = [af; 2*pi + af(1)];
        dacirc = diff(acirc);
        con_order = dacirc >= Theta;

        if any(~con_order)
            outcurr = result_std.out{jj, i};
            polish_out{jj, i} = opp_polish_RL(outcurr);
    
            sprintf('polish: k=%d, M=%0.2f', klist(jj), modulationlist(i));
        end
    % end
        save('order_sweep_5_load_mod_polish.mat', 'polish_out');
    end
end


%% analyze
en_polish = cellfun(@(p) p.warm.objective, polish_out);
% en_out = cellfun(@(p) p.warm.objective, polish_out);
% en_gap = en_polish - result_std.energy(1, :)
en_gap = en_polish - result_std.energy;