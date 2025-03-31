
%% sweep over standard
load('experiment_N_5_sweep_3_std.mat');

polish_std_4 = cell(11, 1);
for mi = 1:11
    for ki = 1:length(1)
        osc = out_std_3{mi, ki};
        % osc = out_resolve_3{mi, ki};
        a0 = pi/2-0.1;
        
        if ~isempty(osc)
            osc.pattern.alpha = [a0; pi-a0; pi+a0; 2*pi-a0];
            polish_std_4{mi, ki} = opp_polish_qw(osc);           
        end
        disp([mi, ki])
    end
end