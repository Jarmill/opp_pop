%https://www.techrxiv.org/users/692127/articles/682420-optimized-pulse-patterns-with-bounded-semiconductor-losses

%Table 1
%3-level inverter: [-1, 0, 1]
%up/down: voltage step
%pos/neg: polarity of current
%on/off: diode turns on or off


dispatch = opp_power_dispatch();

% %conduction constants
% cond_gct = [0.97, 0.245e-3];
% cond_diode = [1.19, 0.395e-3];
% 
% %other constants
% Vdc = 5000;
% 
%   %other constants
%                 Vdc = 5000;
% 
% 
% on_gct = [0, 0.095e-3];
% off_gct = [0, 2.6e-3];
% on_diode = 0;
% off_diode = [15.2, 0, 0, 0]; %reverse recovery, make this fancier later
% 
% 
% 
% topology = cell(10, 1);
% %figure 4
% for i = 1:4
%     topology{i} = opp_component(i, true, cond_gct, on_gct, off_gct);
% end
% for i = 5:10
%     topology{i} = opp_component(i, false, cond_diode, on_diode, off_diode);
% end
% 
% 
% %which components turn on/off when
% 
% conduction_pos = {[7, 8], [2, 9], [1, 2]};
% conduction_neg = {[3, 4], [3, 10], [5, 6]};
% 
% 
% jump_up_pos_on = {[2],[1]};
% jump_up_pos_off = {[8],[9]};
% 
% jump_up_neg_on = {[], []};
% jump_up_neg_off = {[4], [3]};
% 
% jump_down_pos_on = {[],[]};
% jump_down_pos_off = {[2],[1]};
% 
% jump_down_neg_on = {[4], [3]};
% jump_down_neg_off = {[10], [5]};
% 
% 
