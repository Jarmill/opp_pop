% get the signal
%original signal
SYM = 2;

if SYM==0
    af = 2*pi*[0.1; 0.15; 0.3; 0.4; 0.6; 0.8; 0.85; 0.9];
    uf = [0; 1; 0; 1; 0; -1; 0; 1; 0];
elseif SYM==1
    ah = [0.5; 0.55;     0.85;     1.2;     1.3;     1.5];
    uvh = [0 1 0 -1 0 1 0]';
    af = [ah; pi+ah];
    uf = [uvh(1:end-1); -uvh];
elseif SYM == 2
    aq= [0.5764188652596557;
     1.4940129440240595;
     1.5098427393429221;
     1.5256472497294145;
     1.562892065118792];
      uq= [0 1 0 1 0 1]';
    af = [aq; pi - (aq(end:-1:1)); pi + aq; 2*pi - aq(end:-1:1)];
    uf = [uq; uq(end-1:-1:2); -uq; -uq(end-1:-1:1)];
end

 % 
 % 
 %    aa_all = [aa_quad; pi - (aa_quad(end:-1:1)); pi + aa_quad; 2*pi - aa_quad(end:-1:1)];
 %    u_all = [u_quad; u_quad(end:-1:1); -u_quad; -u_quad(end:-1:1)];
 % 
 %    af = [0; aa_all; 2*pi];
 %    % uf = [ u_all]
 %    uf = u_all;

%voltage signal
a_full = [0; kron(af, [1; 1]); 2*pi];
u_full = kron(uf, [1; 1]);

%pad the signal
ah = [0; af; 2*pi];
da = diff(ah);

energy_v = sum((uf.^2) .* da)/(2*pi);

du = diff(uf);

%compute the current
% I0 = 0;
I0 = 0;
I_step = uf.*da;
I0_val = cumsum([0; I_step]);
mean_I = mean(I0_val);
I_val = I0_val - mean_I;

%linear interpolation of current
theta = linspace(0, 2*pi, 1000);
I_eval = interp1(ah,I_val,theta);

energy_I_trap = trapz(theta, I_eval.^2)/(2*pi);

% ar1 = ah(2:end) - ah(1:end-1);
% ar2 = ah(2:end).^2 - ah(1:end-1).^2;
% ar3 = ah(2:end).^3 - ah(1:end-1).^3;

% c = I_val(1:end-1);
% k = uf;
energy_raw = 0;
for i = 1:length(da)
    slope = uf(i);
    offset = I_val(i);
    prev = ah(i);

    if slope == 0
        energy_curr = offset.^2 * (ah(i+1)-ah(i));
    else
        pt_end = (offset + slope*(ah(i+1)-prev))^3/(3*slope);
        pt_start = (offset)^3/(3*slope);
        energy_curr = pt_end - pt_start;
    end
    energy_raw = energy_raw + energy_curr/(2*pi);
end

% energy_raw = (ar1'* (k).^2 + ar2' * (k .* c) + ar3' * (k.^2)./3)/(2*pi);
% energy_fund = cos_coeff[2]^2 + sin_coeff[1]^2
% energy = energy_raw - energy_fund



%% plot the signals
c = linspecer(3);
figure(1)
clf
hold on
plot(a_full, u_full, 'linewidth', 3, 'color', c(1, :))
plot(theta, sin(theta), 'k', 'linewidth', 3)
xlabel('\theta')
ylabel('V(\theta)')

figure(3)
hold on
xlabel('\theta')
ylabel('I(\theta)')
plot(ah, I_val, 'linewidth', 3, 'color', c(2, :));
plot(theta, -cos(theta), 'k', 'linewidth', 3)

figure(2)
clf
plot(theta, I_eval.^2);


% u0 = uf(1);
% du = diff(uf);