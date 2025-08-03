function [out] = pulse_current_voltage_ext(u, alpha, sym, I0, uext)
%Find the energy of the pulse function.
%
%energy = integrate x(th)^2 dth for th in [0, 2pi]

%replicate the voltages and angles
if sym == 0
    %full-wave symmetry
    uf = u;
    af = alpha;
elseif sym==1
    %half-wave symmetry
    uf = [u; -u(2:end)];
    af = [alpha; alpha+pi];
else
    %quarter-wave symmetry
    urev = u(end-1:-1:1);
    uf = [u; urev; -u(2:end); -urev];
    af= [alpha; pi - (alpha(end:-1:1)); pi + alpha; 2*pi - alpha(end:-1:1)];
end

ah = [0; af; 2*pi];
da = diff(ah);
Na = (length(ah)-1);

out = struct;
out.alpha = ah;
out.voltage = uf;

%compute the current of the pulse pattern
% I0 = 0;
I_step = uf.*da;
I0_s= cumsum([0; I_step]);
if nargin < 4
    % mean_I = mean(I0_val);
    I_s = I0_s - max(I0_s)/2;
else
    I_s = I0_s + I0;
end

% I_s = I_val;
%compute the current of the external voltage
Aext = abs(uext);
phiext = angle(uext);
kappa = 0;
% uext = @(a) -Aext/(kappa^2+1) * (kappa*sin(a + phiext) - cos(a + phiext));
uphase = @(a) Aext/sqrt(kappa^2+1) * (cos(a + phiext + atan(kappa)));

alpha_val = linspace(0, 2*pi, 10000);

I_s_ext = uphase(ah);
% I_val_ext = uext(alpha_val);
I_val_phase = uphase(alpha_val);

I_s_orig = I_s;
I_val_orig = I_val;

I_s = I_s + I_s_ext;
I_val = I_val + I_val_ext;


figure(2)
clf
hold on
plot(alpha_val, I_val_ext, 'b')
plot(alpha_val, I_val, 'k')
plot(alpha_val, I_val_orig, 'r')
plot(alpha_val, I_val_ext.*I_val_orig, 'g')


out.current = I_s;
out.u = uf;

nalpha = length(alpha);


%compute the voltage energy
energy = 0;
for i = 1:nalpha
    da = (ah(i+1) - ah(i));
    dx = u(i)^2;

    energy = energy + dx*da;        
end

out.energy_V = energy;

[na, nb] = pulse_harmonics(50, uf', af');
out.harmonics.a = na;
out.harmonics.b = nb;

%current energy of the pulse pattern (pure inductance)    
    energy_pulse = 0;
    for i = 1:length(ah)-1
        slope = uf(i);
        offset = I_s(i);
        prev = ah(i);
    
        if slope == 0
            energy_curr = offset.^2 * (ah(i+1)-ah(i));
        else
            pt_end = (offset + slope*(ah(i+1)-prev))^3/(3*slope);
            pt_start = (offset)^3/(3*slope);
            energy_curr = pt_end - pt_start;
        end
        % energy_pulse = energy_pulse + energy_curr/(2*pi);
        energy_pulse = energy_pulse + energy_curr;
    end

%energy of the mix Ip*Ie
    %integrate cos(t+z)*(u*(1-e^(-k*(t)))/k + c*e^(-k*(t))) dt
    E_s_mix = zeros(Na, 1);
    for i = 1:Na
        ucurr = uf(i);
        Iprev = I_s_orig(i);
        dt = da(i);
        acurr = ah(i);

        gain = (Aext);
       
       Esin = Iprev + ucurr*(acurr + dt);
       Ecos = ucurr;


        E_denom =1;

       E_dt = (Esin*sin(dt+acurr + phiext)...
           + Ecos*cos(dt+acurr+phiext));
        
    
        Esin0 = (Iprev* + ucurr*acurr);
        Ecos0 = Ecos;
        E0 = (sin(phiext+acurr)*Esin0 +...
            Ecos0*cos(phiext+acurr));
       
        E_s_mix(i) =(E_dt - E0)*(gain);
    end
    
    % disp(E_s_mix')
    energy_pulse = sum(E_s);
    
    energy_mix = sum(E_s_mix);
    
    %energy of only the trigonometric external component
    energy_ext = pi*Aext^2;
    


    out.energy_I = energy_pulse + energy_ext + 2*energy_mix;




end



