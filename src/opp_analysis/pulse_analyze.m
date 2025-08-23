function [out] = pulse_analyze(u, alpha, sym, N, tau, I0, uext)
%Find the energy of the pulse function.
%Inputs:
%u:     voltage levels
%alpha: switching angles
%sym:   symmetry
%N:     number of interpolation points (for visualization)
%tau:   load ratio R/L
%I0:    initial current
%uext:  externally applied sinusoidal voltage signal

%output:
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

%% compute the current of the pulse pattern
% I_step = uf.*da;
% I0_s= cumsum([0; I_step]);
% if nargin < 4
%     % mean_I = mean(I0_val);
%     I_s = I0_s - max(I0_s)/2;
% else
%     I_s = I0_s + I0;
% end

Aext = abs(uext);
phiext = angle(uext);
% tau = 0;
% uext = @(a) -Aext/(tau^2+1) * (tau*sin(a + phiext) - cos(a + phiext));
uphase = @(a) Aext/sqrt(tau^2+1) * (cos(a + phiext + atan(tau)));

%initial current for the pulse response only
I0 = I0 - uphase(0);

I_s = I0*ones(size(ah));

alpha_val = linspace(0, 2*pi, N);
I_val = zeros(1, N);
I_val(1) = I0;

Na = (length(I_s)-1);

I_s = I0*ones(size(ah));

%compute the current from the pulse pattern
for i = 1:Na
    ucurr = uf(i);
    Iprev = I_s(i);
    dt = da(i);
        a_range = (alpha_val >= ah(i)) & (alpha_val <= ah(i+1));
    dt_range = alpha_val(a_range) - ah(i);
    if tau==0
        I_s(i+1) = Iprev + ucurr*dt;
        I_val(a_range) = Iprev +  ucurr*(dt_range);
        % I0_s= cumsum([0; I_step]);
    else
        I_s(i+1) = ucurr*(1-exp(-tau*dt))/tau + Iprev*exp(-tau*dt);        
        I_val(a_range) = ucurr*(1-exp(-tau*dt_range))/tau + Iprev*exp(-tau*dt_range);
    end
end

% I_s = I_val;
%compute the current of the external voltage


I_s_ext = uphase(ah);
I_val_ext = uphase(alpha_val);

I_s_orig = I_s;
I_val_orig = I_val;

I_s = I_s + I_s_ext;
I_val = I_val + I_val_ext;

%DEBUG
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


%compute harmonics
[na, nb] = pulse_harmonics(50, uf', af');
out.harmonics.a = na;
out.harmonics.b = nb;

%% current energy of the pulse pattern (pure inductance)   

E_s = zeros(Na, 1);
for i = 1:Na
    ucurr = uf(i);
    Iprev = I_s_orig(i);
    dt = da(i);

    if tau ==0
        %reactance
        if ucurr == 0
            E_s(i) = Iprev.^2 * (dt);
        else
            pt_end = (Iprev + ucurr*(dt))^3;
            pt_start = (Iprev)^3;
            E_s(i) = (pt_end - pt_start)/(3*ucurr);
        end

    else
        %resistance
        E0 = (ucurr-tau*Iprev)*(3*ucurr + tau*Iprev);
        E_denom = 2*tau^3;
        E_num_1 = 2*ucurr^2*tau*dt;
        E_num_2 = exp(-2*tau*dt)*(ucurr - tau*Iprev) * ...
            (tau*Iprev + ucurr*(4*exp(tau*dt)-1));
        
        E_s(i) = (E_num_1 + E_num_2 - E0)/E_denom;
    end
end

%% energy of the mix Ip*Ie
    %integrate cos(t+z)*(u*(1-e^(-k*(t)))/k + c*e^(-k*(t))) dt
    E_s_mix = zeros(Na, 1);
    for i = 1:Na
        ucurr = uf(i);
        Iprev = I_s_orig(i);
        dt = da(i);
        acurr = ah(i);

        gain = (Aext/sqrt(tau^2+1));
       
        if tau == 0
            %pure reactance
            %integrate (q+t*u)*cos(t+z) dt 
            Esin = Iprev + ucurr*(dt);
            Ecos = ucurr;
            E_dt = (Esin*sin(dt+acurr + phiext)...
               + Ecos*cos(dt+acurr+phiext));


            Esin0 = (Iprev);
            Ecos0 = Ecos;
            E0 = (Esin0*sin(phiext+acurr) +...
                Ecos0*cos(phiext+acurr));
            
            E_s_mix(i) =(E_dt - E0)*(gain);
            E_s_curr = E_s_mix(i);
            %test 
            a_range = (alpha_val >= ah(i)) & (alpha_val <= ah(i+1));
            e_curr = trapz(alpha_val(a_range), I_val_ext(a_range).*I_val_orig(a_range));

            e_diff_curr = e_curr - E_s_curr;

            % disp(e_diff_curr)
        else
            %integrate cos(t+z)*(u*(1-e^(-k*(t)))/k + c*e^(-k*(t))) dt
            z = atan(tau) + phiext;
            E_denom = tau^3+tau;
            
            
            Esin = Iprev*tau + ...
            ucurr*((tau^2+1)*exp(tau*dt)-1);
            Ecos = -tau*(Iprev*tau - ucurr);
            E_dt =  (exp(-tau*dt))*(Esin*sin(dt+z+acurr)...
            + Ecos*cos(dt+z+acurr));
            
            
            Esin0 = (Iprev*tau + ucurr*(tau^2));
            Ecos0 = Ecos;
            E0 = (sin(z+acurr)*Esin0 +...
            Ecos0*cos(z+acurr));
            
            
            E_s_mix(i) =(E_dt - E0)*(gain/E_denom);

        end
    end
    
    % disp(E_s_mix')
    energy_pulse = sum(E_s);    
    energy_mix = sum(E_s_mix);
    
    %energy of only the trigonometric external component
    energy_ext = pi*Aext^2/(1+tau^2);
    
energy = energy_pulse + energy_ext + 2*energy_mix;

%compare the result against numerical integration
energy_trapz = trapz(alpha_val, I_val.^2);

ediff  = energy - energy_trapz;


    out.energy = energy;

out.alpha = ah;
out.I = I_s;
out.I_val = I_val;
out.alpha_val = alpha_val;
out.u = uf;


end



