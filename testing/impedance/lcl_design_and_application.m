%% CALCULA FILTRO LCL
% based on: REZNIK et al.: LCL FILTER DESIGN AND PERFORMANCE ANALYSIS FOR GRID-INTERCONNECTED SYSTEMS
% Rafael dos Santos
%https://www.researchgate.net/profile/Rafael-Dos-Santos-5

%% HOW TO USE THIS SCRIPT
%1) Insert the desired design variables at the input section below
%2) Choose delta or wye connection of LCL filter by defining delta=1,0
%3) At the end, a bodeplot of the LCL filter will be shown, and the
%filter parameter values will be shown on prompt window, where: L1 is the
%inverter side inductor, L2 is the grid side inductor, Cf is the filter
%capacitor. The equivalent series resistance of each element is defined by
%rL1, rL2 and rCf respectively.

%PS: if the filter ressonance frequency is below fg/10 or above fsw/5,
%where fg and fsw are grid and switching frequency, respectively, an error
%message will be prompted. In this case, a modification on the design
%parameters is needed (e.g choosing different values for ka)

clc;
clear all;
% close all;

% k = 

%% 1-INPUT

   %nominal values
   % En=127*sqrt(3); %line-line voltage in RMS
   % Pn=1e3; %nominal power (W)
   % VDC=265; %DC-link voltage (V)
   % fg=50; %grid frequency (Hz)
   % fsw=18e3; %switching frequency (Hz)   
   % x=0.05; %maximum power factor variation seen by grid
   % ka=0.1; %attenuation factor - the capacity of attenuating harmonic content
   % deltaI1=10/100; %ripple current - that is 10% of maximum current
   % ESR=15/100; %inductor ESR impedance percentage
   % delta=0; %1-for wye, 0-for delta connection
   
   En=127*sqrt(3); %line-line voltage in RMS
   Pn=1e3; %nominal power (W)
   VDC=265; %DC-link voltage (V)
   fg=50; %grid frequency (Hz)
   % fsw=18e3; %switching frequency (Hz)   
   % fsw = fg*20;
   fsw = 18000;
   x=0.05; %maximum power factor variation seen by grid
   % x=0.04; %maximum power factor variation seen by grid
   ka=0.1; %attenuation factor - the capacity of attenuating harmonic content
   % deltaI1=10/100; %ripple current - that is 10% of maximum current
   deltaI1 = 5/100;
   ESR=15/100; %inductor ESR impedance percentage  
   delta=1; %1-for wye, 0-for delta connection
   
%% 2- BASE VALUES
   Zb=En^2/Pn; %base impedance
   Cb=1/(2*pi*fg*Zb); %base capacitance
   
%% 3 - INVERTER SIDE INDUCTOR
    Vph=En/sqrt(3); %phase voltage
    Imax=Pn*sqrt(2)/(3*Vph); %maximum inductor current
    L1=VDC/(6*fsw*deltaI1*Imax); %inverter side inductor
    rL1=ESR*(2*pi*fg*L1); %valor do ESR do Lf
    
%% 4 - CAPACITOR CF
    Cf=x*Cb;

%% 5 - INDUCTOR L2
    L2=(sqrt(1/(ka^2))+1)/((deltaI1/x)*Cb*fsw^2);
    rL2=ESR*(2*pi*fg*L2); %valor do ESR do Lf
    
%% 6 - RESSONANCE FREQUENCY
    fres=sqrt((L1+L2)/(L1*L2*Cf))/(2*pi);

    if(fres<10*fg)
        msg = 'Ressonance frequency is lower than 10x grid frequency. Please increase this ressonance frequency.';
        error(msg)
    end
    
    if(fres>0.5*fsw)
        msg = 'Ressonance frequency is higher than 0.5x switching frequency. Please reduce this ressonance frequency.';
        error(msg)
    end
            
%% 7 - RESISTOR RF
    Rf=1/(3*2*pi*fres*Cf);

    if delta==1
        Rf=3*Rf;
        Cf=Cf/3;
    end
    
%% 8 - OUTPUT
format short e;
s=tf('s');
HLCL.num=Cf*Rf*s+1;
HLCL.den=L1*Cf*L2*s^3+Cf*(L1+L2)*Rf*s^2+(L1+L2)*s;
HLCL.tf=HLCL.num/HLCL.den;

%% sinusoidal input
% NP = 6;
NP = 1;
NT = 10000;
t = linspace(0, NP/fg, NT)';
u = sin(2*pi*fg*t);
[y, to, xo] = lsim(HLCL.tf, u, t);


%% pulse input
aq = [0.423323100643366	0.646289267844426	0.806615717344141	1.55083012024646	]';
% aq = [0.720186460835395,1.447102394507780]';
uq = [0 1 0  1 0]';
% uq = [0 1 0]';

[outq] = pulse_current_voltage(uq, aq, 2);
outq.uu = kron(outq.u, [1; 1]);    
outq.aa = [0; kron(outq.alpha(2:end-1), [1; 1]); 2*pi];

up = pulse_func(2*pi*fg*t, outq.u, outq.alpha(2:end-1));
[yp, top, xop] = lsim(HLCL.tf, up, t);


%% plot
figure(2)
clf
tiledlayout(3, 1);
nexttile
hold on
c = linspecer(2);
plot(t, u, 'color', c(1, :), 'linewidth', 3);
plot(t, y, 'color', c(2, :), 'linewidth', 3);


nexttile
hold on
c = linspecer(2);
plot(t, up, 'color', c(1, :), 'linewidth', 3);
plot(t, yp, 'color', c(2, :), 'linewidth', 3);

nexttile
hold on
plot(t, up, 'color', c(1, :), 'linewidth', 3);
tp0 = reshape((0:(NP-1))*2*pi + outq.alpha, [], 1)/(2*pi*fg);
yp0= kron(ones(NP, 1), outq.current) - outq.current(1);
plot(tp0, yp0, ...
    'color', c(2, :), 'linewidth', 3);

yi = interp1(tp0,yp0,t);

figure(3)
% Mf = 
Mp = max(outq.current)*2;
Ma = max(yp);
ysi = yi/Mp;
ysp = yp/Ma;
plot(t, ysp - ysi);
% options = bodeoptions;
% options.FreqUnits = 'Hz'; % or 'rad/second', 'rpm', etc.
% options.Title.FontSize=14;
% options.Title.FontWeight='bold';
% bode(HLCL.tf,options)
% disp("Sucess!- LCL Filter parameters")
% grid on;
% L1
% rL1
% L2
% rL2
% Cf
% Rf
% fres
% 
    
