%modified from 
%https://www.mathworks.com/matlabcentral/fileexchange/74369-newton-raphson-algorithm-for-selective-harmonic-elimination
%credit  Alla Eddine Toubal Maamar

clc;
clear all;
% A=input('The Modulation index, M=');
% A = 4/pi;
A = 1.05;
% k = 24;
k = 12;
d = k/4;
% u = [0.5 1 0.5 1 0.5 1];
% u = [1, -1, 1];
% u = [0, 0.5, 1];
% u = [1, -1, 1, -1];
% u = [0, 1, 0];
u = [0, 0.5, 1, 0.5];

% u = [0, 0.5, 1, 0.5];
du = diff(u);


M=0:0.01:A;
l = (1:d)*2 - 1;
Theta_list = zeros(length(M), d);

Theta0 = pi/180 * (linspace(20, 80, d));


% for ii=1:length(M)
for ii = 40;
    % Theta1=35*pi/180;
    % Theta2=55*pi/180;
    % Theta3=80*pi/180;
    % Theta = [Theta1, Theta2, Theta3];
    T=[M(ii)*pi/4, zeros(1, d-1)]';
    Theta = Theta0;
    for i=1:100        
        % F=[cos(Theta1)-cos(Theta2)+cos(Theta3);
        %    cos(3*Theta1)-cos(3*Theta2)+cos(3*Theta3);
        %    cos(5*Theta1)-cos(5*Theta2)+cos(5*Theta3)];
        F = (du *  cos(l' * Theta)')';
        % dF = diag(u.*l) * (sin(l' * Theta))';
        dF = -diag(l) * sin(l' * Theta) * diag(du);

        % dF=[-sin(Theta1) +sin(Theta2) -sin(Theta3);
        %     -3*sin(3*Theta1) +3*sin(3*Theta2) -3*sin(3*Theta3);
        %     -5*sin(5*Theta1) +5*sin(5*Theta2) -5*sin(5*Theta3)];
        dTheta=(dF \ (T-F));
        
        % Theta=[Theta1;Theta2;Theta3]*180/pi
        
        % dTheta*180/pi;
    % Theta1=Theta1+dTheta(1);
    % Theta2=Theta2+dTheta(2);
    % Theta3=Theta3+dTheta(3);
        Theta = Theta + dTheta';
        if norm(dTheta)>-1e-15 & norm(dTheta)<1e-15
            break;
        end
    end
    Theta = Theta * 180/pi;
    Theta_list(ii, :) = Theta;
    Theta0 = Theta;


end

%% plot the result
figure(100)
clf
plot(M,Theta_list','LineWidth',2);
% hold on;
% plot(M,Theta22,'r','LineWidth',2);
% hold on;
% plot(M,Theta33,'k','LineWidth',2);
ylim([0 90]);grid on;
xlabel('Modulation index M');
ylabel('Switching Angles(°)');
title('Switching Angles with N-R Algorithm');