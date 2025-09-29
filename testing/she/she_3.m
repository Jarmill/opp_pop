clc;
% clear all;
% A=input('The Modulation index, M=');
%https://www.mathworks.com/matlabcentral/fileexchange/74369-newton-raphson-algorithm-for-selective-harmonic-elimination
%credit  Alla Eddine Toubal Maamar
A = 1.05;
M=0:0.01:A;
for ii=1:length(M)
    Theta1=35*pi/180;
    Theta2=55*pi/180;
    Theta3=80*pi/180;
for i=1:100
    T=[M(ii)*pi/4 0 0]';
    F=[cos(Theta1)-cos(Theta2)+cos(Theta3);
       cos(3*Theta1)-cos(3*Theta2)+cos(3*Theta3);
       cos(5*Theta1)-cos(5*Theta2)+cos(5*Theta3)];
    dF=[-sin(Theta1) +sin(Theta2) -sin(Theta3);
        -3*sin(3*Theta1) +3*sin(3*Theta2) -3*sin(3*Theta3);
        -5*sin(5*Theta1) +5*sin(5*Theta2) -5*sin(5*Theta3)];
    dTheta=(inv(dF))*(T-F);
    i
    Theta=[Theta1;Theta2;Theta3]*180/pi
    F
    dTheta*180/pi
Theta1=Theta1+dTheta(1);
Theta2=Theta2+dTheta(2);
Theta3=Theta3+dTheta(3);
        if dTheta>-1e-15 & dTheta<1e-15
        break;
end
end
Theta11(ii)=Theta1*180/pi;
Theta22(ii)=Theta2*180/pi;
Theta33(ii)=Theta3*180/pi;
end
plot(M,Theta11,'b','LineWidth',2);
hold on;
plot(M,Theta22,'r','LineWidth',2);
hold on;
plot(M,Theta33,'k','LineWidth',2);
ylim([0 90]);grid on;
xlabel('Modulation index M');
ylabel('Switching Angles(°)');
title('Switching Angles with N-R Algorithm');