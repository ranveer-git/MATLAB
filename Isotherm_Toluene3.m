clear all
close all
clc
 
% Toluene's critical Temperature and pressure
Tc = 591.79;
Pc = 4109000;
% Universal Gas Constant
R = 8.314;
% b and a for vdw EOS
b = (R*Tc)/(8*Pc);
a = (27*R^2*Tc^2)/(64*Pc);
% Declaring T as needed
T = [550 560 570 580 590];
disp(a);
disp(b);
    VPoints = 0.0003:0.00001:0.0009;
    for i=1:1:numel(T)
        PT=(R*T(i)./(VPoints-b)) - (a./(VPoints.^2));
        u=-(R*T(i).*(log(VPoints-b)))+((R*T(i).*VPoints)./(VPoints-b))-((2*a)./VPoints);
        figure(1)
        plot(VPoints,PT);
        hold on
        figure(2)
        plot(u,PT);
        hold on
    end

   
    % plotting isotherm for T =580K
    %figure(1)
    %plot(VPoints,PT);
    %h=plot(VPoints,PT);
    %set(h,'color',rand(1,3),'linewidth',1);
    
    xlabel('u')
    ylabel('pressure in N/m^2')
    title('Isotherms for Toluene')