clc
clear
clf

% Toluene's critical Temperature and pressure
Tc = 591.79;
Pc = 4109000;
% Universal Gas Constant
R = 8.314;
% b and a for vdw EOS
b = (R*Tc)/(8*Pc);
a = (27*R^2*Tc^2)/(64*Pc);
% Declaring T as needed
T = 580;

    VPoints = 0.0003:0.00001:0.0009;
    for i=1:1:numel(VPoints)
        PT(i)=((R*T)/(VPoints(i)-b)) - (a/(VPoints(i)^2));
    end
    
    % plotting isotherm for T =580K
    figure(1)
    %h=plot(v,P);
    h=plot(VPoints,PT);
    set(h,'color',rand(1,3),'linewidth',1);
    hold on
    xlabel('Volume in L')
    ylabel('pressure in N/m^2')
    title('Isotherms for Toluene')
    