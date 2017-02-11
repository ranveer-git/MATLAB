% CO2's critical Temperature and pressure
Tc = 304.25;
Pc = 73.8;
omega = 0.225;
% Universal Gas Constant
R = 8.314;
% b and Kpr for PR EOS
b = 0.077796*R*Tc/Pc;
Kpr = 0.37464 + 1.54226*omega - 0.26992*omega^2;
% Declaring Tr as needed
Tr= [0.8];
i = 1;
% 1.0 1.1 2.0
%for i=1:1:numel(Tr)
    % molar volume
    v=1:1:2500;
    % temperature
    T(i) = Tr(i) * Tc;
    % reduced temperature
    % Tr = Tr(i);
    % a for PR EOS
    a = 0.45724*(R*Tc)^2/Pc*(1 + Kpr*(1 - sqrt(Tr(i))))^2; 
    % PR EOS
    P=R*T(i)./(v-b) - a./(v.*(v + b)+b*(v - b));
    %Pv = [Pv P'];
    % plotting isotherm for T varying from Tr [0.8,1,1.1,2.0]
    figure(2)
    h=plot(v,P);
    set(h,'color',rand(1,3),'linewidth',2);
    hold on
    axis([0 1600 -40 60])
    xlabel('Volume in cm3/mol')
    ylabel('pressure in bar')
    title('Isotherms for propane')
%end
