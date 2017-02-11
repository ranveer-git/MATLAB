% CO2's critical Temperature and pressure
Tc = 369.9;
Pc = 42.0;
omega = 0.152;
% Universal Gas Constant
R = 83.14;
% b and Kpr for PR EOS
b = 0.077796*R*Tc/Pc;
Kpr = 0.37464 + 1.54226*omega - 0.26992*omega^2;
for i=40:10:90
    % molar volume
    v=0.001:1:2500;
    % temperature
    T(i) = 273.15+i;
    % reduced temperature
    Tr = T(i)/Tc;
    % a for PR EOS
    a = 0.45724*(R*Tc)^2/Pc*(1 + Kpr*(1 - sqrt(Tr)))^2; 
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
end
