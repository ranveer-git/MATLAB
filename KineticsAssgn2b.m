clc
clear all
close all

k1 = 0.004;
k2 = 0.3;
k3 = 0.25;
Ca = 0:0.001:0.1;

Sbx = (k2/k1)*Ca.^0.5;
figure(1) 
plot (Ca, Sbx)
xlabel('Ca(mol/dm3)')
ylabel('Sbx')

hold on

Sby = k2./(k3.*Ca);
figure(2)
plot (Ca,Sby);
xlabel('Ca(mol/dm3)')
ylabel('Sby')

Sbxy=(k2.*Ca)./((k1.*Ca.^0.5)+(k3.*Ca.^2));
figure(3)
plot(Ca,Sbxy)
xlabel('Ca(mol/dm3)')
ylabel('Sbxy')