clear all
clc
T = [250:1:300];
To = 70;
A = 16.96*(10^12);
Fao=80;
To=75;
V=(1/7.484)*500;
UA=16000;
Ta1=60;
K=16.96e12.*exp(-32400./1.987./(T+460));
Fbo=1000;
Fmo=100;
mc=1000;
thetaCp=35 + Fbo/Fao*18 + Fmo/Fao*19.5;
vo=Fao/0.923 + Fbo/3.45 +Fmo/1.54;
Ta2=T-(T-Ta1).*exp(-UA./(18.*mc));
tau=V/vo;
Q=mc*18*(Ta1-Ta2);
Xmb = (tau.*K)./(1+(tau.*K));
Xeb1 = ((thetaCp.*(T-75))+ (Q/Fao))./(36000);
Xeb2 = ((thetaCp.*(T-70))+ (Q/Fao))./(36000);
plot(T,Xmb,'r');
hold on;
plot(T,Xeb1,'g',T,Xeb2,'b');