clc
clear all
close all
k1=1;
k2=100;
k1a=10;
k2a=1.5;
Ca0=1;
X=(Ca0-Ca)/Ca0;
t(1)=0;
Ca(1)=1;
Cu(1)=0;
Cd(1)=0;

t = 0:1:20;
Ca=0:0.1:1;
Cu=0:0.1:1;
Cd=0:0.1:1;
dCadt=-(k1*(Ca-Cd/k1a)+(k2*(Ca-Cu/k2a)));
dCddt=k1*(Ca-Cd/k1a);
dCudt=k2*(Ca-Cu/k2a);
h=1;
for n= 1:1:numel(t);
    dCadt(n)=-(k1.*(Ca(n)-Cd(n)/k1a)+(k2.*(Ca(n)-Cu(n)/k2a)));
    dCddt(n)=k1.*(Ca(n)-Cd(n)/k1a);
    dCudt(n)=k2.*(Ca(n)-Cu(n)/k2a);
    dCadt(n+1)=dCadt(n)+(h*dCadt(n));
    dCbdt(n+1)=dCbdt(n)+(h*dCbdt(n));
    dCudt(n+1)=dCudt(n)+(h*dCudt(n));
    t(n+1)=t(n)+h;
end
plot(t,dCadt)