clear all
close all
clc
y(1)=1;
t(1)=0;

h=0.5;
for n= 1:25;
    dydt(n)=-y(n)^2;
    y(n+1)=y(n)+(h*dydt(n));
    t(n+1)=t(n)+h;
end
plot(t,y)
hold on

z(1)=1;
t1(1)=1;
h=1;

for n=1:25
k1=-z(n)^2;
k2=-(z(n)+(h/2*k1))^2;
k3=-(z(n)+(h/2*k2))^2;
k4=-(z(n)+(h*k3))^2;
z(n+1)=z(n)+(h/6*(k1+2*k2+2*k3+k4));
t1(n+1)=t1(n)+h;
end
plot(t1,z)

%analytical solution
t2=1:100;
x=1./(1+t2);
plot(t2,x)


legend('euler','RK4','analytical')