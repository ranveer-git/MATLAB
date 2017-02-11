clc
clear all
close all

y(1)=1;
t(1)=0;

h=0.5;
for n= 1:25;
    dydt(n)=-y(n)^2;
    y(n+1)=y(n)+(h*dydt(n));
end
plot(t,y)