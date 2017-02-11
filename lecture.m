clc
clear all
close all

y(1)=1;
t(1)=0;
dydt=-y(1)^2;
h=1;
y(2)=y(1)+(h*dydt);
dydt=-y(2)^2;
y(3)=y(2)+(h*dydt);

n=[1:25]
for n=1:25
    dydt(n)=-y(n)^2;
    y(n+1)=y(n)+(h*dydt(n));
end