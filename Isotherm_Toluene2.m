clear all
close all
clc
% Toluene's critical Temperature and pressure
Tc = 591.79;
Pc = 41.09;

% b and a for vdw EOS
b = (R*Tc)/(8*Pc);
a = (27*R^2*Tc^2)/(64*Pc);
% Declaring T as needed
R=0.0821;
Psat=37.8;
T=580;
%T(i)=[550 560 570 580 590 591.7];
%for i=1:1:numel(T(i));
    %(Psat*v(i)^3)-((Psat*b)+(R*T(i)))*v(i)^2+(a*v(i))+a*b=0;
    a0=-((a*b)/Psat);
    a1=a/Psat;
    a2=-((Psat*b+(R*T))/Psat);
    q=((3*(a1))-((a2)^2))/9;
    r=((a2)*((9*(a1))-(2*(a2)^2))-27*(a0))/54;
    Q= acos(r/(sqrt(-q^3)));
        
    v1= 2*sqrt(-q)*(cos(Q/3))-a2/3
    v2=-2*sqrt(-q)*(cos((Q-pi)/3))-a2/3
    v3=-2*sqrt(-q)*(cos((Q+pi)/3))-a2/3