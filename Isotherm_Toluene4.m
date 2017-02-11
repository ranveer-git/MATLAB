clc
clear all
close all

a=24.539;
b=0.149;
v=0.15:0.01:2.5;
r=0.0821;
t=[560 565 570 575 580 585 590 591 591.57];

for i=1:8
l=r*t(i);
m=v-b;
q=l./m;
n=a./(v.*v);
p=q-n;
plot(v,p,'r')
hold on
end