clear all
close all
clc

A=[1 2 -1; 1 0 1; 4 -4 5];
[V,E]=eig(A)

B=[-1;0;0];

C=B\V