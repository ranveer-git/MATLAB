clear all
clc
close all

%initial concentrations
Cao=1;
Cbo=1;
Cco=0;
Cdo=0;
Ceo=0;
Cfo=0;

%invarient combination of concentrations
c1=Cao+Cbo+Cco+Cdo+Ceo+Cfo;  
k12=0.1;k21=0.2;k13=0.05;k31=0.1;k23=0.1;k32=0.1;
k34=0.1;k43=0.2;k45=0.1;k54=0.05;k46=0.2;
k64=0.1;k56=0.1;k65=0.1;

A=[-k21-k31 k12 k13 0 0 0;
    k21 -k12-k32 k23 0 0 0;
    k31 k32 -k13-k23-k43 k34 0 0;
    0 0 k43 -k34-k54-k64 k45 k46;
    0 0 0 k54 -k45-k65 k56;
    0 0 0 k64 k65 -k46-k56];

AP=[A;1 1 1 1 1 1];
B1=[0;0;0;0;0;0;2];
C0=[1;1;0;0;0;0];

Ceq=AP\B1

cond(A);
det(A);

[V,D]=eig(A);
V(:,[3,4])=V(:,[4,3]);
V(:,[4,6])=V(:,[6,4]);
V
D(:,[3,4])=D(:,[4,3]);
D(:,[4,6])=D(:,[6,4]);
D

At=A.';
[U,B]=eig(At);
B(:,[3,4])=B(:,[4,3]);
B(:,[4,6])=B(:,[6,4]);
B
U(:,[3,4])=U(:,[4,3]);
U(:,[4,6])=U(:,[6,4]);
U
Ut=U.'
DiagMatrix=Ut*V
C=C0-Ceq
for i=1:6
    b(i)=Ut(i,:)*C/(Ut(i,:)*V(:,i));
end
b
for t=0:0.01:0.5
a2=b(2)*V(:,2)*exp(D(2,2)*t);
a3=b(3)*V(:,3)*exp(D(4,3)*t);
a4=b(4)*V(:,4)*exp(D(6,4)*t);
a5=b(5)*V(:,5)*exp(D(5,5)*t);
a6=b(6)*V(:,6)*exp(D(3,6)*t);
a=Ceq+a2+a3+a4+a5+a6;
plot(t,a,'*')
hold on
end
