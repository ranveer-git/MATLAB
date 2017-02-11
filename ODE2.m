clc
close all
clear all
%systems  of differential equations

Cfe = 0.6;
Cinner=0.6;
KonCD= 1*10^-3;
KonCDO=100*KonCD;
KonCBM= 19*10^-3;
KonCBMO=100*KonCBM;
KoffCBM=3*10^-4;
KoffCD=3*10^-4;
Kslide=16;
Krev=1.6;
KbCB=600;
KuCB=100;
Kh=4.5;
Kb=0.07;

a1 =[-KonCBM*Cfe-KonCD*Cfe-KonCBM*Cinner,KoffCBM,KoffCBM,KoffCD,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
a2 = [KonCBM*Cinner,-KoffCBM,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
a3 = [KonCBM*Cfe,0,-KoffCBM-KonCDO,0,KoffCD,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
a4 = [KonCD*Cfe,0,0,-KoffCD-KonCBMO-Kslide,KoffCBM,Krev,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
a5 = [0,0,KonCDO,KonCBMO,-KoffCBM-KoffCD-Kslide,0,Krev,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
a6 = [0,0,0,Kslide,0,-Krev-Kslide-KonCBMO,KoffCBM,Krev,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
a7 = [0,0,0,0,Kslide,KonCBMO,-Krev-Kslide-KoffCBM,0,Krev,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
a8 = [0,0,0,0,0,Kslide,0,-Krev-Kslide-KonCBMO,KoffCBM,Krev,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
a9 = [0,0,0,0,0,0,Kslide,KonCBMO,-Krev-Kslide-KoffCBM,0,Krev,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
a10 = [0,0,0,0,0,0,0,Kslide,0,-Krev-Kslide-KonCBMO,KoffCBM,Krev,0,0,0,0,0,0,0,0,0,0,0,0,0];
a11 = [0,0,0,0,0,0,0,0,Kslide,KonCBMO,-Krev-Kslide-KoffCBM,0,Krev,0,0,0,0,0,0,0,0,0,0,0,0];
a12 = [0,0,0,0,0,0,0,0,0,Kslide,0,-Krev-Kslide-KonCBMO,KoffCBM,Krev,0,0,0,0,0,0,0,0,0,0,0];
a13 = [0,0,0,0,0,0,0,0,0,0,Kslide,KonCBMO,-Krev-Kslide-KoffCBM,0,Krev,0,0,0,0,0,0,0,0,0,0];
a14 = [0,0,0,0,0,0,0,0,0,0,0,Kslide,0,-Krev-Kslide-KonCBMO,KoffCBM,Krev,0,0,0,0,0,0,0,0,0];
a15 = [0,0,0,0,0,0,0,0,0,0,0,0,Kslide,KonCBMO,-Krev-Kslide-KoffCBM,0,Krev,0,0,0,0,0,0,0,0];
a16 = [0,0,0,0,0,0,0,0,0,0,0,0,0,Kslide,0,-Krev-Kslide-KonCBMO,KoffCBM,Krev,0,0,0,0,0,0,0];
a17 = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,Kslide,KonCBMO,-Krev-Kslide-KoffCBM,0,Krev,0,0,0,0,0,0];
a18 = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,Kslide,0,-Krev-Kslide-KonCBMO-KbCB,KoffCBM,Krev,0,0,0,KuCB,0];
a19 = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,Kslide,KonCBMO,-Krev-Kslide-KoffCBM-KbCB,0,Krev,0,0,0,KuCB];
a20 = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,Kslide,0,-Krev-Kslide-KonCBMO,KoffCBM,Krev,0,0,0];
a21 = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,Kslide,KonCBMO,-Krev-Kslide-KoffCBM,0,Krev,0,0];
a22 = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,Kslide,0,-Krev-Kh-KonCBMO,KoffCBM,Kb,0];
a23 = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,Kslide,KonCBMO,-Krev-Kh-KoffCBM,0,Kb];
a24 = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,KbCB,0,0,0,Kh,0,-KuCB-Kb-KonCBMO,KoffCBM];
a25 = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,KbCB,0,0,0,Kh,KonCBMO,-KuCB-Kb-KoffCBM];
A = [a1;a2;a3;a4;a5;a6;a7;a8;a9;a10;a11;a12;a13;a14;a15;a16;a17;a18;a19;a20;a21;a22;a23;a24;a25]

[eig_vec,eig_value]=eig(A)

initialValue= [10^-3 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]';

%C = linsolve(V,B)
%C = linsolve(EV,B);
C = eig_vec\initialValue
D=diag(eig_value);

syms x t
x= vpa(eig_vec*(C.*exp(D*t)),4)

for n=0:10:5000
    x2 = subs(x(2,1),t,n);
     x1 = subs(x(1,1),t,n);
     BER((n/10)+1)= 1-((x2+x1)/0.001);
end
time = 0:10:5000;
hold on
figure1 = plot(time,BER)
%for n=1:1:numel(D)
    %digitsOld = digits(5);
%xt = vpa((C(n)*(exp(D(n)*t))*eig_vec(:,n)),7)
%end
%for col=1:numel(eig_vec)
 %   for row=1:numel(eig_vec)
  %      xt(row,col)=vpa((C(n)*(exp(D(n)*t))*eig_vec(:,n)),7)
   % end
%end
%x = sum (xt,2)

%CB = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -KbCB -KbCB 0 0 0 0 KuCB KuCB
for n=0:10:5000
    CB = KuCB*(subs(x(24,1),t,n))+ KuCB*(subs(x(25,1),t,n))- KbCB*(subs(x(18,1),t,n))-KbCB*(subs(x(19,1),t,n));
    CBP((n/10)+1)=CB;
end
hold on
figure2 = plot(time,CBP)


%for i=0:10:5000
   % CB(i+1)=CB(i)+h*subs(CB,t,i);
   % time2 = i+h
%end
%[CB]= 0

