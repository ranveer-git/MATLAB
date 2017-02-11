clc
clear all

A=magic(3)
[V,E]=eig(A)
D=diag(E);
b=[1 2 3]';
c=V\b
syms tempMatrix
syms t 
for n=1:1:3
    x(t) =vpa (c(n).*exp(D(n)*t).*V(:,n),4)
    if(n ==1)
        tempMatrix = vpa([x(t)],4);
    else
        tempMatrix = vpa([tempMatrix x(t)],4)
    end
end
%t=0:1:10;
for row=1:numel(D)
xt = vpa(sym(zeros(1,numel(D))),4);
    for col=1:numel(D)
        xt(row,col)=vpa(tempMatrix(row,col),4);
       xs=vpa(sum(xt),4)
    end
 %xs = vpa(sym(zeros(numel(D),1)),4);
 %plot(t,xs);
%hold on
end
%y=sum([tempMatrix{:}],2)
%xs = symsum (tempMatrix,1)