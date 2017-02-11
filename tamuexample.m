A= [1 5 0 2;5 4 6 6; 3 3 0 5;9 2 8 7]
%B= ones(1,4)
%C=pascal(6)
%D=eye(4,4)
%E=zeros(3,3)
%column1=A(:,1)
%column2=A(:,2)
%sum=column1+column2
[ev,e]=eig(A)
D=diag(e)
%syms t
xt=zeros(size(A));
for n=1:1:numel(A)
xt(:,n) = (exp(D(n)))*ev(:,n)
%sum = sum+xt
end

%column1t=column1*t
%column2t=column2*t
%sumt = column1t+column2t