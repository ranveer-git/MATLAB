close all
clc
clear all
A=magic(5)
[V,D]=eig(A)
B=min(V)

for col=1:5
%temp = B(n);

    for row=1:5
        N(row, col) = round(V(row,col)/(B(col)));
    end

end
N
x


%V1 = V(:,1)/-0.4472 

%C = round(bsxfun(@rdivide, V, B))
%C(isnan(C)) = 0;
A= [1 5 0 2;5 4 6 6; 3 3 0 5;9 2 8 7]
B= ones(1,4)
C=pascal(6)
D=eye(4,4)
E=zeroes(3,3)
