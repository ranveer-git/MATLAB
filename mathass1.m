clc                                            
clear all

h=.1;                                             
x = 0:h:2;                                         
%y = zeros(1,length(x)); 
y=[0,2,0]';  
syms
for j=1:2
    h=[.1,.2];
    x=0:h(j):2;
for i=1:(length(x)-1) 
    if j==1
    k1 = Fxy(x(i),y(:,i));
    k2 = Fxy(x(i)+0.5*h(j),y(:,i)+0.5*h(j)*k1);
    k3 = Fxy((x(i)+0.5*h(j)),(y(:,i)+0.5*h(j)*k2));
    k4 = Fxy((x(i)+h(j)),(y(:,i)+k3*h(j)));
    y(:,i+1) = y(:,i) + (1/6)*(k1+2*k2+2*k3+k4)*h(j);
    end
    if j==2
        y(:,23)=[0,2,0]';
        k1 = Fxy(x(i),y(:,i+22));
    k2 = Fxy(x(i)+0.5*h(j),y(:,i+22)+0.5*h(j)*k1);
    k3 = Fxy((x(i)+0.5*h(j)),(y(:,i+22)+0.5*h(j)*k2));
    k4 = Fxy((x(i)+h(j)),(y(:,i+22)+k3*h(j)));
    y(:,i+23)=y(:,i+22)+(1/6)*(k1+2*k2+2*k3+k4)*h(j);
    end
end
end
x1=1:21;
%plot(x1,y(2,:));
hold on
z=zeros(2,length(x));
h=[.1,.2];
for f=1:2
    x1=0:h(f):2;
for j=1:(length(x1))
    if f==1
        z(1,j)=2*(1-(tanh(x1(j)*sqrt(2)))^2);
    
    
    end
    
    if f==2
        z(2,j)=2*(1-(tanh(x1(j)*sqrt(2)))^2);
    
    end
    
end
if f==1
    error1=z(1,:)-y(2,1:21);
    plot(x1,abs(error1))
    hold on
end
end
error2=z(2,1:11)-y(2,23:33);
    plot(x1,abs(error2),'r')
%x12=1:21;
%plot (x12,z(1,:),'r')
%hold on
%error1=z-y(2,1:21);
%error2=z-y(2,23:33);
%plot(x,error1,error2,'r')
%g=1;

%[x11,y11] = meshgrid(x,y(2,:));
%z=x11+i*y11;
g = 1 + z + 0.5*z.^2+(z.^3)/6+(z.^4)/24;
gm=abs(g);
%contour(x11,y11,gm)
    