clc
close all
clear all

beta = 7*10^-9;
smin=1;
smax=10;
ds=1;
s=smin:ds:smax;
%Diameter of each bin
d= ((6/pi).*s).^(1/3);

dt=1;
tmax=100;
t = 0:dt:tmax;
N = zeros(tmax+1,smax);
N(1,1)=10^6;

rform=zeros(length(t)-1,length(s));
rdep=zeros(length(t)-1,length(s));
x(1,:)=sum(rform(1,:));
y(1,:)=sum(rdep(1,:));
volume (1,:)=sum(s.*(N(1,:)));
ratio(1,:) = x(1)./y(1,:);
number(1,:)=sum(N(1,:));
%d43(1,:) = (sum (num(1,:)))./(sum (den(1,:)));


for j=1:max(t);
    for j=1:max(s);
        if j>1
            rform(j,j)=0.5*beta*sum(N(j,1:j-1).*N(j,j-1:-1:1));
        end
        if j<max(s)
            rdep(j,j)=beta*N(j,j)*sum(N(j,1:max(s)-j));
        end
        %for x=1:ds:smax;
         %   rform(x,j)=0.5*beta*N(i,j)*N(smax-i,j);
          %  rdep(x,j)=beta*N(i,j)*N(smax-x,j);
   % end
        N(j+1,:)=(rform(j,:)-rdep(j,:))*dt+N(j,:);
        x(j+1,:)=sum(rform(j,:));
        y(j+1,:)=sum(rdep(j,:));
        ratio(j+1,:) = x(j+1)./y(j+1,:);
        volume (j+1,:)=sum(s.*(N(j,:)));
        %Total number of particles at each time point
        number(j+1,:)=sum(N(j,:));
        %Calculating d43
        D = repmat(d,tmax+1,1);
        num = (D.^4).*N;
        den = (D.^3).*N;
        d43(j+1,:) = (sum (num(j,:)))./(sum (den(j,:)));
                
    end
    
      %  N(i,j)=0.5*beta*N(i-1,j)*N(i-1,smax-j)-beta*N(i,j)*N(i-1,smax-j);
      %  N(i+1,j+1)=N(i,j)+h*N(i,j);
end

figure (1)
plot (t,ratio)
hold on
%plot (t,y)
figure (2)
plot (t,volume)
figure (3)
plot (t,number)
figure (4)
plot (t,d43)

%%%%%#8 second set of 8 bins
%d1 = zeros(1,8);
%dmin = 1.2;
%dmax = 2.6;
%dt1= (dmax-dmin)/(length(d1)-1);
%d1 = dmin:dt1:dmax;
%%%%D1 = repmat (d1,length(t),1);
%%%%%S = repmat (s,length(t),1);

%%%%%%%%%%%
%transformationMatrix = zeros(2,length(d));
%for dIndex=1:length(d)
%    dCurrent = d(dIndex);
%    for bucketIndex = 1 : length(d1);
%       currentBucket = d1(bucketIndex);
%       if(currentBucket > dCurrent)
%        break
%       end
%       if(currentBucket <= dCurrent)
%          transformationMatrix(1,dIndex) = dIndex;
%          transformationMatrix(2,dIndex) = bucketIndex;
  %     end
%    end
%end

%N1=zeros(tmax+1,length(d1))
%for i = 1:tmax+1
%for j=1:length(d)
 %   N1(i,transformationMatrix(2,j))=N1(i,transformationMatrix(2,j))+N(i,transformationMatrix(1,j));
%end
%end

    
%%%%%%%%%%%

