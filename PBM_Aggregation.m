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
tmax=500;
t = 0:dt:tmax;
N = zeros(length(t),length(s));
N(1,1)=10^6;

rform=zeros(length(t)-1,length(s));
rdep=zeros(length(t)-1,length(s));
x(1,:)=sum(rform(1,:));
y(1,:)=sum(rdep(1,:));
volume (1,:)=sum(s.*(N(1,:)));
ratio(1,:) = x(1)./y(1,:);
number(1,:)=sum(N(1,:));
%d43(1,:) = (sum (num(1,:)))./(sum (den(1,:)));

for i=1:(length (t)-1);
    for j=1:length (s);
        if j>1
            rform(i,j)=0.5*beta*sum(N(i,1:j-1).*N(i,j-1:-1:1));
        end
        if j<max(s)
            rdep(i,j)=beta*N(i,j)*sum(N(i,1:max(s)-j));
        end
        N(i+1,:)=(rform(i,:)-rdep(i,:))*dt+N(i,:);
       
        x(i+1,:)=sum(rform(i,:));
        y(i+1,:)=sum(rdep(i,:));
        ratio(i+1,:) = x(i+1)./y(i+1,:);
        volume (i+1,:)=sum(s.*(N(i,:)));
        %Total number of particles at each time point
        number(i+1,:)=sum(N(i,:));
        %Calculating d43
        D = repmat(d,length (t),1);
        num = (D.^4).*N;
        den = (D.^3).*N;
        d43(i+1,:) = (sum (num(i,:)))./(sum (den(i,:)));   
        
    end
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

%#8 second set of 8 bins
d1 = zeros(1,8);
dmin = 1.2;
dmax = 2.6;
dt1= (dmax-dmin)/(length(d1)-1);
d1 = dmin:dt1:dmax;
transformationMatrix = zeros(2,length(d));
d_new = unique((floor(d./0.2)).*0.2);
s_new = (pi/6).*(d_new.^3);
N1=zeros(length(t),length(d_new));
M1=zeros(length(t),length(d_new));

for dIndex=1:length(d)
    dCurrent = d(dIndex);
    for bucketIndex = 1 : length(d_new);
       currentBucket = d_new(bucketIndex);
       if(currentBucket > dCurrent)
        break
       end
       if(currentBucket <= dCurrent)
          transformationMatrix(1,dIndex) = dIndex;
          transformationMatrix(2,dIndex) = bucketIndex;
       end
    end
    N1(:,transformationMatrix(2,dIndex))=N1(:,transformationMatrix(2,dIndex))+N(:,transformationMatrix(1,dIndex));
    M1(:,transformationMatrix(2,dIndex))=N1(:,transformationMatrix(2,dIndex)).*s_new(:,transformationMatrix(2,dIndex))+N(:,transformationMatrix(1,dIndex)).*s(:,transformationMatrix(1,dIndex));
end
   sumN1=repmat((sum(N1'))',[1,length(d_new)]);
 sumM1=repmat((sum(M1'))',[1,length(d_new)]);
figure (5)
plot (d_new,(N1((t==5),:)),...
    d_new,(N1((t==500),:)))

%normalised mass frequency
NormNumFreq = N1./sumN1;
NormMassFreq = M1./sumM1;
CumuMassFreq= zeros(length(t),length(d_new));
CumuMassFreq(1,1)= NormMassFreq(1,1);

    for j = 1:length(d_new)
        CumuMassFreq(:,j)=(sum (NormMassFreq(:,1:j)'))';
    end

figure (6)
plot (d_new,NormMassFreq(t==5,:),...
   d_new,NormMassFreq(t==500,:))
legend ('5 sec','500 sec')
title('Normalised Mass Frequency')

%normalised mass frequency
figure (7)
plot (d_new,NormNumFreq(t==5,:),...
   d_new,NormNumFreq(t==500,:))
legend ('5 sec','500 sec')
title('Normalised Number Frequency')

%Calculation of d50

