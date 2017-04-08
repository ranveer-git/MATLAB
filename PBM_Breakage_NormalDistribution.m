clc
close all
clear all

Kbreak = 8*10^-3;
sigma = 1;

s=1:1:10;%volume bins
d= ((6/pi).*s).^(1/3);%Diameter of each bin
dt=1;
tmax=600;
t = 0:dt:tmax;
S= repmat(s,length (t),1);
N = zeros(length(t),length(s));
N(1,end)=50000;

rform=zeros(length(t)-1,length(s));
rdep=zeros(length(t)-1,length(s));

x(1,:)=sum(rform(1,:));
y(1,:)=sum(rdep(1,:));
volume (1,:)=sum(s.*(N(1,:)));
ratio(1,:) = x(1)./y(1,:);
number(1,:)=sum(N(1,:));
%d43(1,:) = (sum (num(1,:)))./(sum (den(1,:)));

%the normal distribution function
bij = zeros(length(d));
Actual_bij = zeros(length(d));
for i=2:length(d)
    for j = 1:i-1
        bij(i,j) = exp(-(j-i/2)^2/(2*sigma.^2));
    end
end
for l=1:length(d)
    Actual_bij(:,l)=2*bij(:,l)./sum(bij,2);
end

NormalDistribution = Actual_bij';
%number of particles over time
for m=1:(length (t)-1);
    for n=1:length (s);
        if n<length(s)
           rform(m,n)=Kbreak*sum(N(m,n+1:end).*NormalDistribution(n,n+1:end));
        end
        if n>1
            rdep(m,n)=Kbreak*N(m,n);
        end
        
        N(m+1,:)=(rform(m,:)-rdep(m,:))*dt+N(m,:);
        x(m+1,:)=sum(rform(m,:));
        y(m+1,:)=sum(rdep(m,:));
        ratio(m+1,:) = x(m+1)./y(m+1,:);
        volume (m+1,:)=sum(s.*(N(m,:)));
        number(m+1,:)=sum(N(m,:));   %Total number of particles at each time point

    end
end
%Calculating d43
 D = repmat(d,length (t),1);
       num = (D.^4).*N;
       den = (D.^3).*N;
       d43 = (sum (num'))'./(sum (den'))';
       
figure (1)
plot (t,N)
title('Total number of particles over time')
xlabel('Time')
ylabel('Number of particles of each size')
legend('1','2','3','4','5','6','7','8','9','10')
hold on

figure (2)
plot (t,volume)
title('Total volume of particles over time')

figure (3)
plot (t,ratio)
title('Ratio of Total formation rate to Total depletion rate over time')

figure (4)
plot (t,d43)
title('Average diameter (d4,3) over time')

%Calculation of d50
M = N.*S;
%normalised mass frequency
NormNumFreq = N./repmat (sum(N')',[1,length(d)]);
NormMassFreq = M./repmat(sum(M')',[1,length(d)]);
CumuMassFreq= zeros(length(t),length(d));
CumuMassFreq(:,1)= NormMassFreq(:,1);

    for j = 2:length(d)
        CumuMassFreq(:,j)=sum (NormMassFreq(:,1:j)')';
    end

figure (6)
plot (d,NormMassFreq(t==120,:),...
    d,NormMassFreq(t==600,:))
legend ('600 sec')
title('Normalised Mass Frequency')
xlabel ('d')
ylabel ('Normalised Mass Frequency')

%Particle size distribution for new diameter

%#8 second set of 8 bins

dmin = 0;
dmax = 4;
dt1= 0.5;
d_new = dmin:dt1:dmax;
transformationMatrix = zeros(2,length(d));
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
 number1=sum(N1')';
 
%normalised mass frequency
NormNumFreq1 = N1./sumN1;
NormMassFreq1 = M1./sumM1;
CumuMassFreq1= zeros(length(t),length(d_new));
CumuMassFreq1(:,1)= NormMassFreq1(:,1);

    for k = 2:length(d_new)
        CumuMassFreq1(:,k)=(sum (NormMassFreq1(:,1:k)'))';
    end

figure (7)
plot (d_new,NormMassFreq1)
title('Normalised Mass Frequency')
legend('0 sec','120 s','240 s','360 s','480 s','600 s')
xlabel ('d_new')
ylabel ('Normalised Mass Frequency New')

%normalised mass frequency
figure (8)
plot (d_new,NormNumFreq1)
title('Normalised Number Frequency')
legend('0 sec','120 s','240 s','360 s','480 s','600 s')
xlabel ('d_new')
ylabel ('Normalised Number Frequency New')


