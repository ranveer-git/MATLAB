clc
close all
clear all

k=2.5*10^-7;
sigma = 1;
P1 = 1.2e5;
P2 = 3;
Gshear = 20;
smin=1;
smax=10;
ds=1;
s=smin:ds:smax;
index=find(s);

%Diameter of each bin
d= ((6/pi).*s).^(1/3);
dt=1;
tmax=500;
t = 0:dt:tmax;

%Calculation of beta
beta = zeros(length(d));
for m=1:length(d)-1
    for n = 1:length(d)-m
        beta(m,n) = k*(d(m)/2+d(n)/2)^2*sqrt(1/d(m)^3+d(n)^3);
    end
end

%Calculation of normal distribution function
bij = zeros(length(d));
Actual_bij = zeros(length(d));
for i=1:length(d)
    for j = i+1:length(d)
        bij(i,j) = exp(-(j-i/2)^2/(2*sigma.^2));
    end
end
for i=2:length(d)
    Actual_bij(i,:)=2*bij(i,:)./sum(bij);
end

%Calculation of Breakage Rate
Kbreakx = P1*Gshear.*((d*10^3).^P2)/2;
Kbreakx1 = repmat (Kbreakx, [length(t),1]);

%Initializing 
N = zeros(length(t),length(s));
N(1,1)=10^6;
rform_agg=zeros(length(t)-1,length(s));
rdep_agg=zeros(length(t)-1,length(s));
rform_break=zeros(length(t)-1,length(s));
rdep_break=zeros(length(t)-1,length(s));

ratio_agg(1,:)=(sum(rform_agg(1,:)))/(sum(rdep_agg(1,:)));
ratio_break(1,:) = sum(rform_break(1,:))/sum(rdep_break(1,:));
volume (1,:)=sum(s.*(N(1,:)));
number(1,:)=sum(N(1,:));
%d43(1,:) = (sum (num(1,:)))./(sum (den(1,:)));

for i=1:length (t)-1
    for k = 1:length(s)
            if k>1
                for l=1:k-1
                    rform_agg(i,k)=rform_agg(i,k)+0.5*(N(i,l)*N(i,k-l).*beta(l,k-l));
                end
            end
            if k<max(s)
               rdep_agg(i,k)=N(i,k)*sum((beta(k,1:max(s)-k).*N(i,1:max(s)-k)));
            end
               
            if k>1
               rdep_break(i,k)=Kbreakx1(i,k).*N(i,k);
            end
            if k<max(s)
               rform_break(i,k)=sum (N(i,k+1:end).*Actual_bij(k,k+1:max(s)).*(Kbreakx1(i,k+1:end)));
            end
            
            N(i+1,:)=(rform_agg(i,:)-rdep_agg(i,:)+rform_break(i,:)-rdep_break(i,:))*dt+N(i,:);

            ratio_agg(i+1,:)=(sum(rform_agg(i,:)))/(sum(rdep_agg(i,:)));
            ratio_break(i+1,:) = sum(rform_break(i,:))/sum(rdep_break(i,:));
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
plot (t,number)
title('Total number of particles over time')
hold on

figure (2)
plot (t,volume)
title('Total volume of particles over time')

figure (3)
plot (t,ratio_agg)
title('Ratio of Total formation rate to Total depletion rate for aggregation over time')

figure (4)
plot (t,ratio_break)
title('Ratio of Total formation rate to Total depletion rate for breakage over time')

figure (5)
plot (t,d43)
title('Average diameter (d4,3) over time')

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
figure (6)
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

figure (7)
plot (d_new,NormMassFreq(t==5,:),...
   d_new,NormMassFreq(t==500,:))
legend ('5 sec','500 sec')
title('Normalised Mass Frequency')

%normalised mass frequency
figure (8)
plot (d_new,NormNumFreq(t==5,:),...
   d_new,NormNumFreq(t==500,:))
legend ('5 sec','500 sec')
title('Normalised Number Frequency')

%Calculation of d50

%}