clc
close all
clear all

Kbreak = 7*10^-2;
s=0.5:0.5:25;
index=find(s);
%Diameter of each bin
d= ((6/pi).*s).^(1/3);

%Time
dt=1;
t = 0:dt:2000;
%flowrate in
Fin=zeros(length(t),length(s));
Fin(:,end)=50;
S= repmat(s,length (t),1);
N = zeros(length(t),length(s)); %define number of particles

rform=zeros(length(t)-1,length(s));
rdep=zeros(length(t)-1,length(s));

x(1,:)=sum(rform(1,:));
y(1,:)=sum(rdep(1,:));
volume (1,:)=sum(s.*(N(1,:)));
ratio(1,:) = x(1)./y(1,:);
number(1,:)=sum(N(1,:));

for i=1:(length (t)-1);
    for j=1:length (s);
        if j>1
            rdep(i,j)= Kbreak.*N(i,j);
        end
        if j<length(s)
            rform(i,j)=2*Kbreak.*sum(N(i,j+1:end)./index(j:end-1));
        end
        
        N(i+1,:)=N(i,:)+(rform(i,:)-rdep(i,:)).*dt+(Fin(i,:)-0.01*N(i,:)).*dt;
        x(i+1,:)=sum(rform(i,:));
        y(i+1,:)=sum(rdep(i,:));
        ratio(i+1,:) = x(i+1)./y(i+1,:);
        volume (i+1,:)=sum(s.*(N(i,:)));
        
        %Total number of particles at each time point
        number(i+1,:)=sum(N(i,:));   
    end
end

%Calculating d43
 D = repmat(d,length (t),1);
       num = (D.^4).*N;
       den = (D.^3).*N;
       d43 = (sum (num'))'./(sum (den'))';
       
figure (1)
plot (t,N)
title('Number of particles over each class size v/s time')
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

figure (5)
plot (t,number)
title('Total number of particles over time')

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

d_new = unique((floor(d./0.5)).*0.5);
transformationMatrix = zeros(2,length(d));
N1=zeros(length(t),length(d_new));
s_new = zeros(length(t),length(d_new));
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
    N1(:,transformationMatrix(2,dIndex))=N1(:,transformationMatrix(2,dIndex))+0.01.*N(:,transformationMatrix(1,dIndex));
    s_new(:,transformationMatrix(2,dIndex))=s_new(:,transformationMatrix(2,dIndex))+0.01.*N(:,transformationMatrix(1,dIndex)).*s(:,transformationMatrix(1,dIndex));
    
end
 sumN1=repmat((sum(N1'))',[1,length(d_new)]);
 sums_new=repmat((sum(s_new'))',[1,length(d_new)]);
 number1=sum(N1')';
 
%normalised mass frequency
NormNumFreq1 = N1./sumN1;
NormVolFreq1 = s_new./sums_new;
CumuVolFreq1= zeros(length(t),length(d_new));
CumuVolFreq1(:,1)= NormVolFreq1(:,1);

    for k = 2:length(d_new)
        CumuVolFreq1(:,k)=(sum (NormVolFreq1(:,1:k)'))';
    end

figure (7)
plot (d_new,NormVolFreq1(t==120,:),...
    d_new,NormVolFreq1(t==240,:),...
    d_new,NormVolFreq1(t==360,:),...
    d_new,NormVolFreq1(t==480,:),...
    d_new,NormVolFreq1(t==600,:))
grid on
title('Normalised Vol Frequency')
legend('120 s','240 s','360 s','480 s','600 s')
xlabel ('d_new')
ylabel ('Normalised Mass Frequency New')

%normalised mass frequency
figure (8)
plot (d_new,NormNumFreq1(t==120,:),...
    d_new,NormNumFreq1(t==240,:),...
    d_new,NormNumFreq1(t==360,:),...
    d_new,NormNumFreq1(t==480,:),...
    d_new,NormNumFreq1(t==600,:))
grid on
title('Normalised Number Frequency')
legend('120 s','240 s','360 s','480 s','600 s')
xlabel ('d_new')
ylabel ('Normalised Number Frequency New')
