tic;
clearvars
close all
clc
Break_kernel=8*(10^-3); % Rate constant for breakage
vol_min=0.5; % Least possible volume
vol_max=25; % Largest possible volume
vol_step=0.5;
x_max=(vol_max-vol_min)/(vol_step)+1; % Number of particle sizes to be considered
vol=zeros(x_max,1);
vol(1)=vol_min; % Smallest particle size
for b=1:x_max % Entering other particle sizes
    vol(b)=vol_step*b;
end
diameter=zeros(x_max,1);
diameter(:)=(6.*vol(:)/pi).^(1/3);
time_step=120;
t=0:time_step:600; % Defining time and step size for simulation
N=zeros(length(t),x_max); % Defining number of particles matrix
N(1,x_max)=50000; % Entering initial vlues
Sum_dep=zeros(length(t),1);
Sum_form=zeros(length(t),1);
TotalN=zeros(length(t),1);
for a=2:length(t)
    for k=1:x_max
        if k==1 %Equation for lowest particle size
            Sf=0; 
            for e=k+1:x_max %Term accounting for gain of lowest particle size due to breakage of  others
                Sf=Sf+(N(a-1,e)/(e-1));
            end
            N(a,k)=(Break_kernel*2*Sf)*time_step+N(a-1,k);
            Sum_form(a)=Sum_form(a)+Break_kernel*2*Sf;
        elseif k==x_max %max size
            N(a,k)=(-Break_kernel*N(a-1,k))*time_step+N(a-1,k); %Term accounting for formation of size k from smaller sizes
            Sum_dep(a)=Sum_dep(a)+Break_kernel*N(a-1,k);
        else %all other sizes
            Sf=0;
            for p=k+1:x_max %Term accounting for formation of size k from larger sizes
                Sf=Sf+(N(a-1,p)/(p-1));
            end
            N(a,k)=(Break_kernel*2*Sf-Break_kernel*N(a-1,k))*time_step+N(a-1,k);
            Sum_dep(a)=Sum_dep(a)+Break_kernel*N(a-1,k);
            Sum_form(a)=Sum_form(a)+Break_kernel*2*Sf;
        end
    end
    TotalN(a)=sum(N(a,:));
end
clear a k;
diameter_PSD=unique(floor(diameter./0.5).*0.5);
N_PSD=zeros(length(t),length(diameter_PSD)); % Defining No. of particles matrix
vol_PSD=zeros(length(t),length(diameter_PSD));
for a=1:length(diameter_PSD)
    for k=1:x_max
        if a~=length(diameter_PSD)
            if diameter_PSD(a)<diameter(k) && diameter(k)<diameter_PSD(a+1)
                N_PSD(:,a)=N_PSD(:,a)+N(:,k);
                vol_PSD(:,a)=vol_PSD(:,a)+N(:,k)*vol(k);
            end
        else
            if diameter_PSD(a)<diameter(k)
                N_PSD(:,a)=N_PSD(:,a)+N(:,k);
                vol_PSD(:,a)=vol_PSD(:,a)+N(:,k)*vol(k);
            end
        end
    end
end
figure(1) % Plotting Number of particles of each size class against time
plot(t,N)
title('Plot of Number of Particles of each size class against Time')
xlabel('Time in seconds -->')
ylabel('Number of Particles -->')
figure(2) % Plotting Total of particles of against time
plot(t,N*vol)
title('Plot of Total Volume of Particles against Time')
xlabel('Time in seconds -->')
ylabel('Total Volume of Varticles in mm^3 -->');
figure(3) % Plotting Total of particles against time
plot(t,Sum_form./Sum_dep)
title('Plot of Formation to Depletion ratio against Time')
xlabel('Time in seconds -->')
ylabel('Formation to Depletion ratio -->');
figure(4) % Plotting d43 against time
plot(t,(N*(diameter.^4))./(N*(diameter.^3)))
title('Plot of d43 against Time')
xlabel('Time in seconds -->')
ylabel('d43 in mm -->');
figure(5) % Plotting Number of particles against time
plot(t,TotalN)
title('Plot of Total Number of Particles against Time')
xlabel('Time in seconds -->')
ylabel('Number of Particles -->');
figure(6) % Plotting Number PSD
%{
plot(diameter_PSD,transpose(N_PSD(1,:)./sum(N_PSD(1,:))),...
    diameter_PSD,transpose(N_PSD(121,:)./sum(N_PSD(121,:))),...
    diameter_PSD,transpose(N_PSD(241,:)./sum(N_PSD(241,:))),...
    diameter_PSD,transpose(N_PSD(361,:)./sum(N_PSD(361,:))),...
    diameter_PSD,transpose(N_PSD(481,:)./sum(N_PSD(481,:))),...
    diameter_PSD,transpose(N_PSD(601,:)./sum(N_PSD(601,:))))

title('Plot of Total Number PSD')
xlabel('Diameter in millimeters -->')
xlim([min(diameter_PSD) max(diameter_PSD)])
ylabel('Number fraction of Particles -->');
figure(7) % Plotting Volume PSD
plot(diameter_PSD,transpose(vol_PSD(1,:)./sum(vol_PSD(1,:))),...
    diameter_PSD,transpose(vol_PSD(121,:)./sum(vol_PSD(121,:))),...
    diameter_PSD,transpose(vol_PSD(241,:)./sum(vol_PSD(241,:))),...
    diameter_PSD,transpose(vol_PSD(361,:)./sum(vol_PSD(361,:))),...
    diameter_PSD,transpose(vol_PSD(481,:)./sum(vol_PSD(481,:))),...
    diameter_PSD,transpose(vol_PSD(601,:)./sum(vol_PSD(601,:))))
title('Plot of Total Volume PSD')
xlabel('Diameter in millimeters -->')
xlim([min(diameter_PSD) max(diameter_PSD)])
ylabel('Volume fraction of Particles -->');
toc;
%}