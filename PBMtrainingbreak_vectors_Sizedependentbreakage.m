tic;
clearvars
close all
clc
v=0.5:0.5:25;
x_max=length(v); % Number of particle sizes to be considered
vol=transpose(v);
diameter=zeros(x_max,1);
diameter(:)=(6.*vol(:)/pi).^(1/3);
Break_rate=zeros(x_max,1);
P1=(16*10^-3)/(20*((diameter(vol==10))^1.2));
Break_rate(:)=P1*20.*(diameter(:).^1.2)/2; % Rate for breakage
time_step=1;
t=0:time_step:600; % Defining time and step size for simulation
N=zeros(length(t),x_max); % Defining number of particles matrix
N(1,x_max)=50000; % Entering initial vlues
BreakR_f=zeros(length(t),x_max);
BreakR_d=zeros(length(t),x_max);
for a=2:length(t)
    BreakR_d(a,x_max)=Break_rate(x_max)*(N(a-1,x_max)); % For largest size
    for k=1:x_max-1 % For all other sizes
        if k>1 % All sizes except smallest size
        BreakR_d(a,k)=Break_rate(k)*N(a-1,k);
        end
        Sf=0;
        for p=k+1:x_max %Term accounting for formation of size k from larger sizes
                Sf=Sf+Break_rate(p)*(N(a-1,p)/(p-1));
        end
        BreakR_f(a,k)=2*Sf;
    end
    N(a,:)=(BreakR_f(a,:)-BreakR_d(a,:))*time_step+N(a-1,:);
end
TotalN=transpose(sum(transpose(N)));
Sum_form=transpose(sum(transpose(BreakR_f)));
Sum_dep=transpose(sum(transpose(BreakR_d)));
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
plot(diameter_PSD,transpose(N_PSD(1,:)./sum(N_PSD(1,:))),'r-o',...
    diameter_PSD,transpose(N_PSD(121,:)./sum(N_PSD(121,:))),'y-o',...
    diameter_PSD,transpose(N_PSD(241,:)./sum(N_PSD(241,:))),'g-o',...
    diameter_PSD,transpose(N_PSD(361,:)./sum(N_PSD(361,:))),'c-o',...
    diameter_PSD,transpose(N_PSD(481,:)./sum(N_PSD(481,:))),'b-o',...
    diameter_PSD,transpose(N_PSD(601,:)./sum(N_PSD(601,:))),'m-o')
title('Plot of Number PSD')
xlabel('Diameter in millimeters -->')
xlim([min(diameter_PSD) max(diameter_PSD)])
ylabel('Number fraction of Particles -->');
figure(7) % Plotting Volume PSD
plot(diameter_PSD,transpose(vol_PSD(1,:)./sum(vol_PSD(1,:))),'r-o',...
    diameter_PSD,transpose(vol_PSD(121,:)./sum(vol_PSD(121,:))),'y-o',...
    diameter_PSD,transpose(vol_PSD(241,:)./sum(vol_PSD(241,:))),'g-o',...
    diameter_PSD,transpose(vol_PSD(361,:)./sum(vol_PSD(361,:))),'c-o',...
    diameter_PSD,transpose(vol_PSD(481,:)./sum(vol_PSD(481,:))),'b-o',...
    diameter_PSD,transpose(vol_PSD(601,:)./sum(vol_PSD(601,:))),'m-o')
title('Plot of Volume PSD')
xlabel('Diameter in millimeters -->')
xlim([min(diameter_PSD) max(diameter_PSD)])
ylabel('Volume fraction of Particles -->');
toc;