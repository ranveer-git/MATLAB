%%
tic;
clearvars
close all
clc
vol=transpose(1:1:100);
% vol=transpose(0.5:0.5:25);
x_max=length(vol); % Number of particle sizes to be considered
diameter=(6.*vol/pi).^(1/3);
%%
Agg_kernel_constant=7*(10^-9); % Kernel constant for aggregation
% Agg_kernel_constant=2.5*(10^-7); % Kernel constant for aggregation
Agg_kernel=zeros(x_max);
for c=1:x_max-1;
    for e=1:x_max-c;
        Agg_kernel(c,e)=Agg_kernel_constant; %*2*((diameter(c)+diameter(e))^2)*sqrt((diameter(c)^-3)+(diameter(e)^-3)); % Aggregation EKE kernel
    end
end
%%
Break_rate=zeros(x_max,1);
Break_rate(:)=8*(10^-3)*2.*(diameter(:).^1.2)/2; % Rate for breakage
% Break_rate(:)=1.2*(10^5)*20.*(diameter(:).^3)/2; % Rate for breakage
Break_distribution=zeros(x_max);
b=zeros(x_max);
Break_sigma=1;
for c=1:x_max
    for e=c+1:x_max
        b(c,e)=exp((-1*(c-(e/2))^2)/(2*Break_sigma^2));
    end    
end
for c=1:x_max
    Break_distribution(c,:)=2*b(c,:)./sum(b);
end
%%
time_step=1;
t=0:time_step:1000; % Defining time and step size for simulation
N=zeros(length(t),x_max); % Defining number of particles matrix
N(1,47)=0.02*1E6;
N(1,48)=0.13*1000000; % Entering initial vlues
N(1,49)=0.35*1000000; % Entering initial vlues
N(1,52)=0.02*1E6;
N(1,51)=0.13*1000000; % Entering initial vlues
N(1,50)=0.35*1000000; % Entering initial vlues
AggR_f=zeros(length(t),x_max);
AggR_d=zeros(length(t),x_max);
BreakR_f=zeros(length(t),x_max);
BreakR_d=zeros(length(t),x_max);
%%
for a=2:length(t)
    %%
    for k=1:x_max
        if k<x_max
            AggR_d(a,k)=N(a-1,k)*sum((N(a-1,1:x_max-k).*Agg_kernel(k,1:x_max-k)));
        end
        if k>1
            for q=1:k-1 % Formation rates calculation
                AggR_f(a,k)=AggR_f(a,k)+0.5*(N(a-1,q)*N(a-1,k-q).*Agg_kernel(q,k-q));
            end
        end     
    end
    %%
    for k=1:x_max
        if k<x_max
            BreakR_f(a,k)=(N(a-1,k+1:x_max).*Break_distribution(k,k+1:x_max))*Break_rate(k+1:x_max,1);
        end
        if k>1
            BreakR_d(a,k)=Break_rate(k)*N(a-1,k);
        end
    end
    %%
    N(a,:)=(AggR_f(a,:)-AggR_d(a,:)+BreakR_f(a,:)-BreakR_d(a,:))*time_step+N(a-1,:);
end
%%
clear a k;
diameter_PSD=unique(floor(diameter./1).*1);
N_PSD=zeros(length(t),length(diameter_PSD)); 
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
%%
TotalN=transpose(sum(transpose(N)));
Sum_form_agg=transpose(sum(transpose(AggR_f)));
Sum_dep_agg=transpose(sum(transpose(AggR_d)));
Sum_form_break=transpose(sum(transpose(BreakR_f)));
Sum_dep_break=transpose(sum(transpose(BreakR_d)));
%%
figure(1) % Plotting Number of particles against time
plot(t,N)
title('Plot of Number of Particles against Time')
xlabel('Time in seconds -->')
ylabel('Number of Particles -->');
%%
figure(2) % Plotting Total volume of particles against time
plot(t,N*vol,'b-')
title('Plot of Total Volume of Particles against Time')
xlabel('Time in seconds -->')
ylabel('Total Volume of Varticles in mm^3 -->');
%%
figure(3) % Plotting Aggregation formation to depletion ratio against time
plot(t,Sum_form_agg./Sum_dep_agg,'b')
title('Plot of Aggregation formation to depletion ratio against Time')
xlabel('Time in seconds -->')
ylabel('Aggregation formation to depletion ratio -->');
%%
figure(4) % Plotting Breakage formation to depletion ratio against time
plot(t,Sum_form_break./Sum_dep_break,'b')
title('Plot of Breakage formation to depletion ratio against Time')
xlabel('Time in seconds -->')
ylabel('Breakage formation to depletion ratio -->');
%%
figure(5) % Plotting Aggregation formation to Breakage depletion  ratio against time
plot(t,(Sum_form_agg./Sum_dep_break),'b')
title('Plot of Aggregation formation to Breakage depletion ratio against Time')
xlabel('Time in seconds -->')
ylabel('Aggregation formation to Breakage depletion ratio -->');
%%
figure(6) % Plotting Aggregation depletion to Breakage formation ratio against time
plot(t,Sum_dep_agg./Sum_form_break,'b')
title('Plot of Aggregation depletion to Breakage formation ratio against Time')
xlabel('Time in seconds -->')
ylabel('Aggregation depletion to Breakage formation ratio -->');
%%
figure(7) % Plotting Total Aggregation to Total Breakage ratio against time
plot(t,(Sum_form_agg+Sum_dep_agg)./(Sum_form_break+Sum_form_agg),'b')
title('Plot of Total Aggregation to Total Breakage ratio against Time')
xlabel('Time in seconds -->')
ylabel('Total Aggregation to Total Breakage ratio -->');
%%
figure(8) % Plotting Total formation to Total depletion ratio against time
plot(t,(Sum_form_agg+Sum_form_break)./(Sum_dep_agg+Sum_form_agg),'b')
title('Plot of Total formation to Total depletion ratio against Time')
xlabel('Time in seconds -->')
ylabel('Total formation to Total depletion ratio -->');
%%
figure(9) % Plotting d43 against time
plot(t,(N*(diameter.^4))./(N*(diameter.^3)),'b-')
title('Plot of d43 against Time')
xlabel('Time in seconds -->')
ylabel('d43 in mm -->');
%%
figure(10) % Plotting Number of particles against time
plot(t,TotalN,'b-')
title('Plot of Number of Particles against Time')
xlabel('Time in seconds -->')
ylabel('Number of Particles -->');
%%
figure(11) % Plotting Number PSD
plot(diameter_PSD,transpose(N_PSD(1,:)./sum(N_PSD(1,:))),'r-o',...
    diameter_PSD,transpose(N_PSD(6,:)./sum(N_PSD(6,:))),'g-o',...
    diameter_PSD,transpose(N_PSD(501,:)./sum(N_PSD(501,:))),'b-o')
title('Plot of Number PSD')
xlabel('Diameter in millimeters -->')
xlim([min(diameter_PSD) max(diameter_PSD)])
ylabel('Number fraction of Particles -->');
%%
figure(12) % Plotting Volume PSD
plot(diameter_PSD,transpose(vol_PSD(1,:)./sum(vol_PSD(1,:))),'r-o',...
    diameter_PSD,transpose(vol_PSD(6,:)./sum(vol_PSD(6,:))),'g-o',...
    diameter_PSD,transpose(vol_PSD(501,:)./sum(vol_PSD(501,:))),'b-o')
title('Plot of Volume PSD')
xlabel('Diameter in millimeters -->')
xlim([min(diameter_PSD) max(diameter_PSD)])
ylabel('Volume fraction of Particles -->');
%%
toc;