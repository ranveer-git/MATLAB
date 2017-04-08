clear all
clc
lambda = 0:0.1:2;
Pe1 = 20;
Pe2 = 1;
Pe3 = 0.05;
theta = 0.5;
theta2 = 0.2;
E1 = ((Pe1)^0.5 /(4*pi*(theta)^3))*exp(-Pe1*(1-theta)^2/(4*theta))
E2 = ((Pe2)^0.5 /(4*pi*(theta)^3))*exp(-Pe2*(1-theta)^2/(4*theta))
E3 = ((Pe3)^0.5 /(4*pi*(theta)^3))*exp(-Pe3*(1-theta)^2/(4*theta))
C1 = ((Pe1).^0.5 ./(4.*pi.*(theta))).*exp(-Pe1.*(lambda-theta).^2./(4.*theta));
C2 = ((Pe2).^0.5 ./(4.*pi.*(theta))).*exp(-Pe2.*(lambda-theta).^2./(4.*theta));
C3 = ((Pe3).^0.5 ./(4.*pi.*(theta))).*exp(-Pe3.*(lambda-theta).^2./(4.*theta));
plot(lambda,C1,'r',lambda,C2,'g',lambda,C3,'b')
xlabel('Lambda')
ylabel('C (theta=0.5)')
legend('C1','C2','C3')
h_title = title('HW 9 Q1');
figure
E1 = ((Pe1)^0.5 /(4*pi*(theta2)^3))*exp(-Pe1*(1-theta2)^2/(4*theta2))
E2 = ((Pe2)^0.5 /(4*pi*(theta2)^3))*exp(-Pe2*(1-theta2)^2/(4*theta2))
E3 = ((Pe3)^0.5 /(4*pi*(theta2)^3))*exp(-Pe3*(1-theta2)^2/(4*theta2))
C1 = ((Pe1).^0.5 ./(4.*pi.*(theta2))).*exp(-Pe1.*(lambda-theta2).^2./(4.*theta2));
C2 = ((Pe2).^0.5 ./(4.*pi.*(theta2))).*exp(-Pe2.*(lambda-theta2).^2./(4.*theta2));
C3 = ((Pe3).^0.5 ./(4.*pi.*(theta2))).*exp(-Pe3.*(lambda-theta2).^2./(4.*theta2));
plot(lambda,C1,'r',lambda,C2,'g',lambda,C3,'b')
xlabel('Lambda')
ylabel('C for theta=0.2')
legend('C1','C2','C3')
h_title = title('HW 9 Q1');