%PHSY371 Assignment 3
%Ahmet Burak Catli
%Pendulum string system differential equation solution


%initialize variables

g = 10;
m1 = 1;
m2 = 1;
l1 = 10;
l2 = 10;
kspring = 200;
theta1_init = 0.2; 
theta2_init = 0.1;
omega1_init = -5;
omega2_init = 3;
h = 0.001;
t_end = 10;

%---------------------------%


t = [0:h:t_end];
tprime = [0:h:t_end+h];
theta1 = zeros(1, length(t)+1);
theta2 = zeros(1, length(t)+1);
omega1 = zeros(1, length(t)+1);
omega2 = zeros(1, length(t)+1);

theta1(1) = theta1_init;
theta2(1) = theta2_init;
omega1(1) = omega1_init;
omega2(1) = omega2_init;

dw2_dt = inline('-t2.*(gl2+km)+t1.*km.*l1l2', 't1', 't2','gl1', 'gl2', 'km', 'l1l2');
dw1_dt = inline('-t1.*(gl1+km)+t2.*km.*l2l1', 't1', 't2','gl1', 'gl2', 'km', 'l2l1');

km1 = kspring/m1;
km2 = kspring/m2;
l1l2 = gl2./gl1;
l2l1 = 1./l1l2;
gl1 = g/l1;
gl2 = g/l2;

%---- Velocity Verlet --------------%

for k = 1:length(t)
    theta2(k+1) = theta2(k) + omega2(k).*h + 0.5*dw2_dt(theta1(k), theta2(k), gl1, gl2, km2, l1l2).*(h^2);
    theta1(k+1) = theta1(k) + omega1(k).*h + 0.5*dw1_dt(theta1(k), theta2(k), gl1, gl2, km1, l2l1).*(h^2);
    omega2(k+1) = omega2(k) + 0.5*(dw2_dt(theta1(k), theta2(k), gl1, gl2, km2, l1l2) + dw2_dt(theta1(k+1), theta2(k+1), gl1, gl2, km2, l1l2))*h;
    omega1(k+1) = omega1(k) + 0.5*(dw1_dt(theta1(k), theta2(k), gl1, gl2, km1, l2l1) + dw1_dt(theta1(k+1), theta2(k+1), gl1, gl2, km1, l2l1))*h;
end

deltax = theta2.*l2 - theta1.*l1;
PE = (0.5).*kspring.*((deltax).^2) + m1.*g.*l1.*(1-cos(theta1)) + m2.*g.*l2.*(1-cos(theta2));
KE = (0.5).*m1.*(l1^2).*((omega1).^2) + (0.5).*m2.*(l2^2).*((omega2).^2);
energy = KE + PE;

plot(tprime, theta1); hold on; plot(tprime,theta2);
title('Theta1 and Theta2 vs time');
ylabel('Theta(rad)');
xlabel('t(s)');
legend('Theta 1','Theta 2');

% subplot(2,2,1);
% plot(tprime, theta1); hold on; plot(tprime,theta2);
% title('Theta1 and Theta2 vs time');
% ylabel('Theta(rad)');
% xlabel('t(s)');
% legend('Theta 1','Theta 2');
% 
% subplot(2,2,2);
% plot(tprime, deltax);
% title('Spring elongation vs time');
% ylabel('delta X (m)');
% xlabel('t(s)');
% 
% subplot(2,2,3);
% plot(tprime, KE); hold on; plot(tprime,PE);
% title('KE and PE vs time');
% ylabel('Kinetic and potential energy(J)');
% xlabel('t(s)');
% legend('KE','PE');
% 
% subplot(2,2,4);
% plot(tprime, energy);
% title('Energy vs time');
% ylabel('Energy(J)');
% xlabel('t(s)');
% ylim([0 105]);