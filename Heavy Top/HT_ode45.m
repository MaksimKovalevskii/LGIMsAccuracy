tic;
% Parameters
b00=1;
b10=0;
b20=0;
b30=0;

x0=0;
y0=1;
z0=0;

wx0=0;
wy0=150;
wz0= -4.61538;

m0=15;
l0=1;

% Inertia tensor
Iqq = diag([0.234375; 0.46875; 0.234375]);

A0 = [1-2*b20.^2-2*b30.^2 2*(b10.*b20-b00.*b30)  2*(b10.*b30+b00.*b20) ;
2*(b10.*b20+b00.*b30)  1-2*b10.^2-2*b30.^2 2*(b20.*b30-b00.*b10);
2*(b10.*b30-b00.*b20)  2*(b20.*b30+b00.*b10) 1-2*b20.^2-2*b10.^2];

gw0 = A0*[wx0 wy0 wz0]';

xd0=gw0(2)*z0 - gw0(3)*y0;
yd0=gw0(3)*x0 - gw0(1)*z0;
zd0=gw0(1)*y0 - gw0(2)*x0;

% Initial conditions
initial_conditions = [x0; xd0; y0; yd0;z0; zd0; wx0; wy0; wz0;b00; b10; b20;b30];

t_end = 100;  % End value for time
% Time span
tspan = 0:0.0001:t_end;

function S = skew(v)
    S = [0 -v(3) v(2); v(3) 0 -v(1); -v(2) v(1) 0];
end
sj=[0;-l0;0];
Sj_hat = skew(sj);

r0=A0*sj+[0;1;0];
w_vec0=[wx0;
wy0;
wz0];

    M0 = [m0*eye(3) zeros(3);
        zeros(3) Iqq];
qd0= [xd0 yd0 zd0 wx0 wy0 wz0];
    Ener0 = 0.5 * (qd0 * M0* qd0') +m0*9.81*z0;

% Define the ODE system as a function
function [dydt, second_derivatives, C, dC,ddC, Ener, UC] = odesystem(t, F)
%global wx0 wy0 wz0 xd0 yd0 zd0 l0 m0;
    x = F(1);
    xd = F(2);
    y = F(3);
    yd = F(4);
    z = F(5);
    zd = F(6);
    wx = F(7);
    wy = F(8);
    wz = F(9);
    b0 = F(10);
    b1 = F(11);
    b2 = F(12);
    b3 = F(13);

    l0=1;
    m0=15;
Gint=[-b1 b0  b3 -b2;
-b2 -b3  b0 b1;
-b3 b2  -b1 b0];

A = [1-2*b2^2-2*b3^2 2*(b1*b2-b0*b3)  2*(b1*b3+b0*b2) ;
2*(b1*b2+b0*b3)  1-2*b1^2-2*b3^2 2*(b2*b3-b0*b1);
2*(b1*b3-b0*b2)  2*(b2*b3+b0*b1) 1-2*b2^2-2*b1^2];

w_vec=[wx;
wy;
wz];

b_ans=0.5*Gint'*w_vec;

    Phi_rpi = zeros (3,6);
Phi_r= eye(3);
sj=[0;-l0;0];

Sj_hat = skew(sj);

Phi_pi= -A*Sj_hat;
Phi_rpi = [Phi_r Phi_pi];

Iqq = diag([0.234375; 0.46875; 0.234375]);
    M = [m0*eye(3) zeros(3);
        zeros(3) Iqq];

    Mat = [M Phi_rpi';
       Phi_rpi zeros(3)];

    w_hat=skew(w_vec);

    tempC=[x;y;z]+A*sj;
C=tempC(1:3); %position level
    VelCon = Phi_r*[xd;yd;zd] + Phi_pi*w_vec;
dC=VelCon(1:3);
%gamma = -A*w_hat*w_hat*sj;
gamma = -A*w_hat*w_hat*sj -20*dC-100*C;

% Gyroscopic torques 
Qv_gyro = -w_hat*Iqq*w_vec;

Forces=[0; 0; -m0 *9.81; Qv_gyro; gamma];
NewRes=Mat\Forces;

 dydt = zeros(13,1);
    dydt(1) = xd; % x' = xd
    dydt(2) = NewRes(1); % xd' (x'')
    dydt(3) = yd; % y' = yd
    dydt(4) = NewRes(2);  % yd' (y'')  
    dydt(5) = zd; % z' = zd
    dydt(6) = NewRes(3); % zd' (z'')  
    dydt(7) = NewRes(4); % wx'
    dydt(8) = NewRes(5); % wy' 
    dydt(9) = NewRes(6); % wz' 
    dydt(10) = b_ans(1) ; % b0
    dydt(11) = b_ans(2); % b1 
    dydt(12) = b_ans(3); % b2
    dydt(13) = b_ans(4); % b3

second_derivatives = NewRes(1:6);
tempC=[x;y;z]+A*sj;
C=tempC(1:3); %position level
VelCon = Phi_r*[xd;yd;zd] + Phi_pi*w_vec;
dC=VelCon(1:3);
AccCon = Phi_r*NewRes(1:3) + Phi_pi*NewRes(4:6)+A*w_hat*w_hat*sj;
ddC=AccCon(1:3);

qd= [xd yd zd wx wy wz];
Ener = 0.5 * (qd * M* qd') +m0*9.81*z - 5.435696790865547e+03;

%Unit constraint
UC  = 1 - b0^2  - b1^2 - b2^2 - b3^2;
end

options = odeset('RelTol', 1e-14, 'AbsTol', 1e-14);
% Solve the ODE using ode45 with the anonymous function
%[t, F, second_derivatives, C, dC, ddC,Ener, UC] = ode45(@odesystem, tspan, initial_conditions, options);
[t, F] = ode45(@odesystem, tspan, initial_conditions, options);

n = length(t);
second_derivatives = zeros(n,6);
C  = zeros(n,3);   dC  = zeros(n,3);   ddC  = zeros(n,3);
Ener = zeros(n,1); UC  = zeros(n,1); 

for i = 1:n
    [~, second_derivatives(i,:), C(i,:), dC(i,:), ddC(i,:), Ener(i), UC(i)] = ...
        odesystem(t(i), F(i,:)');
end

% Extract results from Y matrix
x = F(:,1);
xd = F(:,2);
y = F(:,3);
yd = F(:,4);
z = F(:,5);
zd = F(:,6);
wx = F(:,7);
wy = F(:,8);
wz = F(:,9);


b0 = F(:,10);
b1 = F(:,11);
b2 = F(:,12);
b3 = F(:,13);

x_double_prime = second_derivatives(:,1);
y_double_prime = second_derivatives(:,2);
z_double_prime = second_derivatives(:,3);
AngAccX = second_derivatives(:,4);
AngAccY = second_derivatives(:,5);
AngAccZ = second_derivatives(:,6);

%constraint violations
grC=zeros(size(t));
grdC=zeros(size(t));
grddC=zeros(size(t));
n=length(t);
for i=1:n
grC(i)=(C(i,1)+C(i,2)+C(i,3))/3;  
grdC(i)=(dC(i,1)+dC(i,2)+dC(i,3))/3;
grddC(i)=(ddC(i,1)+ddC(i,2)+ddC(i,3))/3;
end

executionTime = toc;
if ~exist('save_filename','var') || isempty(save_filename)
    save_filename = 'HTRef.mat';
end
save(save_filename, '-v7.3');

% Plotting results for x y z
figure;
subplot(3,1,1);
plot(t,x,'b--','LineWidth',1);
hold on;
plot(t,y,'r-','LineWidth',1); 
hold on;
plot(t,z,'g-','LineWidth',1); 
ylim([-1,1]);
title('Rx, Ry, Rz');
ylabel('Displacements, m');
grid on;

subplot(3,1,2);
plot(t,xd,'b--','LineWidth',1); 
hold on;
plot(t,yd,'r-','LineWidth',1); 
hold on;
plot(t,zd,'g-','LineWidth',1); 
ylim([-10,10]);
title('Rx'', Ry'', Rz'' (t)');
ylabel({'Velocities'; '(m/s)'});
grid on;

subplot(3,1,3);
plot(t,x_double_prime,'b--','LineWidth',1); 
hold on;
plot(t,y_double_prime,'r-','LineWidth',1); 
hold on;
plot(t,z_double_prime,'g-','LineWidth',1); 
ylim([-50,45]);
title('Rx'''', Ry'''', Rz''''(t)');
xlabel('Time (s)');
ylabel({'Accelerations'; '(m/s^2)'});
legend ('Rx','Ry','Rz');
grid on;

% Plotting results for x
figure;
subplot(3,1,1);
plot(t,x,'b-','LineWidth',2);
ylim([-1,1]);
title('Rx(t)');
ylabel('Displacement X, m');
grid on;

subplot(3,1,2);
plot(t,xd,'r-','LineWidth',2); 
ylim([-10,10]);
title('Rx''(t)');
ylabel({'Velocity X'; '(m/s)'});
grid on;

subplot(3,1,3);
plot(t,x_double_prime,'g-','LineWidth',2); 
ylim([-50,45]);
title('Rx''''(t)');
xlabel('Time (s)');
ylabel({'Acceleration X'; '(m/s^2)'});
grid on;

% Plotting results for y
figure;
subplot(3,1,1);
plot(t,y,'b-','LineWidth',2); 
ylim([-1,1]);
title('Ry(t)');
ylabel('Displacement Y, m');
grid on;

subplot(3,1,2);
plot(t,yd,'r-','LineWidth',2); 
ylim([-10,10]);
title('Ry''(t)');
ylabel({'Velocity Y'; '(m/s)'});
grid on;

subplot(3,1,3);
plot(t,y_double_prime,'g-','LineWidth',2); 
ylim([-50,50]);
title('Ry''''(t)');
xlabel('Time (s)');
ylabel({'Acceleration Y'; '(m/s^2)'});
grid on;

% Plotting results for z
figure;
subplot(3,1,1);
plot(t,z,'b-','LineWidth',2); 
ylim([-1.5,0.5]);
title('Rz(t)');
ylabel('Displacement Z, m');
grid on;

subplot(3,1,2);
plot(t,zd,'r-','LineWidth',2); 
ylim([-10,10]);
title('Rz''(t)');
ylabel({'Velocity Z'; '(m/s)'});
grid on;

subplot(3,1,3);
plot(t,z_double_prime,'g-','LineWidth',2); 
ylim([-40,50]);
title('Rz''''(t)');
xlabel('Time (s)');
ylabel({'Acceleration Z'; '(m/s^2)'});
grid on;

% Plotting results for Angular velocities
figure;
subplot(3,1,1);
plot(t,wx,'b-','LineWidth',2); 
ylim([-7,7]);
title('Wx');
ylabel({'Ang velocity Wx'; 's-1'});
grid on;

subplot(3,1,2);
plot(t,wy,'r-','LineWidth',2); 
ylim([0,200]);
title('Wy');
ylabel({'Ang velocity Wy'; 's-1'});
grid on;

subplot(3,1,3);
plot(t,wz,'g-','LineWidth',2); 
ylim([-7,7]);
title('Wz');
ylabel({'Ang velocity Wz';'s-1'});
grid on;

%Plotting results for Angular accelerations
figure;
subplot(3,1,1);
plot(t,AngAccX,'b-','LineWidth',2); 
ylim([-1000,1000]);
title('Wx''');
ylabel({'Ang acceleration Wx';' s-2'});
grid on;

subplot(3,1,2);
plot(t,AngAccY,'r-','LineWidth',2); 
ylim([-3,3]);
title('Wy''');
ylabel({'Ang acceleration Wy';' s-2'});
grid on;

subplot(3,1,3);
plot(t,AngAccZ,'g-','LineWidth',2); 
ylim([-1000,1000]);
title('Wz''');
ylabel({'Ang acceleration Wz';' s-2'});
grid on;


% Plotting results for Constraint violation and Energy balance
figure;
subplot(4,1,1);
plot(t,grC,'b-','LineWidth',0.5); 
ylim([-1*10^-12,1*10^-12]);
title('∣∣C∣∣');
ylabel('Violation');
grid on;

subplot(4,1,2);
plot(t,grdC,'g-','LineWidth',1); 
ylim([-1*10^-12,1*10^-12]);
title('∣∣C''∣∣');
ylabel('Violation');
grid on;

subplot(4,1,3);
plot(t,grddC,'r-','LineWidth',1); 
ylim([-5*10^-12,5*10^-12]);
title('∣∣C''''∣∣');
ylabel('Violation');
grid on;

subplot(4,1,4);
plot(t,Ener,'b-','LineWidth',1); 
ylim([-1*10^-9,1*10^-9]);
title('Total Energy');
xlabel('Time (s)');
ylabel('Energy balance');
grid on;

% Plotting results for Unit constraints
figure;
plot(t,UC,'b-','LineWidth',0.5); 
ylim([-1*10^-13,1*10^-13]);
title('Unit costraint - EP Lie');
ylabel('Violation');
grid on;

% Plotting results for Constraint violation and Energy balance
figure;
subplot(4,1,1);
plot(t,abs(grC),'b-','LineWidth',0.5);
set(gca, 'YScale', 'log');
	set(gca, 'YMinorTick', 'on');
title('||C|| (Log scale)');
ylabel('Violation');
grid on;

subplot(4,1,2);
plot(t,abs(grdC),'g-','LineWidth',1);
set(gca, 'YScale', 'log');
	set(gca, 'YMinorTick', 'on');
title('||C''|| (Log scale)');
ylabel('Violation');
grid on;

subplot(4,1,3);
plot(t,abs(grddC),'r-','LineWidth',1);
set(gca, 'YScale', 'log');
	set(gca, 'YMinorTick', 'on');
title('||C''''|| (Log scale)');
ylabel('Violation');
grid on;

subplot(4,1,4);
plot(t,abs(Ener),'b-','LineWidth',1);
set(gca, 'YScale', 'log');
	set(gca, 'YMinorTick', 'on');
title('Total Energy (Log scale)');
xlabel('Time (s)');
ylabel('Energy balance');
grid on;

figure;
subplot(3,1,1);
plot(t,x,'b-','LineWidth',2);
ylim([-1,1]);
title('$R_X(t)$', 'Interpreter', 'latex');
ylabel('Displacement $X$ (m)', 'Interpreter', 'latex');
grid on;
set(gca, 'FontName','Times', 'FontSize',11);

subplot(3,1,2);
plot(t,xd,'r-','LineWidth',2); 
ylim([-10,10]);
title('$\dot{R}_X(t)$', 'Interpreter', 'latex');
ylabel('Velocity $X$ (m/s)', 'Interpreter', 'latex');
grid on;
set(gca, 'FontName','Times', 'FontSize',11);

subplot(3,1,3);
plot(t,x_double_prime,'g-','LineWidth',2); 
ylim([-50,45]);
title('$\ddot{R}_X(t)$', 'Interpreter', 'latex');
xlabel('Time (s)', 'Interpreter','latex');
ylabel('Acceleration $X$ (m/s$^2$)', 'Interpreter', 'latex');
grid on;
set(gca, 'FontName','Times', 'FontSize',11);


figure;

% --- Subplot 1: wy ---
subplot(2,1,1);
plot(t, wy, 'r-', 'LineWidth', 1.5);
xlim([0 2]);
ylim([0 200]);
ylabel({'Angular velocity $\omega_y$ (s$^{-1}$)'}, 'Interpreter','latex');
grid on;
set(gca, 'FontName','Times', 'FontSize',12);

% --- Subplot 2: wx and wz ---
subplot(2,1,2);
plot(t, wx, 'b-', 'LineWidth', 1.5); hold on;
plot(t, wz, 'g-', 'LineWidth', 1.5);
xlim([0 2]);
ylim([-7 7]);
xlabel('Time (s)', 'Interpreter','latex');
ylabel({'Angular velocity $\omega_x,\ \omega_z$ (s$^{-1}$)'}, 'Interpreter','latex');
legend({'$\omega_x$', '$\omega_z$'}, 'Interpreter','latex', ...
       'Location','best', 'Box','off');
grid on;
set(gca, 'FontName','Times', 'FontSize',12);

