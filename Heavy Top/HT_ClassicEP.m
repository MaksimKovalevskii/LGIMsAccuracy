tic;
% Parameters
b00=1;
b10=0;
b20=0;
b30=0;

x0=0;
y0=1;
z0=0;

%wx0=0;
%wy0=0;
%wz0=0;
wx0=0;
wy0=150;
wz0= -4.61538;
w_vec0=[wx0;
wy0;
wz0];

m0=15;
l0=1;

% Inertia tensor
Iqq = diag([0.234375; 0.46875; 0.234375]);

A0 = [1-2*b20.^2-2*b30.^2 2*(b10.*b20-b00.*b30)  2*(b10.*b30+b00.*b20) ;
2*(b10.*b20+b00.*b30)  1-2*b10.^2-2*b30.^2 2*(b20.*b30-b00.*b10);
2*(b10.*b30-b00.*b20)  2*(b20.*b30+b00.*b10) 1-2*b20.^2-2*b10.^2];

Gint0=[-b10 b00  b30 -b20;
-b20 -b30  b00 b10;
-b30 b20  -b10 b00];

gw0 = A0*[wx0 wy0 wz0]';

xd0=gw0(2)*z0 - gw0(3)*y0;
yd0=gw0(3)*x0 - gw0(1)*z0;
zd0=gw0(1)*y0 - gw0(2)*x0;

bd_init=0.5*Gint0'*w_vec0;
%bd_init=0.25*Gint0'*w_vec0;
%bd_init=1*Gint0'*w_vec0;

% Initial conditions
initial_conditions = [x0;  y0; z0; xd0;yd0;zd0;b00; b10; b20;b30;bd_init];

sj=[0;-l0;0];
r0=A0*sj+[0;1;0]

t_end = 20;  % End value for time
% Time step (seconds)
if ~exist('dt','var') || isempty(dt)
    dt = 0.003;
end
% Time span
tspan = 0:dt:t_end;
function S = skew(v)
    S = [0 -v(3) v(2); v(3) 0 -v(1); -v(2) v(1) 0];
end

    M0 = [m0*eye(3) zeros(3);
        zeros(3) Iqq];
qd0= [xd0 yd0 zd0 wx0 wy0 wz0];
    Ener0 = 0.5 * (qd0 * M0* qd0') +m0*9.81*z0;

% Define the ODE system as a function
function [dydt, second_derivatives, C, dC,ddC, Ener, UC, theta,omegaint] = odesystem(t, F)
    x = F(1);
    y = F(2);
    z = F(3);
    xd = F(4);
    yd = F(5);
    zd = F(6);
    b0 = F(7);
    b1 = F(8);
    b2 = F(9);
    b3 = F(10);
    b0d = F(11);
    b1d = F(12);
    b2d = F(13);
    b3d = F(14);

betasd=[b0d b1d b2d b3d];

Gint=GmatrixEP(x, y, z, b0, b1, b2, b3);
%Gint=2.*[-b1 b0  b3 -b2;
%-b2 -b3  b0 b1;
%-b3 b2  -b1 b0];
%disp('GmatrixEP =');  disp(Gtest);
%disp('Ghardcoded =');  disp(Gint);
%difference = Gint - Gtest;
%disp('Difference =');  disp(difference);


%G2=GmatrixEP(x, y, z, b0d, b1d, b2d, b3d);
G2=2.*[-b1d b0d  b3d -b2d;
-b2d -b3d  b0d b1d;
-b3d b2d  -b1d b0d];

A = [1-2*b2^2-2*b3^2 2*(b1*b2-b0*b3)  2*(b1*b3+b0*b2) ;
2*(b1*b2+b0*b3)  1-2*b1^2-2*b3^2 2*(b2*b3-b0*b1);
2*(b1*b3-b0*b2)  2*(b2*b3+b0*b1) 1-2*b2^2-2*b1^2];

l0=1;
m0=15;
Iqq = diag([0.234375; 0.46875; 0.234375]);
sj=[0;-l0;0];

Mqq=Gint'*Iqq*Gint;

M = [m0*eye(3) zeros(3,4);
    zeros(4,3) Mqq];

Cq = [1 0 0 2.*b3 -2.*b2 -2.*b1 2.*b0 ;
    0 1 0 0 4.*b1 0 4.*b3 ;
    0 0 1 -2.*b1 -2.*b0 -2.*b3 -2.*b2 ;
    0 0 0 2.*b0 2.*b1 2.*b2 2.*b3 ];

Mat = [M Cq';
       Cq zeros(4)];

omegaint=Gint*betasd';
%Qv1=cross(omegaint,Iqq*omegaint);
Qv1=cross(omegaint,Iqq*omegaint)+Iqq*G2*betasd';
result = Iqq*G2*betasd';
%fprintf('%g ', result');  % Transpose for proper display
%fprintf('\n');  % Add newline at the end
Qv=-Gint'*Qv1; 

Qc = [-4.*b0d*b3d + 4.*b1d*b2d;
    -4.*b1d*b1d-4.*b3d*b3d;
    4.*b0d*b1d + 4.*b2d*b3d;
    -2.*b1d*b1d-2.*b2d*b2d-2.*b0d*b0d-2.*b3d*b3d ];

Forces=[0; 0; -m0 *9.81; Qv; Qc];

NewRes=Mat\Forces;

    dydt = zeros(14,1);
    dydt(1) = xd; % x' = xd
    dydt(2) = yd; % y' = yd
    dydt(3) = zd; % z' = zd
    dydt(4) = NewRes(1); % xd' (x'')
    dydt(5) = NewRes(2);  % yd' (y'')  
    dydt(6) = NewRes(3); % zd' (z'') 
    dydt(7) = b0d; 
    dydt(8) = b1d; 
    dydt(9) = b2d; 
    dydt(10) = b3d; 
    dydt(11) = NewRes(4); 
    dydt(12) = NewRes(5);  
    dydt(13) = NewRes(6);  
    dydt(14) = NewRes(7); 

second_derivatives = NewRes(1:7);
Consts=[x;y;z]+A*sj;
tempC = [Consts;
    b0.^2+b1.^2+b2.^2+b3.^2-1];
C=tempC(1:4); %position level

qd= [xd yd zd b0d b1d b2d b3d];
VelCon=Cq*qd';
dC=VelCon(1:4);

AccCon=Cq*second_derivatives-Qc;
ddC=AccCon(1:4);

wx=omegaint(1);
wy=omegaint(2);
wz=omegaint(3);
%qd= [xd yd zd wx wy wz];
Ener = 0.5.*(qd * M* qd') +m0*9.81*z - 5.435696790865547e+03;
%Ener = m0*9.81*z - 5.435696790865547e+03;

%Unit constraint
UC  = 1 - b0^2  - b1^2 - b2^2 - b3^2;
theta=2*acos(b0);
if theta > pi                     % Flip sign if needed
    theta = 2 * acos(-b0);
end
end

% Custom Runge-Kutta 4th order method
function [t, F, second_derivatives,C,dC,ddC,Ener,UC,theta,omegaint] = custom_rk4(odefun, tspan, y0)
    n = length(tspan);
    F_temp = zeros(length(y0), n);
    F_temp(:,1) = y0;

    second_derivatives = zeros(7, n);
    C = zeros(4,n);
    dC = zeros(4,n);
    ddC = zeros(4,n);
    Ener=zeros(1,n);
    UC = zeros(1,n);
    theta = zeros(1,n);
    omegaint=zeros(3,n);

    for i = 1:(n-1)
        h = tspan(i+1) - tspan(i);
        ti = tspan(i);
        yi = F_temp(:,i);
        
        [k1, sd1,C1,dC1,ddC1,Ener1, UC1,theta1,omegaint1] = odefun(ti, yi);
        [k2, ~] = odefun(ti + h/2, yi + h*k1/2);
        [k3, ~] = odefun(ti + h/2, yi + h*k2/2);
        [k4, ~] = odefun(ti + h, yi + h*k3);
        
        F_temp(:,i+1) = yi + (h/6)*(k1 + 2*k2 + 2*k3 + k4);
        second_derivatives(:,i) = sd1;  
        C(:,i) = C1;
        dC(:,i) = dC1;
        ddC(:,i) = ddC1;
        Ener(:,i) = Ener1;
        UC(:,i) = UC1;
        theta(:,i) = theta1;
        omegaint(:,i) = omegaint1;
    end
    [~, sd_final,c_final, dc_final,ddc_final,Ener_final,UC_final,theta_final,omegaint_final] = odefun(tspan(end), F_temp(:,end));
    second_derivatives(:,end) = sd_final;
    C(:,end) = c_final;
    dC(:,end) = dc_final;
    ddC(:,end) = ddc_final;
    Ener(:,end) = Ener_final;
    UC(:,end) = UC_final;
    theta(:,end) = theta_final;
    omegaint(:,end) = omegaint_final;

    t = tspan;
    F = F_temp';
    second_derivatives = second_derivatives';
    C=C';
    dC=dC';
    ddC=ddC';
    Ener=Ener';
    UC=UC';    
    theta=theta'; 
    omegaint=omegaint';
end

% Solve the ODE using custom RK4
[t, F, second_derivatives, C, dC, ddC,Ener, UC, theta,omegaint] = custom_rk4(@odesystem, tspan, initial_conditions);

% Extract results from F matrix
    x = F(:,1);
    y = F(:,2);
    z = F(:,3);
    xd = F(:,4);
    yd = F(:,5);
    zd = F(:,6);
    b0 = F(:,7);
    b1 = F(:,8);
    b2 = F(:,9);
    b3 = F(:,10);
    b0d = F(:,11);
    b1d = F(:,12);
    b2d = F(:,13);
    b3d = F(:,14);

x_double_prime = second_derivatives(:,1);
y_double_prime = second_derivatives(:,2);
z_double_prime = second_derivatives(:,3);
b0_double_prime = second_derivatives(:,4);
b1_double_prime = second_derivatives(:,5);
b2_double_prime = second_derivatives(:,6);
b3_double_prime = second_derivatives(:,7);

wx=omegaint(:,1);
wy=omegaint(:,2);
wz=omegaint(:,3);

%constraint violations
grC=zeros(size(t));
grdC=zeros(size(t));
grddC=zeros(size(t));
n=length(t);
for i=1:n
grC(i)=(C(i,1)+C(i,2)+C(i,3)+C(i,4))/4;  
grdC(i)=(dC(i,1)+dC(i,2)+dC(i,3)+dC(i,4))/4;
grddC(i)=(ddC(i,1)+ddC(i,2)+ddC(i,3)+ddC(i,4))/4;
end

executionTime = toc;
if ~exist('save_filename','var') || isempty(save_filename)
    save_filename = sprintf('HTClassicEP_dt_%.2fms.mat', dt*1000);
end
save(save_filename, '-v7.3');


% % Plotting results for Angular velocities
% figure;
% subplot(3,1,1);
% plot(t,wx,'b-','LineWidth',2); 
% ylim([-7,7]);
% title('Wx');
% ylabel({'Ang velocity Wx'; 's-1'});
% grid on;
% 
% subplot(3,1,2);
% plot(t,wy,'r-','LineWidth',2); 
% ylim([0,200]);
% title('Wy');
% ylabel({'Ang velocity Wy'; 's-1'});
% grid on;
% 
% subplot(3,1,3);
% plot(t,wz,'g-','LineWidth',2); 
% ylim([-7,7]);
% title('Wz');
% ylabel({'Ang velocity Wz';'s-1'});
% grid on;
% 
% 
% % Plotting results for x y z
% figure;
% subplot(3,1,1);
% plot(t,x,'b--','LineWidth',1);
% hold on;
% plot(t,y,'r-','LineWidth',1); 
% hold on;
% plot(t,z,'g-','LineWidth',1); 
% ylim([-1,1]);
% title('Rx, Ry, Rz');
% ylabel('Displacements, m');
% grid on;
% 
% subplot(3,1,2);
% plot(t,xd,'b--','LineWidth',1); 
% hold on;
% plot(t,yd,'r-','LineWidth',1); 
% hold on;
% plot(t,zd,'g-','LineWidth',1); 
% ylim([-10,10]);
% title('Rx'', Ry'', Rz'' (t)');
% ylabel({'Velocities'; '(m/s)'});
% grid on;
% 
% subplot(3,1,3);
% plot(t,x_double_prime,'b--','LineWidth',1); 
% hold on;
% plot(t,y_double_prime,'r-','LineWidth',1); 
% hold on;
% plot(t,z_double_prime,'g-','LineWidth',1); 
% ylim([-50,45]);
% title('Rx'''', Ry'''', Rz''''(t)');
% xlabel('Time (s)');
% ylabel({'Accelerations'; '(m/s^2)'});
% legend ('Rx','Ry','Rz');
% grid on;
% 
% % Plotting results for ep
% figure;
% plot(t,b0,'b--','LineWidth',1);
% hold on;
% plot(t,b1,'r-','LineWidth',1); 
% hold on;
% plot(t,b2,'g-','LineWidth',1); 
% hold on;
% plot(t,b3,'c--','LineWidth',1); 
% xlim([0,0.2]);
% ylim([-1,1]);
% title('Euler parameters');
% ylabel('Quaternion components');
% legend ('b0','b1','b2','b3');
% grid on;
% 
% % Plotting results for x
% figure;
% subplot(3,1,1);
% plot(t,x,'b-','LineWidth',2);
% ylim([-1,1]);
% title('Rx(t)');
% ylabel('Displacement X, m');
% grid on;
% 
% subplot(3,1,2);
% plot(t,xd,'r-','LineWidth',2); 
% ylim([-10,10]);
% title('Rx''(t)');
% ylabel({'Velocity X'; '(m/s)'});
% grid on;
% 
% subplot(3,1,3);
% plot(t,x_double_prime,'g-','LineWidth',2); 
% ylim([-50,45]);
% title('Rx''''(t)');
% xlabel('Time (s)');
% ylabel({'Acceleration X'; '(m/s^2)'});
% grid on;
% 
% % Plotting results for y
% figure;
% subplot(3,1,1);
% plot(t,y,'b-','LineWidth',2); 
% ylim([-1,1]);
% title('Ry(t)');
% ylabel('Displacement Y, m');
% grid on;
% 
% subplot(3,1,2);
% plot(t,yd,'r-','LineWidth',2); 
% ylim([-10,10]);
% title('Ry''(t)');
% ylabel({'Velocity Y'; '(m/s)'});
% grid on;
% 
% subplot(3,1,3);
% plot(t,y_double_prime,'g-','LineWidth',2); 
% ylim([-50,50]);
% title('Ry''''(t)');
% xlabel('Time (s)');
% ylabel({'Acceleration Y'; '(m/s^2)'});
% grid on;
% 
% % Plotting results for z
% figure;
% subplot(3,1,1);
% plot(t,z,'b-','LineWidth',2); 
% ylim([-1.5,0.5]);
% title('Rz(t)');
% ylabel('Displacement Z, m');
% grid on;
% 
% subplot(3,1,2);
% plot(t,zd,'r-','LineWidth',2); 
% ylim([-10,10]);
% title('Rz''(t)');
% ylabel({'Velocity Z'; '(m/s)'});
% grid on;
% 
% subplot(3,1,3);
% plot(t,z_double_prime,'g-','LineWidth',2); 
% ylim([-40,50]);
% title('Rz''''(t)');
% xlabel('Time (s)');
% ylabel({'Acceleration Z'; '(m/s^2)'});
% grid on;
% 
% % Plotting results for theta
% figure;
% plot(t,theta,'b-','LineWidth',2); 
% xlim([0,2]);
% ylim([0,4]);
% title('Derived from EP - Angle Theta');
% xlabel('Time (s)');
% ylabel('Norm of rotation vector Theta(t)');
% grid on;
% 
% 
% % Plotting results for Constraint violation and Energy balance
% figure;
% subplot(4,1,1);
% plot(t,grC,'b-','LineWidth',0.5); 
% ylim([-1*10^-4,1*10^-4]);
% title('∣∣C∣∣');
% ylabel('Violation');
% grid on;
% 
% subplot(4,1,2);
% plot(t,grdC,'g-','LineWidth',1); 
% ylim([-1*10^-4,1*10^-4]);
% title('∣∣C''∣∣');
% ylabel('Violation');
% grid on;
% 
% subplot(4,1,3);
% plot(t,grddC,'r-','LineWidth',1); 
% ylim([-1*10^-12,1*10^-12]);
% title('∣∣C''''∣∣');
% ylabel('Violation');
% grid on;
% 
% subplot(4,1,4);
% plot(t,Ener,'b-','LineWidth',1); 
% ylim([-1*10^-4,1*10^-4]);
% title('Total Energy');
% xlabel('Time (s)');
% ylabel('Energy balance');
% grid on;
% 
% % Plotting results for Unit constraints
% figure;
% plot(t,UC,'b-','LineWidth',0.5); 
% ylim([-1*10^-5,1*10^-5]);
% title('Unit costraint - Classic EP');
% ylabel('Violation');
% grid on;
% 
% % Plotting results for Constraint violation and Energy balance
% figure;
% subplot(4,1,1);
% plot(t,abs(grC),'b-','LineWidth',0.5);
% set(gca, 'YScale', 'log');
% 	set(gca, 'YMinorTick', 'on');
% title('||C|| (Log scale)');
% ylabel('Violation');
% grid on;
% 
% subplot(4,1,2);
% plot(t,abs(grdC),'g-','LineWidth',1);
% set(gca, 'YScale', 'log');
% 	set(gca, 'YMinorTick', 'on');
% title('||C''|| (Log scale)');
% ylabel('Violation');
% grid on;
% 
% subplot(4,1,3);
% plot(t,abs(grddC),'r-','LineWidth',1);
% set(gca, 'YScale', 'log');
% 	set(gca, 'YMinorTick', 'on');
% title('||C''''|| (Log scale)');
% ylabel('Violation');
% grid on;
% 
% subplot(4,1,4);
% plot(t,abs(Ener),'b-','LineWidth',1);
% set(gca, 'YScale', 'log');
% 	set(gca, 'YMinorTick', 'on');
% title('Total Energy (Log scale)');
% xlabel('Time (s)');
% ylabel('Energy balance');
% grid on;