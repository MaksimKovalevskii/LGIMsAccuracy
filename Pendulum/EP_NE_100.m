% --- timing: tic/toc now brackets only the custom_rk4 integration call ---
global wx0 wy0 wz0 xd0 yd0 zd0;
% Parameters
b10=0.27059805007;
b20=0.27059805007;
b30=0;
b00=sqrt(1-b10^2-b20^2);

x0=6.*(b10.*b30+b00.*b20);
y0=6.*(b20.*b30-b00.*b10);
z0=3.*(1-2.*b10.^2-2.*b20.^2);

wx0=0;
wy0=0;
wz0=0;

A0 = [1-2*b20.^2-2*b30.^2 2*(b10.*b20-b00.*b30)  2*(b10.*b30+b00.*b20) ;
2*(b10.*b20+b00.*b30)  1-2*b10.^2-2*b30.^2 2*(b20.*b30-b00.*b10);
2*(b10.*b30-b00.*b20)  2*(b20.*b30+b00.*b10) 1-2*b20.^2-2*b10.^2];

gw0 = A0*[wx0 wy0 wz0]';

xd0=gw0(2)*z0 - gw0(3)*y0;
yd0=gw0(3)*x0 - gw0(1)*z0;
zd0=gw0(1)*y0 - gw0(2)*x0;

t_end = 100;  % End value for time
if ~exist('dt', 'var') || isempty(dt)
    dt = 0.1;   % Time step (default 100 ms)
end
Iqq=[3.0025 0 0; 0 3.0025 0;0 0 0.005];

% Initial conditions
initial_conditions = [x0; xd0; y0; yd0;z0; zd0; wx0; wy0; wz0;b00; b10; b20;b30];


% Time span
tspan = 0:dt:t_end;

function S = skew(v)
    S = [0 -v(3) v(2); v(3) 0 -v(1); -v(2) v(1) 0];
end
sj=[0;0;-3];
Sj_hat = skew(sj);

% Define the ODE system as a function
function [dydt, second_derivatives, dC,ddC, Ener, UC, theta] = odesystem(t, F)
global wx0 wy0 wz0 xd0 yd0 zd0;
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
sj=[0;0;-3];
Sj_hat = skew(sj);

Phi_pi= -A*Sj_hat;
Phi_rpi = [Phi_r Phi_pi];

    Iqq=[3.0025 0 0; 0 3.0025 0;0 0 0.005];
    M = [eye(3) zeros(3);
        zeros(3) Iqq];

Mat = [M Phi_rpi';
       Phi_rpi zeros(3)];

w_hat=skew(w_vec);
gamma = -A*w_hat*w_hat*sj;

Qv = -w_hat*Iqq*w_vec;

Forces=[0; 0; -9.8; Qv; gamma];
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
VelCon = Phi_r*[xd;yd;zd] + Phi_pi*w_vec;
dC=VelCon(1:3);
AccCon = Phi_r*NewRes(1:3) + Phi_pi*NewRes(4:6)-gamma;
ddC=AccCon(1:3);

qd= [xd yd zd wx wy wz];
qd0=[xd0 yd0 zd0 wx0 wy0 wz0];
Ener = 0.5 * (qd * M* qd')+9.8*z-9.8*3.*(1-2*0.27059805007.^2-2*0.27059805007.^2) - 0.5*(qd0 * M* qd0');

%Unit constraint
UC  = 1 - b0^2  - b1^2 - b2^2 - b3^2;
theta=2*acos(b0);
if theta > pi                     
    theta = 2 * acos(-b0);
end
end

% Custom Runge-Kutta 4th order method
function [t, F, second_derivatives,dC,ddC,Ener,UC,theta] = custom_rk4(odefun, tspan, y0)
    n = length(tspan);
    F_temp = zeros(length(y0), n);
    F_temp(:,1) = y0;

    second_derivatives = zeros(6, n);
    dC = zeros(3,n);
    ddC = zeros(3,n);
    Ener=zeros(1,n);
    UC = zeros(1,n);
    theta = zeros(1,n);

    for i = 1:(n-1)
        h = tspan(i+1) - tspan(i);
        ti = tspan(i);
        yi = F_temp(:,i);
        
        [k1, sd1,dC1,ddC1,Ener1, UC1,theta1] = odefun(ti, yi);
        [k2, ~] = odefun(ti + h/2, yi + h*k1/2);
        [k3, ~] = odefun(ti + h/2, yi + h*k2/2);
        [k4, ~] = odefun(ti + h, yi + h*k3);
        
        F_temp(:,i+1) = yi + (h/6)*(k1 + 2*k2 + 2*k3 + k4);
        second_derivatives(:,i) = sd1;  
        dC(:,i) = dC1;
        ddC(:,i) = ddC1;
        Ener(:,i) = Ener1;
        UC(:,i) = UC1;
        theta(:,i) = theta1;
    end
    [~, sd_final,dc_final,ddc_final,Ener_final,UC_final,theta_final] = odefun(tspan(end), F_temp(:,end));
    second_derivatives(:,end) = sd_final;
    dC(:,end) = dc_final;
    ddC(:,end) = ddc_final;
    Ener(:,end) = Ener_final;
    UC(:,end) = UC_final;
    theta(:,end) = theta_final;

    t = tspan;
    F = F_temp';
    second_derivatives = second_derivatives';
    dC=dC';
    ddC=ddC';
    Ener=Ener';
    UC=UC';    
    theta=theta'; 
end

% Solve the ODE using custom RK4
% --- timing: integration is run n_timing_repeats times; executionTime is the median ---
% (only the integration is repeated; setup, post-processing and save run once)
if ~exist('n_timing_repeats', 'var') || isempty(n_timing_repeats)
    n_timing_repeats = 1;   % a batch runner may raise this for a timing study
end
timing_samples = zeros(1, n_timing_repeats);
for timing_rep = 1:n_timing_repeats
    tic;
    [t, F, second_derivatives, dC, ddC,Ener, UC, theta] = custom_rk4(@odesystem, tspan, initial_conditions);
    timing_samples(timing_rep) = toc;  % integration only
end
executionTime = median(timing_samples);  % robust to run-to-run noise

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

A = [1-2*b2.^2-2*b3.^2 2*(b1.*b2-b0.*b3)  2*(b1.*b3+b0.*b2) ;
2*(b1.*b2+b0.*b3)  1-2*b1.^2-2*b3.^2 2*(b2.*b3-b0.*b1);
2*(b1.*b3-b0.*b2)  2*(b2.*b3+b0.*b1) 1-2*b2.^2-2*b1.^2];

% Calculate norms of Constraint violations' matrices for plotting
C1=x-6.*(b1.*b3+b0.*b2);
C2=y-6.*(b2.*b3-b0.*b1);
C3=z-3.*(1-2.*b1.^2-2.*b2.^2);
C4=b0.^2+b1.^2+b2.^2+b3.^2-1;
C=(C1+C2+C3+C4)/4; %C

%dC=zeros(size(t));
grdC=zeros(size(t));
grddC=zeros(size(t));
n=length(t);
for i=1:n
grdC(i)=(dC(i,1)+dC(i,2)+dC(i,3))/3;
grddC(i)=(ddC(i,1)+ddC(i,2)+ddC(i,3))/3;
end

UC2 = 1 - b0.^2  - b1.^2 - b2.^2 - b3.^2;

% executionTime captured at the custom_rk4 call above (timing harness)
if ~exist('save_filename', 'var') || isempty(save_filename)
    save_filename = sprintf('EP_NE_dt_%.1fms.mat', dt * 1000);
end
save(save_filename, '-v7.3');

% % Plotting results for x
% figure;
% subplot(3,1,1);
% plot(t,x,'b-','LineWidth',2);
% ylim([-3,3]);
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
% ylim([-20,20]);
% title('Rx''''(t)');
% xlabel('Time (s)');
% ylabel({'Acceleration X'; '(m/s^2)'});
% grid on;
% 
% % Plotting results for y
% figure;
% subplot(3,1,1);
% plot(t,y,'b-','LineWidth',2); 
% ylim([-5,5]);
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
% ylim([-20,20]);
% title('Ry''''(t)');
% xlabel('Time (s)');
% ylabel({'Acceleration Y'; '(m/s^2)'});
% grid on;
% 
% % Plotting results for z
% figure;
% subplot(3,1,1);
% plot(t,z,'b-','LineWidth',2); 
% ylim([-4,4]);
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
% ylim([-15,30]);
% title('Rz''''(t)');
% xlabel('Time (s)');
% ylabel({'Acceleration Z'; '(m/s^2)'});
% grid on;
% 
% % Plotting results for Angular velocities
% figure;
% subplot(3,1,1);
% plot(t,wx,'b-','LineWidth',2); 
% ylim([-3,3]);
% title('Wx');
% ylabel({'Ang velocity Wx'; 's-1'});
% grid on;
% 
% subplot(3,1,2);
% plot(t,wy,'r-','LineWidth',2); 
% ylim([-3,3]);
% title('Wy');
% ylabel({'Ang velocity Wy'; 's-1'});
% grid on;
% 
% subplot(3,1,3);
% plot(t,wz,'g-','LineWidth',2); 
% ylim([-1,1]);
% title('Wz');
% ylabel({'Ang velocity Wz';'s-1'});
% grid on;
% 
% % Plotting results for Angular accelerations
% figure;
% subplot(3,1,1);
% plot(t,AngAccX,'b-','LineWidth',2); 
% ylim([-3,3]);
% title('Wx''');
% ylabel({'Ang acceleration Wx';' s-2'});
% grid on;
% 
% subplot(3,1,2);
% plot(t,AngAccY,'r-','LineWidth',2); 
% ylim([-3,3]);
% title('Wy''');
% ylabel({'Ang acceleration Wy';' s-2'});
% grid on;
% 
% subplot(3,1,3);
% plot(t,AngAccZ,'g-','LineWidth',2); 
% ylim([-1,1]);
% title('Wz''');
% ylabel({'Ang acceleration Wz';' s-2'});
% grid on;
% 
% % Plotting results for theta
% figure;
% plot(t,theta,'b-','LineWidth',2); 
% ylim([0,7]);
% title('Derived from EP - Angle Theta');
% xlabel('Time (s)');
% ylabel('Norm of rotation vector Theta(t)');
% grid on;
% 
% % Plotting results for Constraint violation and Energy balance
% figure;
% subplot(4,1,1);
% plot(t,C,'b-','LineWidth',0.5); 
% ylim([-1*10^-3,1*10^-3]);
% title('âˆ£âˆ£Câˆ£âˆ£');
% ylabel('Violation');
% grid on;
% 
% subplot(4,1,2);
% plot(t,grdC,'g-','LineWidth',1); 
% ylim([-5*10^-4,5*10^-4]);
% title('âˆ£âˆ£C''âˆ£âˆ£');
% ylabel('Violation');
% grid on;
% 
% subplot(4,1,3);
% plot(t,grddC,'r-','LineWidth',1); 
% ylim([-5*10^-14,5*10^-14]);
% title('âˆ£âˆ£C''''âˆ£âˆ£');
% ylabel('Violation');
% grid on;
% 
% subplot(4,1,4);
% plot(t,Ener,'b-','LineWidth',1); 
% ylim([-1*10^-1,1*10^-1]);
% title('Total Energy');
% xlabel('Time (s)');
% ylabel('Energy balance');
% grid on;
% 
% % Plotting results for Unit constraints
% figure;
% subplot(2,1,1);
% plot(t,UC,'b-','LineWidth',0.5); 
% ylim([-1*10^-3,1*10^-3]);
% title('Unit costraint - EP Lie');
% ylabel('Violation');
% grid on;
% 
% subplot(2,1,2);
% plot(t,UC2,'b-','LineWidth',0.5); 
% ylim([-1*10^-3,1*10^-3]);
% title('Unit costraint - EP Lie');
% ylabel('Violation');
% grid on;
% 
% % Plotting results for Constraint violation and Energy balance
% figure;
% subplot(4,1,1);
% plot(t,abs(C),'b-','LineWidth',0.5);
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
