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
initial_conditions = [x0; y0; z0;b00;b10;b20;b30;xd0; yd0; zd0; wx0;wy0;wz0];

t_end = 20;  % End value for time
% Time step (seconds)
if ~exist('dt','var') || isempty(dt)
    dt = 0.01;
end
% Time span
tspan = 0:dt:t_end;

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

  function qf = quat_mult(q1, q2)   
        % Quaternion multiplication
        qf = [q1(1)*q2(1) - dot(q1(2:4), q2(2:4));
             q1(1)*q2(2:4) + q2(1)*q1(2:4) + cross(q1(2:4), q2(2:4))]; 
    end

function q_new = LGT(q_old, u)
phi = 0.5*norm(u);
if phi<1e-2
rot_exp=[cos(phi); 0.5*u*(1-phi^2/6+phi^4/120-phi^6/5040)];
else
rot_exp=[cos(phi); 0.5*u*sin(phi)/phi];
end 
    q_new = quat_mult(q_old, rot_exp);
end
%yi(10:12)
function dexpInv = invdexp(w,k)
 phi = norm(k);

if phi < 1e-2
    dexpInv = w + 0.5*skew(k)*w + ...
                    (1/12 + phi^2/720+phi^4/30240+phi^6/1209600)*((skew(k))^2)*w;
else
dexpInv= w + 0.5*skew(k)*w + ...
                    (1 - (phi/2)*cot(phi/2))/phi^2 * ((skew(k))^2)*w;
end
end  

% Define the ODE system as a function
function [dydt, second_derivatives, C, dC,ddC, Ener, UC] = odesystem(t, F)
    x = F(1);
    y = F(2);
    z = F(3);
    b0 = F(4);
    b1 = F(5);
    b2 = F(6);
    b3 = F(7);    
    xd = F(8);
    yd = F(9);
    zd = F(10);
    wx = F(11);
    wy = F(12);
    wz = F(13);

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
gamma = -A*w_hat*w_hat*sj;

Qv_gyro = -w_hat*Iqq*w_vec;

Forces=[0; 0; -m0 *9.81; Qv_gyro; gamma];
NewRes=Mat\Forces;

 dydt = zeros(13,1);
    dydt(1) = xd; % x' = xd
    dydt(2) = yd; % y' = yd
    dydt(3) = zd; % z' = zd
    dydt(4) = b_ans(1) ; % 
    dydt(5) = b_ans(2); %  
    dydt(6) = b_ans(3); % 
    dydt(7) = b_ans(4);
    dydt(8) = NewRes(1); % xd' (x'')
    dydt(9) = NewRes(2);  % yd' (y'')  
    dydt(10) = NewRes(3); % zd' (z'') 
    dydt(11) = NewRes(4); % wx'
    dydt(12) = NewRes(5); % wy' 
    dydt(13) = NewRes(6); % wz' 

second_derivatives = NewRes(1:6);
tempC=[x;y;z]+A*sj;
C=tempC(1:3); %position level
VelCon = Phi_r*[xd;yd;zd] + Phi_pi*w_vec;
dC=VelCon(1:3);
AccCon = Phi_r*NewRes(1:3) + Phi_pi*NewRes(4:6)-gamma;
ddC=AccCon(1:3);

qd= [xd yd zd wx wy wz];
Ener = 0.5 * (qd * M* qd') +m0*9.81*z - 5.435696790865547e+03;

%Unit constraint
UC  = 1 - b0^2  - b1^2 - b2^2 - b3^2;
end

% Custom Runge-Kutta 4th order method
function [t, F, second_derivatives, C, dC,ddC,Ener,UC] = custom_rk4(odefun, tspan, y0)
    n = length(tspan);
    F_temp = zeros(length(y0), n);
    F_temp(:,1) = y0;

    second_derivatives = zeros(6, n);
    C = zeros(3,n);
    dC = zeros(3,n);
    ddC = zeros(3,n);
    Ener=zeros(1,n);
    UC = zeros(1,n);

    for i = 1:(n-1)
        h = tspan(i+1) - tspan(i);
        ti = tspan(i);
        yi = F_temp(:,i);
        
        % Split into q (positions/orientations) and v (velocities)
        tri= yi(1:3);  % Translational parameters: [x,y,z]
        qi = yi(4:7);  % Rotational parameters: [b0,b1,b2,b3]
        vi = yi(8:13); % Velocities: [xd,yd,zd,wx,wy,wz]
        
% Stage 1
        [k1, sd1,C1, dC1,ddC1,Ener1, UC1] = odefun(ti, yi);
        k1_tr= h*k1(1:3);
        k1_q = h*yi(11:13);
        k1_v = h*k1(8:13);
% Stage 2

tr2 = tri + (1/2)*k1_tr;
q2 = LGT(qi, k1_q/2);
v2 = vi + (1/2)*k1_v;
y2 = [tr2;q2; v2];

[k2, ~] = odefun(ti + h/2, y2);
        k2_tr = h*k2(1:3);
        k2_q =  h*invdexp(v2(4:6),k1_q/2);
        k2_v = h*k2(8:13);

% Stage 3

tr3 = tri + (1/2)*k2_tr;
q3 = LGT(qi, k2_q/2);
v3 = vi + (1/2)*k2_v;
y3 = [tr3;q3; v3];

        [k3,  ~] = odefun(ti + h/2, y3);
        k3_tr = h*k3(1:3);
        k3_q =  h*invdexp(v3(4:6),k2_q/2);
        k3_v = h*k3(8:13);


% Stage 4            
tr4 = tri + k3_tr;
q4 = LGT(qi, k3_q);
v4 = vi + k3_v;
y4 = [tr4;q4; v4];

        [k4, ~] = odefun(ti + h, y4);
        k4_tr = h*k4(1:3);
        k4_q =  h*invdexp(v4(4:6),k3_q);
        k4_v = h*k4(8:13);

Delta_q = (1/6)*(k1_q +2*k2_q +2*k3_q +k4_q);
tr_new = tri + (1/6)*(k1_tr + 2*k2_tr + 2*k3_tr + k4_tr);
q_new = LGT(qi, Delta_q);
v_new = vi + (1/6)*(k1_v + 2*k2_v + 2*k3_v + k4_v);
        
        % Store full state (13 elements)
        F_temp(:,i+1) = [tr_new;q_new; v_new];

        second_derivatives(:,i) = sd1;  
        C(:,i) = C1;
        dC(:,i) = dC1;
        ddC(:,i) = ddC1;
        Ener(:,i) = Ener1;
        UC(:,i) = UC1;
    end
    [~, sd_final,c_final, dc_final,ddc_final,Ener_final,UC_final] = odefun(tspan(end), F_temp(:,end));
    second_derivatives(:,end) = sd_final;
    C(:,end) = c_final;
    dC(:,end) = dc_final;
    ddC(:,end) = ddc_final;
    Ener(:,end) = Ener_final;
    UC(:,end) = UC_final;

    t = tspan;
    F = F_temp';
    second_derivatives = second_derivatives';
    C=C';
    dC=dC';
    ddC=ddC';
    Ener=Ener';
    UC=UC';
end

% Solve the ODE using custom RK4
[t, F, second_derivatives, C, dC, ddC,Ener, UC] = custom_rk4(@odesystem, tspan, initial_conditions);

% Extract results from Y matrix
% Extract results from Y matrix
x = F(:,1);
y = F(:,2);
z = F(:,3);
b0 = F(:,4);
b1 = F(:,5);
b2 = F(:,6);
b3 = F(:,7);
xd = F(:,8);
yd = F(:,9);
zd = F(:,10);
wx= F(:,11);
wy = F(:,12);
wz = F(:,13);

% Accelerations
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
    save_filename = sprintf('HT_EPLGIM_dt_%.2fms.mat', dt*1000);
end
save(save_filename, '-v7.3');

% Plotting results for x y z (commented out)
% figure;
% subplot(3,1,1);
% % plot(t,x,'b--','LineWidth',1);
% % hold on;
% % plot(t,y,'r-','LineWidth',1); 
% % hold on;
% % plot(t,z,'g-','LineWidth',1); 
% % ylim([-1,1]);
% % title('Rx, Ry, Rz');
% % ylabel('Displacements, m');
% % grid on;

% subplot(3,1,2);
% % plot(t,xd,'b--','LineWidth',1); 
% % hold on;
% % plot(t,yd,'r-','LineWidth',1); 
% % hold on;
% % plot(t,zd,'g-','LineWidth',1); 
% % ylim([-10,10]);
% % title('Rx'', Ry'', Rz'' (t)');
% % ylabel({'Velocities'; '(m/s)'});
% % grid on;
%
% subplot(3,1,3);
% % plot(t,x_double_prime,'b--','LineWidth',1); 
% % hold on;
% % plot(t,y_double_prime,'r-','LineWidth',1); 
% % hold on;
% % plot(t,z_double_prime,'g-','LineWidth',1); 
% % ylim([-50,45]);
% % title('Rx'''', Ry'''', Rz''''(t)');
% % xlabel('Time (s)');
% % ylabel({'Accelerations'; '(m/s^2)'});
% % legend ('Rx','Ry','Rz');
% % grid on;
%
% % Plotting results for x
% figure;
% subplot(3,1,1);
% plot(t,x,'b-','LineWidth',2);
% ylim([-1,1]);
% title('Rx(t)');
% ylabel('Displacement X, m');
% grid on;

% subplot(3,1,2);
% plot(t,xd,'r-','LineWidth',2); 
% ylim([-10,10]);
% title('Rx''(t)');
% ylabel({'Velocity X'; '(m/s)'});
% grid on;

% subplot(3,1,3);
% plot(t,x_double_prime,'g-','LineWidth',2); 
% ylim([-50,45]);
% title('Rx''''(t)');
% xlabel('Time (s)');
% ylabel({'Acceleration X'; '(m/s^2)'});
% grid on;

% Plotting results for y
% figure;
% subplot(3,1,1);
% plot(t,y,'b-','LineWidth',2); 
% ylim([-1,1]);
% title('Ry(t)');
% ylabel('Displacement Y, m');
% grid on;

% subplot(3,1,2);
% plot(t,yd,'r-','LineWidth',2); 
% ylim([-10,10]);
% title('Ry''(t)');
% ylabel({'Velocity Y'; '(m/s)'});
% grid on;

% subplot(3,1,3);
% plot(t,y_double_prime,'g-','LineWidth',2); 
% ylim([-50,50]);
% title('Ry''''(t)');
% xlabel('Time (s)');
% ylabel({'Acceleration Y'; '(m/s^2)'});
% grid on;

% Plotting results for z
% figure;
% subplot(3,1,1);
% plot(t,z,'b-','LineWidth',2); 
% ylim([-1.5,0.5]);
% title('Rz(t)');
% ylabel('Displacement Z, m');
% grid on;

% subplot(3,1,2);
% plot(t,zd,'r-','LineWidth',2); 
% ylim([-10,10]);
% title('Rz''(t)');
% ylabel({'Velocity Z'; '(m/s)'});
% grid on;

% subplot(3,1,3);
% plot(t,z_double_prime,'g-','LineWidth',2); 
% ylim([-40,50]);
% title('Rz''''(t)');
% xlabel('Time (s)');
% ylabel({'Acceleration Z'; '(m/s^2)'});
% grid on;

% Plotting results for EP
% figure;
% subplot(4,1,1);
% plot(t,b0,'b-','LineWidth',2); 
% ylim([-1,1]);
% title('b0');
%% ylabel({'b0'});
% grid on;

% subplot(4,1,2);
% plot(t,b1,'r-','LineWidth',2); 
% ylim([-1,1]);
% title('b1');
%% ylabel({'b1'});
% grid on;

% subplot(4,1,3);
% plot(t,b2,'g-','LineWidth',2); 
% ylim([-1,1]);
% title('b2');
%% ylabel({'b2'});
% grid on;

% subplot(4,1,4);
% plot(t,b3,'c-','LineWidth',2); 
% ylim([-1,1]);
% title('b3');
%% ylabel({'b3'});
% grid on;

% Plotting results for Angular velocities
% figure;
% subplot(3,1,1);
% plot(t,wx,'b-','LineWidth',2); 
% ylim([-7,7]);
% title('Wx');
% ylabel({'Ang velocity Wx'; 's-1'});
% grid on;

% subplot(3,1,2);
% plot(t,wy,'r-','LineWidth',2); 
% ylim([0,200]);
% title('Wy');
% ylabel({'Ang velocity Wy'; 's-1'});
% grid on;

% subplot(3,1,3);
% plot(t,wz,'g-','LineWidth',2); 
% ylim([-7,7]);
% title('Wz');
% ylabel({'Ang velocity Wz';'s-1'});
% grid on;

% Plotting results for Angular accelerations
% figure;
% subplot(3,1,1);
% plot(t,AngAccX,'b-','LineWidth',2); 
% ylim([-1000,1000]);
% title('Wx''');
% ylabel({'Ang acceleration Wx';' s-2'});
% grid on;

% subplot(3,1,2);
% plot(t,AngAccY,'r-','LineWidth',2); 
% ylim([-3,3]);
% title('Wy''');
% ylabel({'Ang acceleration Wy';' s-2'});
% grid on;

% subplot(3,1,3);
% plot(t,AngAccZ,'g-','LineWidth',2); 
% ylim([-1000,1000]);
% title('Wz''');
% ylabel({'Ang acceleration Wz';' s-2'});
% grid on;

% Plotting results for Constraint violation and Energy balance
% figure;
% subplot(4,1,1);
% plot(t,grC,'b-','LineWidth',0.5); 
% ylim([-1*10^-3,1*10^-3]);
% title('∣∣C∣∣');
% ylabel('Violation');
% grid on;

% subplot(4,1,2);
% plot(t,grdC,'g-','LineWidth',1); 
% ylim([-1*10^-4,1*10^-4]);
% title('∣∣C''∣∣');
% ylabel('Violation');
% grid on;

% subplot(4,1,3);
% plot(t,grddC,'r-','LineWidth',1); 
% ylim([-1*10^-12,1*10^-12]);
% title('∣∣C''''∣∣');
% ylabel('Violation');
% grid on;

% subplot(4,1,4);
% plot(t,Ener,'b-','LineWidth',1); 
% ylim([-1*10^-0,1*10^-0]);
% title('Total Energy');
% xlabel('Time (s)');
% ylabel('Energy balance');
% grid on;

% Plotting results for Unit constraints
% figure;
% plot(t,UC,'b-','LineWidth',0.5); 
% ylim([-1*10^-13,1*10^-13]);
% title('Unit costraint - EP Lie');
% ylabel('Violation');
% grid on;

% Plotting results for Constraint violation and Energy balance
% figure;
% subplot(4,1,1);
% % plot(t,abs(grC),'b-','LineWidth',0.5);
% set(gca, 'YScale', 'log');
	% set(gca, 'YMinorTick', 'on');
% title('||C|| (Log scale)');
% ylabel('Violation');
% grid on;

% subplot(4,1,2);
% % plot(t,abs(grdC),'g-','LineWidth',1);
% set(gca, 'YScale', 'log');
	% set(gca, 'YMinorTick', 'on');
% title('||C''|| (Log scale)');
% ylabel('Violation');
% grid on;

% subplot(4,1,3);
% % plot(t,abs(grddC),'r-','LineWidth',1);
% set(gca, 'YScale', 'log');
	% set(gca, 'YMinorTick', 'on');
% title('||C''''|| (Log scale)');
% ylabel('Violation');
% grid on;

% subplot(4,1,4);
% % plot(t,abs(Ener),'b-','LineWidth',1);
% set(gca, 'YScale', 'log');
	% set(gca, 'YMinorTick', 'on');
% title('Total Energy (Log scale)');
% xlabel('Time (s)');
% ylabel('Energy balance');
% grid on;