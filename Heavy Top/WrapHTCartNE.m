tic;
% Initial values for psi - rotational part
psi10= 0;
psi20= 0;
psi30= 0;
t_end = 20;  % End value for time

%dt=0.0012;
% Time step (seconds)
if ~exist('dt','var') || isempty(dt)
    dt = 0.1;
end
% Time span
tspan = 0:dt:t_end;

x0=0;
y0=1;
z0=0;

wx0=0;
wy0=150;
wz0= -4.61538;
% wx0=0;
% wy0=0;
% wz0= 0;

xd0=4.615380000000000;
%xd0=0;
yd0=0;
zd0=0;

m0=15;
l0=1;

% Inertia tensor
Iqq = diag([0.234375; 0.46875; 0.234375]);

% Initial conditions
initial_conditions =[x0; xd0; y0; yd0;z0; zd0; wx0; wy0; wz0;psi10;psi20;psi30];

function S = skew(v)
    S = [0 -v(3) v(2); v(3) 0 -v(1); -v(2) v(1) 0];
end

sj=[0;-l0;0];
Sj_hat = skew(sj);

Ener0 =5.435696790865547e+03;

% Define the ODE system as a function
function [dydt, second_derivatives, C, dC,ddC, Ener] = odesystem(t, F)
    x = F(1);
    xd = F(2);
    y = F(3);
    yd = F(4);
    z = F(5);
    zd = F(6);
    wx = F(7);
    wy = F(8);
    wz = F(9);
    psi1 = F(10);
    psi2 = F(11);
    psi3 = F(12);
l0=1;
m0=15;

Iqq = diag([0.234375; 0.46875; 0.234375]);

%  Determining angle phi (theta in report) from psi
Phi = sqrt(psi1^2 + psi2^2 + psi3^2);
p = [psi1, psi2, psi3];
Ps = skew (p);

%wrapping 
%    if Phi > pi
 %           fprintf('Wrapping: Phi = %.3f -> ', Phi)
 %       p = -(2*pi - Phi) * (p / Phi);
%Phi=norm(p);
%psi1=p(1);
%psi2=p(2);
%psi3=p(3);
%Ps = skew (p);
%    end

%  Transform matrix A and A2 (for small phi) from Rodriguez formula
if Phi > 10^-4 
A = eye(3) + Ps*(sin(Phi)/Phi) + 2*Ps*Ps*sin(Phi/2)*sin(Phi/2)/(Phi*Phi);
    else
A = eye(3) + Ps*(1 - Phi^2/6 + Phi^4/120) + 2*Ps*Ps*(1/2 - Phi^2/48 + Phi^4/3840);
end

%Gint=GCartMarch(x, y, z, psi1, psi2, psi3);%Function for matrix G
%  Matrix G from Rodriguez formula ???ETH ZURICH (with Taylor series for small phi)
if abs(Phi) < 10^-4
Emap = eye(3) + Ps*(-0.5+Phi^2/24-Phi^4/720+Phi^6/40320) + Ps*Ps*(1/6-Phi^2/120+Phi^4/5040-Phi^6/362880);
Emap2 = eye(3) + 0.5*Ps + Ps*Ps*(1/12 + Phi^2/720);
else
    Emap = eye(3) - Ps*(1-cos(Phi))/(Phi^2) + Ps*Ps*(Phi-sin(Phi))/(Phi^3);
Emap2 = eye(3) + 0.5*Ps + Ps*Ps*(1-0.5*Phi*sin(Phi)/(1-cos(Phi)))/(Phi^2);
end


w_vec=[wx;
wy;
wz];

%Gint= inv (Emap);
%psi_ans=Gint*w_vec;
psi_ans=Emap2*w_vec;

    Phi_rpi = zeros (3,6);
Phi_r= eye(3);
sj=[0;-l0;0];
Sj_hat = skew(sj);

Phi_pi= -A*Sj_hat;
Phi_rpi = [Phi_r Phi_pi];

    M = [m0*eye(3) zeros(3);
        zeros(3) Iqq];

Mat = [M Phi_rpi';
       Phi_rpi zeros(3)];

w_hat=skew(w_vec);

gamma = -A*w_hat*w_hat*sj;

Qv = -w_hat*Iqq*w_vec;

Forces=[0; 0; -m0 *9.81; Qv; gamma];
NewRes=Mat\Forces;

 dydt = zeros(12,1);
    dydt(1) = xd; % x' = xd
    dydt(2) = NewRes(1); % xd' (x'')
    dydt(3) = yd; % y' = yd
    dydt(4) = NewRes(2);  % yd' (y'')  
    dydt(5) = zd; % z' = zd
    dydt(6) = NewRes(3); % zd' (z'')  
    dydt(7) = NewRes(4); % wx'
    dydt(8) = NewRes(5); % wy' 
    dydt(9) = NewRes(6); % wz' 
    dydt(10) = psi_ans(1) ; % b0
    dydt(11) = psi_ans(2); % b1 
    dydt(12) = psi_ans(3); % b2


second_derivatives = NewRes(1:6);

tempC=[x;y;z]+A*sj;
C=tempC(1:3);
VelCon = Phi_r*[xd;yd;zd] + Phi_pi*w_vec;
dC=VelCon(1:3);
AccCon = Phi_r*NewRes(1:3) + Phi_pi*NewRes(4:6)-gamma;
ddC=AccCon(1:3);

qd= [xd yd zd wx wy wz];
Ener = 0.5 * (qd * M* qd')+m0*9.81*z - 5.435696790865547e+03;
end

% Custom Runge-Kutta 4th order method
function [t, F, second_derivatives,C, dC,ddC,Ener] = custom_rk4(odefun, tspan, y0)
    n = length(tspan);
    F_temp = zeros(length(y0), n);
    F_temp(:,1) = y0;

    second_derivatives = zeros(6, n);
    C = zeros(3,n);
    dC = zeros(3,n);
    ddC = zeros(3,n);
    Ener=zeros(1,n);

    for i = 1:(n-1)
        h = tspan(i+1) - tspan(i);
        ti = tspan(i);
        yi = F_temp(:,i);
    
        % --- WRAP at the BEGINNING of each RK4 step ---
        psi = yi(10:12);
        Phi = norm(psi);
        if Phi > pi
fprintf('Wrapping RK4: Phi = %.3f -> ', Phi)
            psi = -(2*pi - Phi) * (psi / Phi);
            yi(10:12) = psi;  % Update yi with wrapped psi
        end

        [k1, sd1,C1, dC1,ddC1,Ener1] = odefun(ti, yi);
        [k2, ~] = odefun(ti + h/2, yi + h*k1/2);
        [k3, ~] = odefun(ti + h/2, yi + h*k2/2);
        [k4, ~] = odefun(ti + h, yi + h*k3);
        
        F_temp(:,i+1) = yi + (h/6)*(k1 + 2*k2 + 2*k3 + k4);
        second_derivatives(:,i) = sd1;  
        C(:,i) = C1;
        dC(:,i) = dC1;
        ddC(:,i) = ddC1;
        Ener(:,i) = Ener1;
    end
    [~, sd_final,c_final, dc_final,ddc_final,Ener_final] = odefun(tspan(end), F_temp(:,end));
    second_derivatives(:,end) = sd_final;
    C(:,end) = c_final;
    dC(:,end) = dc_final;
    ddC(:,end) = ddc_final;
    Ener(:,end) = Ener_final;

    t = tspan;
    F = F_temp';
    second_derivatives = second_derivatives';
    C=C';
    dC=dC';
    ddC=ddC';
    Ener=Ener';
end

% Solve the ODE using custom RK4
[t, F, second_derivatives, C, dC, ddC,Ener] = custom_rk4(@odesystem, tspan, initial_conditions);

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
psi1 = F(:,10);
psi2 = F(:,11);
psi3 = F(:,12);

x_double_prime = second_derivatives(:,1);
y_double_prime = second_derivatives(:,2);
z_double_prime = second_derivatives(:,3);
AngAccX = second_derivatives(:,4);
AngAccY = second_derivatives(:,5);
AngAccZ = second_derivatives(:,6);

%dC=zeros(size(t));
grC=zeros(size(t));
grdC=zeros(size(t));
grddC=zeros(size(t));
n=length(t);
for i=1:n
grC(i)=(C(i,1)+C(i,2)+C(i,3))/3;  
grdC(i)=(dC(i,1)+dC(i,2)+dC(i,3))/3;
grddC(i)=(ddC(i,1)+ddC(i,2)+ddC(i,3))/3;
end

theta = sqrt(psi1.^2 + psi2.^2 + psi3.^2);

executionTime = toc;
if ~exist('save_filename','var') || isempty(save_filename)
    save_filename = sprintf('HT_CartNE_dt_%.2fms.mat', dt*1000);
end
save(save_filename, '-v7.3');

% %Plotting results for x y z (commented out)
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
% %Plotting results for x
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
% %Plotting results for y
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
% %Plotting results for z
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
% %Plotting results for theta
% figure;
% plot(t,theta,'b-','LineWidth',2); 
% xlim([0,0.2]);
% ylim([0,4]);
% title('Angle Theta');
% xlabel('Time (s)');
% ylabel('Norm of rotation vector Theta(t)');
% grid on;
% 
% %Plotting results for psi
% figure;
% subplot(2,1,1);
% plot(t,psi1,'b-','LineWidth',2); 
% hold on;
% plot(t,psi2,'r--','LineWidth',2); 
% xlim([0,0.2]);
% ylim([-4,4]);
% title('psi1(t), psi2(t)');
% ylabel('psi1(t), psi2(t)');
% grid on;
% lgd = legend('psi1', 'psi2', 'Location', 'northeast', 'FontSize', 8, 'TextColor', 'blue');
% 
% subplot(2,1,2);
% plot(t,psi3,'b-','LineWidth',2); 
% ylim([-4,4]);
% title('psi3(t)');
% ylabel('Psi3');
% grid on;
% 
% %Plotting results for Angular velocities
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
% %Plotting results for Angular accelerations
% figure;
% subplot(3,1,1);
% plot(t,AngAccX,'b-','LineWidth',2); 
% ylim([-1000,1000]);
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
% ylim([-1000,1000]);
% title('Wz''');
% ylabel({'Ang acceleration Wz';' s-2'});
% grid on;
% 
% 
% %Plotting results for Constraint violation and Energy balance
% figure;
% subplot(4,1,1);
% plot(t,grC,'b-','LineWidth',0.5); 
% ylim([-1*10^-10,1*10^-10]);
% title('∣∣C∣∣');
% ylabel('Violation');
% grid on;
% 
% subplot(4,1,2);
% plot(t,grdC,'g-','LineWidth',1); 
% ylim([-1*10^-10,1*10^-10]);
% title('∣∣C''∣∣');
% ylabel('Violation');
% grid on;
% 
% subplot(4,1,3);
% plot(t,grddC,'r-','LineWidth',1); 
% ylim([-5*10^-12,5*10^-12]);
% title('∣∣C''''∣∣');
% ylabel('Violation');
% grid on;
% 
% subplot(4,1,4);
% plot(t,Ener,'b-','LineWidth',1); 
% ylim([-1*10^-7,1*10^-7]);
% title('Total Energy');
% xlabel('Time (s)');
% ylabel('Energy balance');
% grid on;
% 
% %Plotting results for Constraint violation and Energy balance
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