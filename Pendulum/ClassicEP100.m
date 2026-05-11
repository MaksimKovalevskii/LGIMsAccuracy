tic;
b10=0.27059805007;
b20=0.27059805007;
b30=0;
b00=sqrt(1-b10^2-b20^2);
t_end = 100;  % End value for time
if ~exist('dt', 'var') || isempty(dt)
    dt = 0.01;   % Time step (default 10 ms)
end
I = [3.0025 0 0; 0 3.0025 0; 0 0 0.005];

% Initial conditions
initial_conditions = [6.*(b10.*b30+b00.*b20); 6.*(b20.*b30-b00.*b10); 3.*(1-2.*b10.^2-2.*b20.^2);b00;b10;b20;b30;0; 0;0; 0;0;0;0];

% Time span
tspan = 0:dt:t_end;

function [dydt, second_derivatives]= odesystem(t, F)
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
    b0d = F(11);
    b1d = F(12);
    b2d = F(13);
    b3d = F(14);

betasd=[b0d b1d b2d b3d];
I=[3.0025 0 0; 0 3.0025 0;0 0 0.005];

Gint=2.*[-b1 b0  b3 -b2;
-b2 -b3  b0 b1;
-b3 b2  -b1 b0];
G2=2.*[-b1d b0d  b3d -b2d;
-b2d -b3d  b0d b1d;
-b3d b2d  -b1d b0d];

omegaint=Gint*betasd';

Qv1=cross(omegaint,I*omegaint)+I*G2*betasd';
Qv=-Gint'*Qv1; 

Mqq=Gint'*I*Gint;

M = [eye(3) zeros(3,4);
    zeros(4,3) Mqq];
Cq = [1 0 0 -6.*b2 -6.*b3 -6.*b0 -6.*b1 ;
    0 1 0 6.*b1 6.*b0 -6.*b3 -6.*b2 ;
    0 0 1 0 12.*b1 12.*b2 0 ;
    0 0 0 2.*b0 2.*b1 2.*b2 2.*b3 ];

Mat = [M Cq';
       Cq zeros(4)];

Forces=[0; 0; -9.8;Qv;12.*b2d.*b0d+ 12.*b3d.*b1d;-12.*b1d.*b0d + 12.*b2d.*b3d;-12.*b1d.*b1d-12.*b2d.*b2d;-2.*b1d.*b1d-2.*b2d.*b2d-2.*b0d.*b0d-2.*b3d.*b3d ];
NewRes=Mat\Forces;

    dydt = zeros(14,1);
    dydt(1) = xd; % x' = xd
    dydt(2) = yd; % y' = yd
    dydt(3) = zd; % z' = zd
    dydt(4) = b0d; 
    dydt(5) = b1d; 
    dydt(6) = b2d; 
    dydt(7) = b3d; 

    dydt(8) = NewRes(1); % xd' (x'')
    dydt(9) = NewRes(2);  % yd' (y'')  
    dydt(10) = NewRes(3); % zd' (z'') 
    dydt(11) = NewRes(4); 
    dydt(12) = NewRes(5);  
    dydt(13) = NewRes(6);  
    dydt(14) = NewRes(7); 

 second_derivatives = NewRes(1:7);

end

% Custom Runge-Kutta 4th order method
function [t, F, second_derivatives] = custom_rk4(odefun, tspan, y0)
    n = length(tspan);
    F_temp = zeros(length(y0), n);
    F_temp(:,1) = y0;

    second_derivatives = zeros(7, n);
  
    for i = 1:(n-1)
        h = tspan(i+1) - tspan(i);
        ti = tspan(i);
        yi = F_temp(:,i);
        
        [k1, sd1] = odefun(ti, yi);
        [k2, ~] = odefun(ti + h/2, yi + h*k1/2);
        [k3, ~] = odefun(ti + h/2, yi + h*k2/2);
        [k4, ~] = odefun(ti + h, yi + h*k3);
        
        F_temp(:,i+1) = yi + (h/6)*(k1 + 2*k2 + 2*k3 + k4);
        second_derivatives(:,i) = sd1;  
    
    end
    [~, sd_final] = odefun(tspan(end), F_temp(:,end));
    second_derivatives(:,end) = sd_final;

    t = tspan;
    F = F_temp';
    second_derivatives = second_derivatives';

end

% Solve the ODE using custom RK4
[t, F, second_derivatives] = custom_rk4(@odesystem, tspan, initial_conditions);

% Extract results from F matrix
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


% Calculate norms of Constraint violations' matrices for plotting
C1=x-6.*(b1.*b3+b0.*b2);
C2=y-6.*(b2.*b3-b0.*b1);
C3=z-3.*(1-2.*b1.^2-2.*b2.^2);
C4=b0.^2+b1.^2+b2.^2+b3.^2-1;
C=(C1+C2+C3+C4)/4; %C

dC1=xd-6.*b2.*b0d -6.*b3.*b1d -6.*b0.*b2d -6.*b1.*b3d ;
dC2=yd+6.*b1.*b0d +6.*b0.*b1d -6.*b3.*b2d -6.*b2.*b3d;
dC3=zd+12.*b1.*b1d+12.*b2.*b2d ;
dC4=2.*b0.*b0d+2.*b3.*b3d+2.*b1.*b1d+2.*b2.*b2d;
dC=(dC1+dC2+dC3+dC4)/4; %C'

ddC1=x_double_prime-6.*b2.*b0_double_prime -6.*b3.*b1_double_prime -6.*b0.*b2_double_prime -6.*b1.*b3_double_prime-6.*b2d.*b0d -6.*b3d.*b1d -6.*b0d.*b2d -6.*b1d.*b3d;
ddC2=y_double_prime+6.*b1.*b0_double_prime +6.*b0.*b1_double_prime -6.*b3.*b2_double_prime -6.*b2.*b3_double_prime+6.*b1d.*b0d +6.*b0d.*b1d -6.*b3d.*b2d -6.*b2d.*b3d;
ddC3=z_double_prime + 12.*b1.*b1_double_prime+12.*b2.*b2_double_prime +12.*b1d.*b1d+12.*b2d.*b2d;
ddC4=2.*b0.*b0_double_prime+2.*b3.*b3_double_prime+2.*b1.*b1_double_prime+2.*b2.*b2_double_prime +2.*b1d.*b1d+2.*b2d.*b2d+2.*b0d.*b0d+2.*b3d.*b3d;
ddC=(ddC1+ddC2+ddC3+ddC4)/4; %C''

% Calculate Energy balance for plotting
for i = 1:length(t)
qd= [xd(i) yd(i) zd(i) b0d(i) b1d(i) b2d(i) b3d(i)];

Gint=2.*[-b1(i) b0(i)  b3(i) -b2(i);
-b2(i) -b3(i)  b0(i) b1(i);
-b3(i) b2(i)  -b1(i) b0(i)];

Mqq=Gint'*I*Gint;
MEner=[1 0 0 0 0 0 0;
0 1 0 0 0 0 0;
0 0 1 0 0 0 0;
 0 0 0 Mqq(1,1) Mqq(1,2) Mqq(1,3) Mqq(1,4);
    0 0 0 Mqq(2,1) Mqq(2,2) Mqq(2,3) Mqq(2,4);
    0 0 0 Mqq(3,1) Mqq(3,2) Mqq(3,3) Mqq(3,4);
    0 0 0 Mqq(4,1) Mqq(4,2) Mqq(4,3) Mqq(4,4)];

Ener232(i)=0.5 * (qd * MEner * qd')+9.8.*(z(i)-z(1));
end
Enernew2=Ener232';

UC2 = 1 - b0.^2  - b1.^2 - b2.^2 - b3.^2;

executionTime = toc
if ~exist('save_filename', 'var') || isempty(save_filename)
    save_filename = sprintf('EP_dt_%.1fms.mat', dt * 1000);
end
save(save_filename, '-v7.3');

% Plotting results for x
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
% Plotting results for y
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
% Plotting results for z
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
% Plotting results for b0
% figure;
% subplot(3,1,1);
% plot(t,b0,'b-','LineWidth',2); 
% ylim([-1,1]);
% title('b0(t)');
% ylabel('Parameter b0');
% grid on;
% 
% subplot(3,1,2);
% plot(t,b0d,'r-','LineWidth',2); 
% ylim([-3,3]);
% title('b0''(t)');
% ylabel('1st Derivative b0');
% grid on;
% 
% subplot(3,1,3);
% plot(t,b0_double_prime,'g-','LineWidth',2); 
% ylim([-3,3]);
% title('b0''''(t)');
% xlabel('Time (s)');
% ylabel('2nd Derivative b0');
% grid on;
% 
% Plotting results for b1 and b2
% figure;
% subplot(3,1,1);
% plot(t,b1,'b-','LineWidth',2); 
% hold on;
% plot(t,b2,'r--','LineWidth',2); 
% ylim([0,1]);
% title('b1(t), b2(t)');
% ylabel('Parameters b1,b2');
% grid on;
% lgd = legend('b1', 'b2', 'Location', 'northeast', 'FontSize', 8, 'TextColor', 'blue');
% 
% subplot(3,1,2);
% plot(t,b1d,'b-','LineWidth',2); 
% hold on;
% plot(t,b2d,'r--','LineWidth',2); 
% ylim([-1,1]);
% title('b1''(t), b2''(t)');
% ylabel('1st Derivatives');
% grid on;
% 
% subplot(3,1,3);
% plot(t,b1_double_prime,'b-','LineWidth',2); 
% hold on;
% plot(t,b2_double_prime,'r--','LineWidth',2); 
% ylim([-2,1]);
% title('b1''''(t), b2''''(t)');
% xlabel('Time (s)');
% ylabel('2nd Derivatives');
% grid on;
% 
% Plotting results for b3
% figure;
% subplot(3,1,1);
% plot(t,b3,'b-','LineWidth',2); 
% ylim([-1,1]);
% title('b3(t)');
% ylabel('Parameter b3');
% grid on;
% 
% subplot(3,1,2);
% plot(t,b3d,'r-','LineWidth',2); 
% ylim([-1,1]);
% title('b3''(t)');
% ylabel('1st Derivative b3');
% grid on;
% 
% subplot(3,1,3);
% plot(t,b3_double_prime,'g-','LineWidth',2); 
% ylim([-1,1]);
% title('b3''''(t)');
% xlabel('Time (s)');
% ylabel('2nd Derivative b3');
% grid on;
% 
% Plotting results for Constraint violation and Energy balance
% figure;
% subplot(4,1,1);
% plot(t,C,'b-','LineWidth',0.5); 
% ylim([-1*10^-3,1*10^-3]);
% title('∣∣C∣∣');
% ylabel('Violation');
% grid on;
% 
% subplot(4,1,2);
% plot(t,dC,'g-','LineWidth',1); 
% ylim([-5*10^-4,5*10^-4]);
% title('∣∣C''∣∣');
% ylabel('Violation');
% grid on;
% 
% subplot(4,1,3);
% plot(t,ddC,'r-','LineWidth',1); 
% ylim([-5*10^-13,5*10^-13]);
% title('∣∣C''''∣∣');
% ylabel('Violation');
% grid on;
% 
% subplot(4,1,4);
% plot(t,Enernew2,'b-','LineWidth',1); 
% ylim([-1*10^-1,1*10^-1]);
% title('Total Energy');
% xlabel('Time (s)');
% ylabel('Energy balance');
% grid on;
% 
% Plotting results for Unit constraints
% figure;
% plot(t,UC2,'b-','LineWidth',0.5); 
% ylim([-5*10^-1,5*10^-1]);
% title('Unit costraint - EP Classic');
% ylabel('Violation');
% grid on;
% 
% Plotting results for Constraint violation and Energy balance
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
% plot(t,abs(dC),'g-','LineWidth',1);
% set(gca, 'YScale', 'log');
% 	set(gca, 'YMinorTick', 'on');
% title('||C''|| (Log scale)');
% ylabel('Violation');
% grid on;
% 
% subplot(4,1,3);
% plot(t,abs(ddC),'r-','LineWidth',1);
% set(gca, 'YScale', 'log');
% 	set(gca, 'YMinorTick', 'on');
% title('||C''''|| (Log scale)');
% ylabel('Violation');
% grid on;
% 
% subplot(4,1,4);
% plot(t,abs(Enernew2),'b-','LineWidth',1);
% set(gca, 'YScale', 'log');
% 	set(gca, 'YMinorTick', 'on');
% title('Total Energy (Log scale)');
% xlabel('Time (s)');
% ylabel('Energy balance');
% grid on;
% 
% Plotting results for position, velocity, and acceleration with all x,y,z components together
% figure;
% 
% Position
% subplot(3,1,1);
% plot(t,x,'b-','LineWidth',2); hold on;
% plot(t,y,'r-','LineWidth',2);
% plot(t,z,'g-','LineWidth',2);
% legend('R_X(t)', 'R_Y(t)', 'R_Z(t)', 'Location', 'best');
% title('Position Components');
% ylabel('Displacement (m)');
% grid on;
% 
% Velocity
% subplot(3,1,2);
% plot(t,xd,'b-','LineWidth',2); hold on;
% plot(t,yd,'r-','LineWidth',2);
% plot(t,zd,'g-','LineWidth',2);
% legend('$\dot{R}_X(t)$', '$\dot{R}_Y(t)$', '$\dot{R}_Z(t)$', 'Location', 'best', 'Interpreter', 'latex');
% title('Velocity Components');
% ylabel('Velocity (m/s)');
% grid on;
% 
% Acceleration
% subplot(3,1,3);
% plot(t,x_double_prime,'b-','LineWidth',2); hold on;
% plot(t,y_double_prime,'r-','LineWidth',2);
% plot(t,z_double_prime,'g-','LineWidth',2);
% legend('$\ddot{R}_X(t)$', '$\ddot{R}_Y(t)$', '$\ddot{R}_Z(t)$', 'Location', 'best', 'Interpreter', 'latex');
% title('Acceleration Components');
% xlabel('Time (s)');
% ylabel('Acceleration (m/s²)');
% grid on;
% 
% Plotting results for x (paper style)
% paperFont = 'Times New Roman';
% figure;
% subplot(3, 1, 1);
% plot(t, x, 'b-', 'LineWidth', 2);
% ylim([-3, 3]);
% grid on;
% set(gca, 'FontName', paperFont, 'XMinorGrid', 'off', 'YMinorGrid', 'off');
% ylabel('$R_X$ [m]', 'Interpreter', 'latex', 'FontName', paperFont);
% subplot(3, 1, 2);
% plot(t, xd, 'r-', 'LineWidth', 2);
% ylim([-10, 10]);
% grid on;
% set(gca, 'FontName', paperFont, 'XMinorGrid', 'off', 'YMinorGrid', 'off');
% ylabel('$\dot{R}_X$ [m/s]', 'Interpreter', 'latex', 'FontName', paperFont);
% subplot(3, 1, 3);
% plot(t, x_double_prime, 'g-', 'LineWidth', 2);
% ylim([-20, 20]);
% grid on;
% set(gca, 'FontName', paperFont, 'XMinorGrid', 'off', 'YMinorGrid', 'off');
% xlabel('Time [s]', 'FontName', paperFont);
% ylabel('$\ddot{R}_X$ [m/s$^2$]', 'Interpreter', 'latex', 'FontName', paperFont);
