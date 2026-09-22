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

% Inertia tensor (COM)
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

t_end = 20;  % End value for time
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
% unconstrained: Euler about the fixed point, no spherical-joint lambda
function [dydt, second_derivatives, C, dC,ddC, Ener, UC] = odesystem(t, F)
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

sj=[0;-l0;0];
rb=[0;l0;0];

Iqq = diag([0.234375; 0.46875; 0.234375]);
    M = [m0*eye(3) zeros(3);
        zeros(3) Iqq];

    w_hat=skew(w_vec);

% inertia at the fixed point (parallel axis)
IFP = Iqq + m0*((rb'*rb)*eye(3) - rb*rb');
% gravity torque about the fixed point, body frame
TFP = cross(rb, A'*[0; 0; -m0*9.81]);
omegadot = IFP \ (TFP - w_hat*IFP*w_vec);

% COM kinematics of rotation about the origin
v_kin = A*cross(w_vec, rb);
a_kin = A*(cross(omegadot, rb) + w_hat*w_hat*rb);

 dydt = zeros(13,1);
    dydt(1) = xd; % x' = xd
    dydt(2) = a_kin(1); % xd' (x'')
    dydt(3) = yd; % y' = yd
    dydt(4) = a_kin(2);  % yd' (y'')
    dydt(5) = zd; % z' = zd
    dydt(6) = a_kin(3); % zd' (z'')
    dydt(7) = omegadot(1); % wx'
    dydt(8) = omegadot(2); % wy'
    dydt(9) = omegadot(3); % wz'
    dydt(10) = b_ans(1) ; % b0
    dydt(11) = b_ans(2); % b1
    dydt(12) = b_ans(3); % b2
    dydt(13) = b_ans(4); % b3

second_derivatives = [a_kin; omegadot];
tempC=[x;y;z]+A*sj;
C=tempC(1:3); 
dC=[xd;yd;zd] - v_kin;
ddC=[0;0;0];

qd= [xd yd zd wx wy wz];
Ener = 0.5 * (qd * M* qd') +m0*9.81*z - 5.435696790865547e+03;

%Unit constraint
UC  = 1 - b0^2  - b1^2 - b2^2 - b3^2;
end

%options = odeset('RelTol', 1e-14, 'AbsTol', 1e-14);
options = odeset('RelTol', 100*eps, 'AbsTol', 1e-16);
% Solve the ODE using ode45 with the anonymous function
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
grC = zeros(size(t));
grdC = zeros(size(t));
grddC = zeros(size(t));
nt = length(t);
for i = 1:nt
    grC(i) = (C(i, 1) + C(i, 2) + C(i, 3)) / 3;
    grdC(i) = (dC(i, 1) + dC(i, 2) + dC(i, 3)) / 3;
    grddC(i) = (ddC(i, 1) + ddC(i, 2) + ddC(i, 3)) / 3;
end

executionTime = toc;
save_filename = 'UncHTRef.mat';
save(save_filename, '-v7.3');

% --- Paper-style figures (Times New Roman, units in [], no titles, no minor grids)
paperFont = 'Times New Roman';

% Cartesian position, velocity, acceleration ($R_X,R_Y,R_Z$)
figure;
subplot(3, 1, 1);
hold on;
plot(t, x, 'b--', 'LineWidth', 1);
plot(t, y, 'r-', 'LineWidth', 1);
plot(t, z, 'g-', 'LineWidth', 1);
ylim([-1, 1]);
grid on;
set(gca, 'FontName', paperFont, 'XMinorGrid', 'off', 'YMinorGrid', 'off');
ylabel('Displacement [m]', 'FontName', paperFont);
legend({'$R_X$', '$R_Y$', '$R_Z$'}, 'Interpreter', 'latex', 'Location', 'best', 'FontName', paperFont);

subplot(3, 1, 2);
hold on;
plot(t, xd, 'b--', 'LineWidth', 1);
plot(t, yd, 'r-', 'LineWidth', 1);
plot(t, zd, 'g-', 'LineWidth', 1);
ylim([-10, 10]);
grid on;
set(gca, 'FontName', paperFont, 'XMinorGrid', 'off', 'YMinorGrid', 'off');
ylabel('Velocity [m/s]', 'FontName', paperFont);

subplot(3, 1, 3);
hold on;
plot(t, x_double_prime, 'b--', 'LineWidth', 1);
plot(t, y_double_prime, 'r-', 'LineWidth', 1);
plot(t, z_double_prime, 'g-', 'LineWidth', 1);
ylim([-50, 45]);
grid on;
set(gca, 'FontName', paperFont, 'XMinorGrid', 'off', 'YMinorGrid', 'off');
xlabel('Time [s]', 'FontName', paperFont);
ylabel('Acceleration [m/s$^2$]', 'Interpreter', 'latex', 'FontName', paperFont);

% Angular velocity
figure;
subplot(3, 1, 1);
plot(t, wx, 'b-', 'LineWidth', 1.5);
ylim([-7, 7]);
grid on;
set(gca, 'FontName', paperFont, 'XMinorGrid', 'off', 'YMinorGrid', 'off');
ylabel('$\omega_x$ [rad/s]', 'Interpreter', 'latex', 'FontName', paperFont);

subplot(3, 1, 2);
plot(t, wy, 'r-', 'LineWidth', 1.5);
ylim([0, 200]);
grid on;
set(gca, 'FontName', paperFont, 'XMinorGrid', 'off', 'YMinorGrid', 'off');
ylabel('$\omega_y$ [rad/s]', 'Interpreter', 'latex', 'FontName', paperFont);

subplot(3, 1, 3);
plot(t, wz, 'g-', 'LineWidth', 1.5);
ylim([-7, 7]);
grid on;
set(gca, 'FontName', paperFont, 'XMinorGrid', 'off', 'YMinorGrid', 'off');
xlabel('Time [s]', 'FontName', paperFont);
ylabel('$\omega_z$ [rad/s]', 'Interpreter', 'latex', 'FontName', paperFont);

% Angular acceleration
figure;
subplot(3, 1, 1);
plot(t, AngAccX, 'b-', 'LineWidth', 1.5);
ylim([-1000, 1000]);
grid on;
set(gca, 'FontName', paperFont, 'XMinorGrid', 'off', 'YMinorGrid', 'off');
ylabel('$\dot{\omega}_x$ [rad/s$^2$]', 'Interpreter', 'latex', 'FontName', paperFont);

subplot(3, 1, 2);
plot(t, AngAccY, 'r-', 'LineWidth', 1.5);
ylim([-3, 3]);
grid on;
set(gca, 'FontName', paperFont, 'XMinorGrid', 'off', 'YMinorGrid', 'off');
ylabel('$\dot{\omega}_y$ [rad/s$^2$]', 'Interpreter', 'latex', 'FontName', paperFont);

subplot(3, 1, 3);
plot(t, AngAccZ, 'g-', 'LineWidth', 1.5);
ylim([-1000, 1000]);
grid on;
set(gca, 'FontName', paperFont, 'XMinorGrid', 'off', 'YMinorGrid', 'off');
xlabel('Time [s]', 'FontName', paperFont);
ylabel('$\dot{\omega}_z$ [rad/s$^2$]', 'Interpreter', 'latex', 'FontName', paperFont);

% Constraint residuals and energy (log magnitude)
figure;
subplot(4, 1, 1);
plot(t, abs(grC), 'b-', 'LineWidth', 0.75);
set(gca, 'YScale', 'log', 'FontName', paperFont, 'XMinorGrid', 'off', 'YMinorGrid', 'off');
grid on;
ylabel('$||\mathbf{C}||$ [-]', 'Interpreter', 'latex', 'FontName', paperFont);

subplot(4, 1, 2);
plot(t, abs(grdC), 'g-', 'LineWidth', 1);
set(gca, 'YScale', 'log', 'FontName', paperFont, 'XMinorGrid', 'off', 'YMinorGrid', 'off');
grid on;
ylabel('$||\dot{\mathbf{C}}||$ [-]', 'Interpreter', 'latex', 'FontName', paperFont);

subplot(4, 1, 3);
plot(t, abs(grddC), 'r-', 'LineWidth', 1);
set(gca, 'YScale', 'log', 'FontName', paperFont, 'XMinorGrid', 'off', 'YMinorGrid', 'off');
grid on;
ylabel('$||\ddot{\mathbf{C}}||$ [-]', 'Interpreter', 'latex', 'FontName', paperFont);

subplot(4, 1, 4);
plot(t, abs(Ener), 'b-', 'LineWidth', 1);
set(gca, 'YScale', 'log', 'FontName', paperFont, 'XMinorGrid', 'off', 'YMinorGrid', 'off');
grid on;
xlabel('Time [s]', 'FontName', paperFont);
ylabel('Energy deviation magnitude [J]', 'FontName', paperFont);

% EP unit constraint
figure;
plot(t, UC, 'b-', 'LineWidth', 0.75);
ylim([-1e-13, 1e-13]);
grid on;
set(gca, 'FontName', paperFont, 'XMinorGrid', 'off', 'YMinorGrid', 'off');
xlabel('Time [s]', 'FontName', paperFont);
ylabel('EP unit constraint deviation [-]', 'FontName', paperFont);
