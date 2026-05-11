% Parameters
b00=1;
b10=0;
b20=0;
b30=0;

wx0=1e-3;
wy0=1;
wz0=1e-3;

A0 = [1-2*b20.^2-2*b30.^2 2*(b10.*b20-b00.*b30)  2*(b10.*b30+b00.*b20) ;
2*(b10.*b20+b00.*b30)  1-2*b10.^2-2*b30.^2 2*(b20.*b30-b00.*b10);
2*(b10.*b30-b00.*b20)  2*(b20.*b30+b00.*b10) 1-2*b20.^2-2*b10.^2];

t_end = 200;  % End value for time
Iqq=[0.1 0 0; 0 1 0;0 0 10];

% Initial conditions
initial_conditions = [wx0; wy0; wz0;b00; b10; b20;b30];

% Time span
tspan = 0:0.0001:t_end;

function S = skew(v)
    S = [0 -v(3) v(2); v(3) 0 -v(1); -v(2) v(1) 0];
end

% Define the ODE system as a function
function [dydt,UC] = odesystem(t, F)
global wx0 wy0 wz0 xd0 yd0 zd0;
    wx = F(1);
    wy = F(2);
    wz = F(3);
    b0 = F(4);
    b1 = F(5);
    b2 = F(6);
    b3 = F(7);

Gint=[-b1 b0  b3 -b2;
-b2 -b3  b0 b1;
-b3 b2  -b1 b0];

Iqq=[0.1 0 0; 0 1 0;0 0 10];

A = [1-2*b2^2-2*b3^2 2*(b1*b2-b0*b3)  2*(b1*b3+b0*b2) ;
2*(b1*b2+b0*b3)  1-2*b1^2-2*b3^2 2*(b2*b3-b0*b1);
2*(b1*b3-b0*b2)  2*(b2*b3+b0*b1) 1-2*b2^2-2*b1^2];

w_vec=[wx;
wy;
wz];

w_hat=skew(w_vec);
torques= w_hat*Iqq*w_vec;
NewRes=-Iqq\torques;

b_ans=0.5*Gint'*w_vec;

 dydt = zeros(7,1); 
    dydt(1) = NewRes(1); % wx'
    dydt(2) = NewRes(2); % wy' 
    dydt(3) = NewRes(3); % wz' 
    dydt(4) = b_ans(1) ; % b0
    dydt(5) = b_ans(2); % b1 
    dydt(6) = b_ans(3); % b2
    dydt(7) = b_ans(4); % b3

    %Unit constraint
UC  = 1 - b0^2  - b1^2 - b2^2 - b3^2;
end


%options = odeset('RelTol', 1e-14, 'AbsTol', 1e-14);
options = odeset('RelTol', 100*eps, 'AbsTol', 1e-16);
% Solve the ODE using ode45 with the anonymous function
[t, F] = ode45(@odesystem, tspan, initial_conditions, options);


% Extract results from Y matrix
wx = F(:,1);
wy = F(:,2);
wz = F(:,3);
b0 = F(:,4);
b1 = F(:,5);
b2 = F(:,6);
b3 = F(:,7);

UC2 = 1 - b0.^2  - b1.^2 - b2.^2 - b3.^2;

w_vec0=[wx0;
wy0;
wz0];
Ener0=0.5*w_vec0'*Iqq*w_vec0;
H_body0 = [Iqq(1,1)*wx0, Iqq(2,2)*wy0, Iqq(3,3)*wz0];
H_norm0 = sqrt(sum(H_body0.^2, 2));

H_norm = zeros(length(t), 1);
Ener = zeros(length(t), 1);
for i = 1:length(t)
w_vec=[wx(i);
wy(i);
wz(i)];
H_body = [Iqq(1,1)*wx(i), Iqq(2,2)*wy(i), Iqq(3,3)*wz(i)];
H_norm(i) = sqrt(sum(H_body.^2, 2)) - H_norm0; 
Ener(i) = 0.5*w_vec'*Iqq*w_vec-Ener0 ;
end

if ~exist('save_filename', 'var') || isempty(save_filename)
    save_filename = 'TennisRef.mat';
end
save(save_filename, '-v7.3');

paperFont = 'Times New Roman';

% Angular velocity
figure;
subplot(3, 1, 1);
plot(t, wx, 'b-', 'LineWidth', 1.5);
ylim([-3.5, 3.5]);
grid on;
set(gca, 'FontName', paperFont, 'XMinorGrid', 'off', 'YMinorGrid', 'off');
ylabel('$\omega_x$ [rad/s]', 'Interpreter', 'latex', 'FontName', paperFont);

subplot(3, 1, 2);
plot(t, wy, 'r-', 'LineWidth', 1.5);
ylim([-1.5, 1.5]);
grid on;
set(gca, 'FontName', paperFont, 'XMinorGrid', 'off', 'YMinorGrid', 'off');
ylabel('$\omega_y$ [rad/s]', 'Interpreter', 'latex', 'FontName', paperFont);

subplot(3, 1, 3);
plot(t, wz, 'g-', 'LineWidth', 1.5);
ylim([-0.1, 0.2]);
grid on;
set(gca, 'FontName', paperFont, 'XMinorGrid', 'off', 'YMinorGrid', 'off');
xlabel('Time [s]', 'FontName', paperFont);
ylabel('$\omega_z$ [rad/s]', 'Interpreter', 'latex', 'FontName', paperFont);

% Euler parameters (quaternion components)
figure;
subplot(4, 1, 1);
plot(t, b0, 'b-', 'LineWidth', 1.5);
ylim([-1, 1]);
grid on;
set(gca, 'FontName', paperFont, 'XMinorGrid', 'off', 'YMinorGrid', 'off');
ylabel('$b_0$ [-]', 'Interpreter', 'latex', 'FontName', paperFont);

subplot(4, 1, 2);
plot(t, b1, 'r-', 'LineWidth', 1.5);
ylim([-1, 1]);
grid on;
set(gca, 'FontName', paperFont, 'XMinorGrid', 'off', 'YMinorGrid', 'off');
ylabel('$b_1$ [-]', 'Interpreter', 'latex', 'FontName', paperFont);

subplot(4, 1, 3);
plot(t, b2, 'g-', 'LineWidth', 1.5);
ylim([-1, 1]);
grid on;
set(gca, 'FontName', paperFont, 'XMinorGrid', 'off', 'YMinorGrid', 'off');
ylabel('$b_2$ [-]', 'Interpreter', 'latex', 'FontName', paperFont);

subplot(4, 1, 4);
plot(t, b3, 'c-', 'LineWidth', 1.5);
ylim([-1, 1]);
grid on;
set(gca, 'FontName', paperFont, 'XMinorGrid', 'off', 'YMinorGrid', 'off');
xlabel('Time [s]', 'FontName', paperFont);
ylabel('$b_3$ [-]', 'Interpreter', 'latex', 'FontName', paperFont);

% EP unit constraint (algebraic check on trajectory)
figure;
plot(t, UC2, 'b-', 'LineWidth', 0.75);
ylim([-1e-13, 1e-13]);
grid on;
set(gca, 'FontName', paperFont, 'XMinorGrid', 'off', 'YMinorGrid', 'off');
xlabel('Time [s]', 'FontName', paperFont);
ylabel('EP unit constraint deviation [-]', 'FontName', paperFont);

% Energy and angular momentum (deviations from initial invariants)
figure;
subplot(2, 1, 1);
plot(t, Ener, 'b-', 'LineWidth', 1);
ylim([-1e-13, 1e-13]);
grid on;
set(gca, 'FontName', paperFont, 'XMinorGrid', 'off', 'YMinorGrid', 'off');
ylabel('Energy deviation [J]', 'FontName', paperFont);

subplot(2, 1, 2);
plot(t, H_norm, 'r-', 'LineWidth', 1);
ylim([-1e-13, 1e-13]);
grid on;
set(gca, 'FontName', paperFont, 'XMinorGrid', 'off', 'YMinorGrid', 'off');
xlabel('Time [s]', 'FontName', paperFont);
ylabel('Angular momentum norm deviation [J s]', 'FontName', paperFont);

