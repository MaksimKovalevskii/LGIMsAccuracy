% --- timing: tic/toc now brackets only the custom_rk4 integration call ---
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
H_body0 = [Iqq(1,1)*wx0, Iqq(2,2)*wy0, Iqq(3,3)*wz0];
H_norm0 = sqrt(sum(H_body0.^2, 2));

% Initial conditions
initial_conditions = [wx0; wy0; wz0;b00; b10; b20;b30];

% Time step (seconds)
if ~exist('dt','var') || isempty(dt)
    dt = 0.001;
end

% Time span
tspan = 0:dt:t_end;

function S = skew(v)
    S = [0 -v(3) v(2); v(3) 0 -v(1); -v(2) v(1) 0];
end

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
function [dydt,UC, Ener,H_norm, theta] = odesystem(t, F)
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
Ener = 0.5*w_vec'*Iqq*w_vec - 0.500005050000000;
H_body = [Iqq(1,1)*wx, Iqq(2,2)*wy, Iqq(3,3)*wz];
H_norm = sqrt(sum(H_body.^2, 2))-1.000050003749813; 
theta=2*acos(b0);
end

% Custom Runge-Kutta 4th order method
function [t, F, UC,Ener,H_norm,theta] = custom_rk4(odefun, tspan, y0)
    n = length(tspan);
    F_temp = zeros(length(y0), n);
    F_temp(:,1) = y0;

    UC = zeros(1,n);
        Ener = zeros(1,n);
        H_norm=zeros(1,n);
        theta = zeros(1,n);

    for i = 1:(n-1)
        h = tspan(i+1) - tspan(i);
        ti = tspan(i);
        yi = F_temp(:,i);
         
        % Split into q (orientations) and W (ang velocities)
        qi = yi(4:7);  % Rotational parameters: [b0,b1,b2,b3]
        vi = yi(1:3); % Velocities: [wx,wy,wz]
               
% Stage 1
[k1, UC1,Ener1,H_norm1,theta1] = odefun(ti, yi);
        k1_q = h*yi(1:3);
        k1_v = h*k1(1:3);
% Stage 2
q2 = LGT(qi, k1_q/2);
v2 = vi + (1/2)*k1_v;
y2 = [v2;q2];

[k2, ~] = odefun(ti + h/2, y2);
        k2_q =  h*invdexp(v2,k1_q/2);
        k2_v = h*k2(1:3);

        % Stage 3
q3 = LGT(qi, k2_q/2);
v3 = vi + (1/2)*k2_v;
y3 = [v3;q3];  

        [k3,  ~] = odefun(ti + h/2, y3);
        k3_q =  h*invdexp(v3,k2_q/2);
        k3_v = h*k3(1:3);
% Stage 4   
q4 = LGT(qi, k3_q);
v4 = vi + k3_v;
y4 = [v4;q4];

        [k4, ~] = odefun(ti + h, y4);
        k4_q =  h*invdexp(v4,k3_q);
        k4_v = h*k4(1:3);

Delta_q = (1/6)*(k1_q +2*k2_q +2*k3_q +k4_q);
q_new = LGT(qi, Delta_q);
v_new = vi + (1/6)*(k1_v + 2*k2_v + 2*k3_v + k4_v);
        
        % Store full state 
        F_temp(:,i+1) = [v_new;q_new];

       UC(:,i) = UC1;
       Ener(:,i) = Ener1;
       H_norm(:,i) = H_norm1;
       theta(:,i) = theta1;
    end
    [~,UC_final,Ener_final,H_norm_final,theta_final] = odefun(tspan(end), F_temp(:,end));
    UC(:,end) = UC_final;
    Ener(:,end) = Ener_final;
    H_norm(:,end) = H_norm_final;
    theta(:,end) = theta_final;

    t = tspan;
    F = F_temp';
    UC=UC';    
    Ener=Ener';
    H_norm=H_norm';
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
    [t, F, UC,Ener,H_norm, theta] = custom_rk4(@odesystem, tspan, initial_conditions);
    timing_samples(timing_rep) = toc;  % integration only
end
executionTime = median(timing_samples);  % robust to run-to-run noise

% Extract results from Y matrix
wx = F(:,1);
wy = F(:,2);
wz = F(:,3);
b0 = F(:,4);
b1 = F(:,5);
b2 = F(:,6);
b3 = F(:,7);

UC2 = 1 - b0.^2  - b1.^2 - b2.^2 - b3.^2;

% executionTime captured at the custom_rk4 call above (timing harness)

if ~exist('save_filename','var') || isempty(save_filename)
    save_filename = sprintf('EP_LGIM_dt_%0.1fms.mat', dt*1000);
end
save(save_filename, '-v7.3');

% % Plotting results for Angular velocities
% figure;
% subplot(3,1,1);
% plot(t,wx,'b-','LineWidth',2); 
% ylim([-3.5,3.5]);
% title('Wx');
% ylabel({'Ang velocity Wx'; 's-1'});
% grid on;
% 
% subplot(3,1,2);
% plot(t,wy,'r-','LineWidth',2); 
% ylim([-1.5,1.5]);
% title('Wy');
% ylabel({'Ang velocity Wy'; 's-1'});
% grid on;
% 
% subplot(3,1,3);
% plot(t,wz,'g-','LineWidth',2); 
% ylim([-0.1,0.2]);
% title('Wz');
% ylabel({'Ang velocity Wz';'s-1'});
% grid on;
% 
% % Plotting results for EP
% figure;
% subplot(4,1,1);
% plot(t,b0,'b-','LineWidth',2); 
% ylim([-1,1]);
% title('b0');
% %ylabel({'b0'});
% grid on;
% 
% subplot(4,1,2);
% plot(t,b1,'r-','LineWidth',2); 
% ylim([-1,1]);
% title('b1');
% %ylabel({'b1'});
% grid on;
% 
% subplot(4,1,3);
% plot(t,b2,'g-','LineWidth',2); 
% ylim([-1,1]);
% title('b2');
% %ylabel({'b2'});
% grid on;
% 
% subplot(4,1,4);
% plot(t,b3,'c-','LineWidth',2); 
% ylim([-1,1]);
% title('b3');
% %ylabel({'b3'});
% grid on;
% 
% % Plotting results for Unit constraints
% figure;
% subplot(2,1,1);
% plot(t,UC,'b-','LineWidth',0.5); 
% ylim([-1*10^-12,1*10^-12]);
% title('Unit costraint - EP Lie');
% ylabel('Violation');
% grid on;
% 
% subplot(2,1,2);
% plot(t,UC2,'b-','LineWidth',0.5); 
% ylim([-1*10^-12,1*10^-12]);
% title('Unit costraint - EP Lie');
% ylabel('Violation');
% grid on;
% 
% % Plotting results for theta
% figure;
% plot(t,theta,'b-','LineWidth',2); 
% ylim([-1,7]);
% title('Derived from EP - Angle Theta');
% xlabel('Time (s)');
% ylabel('Norm of rotation vector Theta(t)');
% grid on;
% 
% 
% % Plotting results for Energy and Momentum
% figure;
% subplot(2,1,1);
% plot(t,Ener,'b-','LineWidth',1); 
% ylim([-5*10^-2,5*10^-2]);
% title('Total Energy');
% ylabel('Energy balance');
% grid on;
% 
% subplot(2,1,2);
% plot(t,H_norm,'r-','LineWidth',1); 
% ylim([-5*10^-2,5*10^-2]);
% title('Angular momentum(norm)');
% xlabel('Time (s)');
% ylabel('Ang moment');
% grid on;
