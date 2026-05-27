% --- timing: tic/toc now brackets only the custom_rk4 integration call ---
% Parameters
b00=1;
b10=0;
b20=0;
b30=0;

wx0=1e-3;
wy0=1;
wz0=1e-3;
w_vec0=[wx0;
wy0;
wz0];

A0 = [1-2*b20.^2-2*b30.^2 2*(b10.*b20-b00.*b30)  2*(b10.*b30+b00.*b20) ;
2*(b10.*b20+b00.*b30)  1-2*b10.^2-2*b30.^2 2*(b20.*b30-b00.*b10);
2*(b10.*b30-b00.*b20)  2*(b20.*b30+b00.*b10) 1-2*b20.^2-2*b10.^2];

Gint0=[-b10 b00  b30 -b20;
-b20 -b30  b00 b10;
-b30 b20  -b10 b00];

bd_init=0.5*Gint0'*w_vec0;

t_end = 200;  % End value for time
Iqq=[0.1 0 0; 0 1 0;0 0 10];
H_body0 = [Iqq(1,1)*wx0, Iqq(2,2)*wy0, Iqq(3,3)*wz0];
H_norm0 = sqrt(sum(H_body0.^2, 2));

% Initial conditions
initial_conditions = [b00; b10; b20;b30;bd_init];

% Time step (seconds); batch sets dt and save_filename before run()
if ~exist('dt', 'var') || isempty(dt)
    dt = 0.01;  
end
tspan = 0:dt:t_end;

function S = skew(v)
    S = [0 -v(3) v(2); v(3) 0 -v(1); -v(2) v(1) 0];
end

% Define the ODE system as a function
function [dydt,UC, Ener,H_norm, theta,omegaint] = odesystem(t, F)
global wx0 wy0 wz0 H_norm0;
    b0 = F(1);
    b1 = F(2);
    b2 = F(3);
    b3 = F(4);
    bd0 = F(5);
    bd1 = F(6);
    bd2 = F(7);
    bd3 = F(8);

betasd=[bd0 bd1 bd2 bd3];

Gint=2*[-b1 b0  b3 -b2;
-b2 -b3  b0 b1;
-b3 b2  -b1 b0];
G2=2*[-bd1 bd0  bd3 -bd2;
-bd2 -bd3  bd0 bd1;
-bd3 bd2  -bd1 bd0];


Iqq=[0.1 0 0; 0 1 0;0 0 10];

A = [1-2*b2^2-2*b3^2 2*(b1*b2-b0*b3)  2*(b1*b3+b0*b2) ;
2*(b1*b2+b0*b3)  1-2*b1^2-2*b3^2 2*(b2*b3-b0*b1);
2*(b1*b3-b0*b2)  2*(b2*b3+b0*b1) 1-2*b2^2-2*b1^2];

omegaint=Gint*betasd';
Qv1=cross(omegaint,Iqq*omegaint)+Iqq*G2*betasd';
Qv=-Gint'*Qv1; 

Mqq=Gint'*Iqq*Gint;
Cq = [2.*b0 2.*b1 2.*b2 2.*b3 ];

Mat = [Mqq Cq';
       Cq 0];

Forces=[Qv;-2.*bd1.*bd1-2.*bd2.*bd2-2.*bd0.*bd0-2.*bd3.*bd3];
NewRes=Mat\Forces;

 dydt = zeros(8,1); 
    dydt(1) = bd0; % b0'=b0d
    dydt(2) = bd1; % b1'=b1d 
    dydt(3) = bd2; % b2'=b2d 
    dydt(4) = bd3; % b3'=b3d
    dydt(5) = NewRes(1); % b0d
    dydt(6) = NewRes(2); % b1d 
    dydt(7) = NewRes(3); % b2d
    dydt(8) = NewRes(4); % b3d

    %Unit constraint
UC  = 1 - b0^2  - b1^2 - b2^2 - b3^2;
Ener = 0.5*omegaint'*Iqq*omegaint - 0.500005050000000;
H_body = [Iqq(1,1)*omegaint(1), Iqq(2,2)*omegaint(2), Iqq(3,3)*omegaint(3)];
H_norm = sqrt(sum(H_body.^2, 2))-1.000050003749813;
%just cheking crv based on EP - to validate
theta=2*acos(b0);
if theta > pi                     
    theta = 2 * acos(-b0);
end
end

% Custom Runge-Kutta 4th order method
function [t, F, UC,Ener,H_norm,theta,omegaint] = custom_rk4(odefun, tspan, y0)
    n = length(tspan);
    F_temp = zeros(length(y0), n);
    F_temp(:,1) = y0;

    UC = zeros(1,n);
    Ener = zeros(1,n);
    H_norm=zeros(1,n);
    theta = zeros(1,n);
    omegaint=zeros(3,n);

    for i = 1:(n-1)
        h = tspan(i+1) - tspan(i);
        ti = tspan(i);
        yi = F_temp(:,i);
        
        [k1, UC1,Ener1,H_norm1,theta1,omegaint1] = odefun(ti, yi);
        [k2, ~] = odefun(ti + h/2, yi + h*k1/2);
        [k3, ~] = odefun(ti + h/2, yi + h*k2/2);
        [k4, ~] = odefun(ti + h, yi + h*k3);
        
        F_temp(:,i+1) = yi + (h/6)*(k1 + 2*k2 + 2*k3 + k4);
        UC(:,i) = UC1;
        Ener(:,i) = Ener1;
        H_norm(:,i) = H_norm1;
        theta(:,i) = theta1;
        omegaint(:,i) = omegaint1;
    end
    [~,UC_final,Ener_final,H_norm_final,theta_final,omegaint_final] = odefun(tspan(end), F_temp(:,end));
    UC(:,end) = UC_final;
    Ener(:,end) = Ener_final;
    H_norm(:,end) = H_norm_final;
    theta(:,end) = theta_final;
    omegaint(:,end) = omegaint_final;

    t = tspan;
    F = F_temp';
    UC=UC';   
    Ener=Ener';
    H_norm=H_norm';
    theta=theta'; 
    omegaint=omegaint';
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
    [t, F, UC,Ener,H_norm, theta,omegaint] = custom_rk4(@odesystem, tspan, initial_conditions);
    timing_samples(timing_rep) = toc;  % integration only
end
executionTime = median(timing_samples);  % robust to run-to-run noise

% Extract results from Y matrix
b0 = F(:,1);
b1 = F(:,2);
b2 = F(:,3);
b3 = F(:,4);
bd0 = F(:,5);
bd1 = F(:,6);
bd2 = F(:,7);
bd3 = F(:,8);

wx=omegaint(:,1);
wy=omegaint(:,2);
wz=omegaint(:,3);

UC2 = 1 - b0.^2  - b1.^2 - b2.^2 - b3.^2;

% executionTime captured at the custom_rk4 call above (timing harness)
if ~exist('save_filename', 'var') || isempty(save_filename)
    save_filename = sprintf('ClassicEP_dt_%0.1fms.mat', dt * 1000);
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
% plot(t,UC,'b-','LineWidth',0.5); 
% ylim([-1*10^-1,1*10^-1]);
% title('Unit costraint - Classic EP');
% ylabel('Violation');
% grid on;
% 
% 
% % Plotting results for theta
% figure;
% plot(t,theta,'b-','LineWidth',2); 
% ylim([-0.5,3.5]);
% title('Wrapping - Derived from EP - Theta');
% xlabel('Time (s)');
% ylabel('Norm of rotation vector Theta(t)');
% grid on;
% 
% 
% % Plotting results for Energy and Momentum
% figure;
% subplot(2,1,1);
% plot(t,Ener,'b-','LineWidth',1); 
% ylim([-1*10^-3,1*10^-3]);
% title('Total Energy');
% ylabel('Energy balance');
% grid on;
% 
% subplot(2,1,2);
% plot(t,H_norm,'r-','LineWidth',1); 
% ylim([-1*10^-4,1*10^-4]);
% title('Angular momentum(norm)');
% xlabel('Time (s)');
% ylabel('Ang moment');
% grid on;
