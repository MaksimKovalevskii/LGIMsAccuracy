% tic;
% Parameters
psi10=0;
psi20=0;
psi30=0;

wx0=1e-3;
wy0=1;
wz0=1e-3;

t_end = 200;  % End value for time
Iqq=[0.1 0 0; 0 1 0;0 0 10];
H_body0 = [Iqq(1,1)*wx0, Iqq(2,2)*wy0, Iqq(3,3)*wz0];
H_norm0 = sqrt(sum(H_body0.^2, 2));


% Initial conditions
initial_conditions = [wx0; wy0; wz0;psi10; psi20; psi30];

% Time step (seconds)
if ~exist('dt','var') || isempty(dt)
    dt = 0.001;
end
dt = 0.001;
% Time span
tspan = 0:dt:t_end;

function S = skew(v)
    S = [0 -v(3) v(2); v(3) 0 -v(1); -v(2) v(1) 0];
end

% Define the ODE system as a function
function [dydt,theta, Ener,H_norm,OrtDet,OrtIdent] = odesystem(t, F)
global wx0 wy0 wz0;
    wx = F(1);
    wy = F(2);
    wz = F(3);
    psi1 = F(4);
    psi2 = F(5);
    psi3 = F(6);

Iqq=[0.1 0 0; 0 1 0;0 0 10];

psi = [psi1 psi2 psi3];
PsiHat = skew(psi);
%  Determining angle phi (theta in report) from psi
phi=sqrt(psi1^2+psi2^2+psi3^2);

%  Matrix G (Emap here) from Rodriguez formula (with Taylor series for small phi)
if abs(phi) < 10^-4
Emap = eye(3) + PsiHat*(-0.5+phi^2/24-phi^4/720+phi^6/40320) + PsiHat*PsiHat*(1/6-phi^2/120+phi^4/5040-phi^6/362880);
Emap2 = eye(3) + 0.5*PsiHat + PsiHat*PsiHat*(1/12 + phi^2/720);
else
Emap = eye(3) - PsiHat*(1-cos(phi))/(phi^2) + PsiHat*PsiHat*(phi-sin(phi))/(phi^3);
Emap2 = eye(3) + 0.5*PsiHat + PsiHat*PsiHat*(1-0.5*phi*sin(phi)/(1-cos(phi)))/(phi^2);
end

%  Transform matrix A and A2 (for small phi) from Rodriguez formula
if phi > 10^-6
A = eye(3) + PsiHat*(sin(phi)/phi) + 2*PsiHat*PsiHat*sin(phi/2)*sin(phi/2)/(phi*phi);
else
A = eye(3) + PsiHat*(1 - phi^2/6 + phi^4/120) + 2*PsiHat*PsiHat*(1/2 - phi^2/48 + phi^4/3840);
end

w_vec=[wx;
wy;
wz];

psi_ans=Emap2*w_vec;

w_hat=skew(w_vec);
torques= w_hat*Iqq*w_vec;
NewRes=-Iqq\torques;

 dydt = zeros(6,1); 
    dydt(1) = NewRes(1); % wx'
    dydt(2) = NewRes(2); % wy' 
    dydt(3) = NewRes(3); % wz' 
    dydt(4) = psi_ans(1) ; %psi1
    dydt(5) = psi_ans(2); %psi2 
    dydt(6) = psi_ans(3); %psi3

theta=phi;
Ener = 0.5*w_vec'*Iqq*w_vec - 0.500005050000000;
H_body = [Iqq(1,1)*wx, Iqq(2,2)*wy, Iqq(3,3)*wz];
H_norm = sqrt(sum(H_body.^2, 2))-1.000050003749813; 
OrtDet=abs(det(A) - 1);
OrtIdent=norm(A'*A - eye(3));
end

% Custom Runge-Kutta 4th order method
function [t, F, theta,Ener,H_norm,OrtDet,OrtIdent] = custom_rk4(odefun, tspan, y0)
    n = length(tspan);
    F_temp = zeros(length(y0), n);
    F_temp(:,1) = y0;

    theta = zeros(1,n);
    Ener = zeros(1,n);
    H_norm=zeros(1,n);
    OrtDet=zeros(1,n);
    OrtIdent=zeros(1,n);

    for i = 1:(n-1)
        h = tspan(i+1) - tspan(i);
        ti = tspan(i);
        yi = F_temp(:,i);
  
        % --- WRAP at the BEGINNING of each RK4 step ---
       psi = yi(4:6);
       Phi = norm(psi);
      if Phi > pi
fprintf('Wrapping RK4: Phi = %.3f -> ', Phi)
           psi = -(2*pi - Phi) * (psi / Phi);
        yi(4:6) = psi;  % Update yi with wrapped psi
    end

        [k1, theta1,Ener1,H_norm1,OrtDet1,OrtIdent1] = odefun(ti, yi);
        [k2, ~] = odefun(ti + h/2, yi + h*k1/2);
        [k3, ~] = odefun(ti + h/2, yi + h*k2/2);
        [k4, ~] = odefun(ti + h, yi + h*k3);
        
        F_temp(:,i+1) = yi + (h/6)*(k1 + 2*k2 + 2*k3 + k4);
        theta(:,i) = theta1;
                Ener(:,i) = Ener1;
                H_norm(:,i) = H_norm1;
                OrtDet(:,i) = OrtDet1;
                OrtIdent(:,i) = OrtIdent1;
    end
    [~,theta_final,Ener_final,H_norm_final,OrtDet_final,OrtIdent_final] = odefun(tspan(end), F_temp(:,end));
    theta(:,end) = theta_final;
    Ener(:,end) = Ener_final;
    H_norm(:,end) = H_norm_final;
    OrtDet(:,end) = OrtDet_final;
    OrtIdent(:,end) = OrtIdent_final;

    t = tspan;
    F = F_temp';
    theta=theta'; 
        Ener=Ener';
        H_norm=H_norm';
        OrtDet=OrtDet';
        OrtIdent=OrtIdent';
end

% Solve the ODE using custom RK4
[t, F, theta,Ener,H_norm,OrtDet,OrtIdent] = custom_rk4(@odesystem, tspan, initial_conditions);

% Extract results from Y matrix
wx = F(:,1);
wy = F(:,2);
wz = F(:,3);
psi1 = F(:,4);
psi2 = F(:,5);
psi3 = F(:,6);

% executionTime = toc;
% if ~exist('save_filename','var') || isempty(save_filename)
%     save_filename = sprintf('wrCart_NE_dt_%0.1fms.mat', dt*1000);
% end
% save(save_filename, '-v7.3');

% Plotting results for Angular velocities
figure;
subplot(3,1,1);
plot(t,wx,'b-','LineWidth',2); 
ylim([-3.5,3.5]);
title('Wx');
ylabel({'Ang velocity Wx'; 's-1'});
grid on;

subplot(3,1,2);
plot(t,wy,'r-','LineWidth',2); 
ylim([-1.5,1.5]);
title('Wy');
ylabel({'Ang velocity Wy'; 's-1'});
grid on;

subplot(3,1,3);
plot(t,wz,'g-','LineWidth',2); 
ylim([-0.1,0.2]);
title('Wz');
ylabel({'Ang velocity Wz';'s-1'});
grid on;

% Plotting results for Cartesian rotation vector
figure;
subplot(3,1,1);
plot(t,psi1,'b-','LineWidth',2); 
ylim([-3.5,3.5]);
title('Psi1');
ylabel({'X rot vector'});
grid on;

subplot(3,1,2);
plot(t,psi2,'r-','LineWidth',2); 
ylim([-3.5,3.5]);
title('Psi2');
ylabel({'Y rot vector'});
grid on;

subplot(3,1,3);
plot(t,psi3,'g-','LineWidth',2); 
ylim([-3.5,3.5]);
title('Psi3');
ylabel({'Z rot vector'});
grid on;

% Plotting results for rotation angle theta (norm of Cartesian rotation vector)
figure;
hold on;
plot(t, theta, 'b-', 'LineWidth', 1.5);
yline(pi, 'r--', 'LineWidth', 1.2);
ylim([-0.5, 3.5]);
grid on;
set(gca, 'FontName', 'Times New Roman', 'FontSize', 11, ...
    'XMinorGrid', 'off', 'YMinorGrid', 'off', ...
    'YTick', [0, 1, pi/2, 2, pi], ...
    'YTickLabel', {'0', '1', '\pi/2', '2', '\pi'}, ...
    'TickLabelInterpreter', 'tex');
xlabel('Time [s]', 'FontName', 'Times New Roman');
% MATLAB LaTeX does not reliably support \boldsymbol / \|; use \Vert (supported).
ylabel('$\theta = \Vert \psi \Vert$ [rad]', 'Interpreter', 'latex', ...
    'FontName', 'Times New Roman');
legend('$\theta(t)$', '$\theta = \pi$', 'Interpreter', 'latex', ...
    'FontName', 'Times New Roman', 'Location', 'best');
hold off;

% Plotting results for Energy and Momentum
figure;
subplot(2,1,1);
plot(t,Ener,'b-','LineWidth',1); 
ylim([-5*10^-1,5*10^-1]);
title('Total Energy');
ylabel('Energy balance');
grid on;

subplot(2,1,2);
plot(t,H_norm,'r-','LineWidth',1); 
ylim([-5*10^-1,5*10^-1]);
title('Angular momentum(norm)');
xlabel('Time (s)');
ylabel('Ang moment');
grid on;

% Plotting results for orthogonality
figure;
subplot(2,1,1);
plot(t,OrtDet,'b-','LineWidth',1); 
ylim([-1*10^-13,1*10^-13]);
title('Rotation matrix Orthogonality - determinant');
ylabel('Num err');
grid on;

subplot(2,1,2);
plot(t,OrtIdent,'r-','LineWidth',1); 
ylim([-1*10^-13,1*10^-13]);
title('Rotation matrix Orthogonality - identity');
xlabel('Time (s)');
ylabel('Num err');
grid on;