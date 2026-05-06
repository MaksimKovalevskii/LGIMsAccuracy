tic;
% Initial values for psi - rotational part
psi10=0.55536036727;
psi20=0.55536036727;
psi30=0;
t_end = 100;  % End value for time
if ~exist('dt', 'var') || isempty(dt)
    dt = 0.001;   % Time step (default 1 ms)
end

% Inertia tensor
I = [3.0025 0 0; 0 3.0025 0; 0 0 0.005];

% Initial conditions
initial_conditions = [1.5; -1.5; 2.12132034356;psi10;psi20;psi30;0;0;0;0;0;0];

% Time span
tspan = 0:dt:t_end;

function S = skew(v)
    S = [0 -v(3) v(2); v(3) 0 -v(1); -v(2) v(1) 0];
end

function [dydt, second_derivatives, C_val, Cq_qd_val, C_ddot, theta]= odesystem(t, F)
    x = F(1);
    y = F(2);
    z = F(3);
    psi1 = F(4);
    psi2 = F(5);
    psi3 = F(6);
    xd = F(7);
    yd = F(8);
    zd = F(9);
    psi1d = F(10);
    psi2d = F(11);
    psi3d = F(12);

theta = sqrt(psi1^2+psi2^2+psi3^2);

q=[x y z psi1 psi2 psi3]; % Vector of generalised coord q
qd=[xd yd zd psi1d psi2d psi3d]; % Vector of generalised coord q - first derivative
Psid=[psi1d psi2d psi3d];% Cartesian rotational vector
PsiHatDot = skew(Psid);
psi = [psi1 psi2 psi3];
PsiHat = skew(psi);

I=[3.0025 0 0; 0 3.0025 0;0 0 0.005]; %Moment of Inertia

Gint=GmatrixCart(x, y, z, psi1, psi2, psi3);%Function for matrix G

phi=theta;
 if abs(phi) < 1e-8
        % Small angle: G ≈ I - 0.5*PsiHat + (1/12)*PsiHat^2
        GhatDot = -0.5*PsiHatDot + (1/12)*(PsiHatDot*PsiHat + PsiHat*PsiHatDot);
    else
        cos_phi = cos(phi);
        sin_phi = sin(phi);
        phi_dot = (psi' * Psid) / phi;
        
        % Coefficients
        c1 = (1 - cos_phi) / phi^2;
        c2 = (phi - sin_phi) / phi^3;
        
        % Time derivatives of coefficients
        c1_dot = (phi * sin_phi - 2 * (1 - cos_phi)) * phi_dot / phi^3;
        c2_dot = (1 - cos_phi - phi * sin_phi + 2 * (phi - sin_phi)) * phi_dot / phi^4;
        
        GhatDot = -PsiHatDot * c1 - PsiHat * c1_dot + ...
                  (PsiHatDot * PsiHat + PsiHat * PsiHatDot) * c2 + ...
                  PsiHat * PsiHat * c2_dot;
    end

omegaint=Gint*Psid';
temp1=I*GhatDot*Psid';
Qv1=cross(omegaint,I*omegaint)+temp1;
Qv=-Gint'*Qv1; %force vector associated with the orientation coordinates

Mqq=Gint'*I*Gint;
M = [eye(3) zeros(3,3);
    zeros(3,3) Mqq]; %Mass Matrix

%Jacobian matrix calc - separate for phi close to 0, and for phi>10^-3
if sqrt(psi1^2+psi2^2+psi3^2)<10^-3 
    Cq = Cq_zero(q); %symbolic function to decrease error
    Qc = Qc_zero(q,qd);
    C_val = C_zero(q); 
    Cq_qd_val = Cq_qd_zero(q, qd);
else
    Cq = Cq_func(q); %symbolic function to decrease error
    Qc = Qc_func(q,qd);%symbolic function to decrease error
    C_val = C_func(q);
    Cq_qd_val = Cq_qd_func(q, qd);
end

%Full Matrix from 4.65 Aki's notes
Mat = [M Cq';
       Cq zeros(3)];

Forces=[0; 0; -9.8;Qv;Qc];
NewRes=Mat\Forces; %Expression for q'' and lagrange multipliers

    dydt = zeros(12,1);
    dydt(1) = xd; % x' = xd
    dydt(2) = yd; % y' = yd
    dydt(3) = zd; % z' = zd
    dydt(4) = psi1d; 
    dydt(5) = psi2d; 
    dydt(6) = psi3d; 
    dydt(7) = NewRes(1); % xd' (x'')
    dydt(8) = NewRes(2);  % yd' (y'')  
    dydt(9) = NewRes(3); % zd' (z'') 
    dydt(10) = NewRes(4); 
    dydt(11) = NewRes(5);  
    dydt(12) = NewRes(6);  

second_derivatives = NewRes(1:6); % determine 2nd derivatives right here
C_ddot = Cq * second_derivatives - Qc; % Cq*q_ddot + (Cq*qd)_q*qd = Cq*q_ddot - Qc
end

% Custom Runge-Kutta 4th order method
function [t, F, second_derivatives,C_history, Cq_qd_history, C_ddot_history,theta_history] = custom_rk4(odefun, tspan, y0)
    n = length(tspan);
    F_temp = zeros(length(y0), n);
    F_temp(:,1) = y0;

    % determine 2nd derivatives, constraints C, C',C''
    second_derivatives = zeros(6, n);
    C_history = zeros(3, n);
    Cq_qd_history = zeros(3, n);
    C_ddot_history = zeros(3, n);
    theta_history = zeros (1,n);
    

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

        [k1, sd1, C1, Cqd1, C_ddot1,theta1] = odefun(ti, yi);
        [k2, ~] = odefun(ti + h/2, yi + h*k1/2);
        [k3, ~] = odefun(ti + h/2, yi + h*k2/2);
        [k4, ~] = odefun(ti + h, yi + h*k3);
        
        F_temp(:,i+1) = yi + (h/6)*(k1 + 2*k2 + 2*k3 + k4);
        second_derivatives(:,i) = sd1;  
        C_history(:,i) = C1;
        Cq_qd_history(:,i) = Cqd1;
        C_ddot_history(:,i) = C_ddot1;
        theta_history(:,i) = theta1;
    
    end
     % just for end point
    [~, sd_final,C_final,Cqd_final,C_ddot_final,theta_final] = odefun(tspan(end), F_temp(:,end));
    second_derivatives(:,end) = sd_final;
    C_history(:,end) = C_final;
    Cq_qd_history(:,end) = Cqd_final;
    C_ddot_history(:,end) = C_ddot_final;
    theta_history(:,end) = theta_final;

    t = tspan;
    F = F_temp';
    second_derivatives = second_derivatives';
    C_history = C_history';
    Cq_qd_history = Cq_qd_history';
    C_ddot_history = C_ddot_history';
    theta_history = theta_history';

end

% Solve the ODE using custom RK4
[t, F, second_derivatives, C_history, Cq_qd_history, C_ddot_history,theta_history] = custom_rk4(@odesystem, tspan, initial_conditions);
theta_history (1,1) = sqrt(psi10^2+psi20^2+psi30^2);

% Extract results from F matrix
    x = F(:,1);
    y = F(:,2);
    z = F(:,3);
    psi1 = F(:,4);
    psi2 = F(:,5);
    psi3 = F(:,6);

    xd = F(:,7);
    yd = F(:,8);
    zd = F(:,9);
    psi1d = F(:,10);
    psi2d = F(:,11);
    psi3d= F(:,12);

x_double_prime = second_derivatives(:,1);
y_double_prime = second_derivatives(:,2);
z_double_prime = second_derivatives(:,3);
psi1_double_prime = second_derivatives(:,4);
psi2_double_prime = second_derivatives(:,5);
psi3_double_prime = second_derivatives(:,6);

% Calculate norms of Constraint violations' matrices for plotting
C = max(abs(C_history), [], 2);  
dC = max(abs(Cq_qd_history), [], 2); 
ddC = max(abs(C_ddot_history), [], 2);

% Calculate Energy balance for plotting
for i = 1:length(t)
q = [x(i) y(i) z(i) psi1(i) psi2(i) psi3(i)];

    qd= [xd(i) yd(i) zd(i) psi1d(i) psi2d(i) psi3d(i)];

Gint = GmatrixCart(x(i), y(i), z(i), psi1(i), psi2(i), psi3(i));

Mqq=Gint'*I*Gint;
M = [eye(3) zeros(3,3);
    zeros(3,3) Mqq]; %Mass Matrix

Ener2(i)=0.5 * (qd * M * qd')+9.8.*(z(i)-z(1));
end
Enernew=Ener2';

executionTime = toc
if ~exist('save_filename', 'var') || isempty(save_filename)
    save_filename = sprintf('Cart_dt_%.1fms.mat', dt * 1000);
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
% % Plotting results for psi1 and psi2
% figure;
% subplot(3,1,1);
% plot(t,psi1,'b-','LineWidth',2); 
% hold on;
% plot(t,psi2,'r--','LineWidth',2); 
% ylim([-3,3]);
% title('\psi_1(t), \psi_2(t)');
% ylabel('\psi_1(t), \psi_2(t)');
% grid on;
% lgd = legend('\psi_1', '\psi_2', 'Location', 'northeast', 'FontSize', 8, 'TextColor', 'blue');
% 
% subplot(3,1,2);
% plot(t,psi1d,'b-','LineWidth',2); 
% hold on;
% plot(t,psi2d,'r--','LineWidth',2); 
% ylim([-3,3]);
% title('\psi_1''(t), \psi_2''(t)');
% ylabel('1st Derivatives');
% grid on;
% 
% subplot(3,1,3);
% plot(t,psi1_double_prime,'b-','LineWidth',2); 
% hold on;
% plot(t,psi2_double_prime,'r--','LineWidth',2); 
% ylim([-3,3]);
% title('\psi_1''''(t), \psi_2''''(t)');
% xlabel('Time (s)');
% ylabel('2nd Derivatives');
% grid on;
% 
% % Plotting results for theta
% figure;
% plot(t,theta_history,'b-','LineWidth',2); 
% ylim([0,3.5]);
% title('Angle \theta - norm of rotation vector \psi');
% xlabel('Time (s)');
% ylabel('\theta(t)');
% grid on;
% 
% % Plotting results for psi3
% figure;
% subplot(3,1,1);
% plot(t,psi3,'b-','LineWidth',2); 
% ylim([-1,1]);
% title('psi3(t)');
% ylabel('Psi3');
% grid on;
% 
% subplot(3,1,2);
% plot(t,psi3d,'r-','LineWidth',2); 
% ylim([-3,3]);
% title('Psi3''(t)');
% ylabel('1st Derivative Psi3');
% grid on;
% 
% subplot(3,1,3);
% plot(t,psi3_double_prime,'g-','LineWidth',2); 
% ylim([-3,3]);
% title('Psi3''''(t)');
% xlabel('Time (s)');
% ylabel('2nd Derivative Psi3');
% grid on;
% 
% % Plotting results for Constraint violation and Energy balance
% figure;
% subplot(4,1,1);
% plot(t,C,'b-','LineWidth',0.5); 
% ylim([-2*10^-2,2*10^-2]);
% title('∣∣C∣∣');
% ylabel('Violation');
% grid on;
% 
% subplot(4,1,2);
% plot(t,dC,'g-','LineWidth',1); 
% ylim([-2*10^-2,2*10^-2]);
% title('∣∣C''∣∣');
% ylabel('Violation');
% grid on;
% 
% subplot(4,1,3);
% plot(t,ddC,'r-','LineWidth',1); 
% ylim([-5*10^-14,5*10^-14]);
% title('∣∣C''''∣∣');
% ylabel('Violation');
% grid on;
% 
% subplot(4,1,4);
% plot(t,Enernew,'b-','LineWidth',1); 
% ylim([-5*10^-2,5*10^-2]);
% title('Total Energy');
% xlabel('Time (s)');
% ylabel('Energy balance');
% grid on;
% 
