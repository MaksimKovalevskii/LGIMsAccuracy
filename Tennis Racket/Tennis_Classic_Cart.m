tic;
% Parameters
psi10=0;
psi20=0;
psi30=0;

wx0=1e-3;
wy0=1;
wz0=1e-3;
w_vec0=[wx0;
wy0;
wz0];

function S = skew(v)
    S = [0 -v(3) v(2); v(3) 0 -v(1); -v(2) v(1) 0];
end

% find initial psi' vector psid0
psi0=[psi10; psi20; psi30];
PsiHat0 = skew(psi0);
phi0 = sqrt(psi10^2+psi20^2+psi30^2);
Emap20 = eye(3) + 0.5*PsiHat0 + PsiHat0*PsiHat0*(1/12 + phi0^2/720);
psid0=Emap20*w_vec0;

Ghat0 = eye(3);  % since psi0 = [0;0;0]
omega0_calc = Ghat0 * psid0;
%disp([omega0_calc, w_vec0]);

t_end = 200;  % End value for time
Iqq=[0.1 0 0; 0 1 0;0 0 10];
H_body0 = [Iqq(1,1)*wx0, Iqq(2,2)*wy0, Iqq(3,3)*wz0];
H_norm0 = sqrt(sum(H_body0.^2, 2));

% Initial conditions
initial_conditions = [psi10; psi20; psi30;psid0];

% Time step (seconds); batch sets dt and save_filename before run()
if ~exist('dt', 'var') || isempty(dt)
    dt = 1e-4; 
end
tspan = 0:dt:t_end;

function G = computeGhat(psi)
    phi = norm(psi);
    S = skew(psi);
    if abs(phi) < 1e-5
        G = eye(3) + S*(-0.5+phi^2/24-phi^4/720+phi^6/40320) ...
                + S*S*(1/6-phi^2/120+phi^4/5040-phi^6/362880);
    else
        G = eye(3) - S*(1-cos(phi))/phi^2 + S*S*(phi-sin(phi))/phi^3;
    end
end


% Define the ODE system as a function
function [dydt, Ener,H_norm, theta,omegaint, conditioning_data] = odesystem(t, F)

global wx0 wy0 wz0 H_norm0;
    psi1= F(1);
    psi2 = F(2);
    psi3 = F(3);
    psid1 = F(4);
    psid2 = F(5);
    psid3 = F(6);
phi = sqrt(psi1^2+psi2^2+psi3^2);
psi = [psi1 psi2 psi3];
PsiHat = skew(psi);
psid=[psid1 psid2 psid3];
PsiHatDot = skew(psid);
Iqq=[0.1 0 0; 0 1 0;0 0 10];
theta=phi;

 % Compute exact Ghat
    Ghat = computeGhat(psi);

 if abs(phi) < 1e-4
        % Small angle: G ≈ I - 0.5*PsiHat + (1/12)*PsiHat^2
        GhatDot = -0.5*PsiHatDot + (1/12)*(PsiHatDot*PsiHat + PsiHat*PsiHatDot);
    else
       phi_dot = (psi * psid') / phi;  % psi^T * psid / phi

    % Coefficients
    c1 = (1 - cos(phi)) / phi^2;
    c2 = (phi - sin(phi)) / phi^3;

    % TIME DERIVATIVES OF COEFFICIENTS
    dc1_dphi = (phi*sin(phi) + 2*cos(phi) - 2) / phi^3;
    dc2_dphi = (-phi*cos(phi) - 2*phi + 3*sin(phi)) / phi^4;

    GhatDot = -PsiHatDot * c1 - PsiHat * (dc1_dphi * phi_dot) ...
              + (PsiHatDot * PsiHat + PsiHat * PsiHatDot) * c2 ...
              + PsiHat * PsiHat * (dc2_dphi * phi_dot);
 end

omegaint=Ghat*psid';
Qv1=cross(omegaint,Iqq*omegaint)+Iqq*GhatDot*psid';
Qv=-Ghat'*Qv1; 

Mqq=Ghat'*Iqq*Ghat;

% === CONDITIONING ANALYSIS ===
cond_Mqq = cond(Mqq);
det_Mqq = det(Mqq);
norm_Mqq = norm(Mqq);
eigs_Mqq = eig(Mqq);
min_eig = min(eigs_Mqq);
max_eig = max(eigs_Mqq);
eig_ratio = max_eig/min_eig;
conditioning_data = [cond_Mqq; det_Mqq; norm_Mqq; min_eig; max_eig; eig_ratio];

NewRes=Mqq\Qv;

 dydt = zeros(6,1); 
    dydt(1) = psid1;
    dydt(2) = psid2; 
    dydt(3) = psid3; 
    dydt(4) = NewRes(1);
    dydt(5) = NewRes(2); 
    dydt(6) = NewRes(3); 

Ener = 0.5*omegaint'*Iqq*omegaint - 0.500005050000000;
H_body = [Iqq(1,1)*omegaint(1), Iqq(2,2)*omegaint(2), Iqq(3,3)*omegaint(3)];
H_norm = sqrt(sum(H_body.^2, 2))-1.000050003749813;

end

% Custom Runge-Kutta 4th order method
function [t, F, Ener,H_norm, theta,omegaint,cond_history] = custom_rk4(odefun, tspan, y0)
    n = length(tspan);
    F_temp = zeros(length(y0), n);
    F_temp(:,1) = y0;
    
    Ener = zeros(1,n);
    H_norm=zeros(1,n);
    theta = zeros(1,n);
    omegaint=zeros(3,n);
cond_history = zeros(6,n);  % [cond, det, norm, min_eig, max_eig, eig_ratio]

    for i = 1:(n-1)
        h = tspan(i+1) - tspan(i);
        ti = tspan(i);
        yi = F_temp(:,i);
  
        [k1,Ener1,H_norm1,theta1,omegaint1,cond1] = odefun(ti, yi);
        [k2, ~] = odefun(ti + h/2, yi + h*k1/2);
        [k3, ~] = odefun(ti + h/2, yi + h*k2/2);
        [k4, ~] = odefun(ti + h, yi + h*k3);
        
        F_temp(:,i+1) = yi + (h/6)*(k1 + 2*k2 + 2*k3 + k4);
        Ener(:,i) = Ener1;
        H_norm(:,i) = H_norm1;
        theta(:,i) = theta1;
        omegaint(:,i) = omegaint1;
        cond_history(:,i) = cond1;
    end
    [~,Ener_final,H_norm_final,theta_final,omegaint_final,cond_final] = odefun(tspan(end), F_temp(:,end));
    theta(:,end) = theta_final;
    Ener(:,end) = Ener_final;
    H_norm(:,end) = H_norm_final;
    theta(:,end) = theta_final;
    omegaint(:,end) = omegaint_final;
    cond_history(:,end) = cond_final;

    t = tspan;
    F = F_temp';
    Ener=Ener';
    H_norm=H_norm';
    theta=theta'; 
    omegaint=omegaint';
    cond_history = cond_history';
end

% Solve the ODE using custom RK4
[t, F, Ener,H_norm, theta,omegaint,cond_history] = custom_rk4(@odesystem, tspan, initial_conditions);

% Extract results from Y matrix
psi1 = F(:,1);
psi2 = F(:,2);
psi3 = F(:,3);
psid1 = F(:,4);
psid2 = F(:,5);
psid3 = F(:,6);

wx=omegaint(:,1);
wy=omegaint(:,2);
wz=omegaint(:,3);

executionTime = toc
if ~exist('save_filename', 'var') || isempty(save_filename)
    % %g keeps 0.1 and 0.125 distinct in the filename
    save_filename = sprintf('ClassicCart_dt_%gms.mat', dt * 1000);
end
save(save_filename, '-v7.3');

%% === CONDITIONING ANALYSIS PLOTS ===
% Extract conditioning data
cond_Mqq = cond_history(:,1);
det_Mqq = cond_history(:,2);
norm_Mqq = cond_history(:,3);
min_eig = cond_history(:,4);
max_eig = cond_history(:,5);
eig_ratio = cond_history(:,6);

% % Plot 1: Condition Number vs Time
% figure;
% subplot(2,1,1);
% semilogy(t, cond_Mqq, 'r-', 'LineWidth', 2);
% title('Matrix Conditioning: Mqq = G^T*I*G');
% ylabel('Condition Number');
% xlabel('Time (s)');
% grid on;
% 
% subplot(2,1,2);
% semilogy(t, abs(det_Mqq), 'b-', 'LineWidth', 2);
% ylabel('|det(Mqq)|');
% xlabel('Time (s)');
% grid on;
% 
% % Plot 2: Condition Number vs Rotation Magnitude
% figure;
% loglog(theta, cond_Mqq, 'r-', 'LineWidth', 2);
% title('Condition Number vs Rotation Magnitude');
% xlabel('θ = ||ψ|| (rad)');
% ylabel('cond(Mqq)');
% grid on;
% 
% % Plot 3: Eigenvalue Evolution
% figure;
% subplot(2,1,1);
% semilogy(t, max_eig, 'r-', 'LineWidth', 2);
% hold on;
% semilogy(t, min_eig, 'b-', 'LineWidth', 2);
% title('Eigenvalue Evolution');
% ylabel('Eigenvalue Magnitude');
% legend('λ_{max}', 'λ_{min}', 'Location', 'best');
% grid on;
% 
% subplot(2,1,2);
% semilogy(t, eig_ratio, 'g-', 'LineWidth', 2);
% ylabel('λ_{max}/λ_{min}');
% xlabel('Time (s)');
% grid on;
% 
% % Plot 4: Correlation Analysis
% figure;
% subplot(2,1,1);
% loglog(cond_Mqq, abs(Ener), 'ro', 'MarkerSize', 3);
% title('Energy Error vs Matrix Conditioning');
% xlabel('cond(Mqq)');
% ylabel('|Energy Error|');
% grid on;
% 
% subplot(2,1,2);
% loglog(cond_Mqq, abs(H_norm), 'bo', 'MarkerSize', 3);
% xlabel('cond(Mqq)');
% ylabel('|Momentum Error|');
% grid on;
% 
% fprintf('Maximum condition number: %.2e\n', max(cond_Mqq));
% fprintf('Final condition number: %.2e\n', cond_Mqq(end));
% 
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
% % Plotting results for Cartesian rotation vector
% figure;
% subplot(3,1,1);
% plot(t,psi1,'b-','LineWidth',2); 
% ylim([-3.5,3.5]);
% title('Psi1');
% ylabel({'X rot vector'});
% grid on;
% 
% subplot(3,1,2);
% plot(t,psi2,'r-','LineWidth',2); 
% ylim([-3.5,3.5]);
% title('Psi2');
% ylabel({'Y rot vector'});
% grid on;
% 
% subplot(3,1,3);
% plot(t,psi3,'g-','LineWidth',2); 
% ylim([-3.5,3.5]);
% title('Psi3');
% ylabel({'Z rot vector'});
% grid on;
% 
% % Plotting results for theta
% figure;
% plot(t,theta,'b-','LineWidth',2); 
% ylim([-0.5,7.5]);
% title('Wrapping - Theta - norm of rotation vector');
% xlabel('Time (s)');
% ylabel('Theta(t)');
% grid on;
% 
% % Plotting results for Energy and Momentum
% figure;
% subplot(2,1,1);
% plot(t,Ener,'b-','LineWidth',1); 
% ylim([-5*10^-1,5*10^-1]);
% title('Total Energy');
% ylabel('Energy balance');
% grid on;
% 
% subplot(2,1,2);
% plot(t,H_norm,'r-','LineWidth',1); 
% ylim([-5*10^-1,5*10^-1]);
% title('Angular momentum(norm)');
% xlabel('Time (s)');
% ylabel('Ang moment');
% grid on;