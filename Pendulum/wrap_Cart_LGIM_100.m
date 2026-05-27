% --- timing: tic/toc now brackets only the custom_rk4 integration call ---
% Initial values for psi - rotational part
psi10=0.55536036727;
psi20=0.55536036727;
psi30=0;
t_end = 100;  % End value for time
if ~exist('dt', 'var') || isempty(dt)
    dt = 0.001;   % Time step (default 1 ms)
end

% Inertia tensor
Iqq = [3.0025 0 0; 0 3.0025 0; 0 0 0.005];

% Initial conditions
initial_conditions = [1.5; -1.5; 2.12132034356;psi10;psi20;psi30;0; 0; 0; 0;0;0];

% Time span
tspan = 0:dt:t_end;

% Skew symmetric matrix from vector
function S = skew(v)
    S = [0 -v(3) v(2); v(3) 0 -v(1); -v(2) v(1) 0];
end

% Skew symmetric matrix - spherical joint location
sj=[0;0;-3];
Sj_hat = skew(sj);

% Localâ€“global transition - Algorithm 1, function UAC from article
function q_new = LGT(q_old, Delta_q)
    % Split into translation and rotation parts
    Delta_x = Delta_q(1:3);
    Delta_psi = Delta_q(4:6);
    
    % Translation update (simple addition)
    x_new = q_old(1:3) + Delta_x;
    
    % Rotation update using composition operation â‹„
    psi_old = q_old(4:6);
    psi_new = rotation_vector_composition(psi_old, Delta_psi);
    
    q_new = [x_new; psi_new];
end

% Appendix A.3 from article
    function psi = rotation_vector_composition(psi1, psi2)
        % Implementation of Eq. (A10)
        phi1 = norm(psi1);
        phi2 = norm(psi2);
        
        if phi1 < 1e-16
            psi = psi2;
            return;
        end
        if phi2 < 1e-16
            psi = psi1;
            return;
        end
        
        % Quaternion intermediate representation
        q1 = [cos(phi1/2); (psi1/phi1)*sin(phi1/2)];
        q2 = [cos(phi2/2); (psi2/phi2)*sin(phi2/2)];
        
        % Quaternion multiplication
        q = [q1(1)*q2(1) - dot(q1(2:4), q2(2:4));
             q1(1)*q2(2:4) + q2(1)*q1(2:4) + cross(q1(2:4), q2(2:4))];
        
        % Convert back to rotation vector
        phi = 2*acos(q(1));
        if abs(phi) < 1e-8
            psi = zeros(3,1);
        else
            psi = (phi/sin(phi/2))*q(2:4);
        end
    end

    %Appendix A.2 Tangent operator

    function invT = invT_exp(DeltaPhi)
    
    DeltaPsi=DeltaPhi(4:6);
    phi = norm(DeltaPsi);

            if phi < 1e-3
                invTemp = eye(3) + 0.5*skew(DeltaPsi) + (1/12 + phi^2/720+phi^4/30240+phi^6/1209600)*(skew(DeltaPsi))^2;
            else
                invTemp = eye(3) + 0.5*skew(DeltaPsi) + ...
                    (1 - (phi/2)*cot(phi/2))/phi^2 * (skew(DeltaPsi))^2;
            end
invT = [eye(3) zeros(3);
    zeros(3) invTemp];
    end

% Define the ODE system as a function
function [dydt, second_derivatives, C, dC,ddC, Ener] = odesystem(t, F)
    x = F(1);
    y = F(2);
    z = F(3);
    psi1 = F(4);
    psi2 = F(5);
    psi3 = F(6);
    xd = F(7);
    yd = F(8);
    zd = F(9);
    wx = F(10);
    wy = F(11);
    wz = F(12);

Iqq = [3.0025 0 0; 0 3.0025 0; 0 0 0.005];

Gint=LieGmatrixCart(x, y, z, psi1, psi2, psi3);%Function for matrix G

%  Determining angle phi (theta in report) from psi
Phi = sqrt(psi1^2 + psi2^2 + psi3^2);
p = [psi1, psi2, psi3];
Ps = skew (p);
%  Transform matrix A and A2 (for small phi) from Rodriguez formula
if Phi > 10^-3 
A = eye(3) + Ps*(sin(Phi)/Phi) + 2*Ps*Ps*sin(Phi/2)*sin(Phi/2)/(Phi*Phi);
else
A = eye(3) + Ps*(1 - Phi^2/6 + Phi^4/120) + 2*Ps*Ps*(1/2 - Phi^2/48 + Phi^4/3840);
end

% Calculating psi from angular velocities
w_vec=[wx;
wy;
wz];

InvG= inv (Gint);
psi_ans=InvG*w_vec;

% Kinematic constraints and Jacobian matrices
Phi_rpi = zeros (3,6);
Phi_r= eye(3);
sj=[0;0;-3];
Sj_hat = skew(sj);

Phi_pi= -A*Sj_hat;
Phi_rpi = [Phi_r Phi_pi];

% Inertia tensor
Iqq=[3.0025 0 0; 0 3.0025 0;0 0 0.005];

% Augmented mass matrix
M = [eye(3) zeros(3);
        zeros(3) Iqq];

% Left side matrix from EOM
Mat = [M Phi_rpi';
       Phi_rpi zeros(3)];

% Right side EOM - forces
w_hat=skew(w_vec);
gamma = -A*w_hat*w_hat*sj;
Qv = -w_hat*Iqq*w_vec;
Forces=[0; 0; -9.8; Qv; gamma];

%Solving EOM
NewRes=Mat\Forces;

 dydt = zeros(12,1);
    dydt(1) = xd; % x' = xd
    dydt(2) = yd; % y' = yd
    dydt(3) = zd; % z' = zd
    dydt(4) = psi_ans(1) ; % psi1
    dydt(5) = psi_ans(2); % psi2 
    dydt(6) = psi_ans(3); % psi3
    dydt(7) = NewRes(1); % xd' (x'')
    dydt(8) = NewRes(2);  % yd' (y'')  
    dydt(9) = NewRes(3); % zd' (z'') 
    dydt(10) = NewRes(4); % wx'
    dydt(11) = NewRes(5); % wy' 
    dydt(12) = NewRes(6); % wz' 

% accelerations 
second_derivatives = NewRes(1:6);

%constraints violations
tempC=[x;y;z]+A*sj;
C=tempC(1:3); %position level
VelCon = Phi_r*[xd;yd;zd] + Phi_pi*w_vec;
dC=VelCon(1:3); %velocity level
AccCon = Phi_r*NewRes(1:3) + Phi_pi*NewRes(4:6)-gamma;
ddC=AccCon(1:3); %acceleration level

%energy balance
qd= [xd yd zd wx wy wz];
Ener = 0.5 * (qd * M* qd')+9.8*z-9.8*3.*(1-2*0.27059805007.^2-2*0.27059805007.^2);
end

% Custom Runge-Kutta 4th order method
function [t, F, second_derivatives,C, dC,ddC,Ener] = custom_rk4(odefun, tspan, y0)
    n = length(tspan);
    F_temp = zeros(length(y0), n);
    F_temp(:,1) = y0;

    %for accelerations, constraint violations and energy balance
    second_derivatives = zeros(6, n);
    C = zeros(3,n);
    dC = zeros(3,n);
    ddC = zeros(3,n);
    Ener=zeros(1,n);

% Algortihm 3 from article - Runge Kutta Munthe Kaas
 for i = 1:(n-1)
        h = tspan(i+1) - tspan(i);
        ti = tspan(i);
        yi = F_temp(:,i);

        % --- WRAP at the BEGINNING of each RK4 step ---
        psi = yi(4:6);
        Phi = norm(psi);
        if Phi > pi
%fprintf('Wrapping RK4: Phi = %.3f -> ', Phi)
            psi = -(2*pi - Phi) * (psi / Phi);
            yi(4:6) = psi;  % Update yi with wrapped psi
        end


        % Split into q (positions/orientations) and v (velocities)
        qi = yi(1:6);  % First 6 elements: [x,y,z,psi1,psi2,psi3]
        vi = yi(7:12); % Next 6 elements: [xd,yd,zd,wx,wy,wz]
       
% Stage 1
        [k1, sd1,C1, dC1,ddC1,Ener1] = odefun(ti, yi);
        k1_q = k1(1:6);
        k1_q(4:6) = yi(10:12);
        k1_v = k1(7:12);
     
% Stage 2
        Delta_q2 = k1_q * h/2;
        q2 = LGT(qi, Delta_q2);
        v2 = vi + (h/2)*k1_v;
        y2 = [q2; v2];

        [k2, ~] = odefun(ti + h/2, y2);
        k2_q = invT_exp(Delta_q2)*v2 ;
        k2_v = k2(7:12);

% Stage 3
        Delta_q3 = k2_q * h/2;
        q3 = LGT(qi, Delta_q3);
        v3 = vi + (h/2)*k2_v;
        y3 = [q3; v3];
  
        [k3,  ~] = odefun(ti + h/2, y3);
        k3_q = invT_exp(Delta_q3)*v3 ;
        k3_v = k3(7:12);

% Stage 4        
        Delta_q4 = k3_q * h;
        q4 = LGT(qi, Delta_q4);
        v4 = vi + h*k3_v;
        y4 = [q4; v4];

        [k4, ~] = odefun(ti + h, y4);
        k4_q = invT_exp(Delta_q4)*v4 ;
        k4_v = k4(7:12);

%Final delta q, q and v
Delta_q = (h/6)*(k1_q +2*k2_q +2*k3_q +k4_q);
q_new = LGT(qi, Delta_q);
v_new = vi + (h/6)*(k1_v + 2*k2_v + 2*k3_v + k4_v);
        
        % Store full state (12 elements)
        F_temp(:,i+1) = [q_new; v_new];

%for accelerations, constraint violations and energy balance        
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
% --- timing: integration is run n_timing_repeats times; executionTime is the median ---
% (only the integration is repeated; setup, post-processing and save run once)
if ~exist('n_timing_repeats', 'var') || isempty(n_timing_repeats)
    n_timing_repeats = 1;   % a batch runner may raise this for a timing study
end
timing_samples = zeros(1, n_timing_repeats);
for timing_rep = 1:n_timing_repeats
    tic;
    [t, F, second_derivatives, C, dC, ddC,Ener] = custom_rk4(@odesystem, tspan, initial_conditions);
    timing_samples(timing_rep) = toc;  % integration only
end
executionTime = median(timing_samples);  % robust to run-to-run noise

% Extract results from Y matrix
x = F(:,1);
y = F(:,2);
z = F(:,3);
psi1 = F(:,4);
psi2 = F(:,5);
psi3 = F(:,6);
xd = F(:,7);
yd = F(:,8);
zd = F(:,9);
wx= F(:,10);
wy = F(:,11);
wz = F(:,12);

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

% x(end)
% xd(end)
% x_double_prime(end)

theta = sqrt(psi1.^2 + psi2.^2 + psi3.^2);

% executionTime captured at the custom_rk4 call above (timing harness)

if ~exist('save_filename', 'var') || isempty(save_filename)
    save_filename = sprintf('wrLGIM_dt_%.1fms.mat', dt * 1000);
end
save(save_filename, '-v7.3');  

% % Plotting results for x
% figure;
% subplot(3,1,1);
% plot(t,x,'b-','LineWidth',2);
% ylim([-6,12]);
% title('Rx(t)');
% ylabel('Displacement X, m');
% grid on;
% 
% subplot(3,1,2);
% plot(t,xd,'r-','LineWidth',2); 
% ylim([-15,15]);
% title('Rx''(t)');
% ylabel({'Velocity X'; '(m/s)'});
% grid on;
% 
% subplot(3,1,3);
% plot(t,x_double_prime,'g-','LineWidth',2); 
% ylim([-25,25]);
% title('Rx''''(t)');
% xlabel('Time (s)');
% ylabel({'Acceleration X'; '(m/s^2)'});
% grid on;
% 
% % Plotting results for y
% figure;
% subplot(3,1,1);
% plot(t,y,'b-','LineWidth',2); 
% ylim([-5,20]);
% title('Ry(t)');
% ylabel('Displacement Y, m');
% grid on;
% 
% subplot(3,1,2);
% plot(t,yd,'r-','LineWidth',2); 
% ylim([-10,20]);
% title('Ry''(t)');
% ylabel({'Velocity Y'; '(m/s)'});
% grid on;
% 
% subplot(3,1,3);
% plot(t,y_double_prime,'g-','LineWidth',2); 
% ylim([-20,40]);
% title('Ry''''(t)');
% xlabel('Time (s)');
% ylabel({'Acceleration Y'; '(m/s^2)'});
% grid on;
% 
% % Plotting results for z
% figure;
% subplot(3,1,1);
% plot(t,z,'b-','LineWidth',2); 
% ylim([-4,15]);
% title('Rz(t)');
% ylabel('Displacement Z, m');
% grid on;
% 
% subplot(3,1,2);
% plot(t,zd,'r-','LineWidth',2); 
% ylim([-10,20]);
% title('Rz''(t)');
% ylabel({'Velocity Z'; '(m/s)'});
% grid on;
% 
% subplot(3,1,3);
% plot(t,z_double_prime,'g-','LineWidth',2); 
% ylim([-15,40]);
% title('Rz''''(t)');
% xlabel('Time (s)');
% ylabel({'Acceleration Z'; '(m/s^2)'});
% grid on;
% 
% % Plotting results for Angular velocities
% figure;
% subplot(3,1,1);
% plot(t,wx,'b-','LineWidth',2); 
% ylim([-6,12]);
% title('Wx');
% ylabel({'Ang velocity Wx'; 's-1'});
% grid on;
% 
% subplot(3,1,2);
% plot(t,wy,'r-','LineWidth',2); 
% ylim([-6,12]);
% title('Wy');
% ylabel({'Ang velocity Wy'; 's-1'});
% grid on;
% 
% subplot(3,1,3);
% plot(t,wz,'g-','LineWidth',2); 
% ylim([-2,2]);
% title('Wz');
% ylabel({'Ang velocity Wz';'s-1'});
% grid on;
% 
% % Plotting results for Angular accelerations
% figure;
% subplot(3,1,1);
% plot(t,AngAccX,'b-','LineWidth',2); 
% ylim([-3,3]);
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
% ylim([-1,1]);
% title('Wz''');
% ylabel({'Ang acceleration Wz';' s-2'});
% grid on;
% 
% % Plotting results for theta
% figure;
% plot(t,theta,'b-','LineWidth',2); 
% ylim([0,3.5]);
% title('Angle Theta');
% xlabel('Time (s)');
% ylabel('Norm of rotation vector Theta(t)');
% grid on;
% 
% % Plotting results for psi
% figure;
% subplot(2,1,1);
% plot(t,psi1,'b-','LineWidth',2); 
% hold on;
% plot(t,psi2,'r--','LineWidth',2); 
% ylim([-3.5,3.5]);
% title('psi1(t), psi2(t)');
% ylabel('psi1(t), psi2(t)');
% grid on;
% lgd = legend('psi1', 'psi2', 'Location', 'northeast', 'FontSize', 8, 'TextColor', 'blue');
% 
% subplot(2,1,2);
% plot(t,psi3,'b-','LineWidth',2); 
% ylim([-1,1]);
% title('psi3(t)');
% ylabel('Psi3');
% grid on;
% 
% % Plotting results for Constraint violation and Energy balance
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
