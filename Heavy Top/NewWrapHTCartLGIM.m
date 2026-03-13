tic;
% Initial values for psi - rotational part
psi10= 0;
psi20= 0;
psi30= 0;
t_end = 20;  % End value for time
% Time step (seconds)
if ~exist('dt','var') || isempty(dt)
    dt = 0.00001;
end
% Time span
tspan = 0:dt:t_end;

x0=0;
y0=1;
z0=0;

wx0=0;
wy0=150;
wz0= -4.61538;
%wx0=0;
%wy0=0;
%wz0= 0;

xd0=4.615380000000000;
%xd0=0;
yd0=0;
zd0=0;

m0=15;
l0=1;

% Inertia tensor
Iqq = diag([0.234375; 0.46875; 0.234375]);

% Initial conditions
initial_conditions =[x0; y0;z0;psi10;psi20;psi30; xd0;yd0;zd0; wx0; wy0; wz0;];

function S = skew(v)
    S = [0 -v(3) v(2); v(3) 0 -v(1); -v(2) v(1) 0];
end

sj=[0;-l0;0];
Sj_hat = skew(sj);
A0 = eye(3);
C0 = [x0;y0;z0] + A0*sj;
%fprintf('Initial constraint violation: %.15e\n', norm(C0));
gw0 = A0*[wx0 wy0 wz0]';
xd0_calc = gw0(2)*z0 - gw0(3)*y0;  % Should be 4.615380000000000
%fprintf('xd0 calculated: %.15f\n', xd0_calc);

Ener0 =5.435696790865547e+03;

% Local–global transition - Algorithm 1, function UAC from article
function q_new = NewLGT(q_old, Delta_q)
    Delta_psi = Delta_q;
  
    % Rotation update using composition operation ⋄
    psi_old = q_old;
    psi_new = rotation_vector_composition(psi_old, Delta_psi);
    
    q_new = psi_new;
end

% Appendix A.3 from article
    function psi = rotation_vector_composition(psi1, psi2)
        % Implementation of Eq. (A10)
        phi1 = norm(psi1);
        phi2 = norm(psi2);
        
        if phi1 < 1e-16
            psi = psi2;
            %fprintf('Conv: phi1 = %.3f -> ', phi1)
            return;
        end
        if phi2 < 1e-16
            psi = psi1;
            %fprintf('Conv: phi2 = %.3f -> ', phi2)
            return;
        end
        
        % Quaternion intermediate representation
        q1 = [cos(phi1/2); (psi1/phi1)*sin(phi1/2)];
        q2 = [cos(phi2/2); (psi2/phi2)*sin(phi2/2)];
        
        % Quaternion multiplication
       % q = [q1(1)*q2(1) - dot(q1(2:4), q2(2:4));
        %     q1(1)*q2(2:4) + q2(1)*q1(2:4) + cross(q2(2:4), q1(2:4))];

        q = [q1(1)*q2(1) - q1(2:4)'*q2(2:4);
             q1(1)*q2(2:4) + q2(1)*q1(2:4) + cross(q1(2:4), q2(2:4))];
        
        % Convert back to rotation vector
        phi = 2*acos(q(1));
        if abs(phi) < 1e-8
            psi = q(2:4)/2;
            %fprintf('Conv: phi = %.3f -> ', phi)
        else
            psi = (phi/sin(phi/2))*q(2:4);
        end
    end

    %Appendix A.2 Tangent operator

    function invT = NewinvT_exp(DeltaPhi)
    phi = norm(DeltaPhi);
            if phi < 1e-3
                invTemp = eye(3) + 0.5*skew(DeltaPhi) + (1/12 + phi^2/720+phi^4/30240+phi^6/1209600)*(skew(DeltaPhi))^2;
            else
                invTemp = eye(3) + 0.5*skew(DeltaPhi) + ...
                    (1 - (phi/2)*cot(phi/2))/phi^2 * (skew(DeltaPhi))^2;
            end
invT = invTemp;
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
    dydt(2) = yd; % y' = yd
    dydt(3) = zd; % z' = zd
    dydt(4) = psi_ans(1) ; % psi1
    dydt(5) = psi_ans(2); % psi2 
    dydt(6) = psi_ans(3); % psi3
    %dydt(4) = w_vec(1) ; % psi1
    %dydt(5) = w_vec(2); % psi2 
    %dydt(6) = w_vec(3); % psi3    
    dydt(7) = NewRes(1); % xd' (x'')
    dydt(8) = NewRes(2);  % yd' (y'')  
    dydt(9) = NewRes(3); % zd' (z'') 
    dydt(10) = NewRes(4); % wx'
    dydt(11) = NewRes(5); % wy' 
    dydt(12) = NewRes(6); % wz' 


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
fprintf('Wrapping RK4: Phi = %.3f -> ', Phi)
            psi = -(2*pi - Phi) * (psi / Phi);
            yi(4:6) = psi;  % Update yi with wrapped psi
        end

        % Split into q (positions/orientations) and v (velocities)
        tri= yi(1:3);  % Translational parameters: [x,y,z]
        qi = yi(4:6);  % Rotational parameters: [psi1,psi2,psi3]
        vi = yi(7:12); % Velocities: [xd,yd,zd,wx,wy,wz]
       
% Stage 1
        [k1, sd1,C1, dC1,ddC1,Ener1] = odefun(ti, yi);
        k1_tr= h*k1(1:3);
        k1_q = h*yi(10:12);
        %k1_q = h*k1(4:6);
        k1_v = h*k1(7:12);
     
% Stage 2
tr2 = tri + (1/2)*k1_tr;
        Delta_q2 = k1_q * 1/2;
        q2 = rotation_vector_composition(qi, Delta_q2);
        %q2 = qi + Delta_q2;
        v2 = vi + (1/2)*k1_v;
        y2 = [tr2;q2; v2];

        [k2, ~] = odefun(ti + h/2, y2);
        k2_tr = h*k2(1:3);
        
        %k2_q = h*k2(4:6) ;
        k2_q = h*invdexp(v2(4:6),Delta_q2);
        %k2_q = h*NewinvT_exp(Delta_q2)*v2(4:6);
        %k2_q = h*NewinvT_exp(Delta_q2)*k2(4:6) ;
        k2_v = h*k2(7:12);

% Stage 3
tr3 = tri + (1/2)*k2_tr;
        Delta_q3 = k2_q * 1/2;
        q3 = rotation_vector_composition(qi, Delta_q3);
%q3 = qi + Delta_q3;
        v3 = vi + (1/2)*k2_v;
        y3 = [tr3;q3; v3];
  
        [k3,  ~] = odefun(ti + h/2, y3);
        k3_tr = h*k3(1:3);
        %k3_q = h*k3(4:6) ;
        k3_q = h*invdexp(v3(4:6),Delta_q3);
        %k3_q = h*NewinvT_exp(Delta_q3)*v3(4:6) ;
        %k3_q = h*NewinvT_exp(Delta_q3)*k3(4:6) ;
        k3_v = h*k3(7:12);

% Stage 4        
tr4 = tri + k3_tr;
        Delta_q4 = k3_q;
        q4 = rotation_vector_composition(qi, Delta_q4);
        %q4 = qi + Delta_q4;
        v4 = vi + k3_v;
        y4 = [tr4;q4; v4];

        [k4, ~] = odefun(ti + h, y4);
        k4_tr = h*k4(1:3);
        %k4_q = h*k4(4:6) ;
        k4_q = h*invdexp(v4(4:6),Delta_q4);
        %k4_q = h*NewinvT_exp(Delta_q4)*v4(4:6) ;
        %k4_q = h*NewinvT_exp(Delta_q4)*k4(4:6) ;
        k4_v = h*k4(7:12);

%Final delta q, q and v
tr_new = tri + (1/6)*(k1_tr + 2*k2_tr + 2*k3_tr + k4_tr);
Delta_q = (1/6)*(k1_q +2*k2_q +2*k3_q +k4_q);
q_new = rotation_vector_composition(qi, Delta_q);
%q_new = qi + Delta_q;
v_new = vi + (1/6)*(k1_v + 2*k2_v + 2*k3_v + k4_v);
        
        % Store full state (12 elements)
        F_temp(:,i+1) = [tr_new;q_new; v_new];

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
[t, F, second_derivatives, C, dC, ddC,Ener] = custom_rk4(@odesystem, tspan, initial_conditions);

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
    save_filename = sprintf('HT_CartLGIM_dt_%.2fms.mat', dt*1000);
end
save(save_filename, '-v7.3');

% Plotting results for x y z (commented out)
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
%legend ('Rx','Ry','Rz');
% grid on;

% Plotting results for x
% figure;
% subplot(3,1,1);
% plot(t,x,'b-','LineWidth',2);
% ylim([-1,1]);
% title('Rx(t)');
% ylabel('Displacement X, m');
% grid on;

% subplot(3,1,2);
% plot(t,xd,'r-','LineWidth',2); 
% ylim([-10,10]);
% title('Rx''(t)');
% ylabel({'Velocity X'; '(m/s)'});
% grid on;

% subplot(3,1,3);
% plot(t,x_double_prime,'g-','LineWidth',2); 
% ylim([-50,45]);
% title('Rx''''(t)');
% xlabel('Time (s)');
% ylabel({'Acceleration X'; '(m/s^2)'});
% grid on;

% Plotting results for y
% figure;
% subplot(3,1,1);
% plot(t,y,'b-','LineWidth',2); 
% ylim([-1,1]);
% title('Ry(t)');
% ylabel('Displacement Y, m');
% grid on;

% subplot(3,1,2);
% plot(t,yd,'r-','LineWidth',2); 
% ylim([-10,10]);
% title('Ry''(t)');
% ylabel({'Velocity Y'; '(m/s)'});
% grid on;

% subplot(3,1,3);
% plot(t,y_double_prime,'g-','LineWidth',2); 
% ylim([-50,50]);
% title('Ry''''(t)');
% xlabel('Time (s)');
% ylabel({'Acceleration Y'; '(m/s^2)'});
% grid on;

% Plotting results for z
% figure;
% subplot(3,1,1);
% plot(t,z,'b-','LineWidth',2); 
% ylim([-1.5,0.5]);
% title('Rz(t)');
% ylabel('Displacement Z, m');
% grid on;

% subplot(3,1,2);
% plot(t,zd,'r-','LineWidth',2); 
% ylim([-10,10]);
% title('Rz''(t)');
% ylabel({'Velocity Z'; '(m/s)'});
% grid on;

% subplot(3,1,3);
% plot(t,z_double_prime,'g-','LineWidth',2); 
% ylim([-40,50]);
% title('Rz''''(t)');
% xlabel('Time (s)');
% ylabel({'Acceleration Z'; '(m/s^2)'});
% grid on;

% Plotting results for theta
% figure;
% plot(t,theta,'b-','LineWidth',2); 
%xlim([0,0.2]);
% ylim([0,4]);
% title('Angle Theta');
% xlabel('Time (s)');
% ylabel('Norm of rotation vector Theta(t)');
% grid on;

% Plotting results for psi
% figure;
% subplot(2,1,1);
% plot(t,psi1,'b-','LineWidth',2); 
% hold on;
% plot(t,psi2,'r--','LineWidth',2); 
%xlim([0,0.2]);
% ylim([-4,4]);
% title('psi1(t), psi2(t)');
% ylabel('psi1(t), psi2(t)');
% grid on;
% lgd = legend('psi1', 'psi2', 'Location', 'northeast', 'FontSize', 8, 'TextColor', 'blue');

% subplot(2,1,2);
% plot(t,psi3,'b-','LineWidth',2); 
% ylim([-4,4]);
% title('psi3(t)');
% ylabel('Psi3');
% grid on;

% Plotting results for Angular velocities
% figure;
% subplot(3,1,1);
% plot(t,wx,'b-','LineWidth',2); 
% ylim([-7,7]);
% title('Wx');
% ylabel({'Ang velocity Wx'; 's-1'});
% grid on;

% subplot(3,1,2);
% plot(t,wy,'r-','LineWidth',2); 
% ylim([0,200]);
% title('Wy');
% ylabel({'Ang velocity Wy'; 's-1'});
% grid on;

% subplot(3,1,3);
% plot(t,wz,'g-','LineWidth',2); 
% ylim([-7,7]);
% title('Wz');
% ylabel({'Ang velocity Wz';'s-1'});
% grid on;

%Plotting results for Angular accelerations
% figure;
% subplot(3,1,1);
% plot(t,AngAccX,'b-','LineWidth',2); 
% ylim([-1000,1000]);
% title('Wx''');
% ylabel({'Ang acceleration Wx';' s-2'});
% grid on;

% subplot(3,1,2);
% plot(t,AngAccY,'r-','LineWidth',2); 
% ylim([-3,3]);
% title('Wy''');
% ylabel({'Ang acceleration Wy';' s-2'});
% grid on;

% subplot(3,1,3);
% plot(t,AngAccZ,'g-','LineWidth',2); 
% ylim([-1000,1000]);
% title('Wz''');
% ylabel({'Ang acceleration Wz';' s-2'});
% grid on;


% Plotting results for Constraint violation and Energy balance
% figure;
% subplot(4,1,1);
% plot(t,grC,'b-','LineWidth',0.5); 
% ylim([-1*10^-3,1*10^-3]);
% title('∣∣C∣∣');
% ylabel('Violation');
% grid on;

% subplot(4,1,2);
% plot(t,grdC,'g-','LineWidth',1); 
% ylim([-1*10^-4,1*10^-4]);
% title('∣∣C''∣∣');
% ylabel('Violation');
% grid on;

% subplot(4,1,3);
% plot(t,grddC,'r-','LineWidth',1); 
% ylim([-5*10^-12,5*10^-12]);
% title('∣∣C''''∣∣');
% ylabel('Violation');
% grid on;

% subplot(4,1,4);
% plot(t,Ener,'b-','LineWidth',1); 
% ylim([-1*10^-0,1*10^-0]);
% title('Total Energy');
% xlabel('Time (s)');
% ylabel('Energy balance');
% grid on;

% Plotting results for Constraint violation and Energy balance
% Version with at least 3 ticks and sparser grid
% figure;
% subplot(4,1,1);
% % plot(t,abs(grC),'b-','LineWidth',0.5);
% set(gca, 'YScale', 'log');
% set(gca, 'YMinorTick', 'on');
% Set at least 3 ticks, including the highest error level
% ylims = ylim;
% ymax = max(abs(grC));
% ymin = min(abs(grC(abs(grC) > 0))); % minimum non-zero value
% Create at least 3 ticks, with one close to the maximum
% tick_range = floor(log10(ymin)):ceil(log10(ymax));
% if length(tick_range) < 3
%     % Extend range to have at least 3 ticks
%     tick_range = (ceil(log10(ymax))-2):ceil(log10(ymax));
% end
% ytick_values = 10.^tick_range;
% set(gca, 'YTick', ytick_values);
% set(gca, 'YTickLabelMode', 'auto');
% title('||C|| (Log scale)');
% ylabel('Violation');
% grid on;
% set(gca, 'GridAlpha', 0.3); % Make grid less prominent
% set(gca, 'MinorGridAlpha', 0.1); % Make minor grid very subtle

% subplot(4,1,2);
% % plot(t,abs(grdC),'g-','LineWidth',1);
% set(gca, 'YScale', 'log');
% set(gca, 'YMinorTick', 'on');
% ylims = ylim;
% ymax = max(abs(grdC));
% ymin = min(abs(grdC(abs(grdC) > 0)));
% tick_range = floor(log10(ymin)):ceil(log10(ymax));
% if length(tick_range) < 3
%     tick_range = (ceil(log10(ymax))-2):ceil(log10(ymax));
% end
% ytick_values = 10.^tick_range;
% set(gca, 'YTick', ytick_values);
% set(gca, 'YTickLabelMode', 'auto');
% title('||C'' || (Log scale)');
% ylabel('Violation');
% grid on;
% set(gca, 'GridAlpha', 0.3);
% set(gca, 'MinorGridAlpha', 0.1);

% subplot(4,1,3);
% % plot(t,abs(grddC),'r-','LineWidth',1);
% set(gca, 'YScale', 'log');
% set(gca, 'YMinorTick', 'on');
% ylims = ylim;
% ymax = max(abs(grddC));
% ymin = min(abs(grddC(abs(grddC) > 0)));
% tick_range = floor(log10(ymin)):ceil(log10(ymax));
% if length(tick_range) < 3
%     tick_range = (ceil(log10(ymax))-2):ceil(log10(ymax));
% end
% ytick_values = 10.^tick_range;
% set(gca, 'YTick', ytick_values);
% set(gca, 'YTickLabelMode', 'auto');
% title('||C'''' || (Log scale)');
% ylabel('Violation');
% grid on;
% set(gca, 'GridAlpha', 0.3);
% set(gca, 'MinorGridAlpha', 0.1);

% subplot(4,1,4);
% % plot(t,abs(Ener),'b-','LineWidth',1);
% set(gca, 'YScale', 'log');
% set(gca, 'YMinorTick', 'on');
% ylims = ylim;
% ymax = max(abs(Ener));
% ymin = min(abs(Ener(abs(Ener) > 0)));
% tick_range = floor(log10(ymin)):ceil(log10(ymax));
% if length(tick_range) < 3
%     tick_range = (ceil(log10(ymax))-2):ceil(log10(ymax));
% end
% ytick_values = 10.^tick_range;
% set(gca, 'YTick', ytick_values);
% set(gca, 'YTickLabelMode', 'auto');
% title('Total Energy (Log scale)');
% xlabel('Time (s)');
% ylabel('Energy balance');
% grid on;
% set(gca, 'GridAlpha', 0.3);
% set(gca, 'MinorGridAlpha', 0.1);
