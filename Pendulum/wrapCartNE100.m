tic;
% Initial values for psi - rotational part
psi10=0.55536036727;
psi20=0.55536036727;
psi30=0;
t_end = 100;  % End value for time
if ~exist('dt', 'var') || isempty(dt)
    dt = 0.1;   % Time step (default 100 ms)
end

% Inertia tensor
Iqq = [3.0025 0 0; 0 3.0025 0; 0 0 0.005];

% Initial conditions
initial_conditions = [1.5;0; -1.5;0; 2.12132034356;0; 0; 0; 0;psi10;psi20;psi30];

% Time span
tspan = 0:dt:t_end;

function S = skew(v)
    S = [0 -v(3) v(2); v(3) 0 -v(1); -v(2) v(1) 0];
end
sj=[0;0;-3];
Sj_hat = skew(sj);

% Define the ODE system as a function
function [dydt, second_derivatives, C, dC,ddC, Ener] = odesystem(t, F)
    x = F(1);
    xd = F(2);
    y = F(3);
    yd = F(4);
    z = F(5);
    zd = F(6);
    wx = F(7);
    wy = F(8);
    wz = F(9);
    psi1 = F(10);
    psi2 = F(11);
    psi3 = F(12);

Iqq = [3.0025 0 0; 0 3.0025 0; 0 0 0.005];

%  Determining angle phi (theta in report) from psi
Phi = sqrt(psi1^2 + psi2^2 + psi3^2);
p = [psi1, psi2, psi3];
Ps = skew (p);
%  Transform matrix A and A2 (for small phi) from Rodriguez formula
if Phi > 10^-2 
A = eye(3) + Ps*(sin(Phi)/Phi) + 2*Ps*Ps*sin(Phi/2)*sin(Phi/2)/(Phi*Phi);
else
A = eye(3) + Ps*(1 - Phi^2/6 + Phi^4/120) + 2*Ps*Ps*(1/2 - Phi^2/48 + Phi^4/3840);
end

%  Matrix G from Rodriguez formula 
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

psi_ans=Emap2*w_vec;

    Phi_rpi = zeros (3,6);
Phi_r= eye(3);
sj=[0;0;-3];
Sj_hat = skew(sj);

Phi_pi= -A*Sj_hat;
Phi_rpi = [Phi_r Phi_pi];


    Iqq=[3.0025 0 0; 0 3.0025 0;0 0 0.005];
    M = [eye(3) zeros(3);
        zeros(3) Iqq];


Mat = [M Phi_rpi';
       Phi_rpi zeros(3)];

w_hat=skew(w_vec);

gamma = -A*w_hat*w_hat*sj;

Qv = -w_hat*Iqq*w_vec;

Forces=[0; 0; -9.8; Qv; gamma];
NewRes=Mat\Forces;

 dydt = zeros(12,1);
    dydt(1) = xd; % x' = xd
    dydt(2) = NewRes(1); % xd' (x'')
    dydt(3) = yd; % y' = yd
    dydt(4) = NewRes(2);  % yd' (y'')  
    dydt(5) = zd; % z' = zd
    dydt(6) = NewRes(3); % zd' (z'')  
    dydt(7) = NewRes(4); % wx'
    dydt(8) = NewRes(5); % wy' 
    dydt(9) = NewRes(6); % wz' 
    dydt(10) = psi_ans(1) ; 
    dydt(11) = psi_ans(2); 
    dydt(12) = psi_ans(3); 


second_derivatives = NewRes(1:6);

tempC=[x;y;z]+A*sj;
C=tempC(1:3);
VelCon = Phi_r*[xd;yd;zd] + Phi_pi*w_vec;
dC=VelCon(1:3);
AccCon = Phi_r*NewRes(1:3) + Phi_pi*NewRes(4:6)-gamma;
ddC=AccCon(1:3);

qd= [xd yd zd wx wy wz];
Ener = 0.5 * (qd * M* qd')+9.8*z-9.8*3.*(1-2*0.27059805007.^2-2*0.27059805007.^2);
end

% Custom Runge-Kutta 4th order method
function [t, F, second_derivatives,C, dC,ddC,Ener] = custom_rk4(odefun, tspan, y0)
    n = length(tspan);
    F_temp = zeros(length(y0), n);
    F_temp(:,1) = y0;

    second_derivatives = zeros(6, n);
    C = zeros(3,n);
    dC = zeros(3,n);
    ddC = zeros(3,n);
    Ener=zeros(1,n);

    for i = 1:(n-1)
        h = tspan(i+1) - tspan(i);
        ti = tspan(i);
        yi = F_temp(:,i);

       % --- WRAP at the BEGINNING of each RK4 step ---
        psi = yi(10:12);
        Phi = norm(psi);
        if Phi > pi
fprintf('Wrapping RK4: Phi = %.3f -> ', Phi)
            psi = -(2*pi - Phi) * (psi / Phi);
            yi(10:12) = psi;  % Update yi with wrapped psi
        end

        [k1, sd1,C1, dC1,ddC1,Ener1] = odefun(ti, yi);
        [k2, ~] = odefun(ti + h/2, yi + h*k1/2);
        [k3, ~] = odefun(ti + h/2, yi + h*k2/2);
        [k4, ~] = odefun(ti + h, yi + h*k3);
        
        F_temp(:,i+1) = yi + (h/6)*(k1 + 2*k2 + 2*k3 + k4);
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
xd = F(:,2);
y = F(:,3);
yd = F(:,4);
z = F(:,5);
zd = F(:,6);
wx = F(:,7);
wy = F(:,8);
wz = F(:,9);
psi1 = F(:,10);
psi2 = F(:,11);
psi3 = F(:,12);

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

% function pendulum_slider_xyz(t, x, y, z)
%     % pendulum_slider_xyz  Animate pendulum using direct x,y,z coordinates
%     %
%     % Inputs:
%     %   t – time vector (N×1)
%     %   x, y, z – position vectors (N×1) of the mass point
% 
%     N = numel(t);
%     L = 6;  % rod length for setting axis limits
% 
%     % Figure and static world axes
%     fig = figure('Name','Pendulum from XYZ','NumberTitle','off');
%     ax = axes(fig);
%     hold(ax,'on'), grid(ax,'on'), view(3), axis(ax,'equal');
%     xlabel('X'), ylabel('Y'), zlabel('Z');
%     title('Large rotations of thin rod');
%     lim = L * 1.1;
%     xlim([-L, L]), ylim([-L, L]), zlim([-lim, 0.7*L]);
%     ax.XLimMode='manual'; ax.YLimMode='manual'; ax.ZLimMode='manual';
%     ax.CameraViewAngleMode='manual';
%     %--- After setting azimuth/elevation or view ---%
%     view(ax, [-45, 30]);         % initial view (az,el)
%     ax.Projection = 'perspective';
% 
%     % Then rotate 135° about Z
%     camorbit(ax, 155, 0, 'data', [0 0 1]);
%     ax.CameraViewAngleMode = 'manual';
%         % Set camera so X axis points toward us, Y to the right, Z up
%    % ax.CameraPosition = [1 0 0];      
%   %  ax.CameraTarget   = [0 0 0];
%   %  ax.CameraUpVector = [0 0 1];
%   %  ax.Projection     = 'perspective';
%     % Draw world axes at origin
%     quiver3(0,0,0,1,0,0,'r','LineWidth',2,'MaxHeadSize',0.5);
%     quiver3(0,0,0,0,1,0,'g','LineWidth',2,'MaxHeadSize',0.5);
%     quiver3(0,0,0,0,0,1,'b','LineWidth',2,'MaxHeadSize',0.5);
% 
%     % Plane of rotation: original YZ (x=0) rotated about Z by 45°
% thetaPlane = deg2rad(-45);
% % Normal of original YZ is [1,0,0]; rotate this normal about Z
% n = [1, 0, 0] * [cos(thetaPlane), -sin(thetaPlane), 0;
%                  sin(thetaPlane),  cos(thetaPlane), 0;
%                                   0,               0, 1];
% % Create a grid in Y–Z coordinates spanning the rod swing
% [YY, ZZ] = meshgrid(linspace(-L, L, 50), linspace(-L, L, 50));
% % For each (Y,Z), find X so that dot([X,Y,Z],n)==0 => X = -(n2*Y + n3*Z)/n1
% XX = -(n(2)*YY + n(3)*ZZ) / n(1);
% % Plot the plane
% surf(XX, YY, ZZ, 'FaceColor', [0.8 0.8 0.8], ...
%      'FaceAlpha', 0.2, 'EdgeColor', 'none');
% 
%     % Draw global axes from origin with ticks
%     axisLength = L * 1.0;
%     plot3([0,axisLength],[0,0],[0,0],'r-','LineWidth',2);  % X-axis
%     plot3([0,0],[0,axisLength],[0,0],'g-','LineWidth',2);  % Y-axis
%     plot3([0,0],[0,0],[0,axisLength],'b-','LineWidth',2);  % Z-axis
%     % Negative dotted extensions
% plot3([-axisLength,0],[0,0],[0,0],'r--','LineWidth',1);
% plot3([0,0],[-axisLength,0],[0,0],'g--','LineWidth',1);
% plot3([0,0],[0,0],[-axisLength,0],'b--','LineWidth',1);
%     % Add axis ticks
%     for k = 1:axisLength
%         plot3([k,k],[0,0],[0,0], 'r.');  % X tick
%         plot3([0,0],[k,k],[0,0], 'g.');  % Y tick
%         plot3([0,0],[0,0],[k,k], 'b.');  % Z tick
%     end
% 
% 
%     % Pendulum line (initially to first point)
%     hLine = plot3(ax, [0, 2*x(1)], [0, 2*y(1)], [0, 2*z(1)], ...
%                   'k-','LineWidth',3);
% 
% 
% 
%     % Pivot marker
%     plot3(0,0,0,'ko','MarkerFaceColor','k','MarkerSize',6);
% 
%     % Time text
%     hText = text(0,0,lim*0.95,'','FontSize',12,'FontWeight','bold', ...
%                  'HorizontalAlignment','center','Parent',ax);
% 
%     % Before slider and updateFrame, create the COM marker handle:
%     hCOM = plot3(ax, x(1), y(1), z(1), 'ko', ...
%                  'MarkerFaceColor', 'c', 'MarkerSize', 6);
% 
%     % Slider UI
%     hSlider = uicontrol('Style','slider','Min',1,'Max',N,'Value',1, ...
%         'SliderStep',[1/(N-1),10/(N-1)],'Units','normalized', ...
%         'Position',[0.15 0.02 0.7 0.04],'Callback',@updateFrame);
%     addlistener(hSlider,'Value','PostSet',@(~,~) updateFrame());
% 
%     % Initial frame
%     updateFrame();
% 
%     function updateFrame(~,~)
%         idx = round(hSlider.Value);
%         % Tip position = 2 * COM
%         xt = 2*x(idx);
%         yt = 2*y(idx);
%         zt = 2*z(idx);
%         set(hLine,'XData',[0, xt], ...
%                   'YData',[0, yt], ...
%                   'ZData',[0, zt]);
%         % Move COM marker
%         set(hCOM, 'XData', x(idx), 'YData', y(idx), 'ZData', z(idx));
%         hText.String = sprintf('t = %.2f s', t(idx));
%         drawnow;
%     end
% end
% 
% function exportPendulumGif(t, x, y, z, gifFilename, skip)
%     % exportPendulumGif Save pendulum animation as GIF matching slider view
%     %
%     % t, x, y, z       – vectors (N×1)
%     % gifFilename      – string for output, e.g. 'pendulum.gif'
%     % skip             – integer frame subsample interval
% 
%     N = numel(t);
%     L = 6;
%     fig = figure('Visible','off');
%     ax = axes('Parent',fig);
%     hold(ax,'on');
% 
%     % Camera and view (same as slider)
%     %view(ax,[-45,30]);
%     view(ax, [120, 20]);
%     ax.Projection = 'perspective';
%     %camorbit(ax,155,0,'data',[0 0 1]);
%    % ax.CameraViewAngleMode = 'manual';
% 
%     % Axis limits and box
%     lim = L * 1.1;
%     xlim(ax,[-lim,lim]); ylim(ax,[-lim,lim]); zlim(ax,[-lim,lim]);
%     axis(ax,'equal');
%     grid(ax,'on');
%     xlabel('X'); ylabel('Y'); zlabel('Z');
%     title('Large rotations of thin rod');
% 
%     % Plane of rotation: YZ plane rotated about Z by –45°
%     thetaPlane = deg2rad(45);
%     n = [cos(thetaPlane), sin(thetaPlane), 0];
%     [YY,ZZ] = meshgrid(linspace(-L,L,50), linspace(-L,L,50));
%     XX = -(n(2)*YY + n(3)*ZZ)/n(1);
%     surf(ax,XX,YY,ZZ,'FaceColor',[0.8 0.8 0.8], ...
%          'FaceAlpha',0.2,'EdgeColor','none');
% 
%     % Global axes with dotted negatives
%     axisLength = L;
%     plot3(ax,[0,axisLength],[0,0],[0,0],'r-','LineWidth',2);
%     plot3(ax,[-axisLength,0],[0,0],[0,0],'r--');
%     plot3(ax,[0,0],[0,axisLength],[0,0],'g-','LineWidth',2);
%     plot3(ax,[0,0],[-axisLength,0],[0,0],'g--');
%     plot3(ax,[0,0],[0,0],[0,axisLength],'b-','LineWidth',2);
%     plot3(ax,[0,0],[0,0],[-axisLength,0],'b--');
%     % Ticks
%     for k = 1:axisLength
%         plot3(ax,[k,k],[0,0],[0,0],'r.');
%         plot3(ax,[0,0],[k,k],[0,0],'g.');
%         plot3(ax,[0,0],[0,0],[k,k],'b.');
%     end
% 
%     % Pendulum line and COM marker
%     hLine = plot3(ax,[0,2*x(1)],[0,2*y(1)],[0,2*z(1)],'k-','LineWidth',3);
%     hCOM  = plot3(ax,x(1),y(1),z(1),'ko','MarkerFaceColor','c','MarkerSize',6);
%     hText = text(ax,0,0,lim*0.95,'','FontSize',12,'FontWeight','bold', ...
%                  'HorizontalAlignment','center');
% 
%     first = true;
%     dt = t(2)-t(1);
%     for k = 1:skip:N
%         % Update pendulum line tip and COM marker
%         xt = 2*x(k); yt = 2*y(k); zt = 2*z(k);
%         set(hLine,'XData',[0,xt],'YData',[0,yt],'ZData',[0,zt]);
%         set(hCOM,'XData',x(k),'YData',y(k),'ZData',z(k));
%         hText.String = sprintf('t = %.2f s', t(k));
%         drawnow;
% 
%         % Capture and write GIF frame
%         frame = getframe(fig);
%         im = frame2im(frame);
%         [imind,cm] = rgb2ind(im,256);
%         if first
%             imwrite(imind,cm,gifFilename,'gif','Loopcount',inf,'DelayTime',skip*dt);
%             first = false;
%         else
%             imwrite(imind,cm,gifFilename,'gif','WriteMode','append','DelayTime',skip*dt);
%         end
%     end
% 
%     close(fig);
%     fprintf('Animated GIF saved to %s with skip=%d\n', gifFilename, skip);
% end
% 
% 
% exportPendulumGif(t, x, y, z, 'pendulumlong.gif', 100);
% 
% pendulum_slider_xyz(t, x, y, z);

executionTime = toc
if ~exist('save_filename', 'var') || isempty(save_filename)
    save_filename = sprintf('wrCart_NE_dt_%.1fms.mat', dt * 1000);
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
% % Plotting results for Angular velocities
% figure;
% subplot(3,1,1);
% plot(t,wx,'b-','LineWidth',2); 
% ylim([-3,3]);
% title('Wx');
% ylabel({'Ang velocity Wx'; 's-1'});
% grid on;
% 
% subplot(3,1,2);
% plot(t,wy,'r-','LineWidth',2); 
% ylim([-3,3]);
% title('Wy');
% ylabel({'Ang velocity Wy'; 's-1'});
% grid on;
% 
% subplot(3,1,3);
% plot(t,wz,'g-','LineWidth',2); 
% ylim([-1,1]);
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
% ylim([0,4]);
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
% ylim([-4,4]);
% title('psi1(t), psi2(t)');
% ylabel('psi1(t), psi2(t)');
% grid on;
% lgd = legend('psi1', 'psi2', 'Location', 'northeast', 'FontSize', 8, 'TextColor', 'blue');
% 
% subplot(2,1,2);
% plot(t,psi3,'b-','LineWidth',2); 
% ylim([-4,4]);
% title('psi3(t)');
% ylabel('Psi3');
% grid on;
% 
% % Plotting results for Constraint violation and Energy balance
% figure;
% subplot(4,1,1);
% plot(t,grC,'b-','LineWidth',0.5); 
% ylim([-1*10^-1,1*10^-1]);
% title('∣∣C∣∣');
% ylabel('Violation');
% grid on;
% 
% subplot(4,1,2);
% plot(t,grdC,'g-','LineWidth',1); 
% ylim([-5*10^-2,5*10^-2]);
% title('∣∣C''∣∣');
% ylabel('Violation');
% grid on;
% 
% subplot(4,1,3);
% plot(t,grddC,'r-','LineWidth',1); 
% ylim([-5*10^-13,5*10^-13]);
% title('∣∣C''''∣∣');
% ylabel('Violation');
% grid on;
% 
% subplot(4,1,4);
% plot(t,Ener,'b-','LineWidth',1); 
% ylim([-1*10^1,1*10^1]);
% title('Total Energy');
% xlabel('Time (s)');
% ylabel('Energy balance');
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