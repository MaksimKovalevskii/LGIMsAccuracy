% Compare unconstrained ode45 with constrained HTRef (pos, vel, acc)
unconst = load('UncHTRef.mat');
ref = load('..\Heavy Top\HTRef.mat');

t = unconst.t;

x_ref = interp1(ref.t, ref.x, t, 'spline');
y_ref = interp1(ref.t, ref.y, t, 'spline');
z_ref = interp1(ref.t, ref.z, t, 'spline');

xd_ref = interp1(ref.t, ref.xd, t, 'spline');
yd_ref = interp1(ref.t, ref.yd, t, 'spline');
zd_ref = interp1(ref.t, ref.zd, t, 'spline');

xdd_ref = interp1(ref.t, ref.x_double_prime, t, 'spline');
ydd_ref = interp1(ref.t, ref.y_double_prime, t, 'spline');
zdd_ref = interp1(ref.t, ref.z_double_prime, t, 'spline');

dx = unconst.x(:) - x_ref(:);
dy = unconst.y(:) - y_ref(:);
dz = unconst.z(:) - z_ref(:);

dvx = unconst.xd(:) - xd_ref(:);
dvy = unconst.yd(:) - yd_ref(:);
dvz = unconst.zd(:) - zd_ref(:);

dax = unconst.x_double_prime(:) - xdd_ref(:);
day = unconst.y_double_prime(:) - ydd_ref(:);
daz = unconst.z_double_prime(:) - zdd_ref(:);

disp('RMS |unconst - HTRef|')
disp(['pos  x y z: ', num2str([rms(dx) rms(dy) rms(dz)])])
disp(['vel  x y z: ', num2str([rms(dvx) rms(dvy) rms(dvz)])])
disp(['acc  x y z: ', num2str([rms(dax) rms(day) rms(daz)])])

paperFont = 'Times New Roman';

figure;
subplot(3,1,1);
plot(t, dx, 'b-', 'LineWidth', 1);
grid on;
set(gca, 'FontName', paperFont, 'XMinorGrid', 'off', 'YMinorGrid', 'off');
ylabel('$\Delta R_X$ [m]', 'Interpreter', 'latex', 'FontName', paperFont);

subplot(3,1,2);
plot(t, dy, 'r-', 'LineWidth', 1);
grid on;
set(gca, 'FontName', paperFont, 'XMinorGrid', 'off', 'YMinorGrid', 'off');
ylabel('$\Delta R_Y$ [m]', 'Interpreter', 'latex', 'FontName', paperFont);

subplot(3,1,3);
plot(t, dz, 'g-', 'LineWidth', 1);
grid on;
set(gca, 'FontName', paperFont, 'XMinorGrid', 'off', 'YMinorGrid', 'off');
xlabel('Time [s]', 'FontName', paperFont);
ylabel('$\Delta R_Z$ [m]', 'Interpreter', 'latex', 'FontName', paperFont);

figure;
subplot(3,1,1);
plot(t, dvx, 'b-', 'LineWidth', 1);
grid on;
set(gca, 'FontName', paperFont, 'XMinorGrid', 'off', 'YMinorGrid', 'off');
ylabel('$\Delta \dot{R}_X$ [m/s]', 'Interpreter', 'latex', 'FontName', paperFont);

subplot(3,1,2);
plot(t, dvy, 'r-', 'LineWidth', 1);
grid on;
set(gca, 'FontName', paperFont, 'XMinorGrid', 'off', 'YMinorGrid', 'off');
ylabel('$\Delta \dot{R}_Y$ [m/s]', 'Interpreter', 'latex', 'FontName', paperFont);

subplot(3,1,3);
plot(t, dvz, 'g-', 'LineWidth', 1);
grid on;
set(gca, 'FontName', paperFont, 'XMinorGrid', 'off', 'YMinorGrid', 'off');
xlabel('Time [s]', 'FontName', paperFont);
ylabel('$\Delta \dot{R}_Z$ [m/s]', 'Interpreter', 'latex', 'FontName', paperFont);

figure;
subplot(3,1,1);
plot(t, dax, 'b-', 'LineWidth', 1);
grid on;
set(gca, 'FontName', paperFont, 'XMinorGrid', 'off', 'YMinorGrid', 'off');
ylabel('$\Delta \ddot{R}_X$ [m/s$^2$]', 'Interpreter', 'latex', 'FontName', paperFont);

subplot(3,1,2);
plot(t, day, 'r-', 'LineWidth', 1);
grid on;
set(gca, 'FontName', paperFont, 'XMinorGrid', 'off', 'YMinorGrid', 'off');
ylabel('$\Delta \ddot{R}_Y$ [m/s$^2$]', 'Interpreter', 'latex', 'FontName', paperFont);

subplot(3,1,3);
plot(t, daz, 'g-', 'LineWidth', 1);
grid on;
set(gca, 'FontName', paperFont, 'XMinorGrid', 'off', 'YMinorGrid', 'off');
xlabel('Time [s]', 'FontName', paperFont);
ylabel('$\Delta \ddot{R}_Z$ [m/s$^2$]', 'Interpreter', 'latex', 'FontName', paperFont);
