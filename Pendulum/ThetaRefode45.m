% High-accuracy reference trajectory for the pendulum (ode45).
% Optional workspace variables (set by Run_Pendulum_Batch or before run):
%   save_filename — .mat file to write (default: ExactRef_HighAccuracy.mat)

[t, theta, x, y, z, xd, yd, zd, x_double_prime, y_double_prime, z_double_prime] = exact_pendulum_solution();

if ~exist('save_filename', 'var') || isempty(save_filename)
    save_filename = 'ExactRef_HighAccuracy.mat';
end
save(save_filename, 't', 'theta', 'x', 'y', 'z', ...
    'xd', 'yd', 'zd', 'x_double_prime', 'y_double_prime', 'z_double_prime', '-v7.3');

function [t, theta, x, y, z, xd, yd, zd, xdd, ydd, zdd] = exact_pendulum_solution()
    % Solve the exact nonlinear pendulum equation θ̈ + (mgd/I_eff)sin(θ) = 0
    % using high-accuracy numerical integration

    % Parameters
    d = 3.0;          % COM distance (m)
    m = 1.0;          % mass (kg)
    g = 9.8;          % gravity (m/s²)

    % Correct inertia calculation
    I_cm_trans = 3.0025;     % Transverse inertia about COM
    I_eff = I_cm_trans + m * d^2;  % Parallel axis theorem = 12.0025

    % Initial conditions
    theta0 = 135 * pi/180;   % Initial angle from vertical (radians)
    theta_dot0 = 0.0;        % Initial angular velocity

    fprintf('=== EXACT PENDULUM SOLUTION ===\n');
    fprintf('Parameters:\n');
    fprintf('- Mass: %.1f kg\n', m);
    fprintf('- Distance to COM: %.1f m\n', d);
    fprintf('- I_eff: %.4f kg⋅m²\n', I_eff);
    fprintf('- Initial angle: %.1f°\n', theta0*180/pi);
    fprintf('\n');

    % Define the ODE system:
    function dydt = pendulum_ode(~, y)
        theta = y(1);
        theta_dot = y(2);
        theta_ddot = -(m*g*d/I_eff) * sin(theta);
        dydt = [theta_dot; theta_ddot];
    end

    % Time span and evaluation points
    t_eval = 0:0.0001:100;

    % Initial conditions: [θ(0), θ̇(0)]
    y0 = [theta0; theta_dot0];

    % Solve with high accuracy
    fprintf('Solving nonlinear pendulum ODE with high accuracy...\n');
    % ode45 enforces RelTol >= 100*eps; match that explicitly (avoids odeset warning).
    options = odeset('RelTol', 100*eps, 'AbsTol', 1e-16);
    [t, y] = ode45(@pendulum_ode, t_eval, y0, options);

    theta = y(:, 1);
    theta_dot = y(:, 2);

    % Calculate acceleration
    theta_ddot = -(m*g*d/I_eff) * sin(theta);

    fprintf('Integration successful! %d time points computed.\n', length(t));

    % Convert to Cartesian coordinates
    % Motion plane coordinates (pendulum swings in 45° plane)
    u_coord = d * sin(theta);      % u-axis component
    v_coord = -d * cos(theta);     % v-axis component (hanging down)

    % Transform to original (x, y, z) coordinates
    % Motion is in 45° plane: x = u/√2, y = -u/√2, z = v
    x = u_coord / sqrt(2);
    y = -u_coord / sqrt(2);
    z = v_coord;

    % Calculate velocities
    ud_coord = d * cos(theta) .* theta_dot;
    vd_coord = d * sin(theta) .* theta_dot;

    xd = ud_coord / sqrt(2);
    yd = -ud_coord / sqrt(2);
    zd = vd_coord;

    % Calculate accelerations
    u_ddot = d * (cos(theta) .* theta_ddot - sin(theta) .* theta_dot.^2);
    v_ddot = d * (sin(theta) .* theta_ddot + cos(theta) .* theta_dot.^2);

    xdd = u_ddot / sqrt(2);
    ydd = -u_ddot / sqrt(2);  % Note: same as xdd but negative
    zdd = v_ddot;

    fprintf('Cartesian coordinates computed.\n');
    fprintf('Position ranges:\n');
    fprintf('- X: [%.3f, %.3f] m\n', min(x), max(x));
    fprintf('- Y: [%.3f, %.3f] m\n', min(y), max(y));
    fprintf('- Z: [%.3f, %.3f] m\n', min(z), max(z));

    % Plot results
    figure;
    subplot(3,1,1);
    plot(t, x, 'r-', 'LineWidth', 2);
    title('X Position - Exact Solution (High-Accuracy ODE)');
    ylabel('x (m)');
    grid on;

    subplot(3,1,2);
    plot(t, y, 'g-', 'LineWidth', 2);
    title('Y Position - Exact Solution (High-Accuracy ODE)');
    ylabel('y (m)');
    grid on;

    subplot(3,1,3);
    plot(t, z, 'b-', 'LineWidth', 2);
    title('Z Position - Exact Solution (High-Accuracy ODE)');
    ylabel('z (m)');
    xlabel('Time (s)');
    grid on;

    sgtitle('Exact Analytical Solution for Pendulum - Positions');

    figure;
    subplot(3,1,1);
    plot(t, xd, 'r-', 'LineWidth', 2);
    title('X Velocity - Exact Solution');
    ylabel('xd (m/s)');
    grid on;

    subplot(3,1,2);
    plot(t, yd, 'g-', 'LineWidth', 2);
    title('Y Velocity - Exact Solution');
    ylabel('yd (m/s)');
    grid on;

    subplot(3,1,3);
    plot(t, zd, 'b-', 'LineWidth', 2);
    title('Z Velocity - Exact Solution');
    ylabel('zd (m/s)');
    xlabel('Time (s)');
    grid on;

    sgtitle('Exact Analytical Solution for Pendulum - Velocities');

    figure;
    subplot(3,1,1);
    plot(t, xdd, 'r-', 'LineWidth', 2);
    title('X Acceleration - Exact Solution');
    ylabel('xdd (m/s2)');
    grid on;

    subplot(3,1,2);
    plot(t, ydd, 'g-', 'LineWidth', 2);
    title('Y Acceleration - Exact Solution');
    ylabel('ydd (m/s2)');
    grid on;

    subplot(3,1,3);
    plot(t, zdd, 'b-', 'LineWidth', 2);
    title('Z Acceleration - Exact Solution');
    ylabel('zdd (m/s2)');
    xlabel('Time (s)');
    grid on;

    sgtitle('Exact Analytical Solution for Pendulum - Accelerations');
end
