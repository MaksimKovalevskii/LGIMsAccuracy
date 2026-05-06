% Pendulum: load method outputs vs ExactRef_HighAccuracy, then plot comparisons.

time_steps = [0.2, 0.5, 1, 2, 5, 10, 20, 50, 100];
method_names = {'Conv-Quat', 'Conv-CRV', 'NE-Quat', 'NE-CRV', 'CRV-LGIM', 'Quat-LGIM'};
file_patterns = {'EP_dt_%.1fms.mat', 'Cart_dt_%.1fms.mat', 'EP_NE_dt_%.1fms.mat', ...
    'wrCart_NE_dt_%.1fms.mat', 'wrLGIM_dt_%.1fms.mat', 'EP_LGIM_dt_%.1fms.mat'};

pendulumFont = 'Times New Roman';

reference_data = load('ExactRef_HighAccuracy.mat');
t_ref = reference_data.t;

num_methods = length(method_names);
position_errors = zeros(num_methods, length(time_steps));
velocity_errors = zeros(num_methods, length(time_steps));
accel_errors = zeros(num_methods, length(time_steps));
position_errors_rms = nan(num_methods, length(time_steps));
velocity_errors_rms = nan(num_methods, length(time_steps));
accel_errors_rms = nan(num_methods, length(time_steps));
position_errors_rms_y = nan(num_methods, length(time_steps));
velocity_errors_rms_y = nan(num_methods, length(time_steps));
accel_errors_rms_y = nan(num_methods, length(time_steps));
position_errors_rms_z = nan(num_methods, length(time_steps));
velocity_errors_rms_z = nan(num_methods, length(time_steps));
accel_errors_rms_z = nan(num_methods, length(time_steps));
position_errors_max = nan(num_methods, length(time_steps));
velocity_errors_max = nan(num_methods, length(time_steps));
accel_errors_max = nan(num_methods, length(time_steps));
computation_times = zeros(num_methods, length(time_steps));
energy_last_step = nan(num_methods, length(time_steps));
position_errors_rms_norm = nan(num_methods, length(time_steps));
velocity_errors_rms_norm = nan(num_methods, length(time_steps));
accel_errors_rms_norm = nan(num_methods, length(time_steps));
energy_errors_rms = nan(num_methods, length(time_steps));
uc_last_step = nan(num_methods, length(time_steps));

for i = 1:length(time_steps)
    dt = time_steps(i);

    for m = 1:num_methods
        filename = sprintf(file_patterns{m}, dt);

        try
            if ~exist(filename, 'file')
                error('File not found');
            end
            data = load(filename);
            E_meth = [];

            switch m
                case {1}
                    if isfield(data, 'Enernew2')
                        energy_last_step(m, i) = data.Enernew2(end);
                        E_meth = data.Enernew2;
                    end
                    if isfield(data, 'UC2')
                        uc_last_step(m, i) = data.UC2(end);
                    end
                case {2}
                    if isfield(data, 'Enernew')
                        energy_last_step(m, i) = data.Enernew(end);
                        E_meth = data.Enernew;
                    end
                case {3, 4, 5, 6}
                    if isfield(data, 'Ener')
                        energy_last_step(m, i) = data.Ener(end);
                        E_meth = data.Ener;
                    end
                    if isfield(data, 'UC')
                        uc_last_step(m, i) = data.UC(end);
                    end
                otherwise
            end

            t_method = data.t;

            x_ref_interp = interp1(t_ref, reference_data.x, t_method, 'spline');
            xd_ref_interp = interp1(t_ref, reference_data.xd, t_method, 'spline');
            xdd_ref_interp = interp1(t_ref, reference_data.x_double_prime, t_method, 'spline');

            y_ref_interp = interp1(t_ref, reference_data.y, t_method, 'spline');
            yd_ref_interp = interp1(t_ref, reference_data.yd, t_method, 'spline');
            ydd_ref_interp = interp1(t_ref, reference_data.y_double_prime, t_method, 'spline');

            z_ref_interp = interp1(t_ref, reference_data.z, t_method, 'spline');
            zd_ref_interp = interp1(t_ref, reference_data.zd, t_method, 'spline');
            zdd_ref_interp = interp1(t_ref, reference_data.z_double_prime, t_method, 'spline');

            position_errors_rms(m, i) = rms(abs(x_ref_interp' - data.x));
            velocity_errors_rms(m, i) = rms(abs(xd_ref_interp' - data.xd));
            accel_errors_rms(m, i) = rms(abs(xdd_ref_interp' - data.x_double_prime));

            position_errors_rms_y(m, i) = rms(abs(y_ref_interp' - data.y));
            velocity_errors_rms_y(m, i) = rms(abs(yd_ref_interp' - data.yd));
            accel_errors_rms_y(m, i) = rms(abs(ydd_ref_interp' - data.y_double_prime));

            position_errors_rms_z(m, i) = rms(abs(z_ref_interp' - data.z));
            velocity_errors_rms_z(m, i) = rms(abs(zd_ref_interp' - data.zd));
            accel_errors_rms_z(m, i) = rms(abs(zdd_ref_interp' - data.z_double_prime));

            if ~isempty(E_meth)
                energy_errors_rms(m, i) = rms(E_meth);
            end

            dx = x_ref_interp' - data.x;
            dy = y_ref_interp' - data.y;
            dz = z_ref_interp' - data.z;
            err_mag = sqrt(dx.^2 + dy.^2 + dz.^2);
            position_errors_rms_norm(m, i) = rms(err_mag);

            vdx = xd_ref_interp' - data.xd;
            vdy = yd_ref_interp' - data.yd;
            vdz = zd_ref_interp' - data.zd;
            vel_err_mag = sqrt(vdx.^2 + vdy.^2 + vdz.^2);
            velocity_errors_rms_norm(m, i) = rms(vel_err_mag);

            adx = xdd_ref_interp' - data.x_double_prime;
            ady = ydd_ref_interp' - data.y_double_prime;
            adz = zdd_ref_interp' - data.z_double_prime;
            acc_err_mag = sqrt(adx.^2 + ady.^2 + adz.^2);
            accel_errors_rms_norm(m, i) = rms(acc_err_mag);

            position_errors_max(m, i) = max(abs(x_ref_interp' - data.x));
            velocity_errors_max(m, i) = max(abs(xd_ref_interp' - data.xd));
            accel_errors_max(m, i) = max(abs(xdd_ref_interp' - data.x_double_prime));

            position_errors(m, i) = abs(reference_data.x(end) - data.x(end));
            velocity_errors(m, i) = abs(reference_data.xd(end) - data.xd(end));
            accel_errors(m, i) = abs(reference_data.x_double_prime(end) - data.x_double_prime(end));

            if isfield(data, 'executionTime')
                computation_times(m, i) = data.executionTime;
            end
        catch
            warning('Missing data: %s', filename);
            position_errors_rms(m, i) = NaN;
            velocity_errors_rms(m, i) = NaN;
            accel_errors_rms(m, i) = NaN;
            position_errors_rms_y(m, i) = NaN;
            velocity_errors_rms_y(m, i) = NaN;
            accel_errors_rms_y(m, i) = NaN;
            position_errors_rms_z(m, i) = NaN;
            velocity_errors_rms_z(m, i) = NaN;
            accel_errors_rms_z(m, i) = NaN;
            position_errors_rms_norm(m, i) = NaN;
            velocity_errors_rms_norm(m, i) = NaN;
            accel_errors_rms_norm(m, i) = NaN;
            uc_last_step(m, i) = NaN;
        end
    end
end

line_specs = {
    {'-.o', 'Color', [0 0.4470 0.7410], 'MarkerSize', 8, 'LineWidth', 1.5}; ...
    {'--s', 'Color', [0.8500 0.3250 0.0980], 'MarkerSize', 10, 'LineWidth', 1.3}; ...
    {'k--', 'MarkerSize', 9, 'LineWidth', 1.4}; ...
    {'-.d', 'Color', [0.4940 0.1840 0.5560], 'MarkerSize', 7, 'LineWidth', 1.5}; ...
    {'--', 'Color', [0.4660 0.6740 0.1880], 'MarkerSize', 8, 'LineWidth', 1.6}; ...
    {'--o', 'Color', [0.85 0.6740 0.556], 'MarkerSize', 7, 'LineWidth', 1.8}
    };

%% RMS error — X, Y, Z components
createComparisonPlot(position_errors_rms, velocity_errors_rms, accel_errors_rms, ...
    time_steps, method_names, line_specs, ...
    'X axis component', 'RMS Error');

createComparisonPlot(position_errors_rms_y, velocity_errors_rms_y, accel_errors_rms_y, ...
    time_steps, method_names, line_specs, ...
    'Y axis component', 'RMS Error');

createComparisonPlot(position_errors_rms_z, velocity_errors_rms_z, accel_errors_rms_z, ...
    time_steps, method_names, line_specs, ...
    'Z axis component', 'RMS Error');

%% Maximum absolute error — X component
createComparisonPlot(position_errors_max, velocity_errors_max, accel_errors_max, ...
    time_steps, method_names, line_specs, ...
    'X axis component', 'MAX ABS Error');

%% RMS error — Euclidean norm (magnitude)
createNormComparisonPlot({position_errors_rms_norm, velocity_errors_rms_norm, accel_errors_rms_norm}, ...
    time_steps, method_names, line_specs, ...
    'The norm (magnitude)', 'RMS Error');

%% Error at final time step — X component
createComparisonPlot(position_errors, velocity_errors, accel_errors, ...
    time_steps, method_names, line_specs, ...
    'X axis component', 'Error at the final time step');

%% Computational time vs time step
figure;
hold on;
current_data = [];
plot_handles = gobjects(num_methods, 1);
for m = 1:num_methods
    valid_points = ~isnan(computation_times(m, :));
    if any(valid_points)
        data_to_plot = computation_times(m, valid_points);
        plot_handles(m) = loglog(time_steps(valid_points), data_to_plot, line_specs{m}{:});
        current_data = [current_data, data_to_plot];
    end
end
grid on;
set(gca, 'XScale', 'log', 'YScale', 'log', ...
    'XMinorGrid', 'off', 'YMinorGrid', 'off', 'FontName', pendulumFont);
setLogYTicks(gca, current_data, 'auto');
set(gca, 'XMinorGrid', 'off', 'YMinorGrid', 'off', 'FontName', pendulumFont);
ylabel('Computational time [s]', 'FontName', pendulumFont);
xlabel('Time step [ms]', 'FontName', pendulumFont);
valid_methods = any(~isnan(computation_times), 2);
legend(plot_handles(valid_methods), method_names(valid_methods), ...
    'Location', 'northwest', 'FontName', pendulumFont);

%% Energy magnitude at final time step
figure;
hold on;
current_data = [];
plot_handles = gobjects(num_methods, 1);
for m = 1:num_methods
    valid_points = ~isnan(energy_last_step(m, :));
    if any(valid_points)
        data_to_plot = abs(energy_last_step(m, valid_points));
        plot_handles(m) = loglog(time_steps(valid_points), data_to_plot, line_specs{m}{:});
        current_data = [current_data, data_to_plot];
    end
end
grid on;
set(gca, 'XScale', 'log', 'YScale', 'log', ...
    'XMinorGrid', 'off', 'YMinorGrid', 'off', 'FontName', pendulumFont);
setLogYTicks(gca, current_data, 'auto');
set(gca, 'XMinorGrid', 'off', 'YMinorGrid', 'off', 'FontName', pendulumFont);
xlabel('Time step [ms]', 'FontName', pendulumFont);
ylabel('Energy magnitude [J]', 'FontName', pendulumFont);
valid_methods = any(~isnan(energy_last_step), 2);
legend(plot_handles(valid_methods), method_names(valid_methods), ...
    'Location', 'northwest', 'FontName', pendulumFont);

%% RMS energy magnitude
figure;
hold on;
current_data = [];
plot_handles = gobjects(num_methods, 1);
for m = 1:num_methods
    valid_points = ~isnan(energy_errors_rms(m, :));
    if any(valid_points)
        data_to_plot = abs(energy_errors_rms(m, valid_points));
        plot_handles(m) = loglog(time_steps(valid_points), data_to_plot, line_specs{m}{:});
        current_data = [current_data, data_to_plot];
    end
end
grid on;
set(gca, 'XScale', 'log', 'YScale', 'log', ...
    'XMinorGrid', 'off', 'YMinorGrid', 'off', 'FontName', pendulumFont);
setLogYTicks(gca, current_data, 'auto');
set(gca, 'XMinorGrid', 'off', 'YMinorGrid', 'off', 'FontName', pendulumFont);
xlabel('Time step [ms]', 'FontName', pendulumFont);
ylabel('RMS energy magnitude [J]', 'FontName', pendulumFont);
valid_methods = any(~isnan(energy_errors_rms), 2);
legend(plot_handles(valid_methods), method_names(valid_methods), ...
    'Location', 'northwest', 'FontName', pendulumFont);

%% RMS position error (X) vs computational time
figure;
hold on;
all_posX = position_errors_rms(~isnan(position_errors_rms) & position_errors_rms > 0);
for m = 1:num_methods
    valid_idx = ~isnan(position_errors_rms(m, :)) & ~isnan(computation_times(m, :));
    if any(valid_idx)
        loglog(position_errors_rms(m, valid_idx), computation_times(m, valid_idx), line_specs{m}{:});
    end
end
grid on;
set(gca, 'XScale', 'log', 'YScale', 'log', ...
    'XMinorGrid', 'off', 'YMinorGrid', 'off', 'FontName', pendulumFont);
setLogYTicks(gca, computation_times(~isnan(computation_times)), 'auto');
set(gca, 'XMinorGrid', 'off', 'YMinorGrid', 'off', 'FontName', pendulumFont);
xlabel('RMS position error, X component [m]', 'FontName', pendulumFont);
ylabel('Computational time [s]', 'FontName', pendulumFont);
legend(method_names, 'Location', 'northwest', 'FontName', pendulumFont);

if ~isempty(all_posX)
    xMin = max(min(all_posX) * 0.9, eps);
else
    xMin = eps;
end
if xMin >= 1e3
    xMin = 1e-6;
end
set(gca, 'XLim', [xMin, 1e3]);

%% RMS position norm error vs computational time
figure;
hold on;
for m = 1:num_methods
    valid_idx = ~isnan(position_errors_rms_norm(m, :)) & ~isnan(computation_times(m, :));
    if any(valid_idx)
        loglog(position_errors_rms_norm(m, valid_idx), computation_times(m, valid_idx), line_specs{m}{:});
    end
end
grid on;
set(gca, 'XScale', 'log', 'YScale', 'log', ...
    'XMinorGrid', 'off', 'YMinorGrid', 'off', 'FontName', pendulumFont);
setLogYTicks(gca, computation_times(~isnan(computation_times)), 'auto');
set(gca, 'XMinorGrid', 'off', 'YMinorGrid', 'off', 'FontName', pendulumFont);
xlabel('RMS position norm error [m]', 'FontName', pendulumFont);
ylabel('Computational time [s]', 'FontName', pendulumFont);
legend(method_names, 'Location', 'northwest', 'FontName', pendulumFont);

%% RMS velocity norm error vs computational time
figure;
hold on;
for m = 1:num_methods
    valid_idx = ~isnan(velocity_errors_rms_norm(m, :)) & ~isnan(computation_times(m, :));
    if any(valid_idx)
        loglog(velocity_errors_rms_norm(m, valid_idx), computation_times(m, valid_idx), line_specs{m}{:});
    end
end
grid on;
set(gca, 'XScale', 'log', 'YScale', 'log', ...
    'XMinorGrid', 'off', 'YMinorGrid', 'off', 'FontName', pendulumFont);
setLogYTicks(gca, computation_times(~isnan(computation_times)), 'auto');
set(gca, 'XMinorGrid', 'off', 'YMinorGrid', 'off', 'FontName', pendulumFont);
xlabel('RMS velocity norm error [m/s]', 'FontName', pendulumFont);
ylabel('Computational time [s]', 'FontName', pendulumFont);
legend(method_names, 'Location', 'northwest', 'FontName', pendulumFont);

%% RMS acceleration norm error vs computational time
figure;
hold on;
for m = 1:num_methods
    valid_idx = ~isnan(accel_errors_rms_norm(m, :)) & ~isnan(computation_times(m, :));
    if any(valid_idx)
        loglog(accel_errors_rms_norm(m, valid_idx), computation_times(m, valid_idx), line_specs{m}{:});
    end
end
grid on;
set(gca, 'XScale', 'log', 'YScale', 'log', ...
    'XMinorGrid', 'off', 'YMinorGrid', 'off', 'FontName', pendulumFont);
setLogYTicks(gca, computation_times(~isnan(computation_times)), 'auto');
set(gca, 'XMinorGrid', 'off', 'YMinorGrid', 'off', 'FontName', pendulumFont);
xlabel('RMS acceleration norm error [m/s^2]', 'FontName', pendulumFont);
ylabel('Computational time [s]', 'FontName', pendulumFont);
legend(method_names, 'Location', 'northwest', 'FontName', pendulumFont);

%% EP unit constraint (EP-based methods only)
ep_method_inds = [1, 3, 6];
ep_method_names = method_names(ep_method_inds);
num_ep = length(ep_method_inds);

figure;
hold on;
plot_handles = gobjects(num_ep, 1);
current_data = [];
for k = 1:num_ep
    m = ep_method_inds(k);
    valid_points = ~isnan(uc_last_step(m, :)) & uc_last_step(m, :) ~= 0;
    if any(valid_points)
        ydata = abs(uc_last_step(m, valid_points));
        plot_handles(k) = loglog(time_steps(valid_points), ydata, line_specs{m}{:}, 'LineWidth', 2.5);
        current_data = [current_data, ydata];
    end
end
grid on;
set(gca, 'XScale', 'log', 'YScale', 'log', ...
    'XMinorGrid', 'off', 'YMinorGrid', 'off', 'FontName', pendulumFont);
setLogYTicks(gca, current_data, 'every3');
set(gca, 'XMinorGrid', 'off', 'YMinorGrid', 'off', 'FontName', pendulumFont);

xlabel('Time step [ms]', 'FontName', pendulumFont);
ylabel('EP unit constraint deviation [-]', 'FontName', pendulumFont);

valid_methods = isgraphics(plot_handles);
legend(plot_handles(valid_methods), ep_method_names(valid_methods), ...
    'Location', 'best', 'FontName', pendulumFont);
