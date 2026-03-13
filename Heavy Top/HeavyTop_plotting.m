% Define time steps and methods

%time_steps = [0.01, 0.02,0.05,0.1,0.2, 0.5, 1,1.5, 2,3,4,10];
%time_steps = [1,1.5, 2,3,4,10];
time_steps = [3,4,10];

method_names = {'Conv-Quat', 'NE-Quat', 'NE-CRV', 'CRV-LGIM', 'Quat-LGIM'};
file_patterns = {'HTClassicEP_dt_%.2fms.mat', 'HTNE_EP_dt_%.2fms.mat', ...
                'HT_CartNE_dt_%.2fms.mat', 'HT_CartLGIM_dt_%.2fms.mat', 'HT_EPLGIM_dt_%.2fms.mat'};
% method_names = {'Conv-Quat', 'NE-Quat','Quat-LGIM'};
% file_patterns = {'HTClassicEP_dt_%.2fms.mat', 'HTNE_EP_dt_%.2fms.mat', ...
%                  'HT_EPLGIM_dt_%.2fms.mat'};

% Load reference data 
reference_data = load('HTRef.mat');
t_ref = reference_data.t;  % Reference time vector

% Initialize error matrices with NaN
num_methods = length(method_names);
position_errors = zeros(num_methods, length(time_steps));
velocity_errors = zeros(num_methods, length(time_steps));
accel_errors = zeros(num_methods, length(time_steps));
position_errors_rms = nan(num_methods, length(time_steps));
velocity_errors_rms = nan(num_methods, length(time_steps));
accel_errors_rms = nan(num_methods, length(time_steps));
position_errors_max = nan(num_methods, length(time_steps));
velocity_errors_max = nan(num_methods, length(time_steps));
accel_errors_max = nan(num_methods, length(time_steps));
computation_times = zeros(num_methods, length(time_steps));
energy_last_step = nan(num_methods, length(time_steps));
position_errors_rms_norm = nan(num_methods, length(time_steps));
velocity_errors_rms_norm = nan(num_methods, length(time_steps));
accel_errors_rms_norm = nan(num_methods, length(time_steps));
energy_errors_rms=nan(num_methods, length(time_steps));
uc_last_step = nan(num_methods, length(time_steps));

% Main processing loop
for i = 1:length(time_steps)
    dt = time_steps(i);
    
    for m = 1:num_methods
        filename = sprintf(file_patterns{m}, dt);
        
        try
                if ~exist(filename, 'file')
        error('File not found');
    end
            data = load(filename);

switch m
    case {1,2,3,4,5,6}
    if isfield(data, 'Ener')
            energy_last_step(m,i) = data.Ener(end);
            E_meth = data.Ener;
    end
    if isfield(data, 'UC')
            uc_last_step(m,i) = data.UC(end);
    end
end
            % Get method's time vector (convert ms to seconds if needed)
            t_method = data.t;  
          
            % Interpolate reference to method's time points
            x_ref_interp = interp1(t_ref, reference_data.x, t_method, 'spline');
            xd_ref_interp = interp1(t_ref, reference_data.xd, t_method, 'spline');
            xdd_ref_interp = interp1(t_ref, reference_data.x_double_prime, t_method, 'spline');

            y_ref_interp = interp1(t_ref, reference_data.y, t_method, 'spline');
            yd_ref_interp = interp1(t_ref, reference_data.yd, t_method, 'spline');
            ydd_ref_interp = interp1(t_ref, reference_data.y_double_prime, t_method, 'spline');

            z_ref_interp = interp1(t_ref, reference_data.z, t_method, 'spline');
            zd_ref_interp = interp1(t_ref, reference_data.zd, t_method, 'spline');
            zdd_ref_interp = interp1(t_ref, reference_data.z_double_prime, t_method, 'spline');
            
          position_errors_rms(m,i) = rms(abs(x_ref_interp' - data.x));
          velocity_errors_rms(m,i) = rms(abs(xd_ref_interp' - data.xd));
          accel_errors_rms(m,i) = rms(abs(xdd_ref_interp' - data.x_double_prime));

          position_errors_rms_y(m,i) = rms(abs(y_ref_interp' - data.y));
          velocity_errors_rms_y(m,i) = rms(abs(yd_ref_interp' - data.yd));
          accel_errors_rms_y(m,i) = rms(abs(ydd_ref_interp' - data.y_double_prime));

          position_errors_rms_z(m,i) = rms(abs(z_ref_interp' - data.z));
          velocity_errors_rms_z(m,i) = rms(abs(zd_ref_interp' - data.zd));
          accel_errors_rms_z(m,i) = rms(abs(zdd_ref_interp' - data.z_double_prime));          

          energy_errors_rms(m,i) = rms(E_meth);
% Norms
% Compute vector errors at each time point
dx = x_ref_interp' - data.x;
dy = y_ref_interp' - data.y;
dz = z_ref_interp' - data.z;

% Magnitude of error at each step
err_mag = sqrt(dx.^2 + dy.^2 + dz.^2);

% True RMS of the norm
position_errors_rms_norm(m,i) = rms(err_mag);
% Velocity
vdx = xd_ref_interp' - data.xd;
vdy = yd_ref_interp' - data.yd;
vdz = zd_ref_interp' - data.zd;
vel_err_mag = sqrt(vdx.^2 + vdy.^2 + vdz.^2);
velocity_errors_rms_norm(m,i) = rms(vel_err_mag);

% Acceleration
adx = xdd_ref_interp' - data.x_double_prime;
ady = ydd_ref_interp' - data.y_double_prime;
adz = zdd_ref_interp' - data.z_double_prime;
acc_err_mag = sqrt(adx.^2 + ady.^2 + adz.^2);
accel_errors_rms_norm(m,i) = rms(acc_err_mag);


% max errors
          position_errors_max(m,i) = max(abs(x_ref_interp' - data.x));
          velocity_errors_max(m,i) = max(abs(xd_ref_interp' - data.xd));
          accel_errors_max(m,i) = max(abs(xdd_ref_interp' - data.x_double_prime));

          % Calculate errors at final time step
            position_errors(m, i) = abs(reference_data.x(end) - data.x(end));
            velocity_errors(m, i) = abs(reference_data.xd(end) - data.xd(end));
            accel_errors(m, i) = abs(reference_data.x_double_prime(end) - data.x_double_prime(end));
            
                        % Get computation time if available
            if isfield(data, 'executionTime')
                computation_times(m, i) = data.executionTime;
            end
        catch
            warning('Missing data: %s', filename);
            position_errors_rms(m, i) = NaN;
            velocity_errors_rms(m, i) = NaN;
            accel_errors_rms(m, i) = NaN;
            uc_last_step(m, i) = NaN;
        end
    end
end



% Plotting setup
line_specs = {
    {'-.o', 'Color', [0 0.4470 0.7410], 'MarkerSize', 8, 'LineWidth', 1.5},  % Conv-Quat
    {'k--', 'MarkerSize', 9, 'LineWidth', 1.4},  % NE-Quat
    {'-.d', 'Color', [0.4940 0.1840 0.5560], 'MarkerSize', 7, 'LineWidth', 1.5},  % NE-CRV
    {'--', 'Color', [0.4660 0.6740 0.1880], 'MarkerSize', 8, 'LineWidth', 1.6},  % CRV-LGIM
    {'--o', 'Color', [0.85 0.6740 0.556], 'MarkerSize', 7, 'LineWidth', 1.8}   % Quat-LGIM
};

% 1 RMS Error - X axis component  
createComparisonPlot(position_errors_rms, velocity_errors_rms, accel_errors_rms, ...
                    time_steps, method_names, line_specs, ...
                    'X axis component', 'RMS Error');

% 2 RMS Error - Y axis component
createComparisonPlot(position_errors_rms_y, velocity_errors_rms_y, accel_errors_rms_y, ...
                    time_steps, method_names, line_specs, ...
                    'Y axis component', 'RMS Error');

% 3 RMS Error - Z axis component
createComparisonPlot(position_errors_rms_z, velocity_errors_rms_z, accel_errors_rms_z, ...
                    time_steps, method_names, line_specs, ...
                    'Z axis component', 'RMS Error');

% 4 MAX ABS Error - X axis component
createComparisonPlot(position_errors_max, velocity_errors_max, accel_errors_max, ...
                    time_steps, method_names, line_specs, ...
                    'X axis component', 'MAX ABS Error');

% 5 RMS Error - The norm (magnitude)
createNormComparisonPlot({position_errors_rms_norm, velocity_errors_rms_norm, accel_errors_rms_norm}, ...
                        time_steps, method_names, line_specs, ...
                        'The norm (magnitude)', 'RMS Error');

% 6. Computation time 
figure;
hold on;
current_data = [];
plot_handles = gobjects(num_methods,1);
for m = 1:num_methods
    valid_points = ~isnan(computation_times(m,:));
    if any(valid_points)
        data_to_plot = computation_times(m,valid_points);
        plot_handles(m) = loglog(time_steps(valid_points), data_to_plot, line_specs{m}{:});
        current_data = [current_data, data_to_plot];
    end
end
grid on;
setLogYTicks(gca, current_data);
%title('Computation time');
ylabel('Time (s)');
xlabel('Time Step (ms)');
set(gca, 'XScale', 'log', 'YScale', 'log');
valid_methods = any(~isnan(computation_times), 2);
legend(plot_handles(valid_methods), method_names(valid_methods), 'Location', 'northwest');

% 7. Energy at Final Time Step 
figure;
hold on;
current_data = [];
plot_handles = gobjects(num_methods,1);
for m = 1:num_methods
    valid_points = ~isnan(energy_last_step(m,:));
    if any(valid_points)
        data_to_plot = abs(energy_last_step(m,valid_points));
        plot_handles(m) = loglog(time_steps(valid_points), data_to_plot, line_specs{m}{:});
        current_data = [current_data, data_to_plot];
    end
end
grid on;
setLogYTicks(gca, current_data);
%title('Energy at Final Time Step vs Time Step Size');
xlabel('Time Step (ms)');
ylabel('|Energy| (J)');
set(gca, 'XScale', 'log', 'YScale', 'log');
valid_methods = any(~isnan(energy_last_step), 2);
legend(plot_handles(valid_methods), method_names(valid_methods), 'Location', 'northwest');

% 8. Energy RMS
figure;
hold on;
current_data = [];
plot_handles = gobjects(num_methods,1);
for m = 1:num_methods
    valid_points = ~isnan(energy_errors_rms(m,:));
    if any(valid_points)
        data_to_plot = abs(energy_errors_rms(m,valid_points));
        plot_handles(m) = loglog(time_steps(valid_points), data_to_plot, line_specs{m}{:});
        current_data = [current_data, data_to_plot];
    end
end
grid on;
setLogYTicks(gca, current_data);
%title('RMS Energy balance violation vs Time Step Size');
xlabel('Time Step (ms)');
ylabel('|Energy| (J)');
set(gca, 'XScale', 'log', 'YScale', 'log');
valid_methods = any(~isnan(energy_errors_rms), 2);
legend(plot_handles(valid_methods), method_names(valid_methods), 'Location', 'northwest');

%RMS position X Accuracy vs Computation time
figure;
hold on;
for m = 1:num_methods
    % Only include time steps that have valid errors and valid computational time
    valid_idx = ~isnan(position_errors_rms(m,:)) & ~isnan(computation_times(m,:));
    if any(valid_idx)
        % X-axis: error, Y-axis: computational time
        loglog(position_errors_rms(m,valid_idx), computation_times(m,valid_idx), line_specs{m}{:});
    end
end
grid on;
xlabel('RMS Position X Error (m)');   % or whatever error metric you've chosen
ylabel('Computational Time (s)');
legend(method_names, 'Location','northwest');
%title('Computational Time vs. Accuracy (RMS Position X Error)');
set(gca, 'XScale','log', 'YScale','log');

%RMS norm position Accuracy vs Computation time
figure;
hold on;
for m = 1:num_methods
    % Only include time steps that have valid errors and valid computational time
    valid_idx = ~isnan(position_errors_rms_norm(m,:)) & ~isnan(computation_times(m,:));
    if any(valid_idx)
        % X-axis: error, Y-axis: computational time
        loglog(position_errors_rms_norm(m,valid_idx), computation_times(m,valid_idx), line_specs{m}{:});
    end
end
grid on;
xlabel('RMS Position Norm Error (m)');   % or whatever error metric you've chosen
ylabel('Computational Time (s)');
legend(method_names, 'Location','northwest');
%title('Computational Time vs. Accuracy (RMS Position Norm Error)');
set(gca, 'XScale','log', 'YScale','log');

%RMS norm VELOCITY Accuracy vs Computation time
figure;
hold on;
for m = 1:num_methods
    % Only include time steps that have valid errors and valid computational time
    valid_idx = ~isnan(velocity_errors_rms_norm(m,:)) & ~isnan(computation_times(m,:));
    if any(valid_idx)
        % X-axis: error, Y-axis: computational time
        loglog(velocity_errors_rms_norm(m,valid_idx), computation_times(m,valid_idx), line_specs{m}{:});
    end
end
grid on;
xlabel('RMS Velocity Norm Error (m)');   % or whatever error metric you've chosen
ylabel('Computational Time (s)');
legend(method_names, 'Location','northwest');
%title('Computational Time vs. Accuracy (RMS Velocity Norm Error)');
set(gca, 'XScale','log', 'YScale','log');

%RMS norm Acceleration Accuracy vs Computation time
figure;
hold on;
for m = 1:num_methods
    % Only include time steps that have valid errors and valid computational time
    valid_idx = ~isnan(accel_errors_rms_norm(m,:)) & ~isnan(computation_times(m,:));
    if any(valid_idx)
        % X-axis: error, Y-axis: computational time
        loglog(accel_errors_rms_norm(m,valid_idx), computation_times(m,valid_idx), line_specs{m}{:});
    end
end
grid on;
xlabel('RMS Acceleration Norm Error (m)');   % or whatever error metric you've chosen
ylabel('Computational Time (s)');
legend(method_names, 'Location','northwest');
%title('Computational Time vs. Accuracy (RMS Acceleration Norm Error)');
set(gca, 'XScale','log', 'YScale','log');

% Only include methods with valid data in the legend
valid_methods = any(~isnan(computation_times), 2);
legend(plot_handles(valid_methods), method_names(valid_methods), 'Location', 'northwest');


%UNIT CONSTRAINT violation (only EP methods: Conv-Quat, NE-Quat, Quat-LGIM)
ep_method_inds = [1, 2, 5];  % indices in method_names for Euler-parameter formulations
ep_method_names = method_names(ep_method_inds);
num_ep = length(ep_method_inds);

figure;
hold on;
plot_handles = gobjects(num_ep, 1);
for k = 1:num_ep
    m = ep_method_inds(k);
    valid_points = ~isnan(uc_last_step(m,:)) & uc_last_step(m,:) ~= 0;
    if any(valid_points)
        plot_handles(k) = loglog(time_steps(valid_points), abs(uc_last_step(m,valid_points)), ...
                                 line_specs{m}{:}, 'LineWidth', 2.5);
    end
end
grid on;

xlabel('Time Step (ms)', 'Interpreter', 'latex', 'FontSize', 11);
ylabel('EP Unit Constraint Deviation', 'Interpreter', 'latex', 'FontSize', 11);
set(gca, 'XScale', 'log', 'YScale', 'log', 'FontSize', 10);

valid_methods = isgraphics(plot_handles);
legend(plot_handles(valid_methods), ep_method_names(valid_methods), ...
       'Interpreter', 'latex', 'FontSize', 10, 'Location', 'best');
