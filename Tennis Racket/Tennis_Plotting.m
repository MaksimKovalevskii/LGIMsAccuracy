% Tennis Racket: load method outputs vs TennisRef, then plot comparisons once each.

time_steps = [0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50, 100, 200, 500];

method_names = {'Conv-Quat', 'Conv-CRV', 'NE-Quat', 'NE-CRV', 'CRV-LGIM', 'Quat-LGIM'};
file_patterns = {'ClassicEP_dt_%.1fms.mat', 'ClassicCart_dt_%.1fms.mat', 'EP_NE_dt_%.1fms.mat', ...
    'wrCart_NE_dt_%.1fms.mat', 'wrLGIM_dt_%.1fms.mat', 'EP_LGIM_dt_%.1fms.mat'};

tennisFont = 'Times New Roman';

reference_data = load('TennisRef.mat');
t_ref = reference_data.t;

num_methods = length(method_names);
wx_errors = zeros(num_methods, length(time_steps));
wy_errors = zeros(num_methods, length(time_steps));
wz_errors = zeros(num_methods, length(time_steps));
wx_errors_rms = nan(num_methods, length(time_steps));
wy_errors_rms = nan(num_methods, length(time_steps));
wz_errors_rms = nan(num_methods, length(time_steps));
wx_errors_max = nan(num_methods, length(time_steps));
wy_errors_max = nan(num_methods, length(time_steps));
wz_errors_max = nan(num_methods, length(time_steps));
w_errors_rms_norm = nan(num_methods, length(time_steps));
computation_times = zeros(num_methods, length(time_steps));
energy_last_step = nan(num_methods, length(time_steps));
angular_momentum_last_step = nan(num_methods, length(time_steps));
uc_last_step = nan(num_methods, length(time_steps));

for i = 1:length(time_steps)
    dt = time_steps(i);
    for m = 1:num_methods
        if m == 2
            filename = sprintf('ClassicCart_dt_%gms.mat', dt);
        else
            filename = sprintf(file_patterns{m}, dt);
        end
        try
            if ~exist(filename, 'file')
                error('File not found');
            end
            data = load(filename);

            if isfield(data, 'Ener')
                energy_last_step(m, i) = data.Ener(end);
            end
            if isfield(data, 'H_norm')
                angular_momentum_last_step(m, i) = data.H_norm(end);
            end
            if isfield(data, 'UC')
                uc_last_step(m, i) = data.UC(end);
            end

            t_method = data.t;

            wx_ref_interp = interp1(t_ref, reference_data.wx, t_method, 'spline');
            wy_ref_interp = interp1(t_ref, reference_data.wy, t_method, 'spline');
            wz_ref_interp = interp1(t_ref, reference_data.wz, t_method, 'spline');

            wx_ref_interp = wx_ref_interp(:);
            wy_ref_interp = wy_ref_interp(:);
            wz_ref_interp = wz_ref_interp(:);
            wx_data = data.wx(:);
            wy_data = data.wy(:);
            wz_data = data.wz(:);

            wx_errors_rms(m, i) = rms(abs(wx_ref_interp - wx_data));
            wy_errors_rms(m, i) = rms(abs(wy_ref_interp - wy_data));
            wz_errors_rms(m, i) = rms(abs(wz_ref_interp - wz_data));

            wx_errors_max(m, i) = max(abs(wx_ref_interp - wx_data));
            wy_errors_max(m, i) = max(abs(wy_ref_interp - wy_data));
            wz_errors_max(m, i) = max(abs(wz_ref_interp - wz_data));

            dx = wx_ref_interp - wx_data;
            dy = wy_ref_interp - wy_data;
            dz = wz_ref_interp - wz_data;
            err_mag = sqrt(dx.^2 + dy.^2 + dz.^2);
            w_errors_rms_norm(m, i) = rms(err_mag);

            wx_errors(m, i) = abs(reference_data.wx(end) - data.wx(end));
            wy_errors(m, i) = abs(reference_data.wy(end) - data.wy(end));
            wz_errors(m, i) = abs(reference_data.wz(end) - data.wz(end));

            if isfield(data, 'executionTime')
                computation_times(m, i) = data.executionTime;
            end
        catch
            warning('Missing data: %s', filename);
            wx_errors_rms(m, i) = NaN;
            wy_errors_rms(m, i) = NaN;
            wz_errors_rms(m, i) = NaN;
            w_errors_rms_norm(m, i) = NaN;
            energy_last_step(m, i) = NaN;
            angular_momentum_last_step(m, i) = NaN;
            uc_last_step(m, i) = NaN;
        end
    end
end

line_specs = {
    {'-.o', 'Color', [0 0.4470 0.7410], 'MarkerSize', 8, 'LineWidth', 1.5}, ...
    {'--s', 'Color', [0.8500 0.3250 0.0980], 'MarkerSize', 10, 'LineWidth', 1.3}, ...
    {'k--', 'MarkerSize', 7, 'LineWidth', 1.4}, ...
    {'-.d', 'Color', [0.4940 0.1840 0.5560], 'MarkerSize', 7, 'LineWidth', 1.5}, ...
    {'--', 'Color', [0.4660 0.6740 0.1880], 'MarkerSize', 8, 'LineWidth', 1.6}, ...
    {'--o', 'Color', [0.85 0.6740 0.556], 'MarkerSize', 7, 'LineWidth', 1.8}
};

%% Three-panel figures: RMS, max, and final-time component errors vs time step
plotDefs(1).errLabel = 'RMS error';
plotDefs(1).mats = {wx_errors_rms, wy_errors_rms, wz_errors_rms};
plotDefs(1).ylabs = {'$\omega_x$ [rad/s]', '$\omega_y$ [rad/s]', '$\omega_z$ [rad/s]'};

plotDefs(2).errLabel = 'Maximum absolute error';
plotDefs(2).mats = {wx_errors_max, wy_errors_max, wz_errors_max};
plotDefs(2).ylabs = {'$\omega_x$ [rad/s]', '$\omega_y$ [rad/s]', '$\omega_z$ [rad/s]'};

plotDefs(3).errLabel = 'Final-time absolute error';
plotDefs(3).mats = {wx_errors, wy_errors, wz_errors};
plotDefs(3).ylabs = {'$\omega_x$ [rad/s]', '$\omega_y$ [rad/s]', '$\omega_z$ [rad/s]'};

for pd = 1:numel(plotDefs)
    figure;
    set(gcf, 'Name', plotDefs(pd).errLabel);
    annotation(gcf, 'textbox', [0.005, 0.5, 0.02, 0.1], ...
        'String', plotDefs(pd).errLabel, 'EdgeColor', 'none', 'Rotation', 90, ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
        'FontName', tennisFont);

    mats = plotDefs(pd).mats;
    ylabs = plotDefs(pd).ylabs;

    for sp = 1:3
        subplot(3, 1, sp);
        hold on;
        grid on;
        set(gca, 'XScale', 'log', 'YScale', 'log', ...
            'XMinorGrid', 'off', 'YMinorGrid', 'off', ...
            'FontName', tennisFont);

        plot_handles = gobjects(0);
        legend_entries = {};
        current_data = [];

        for m = 1:num_methods
            idx = find(~isnan(mats{sp}(m, :)));
            idx = idx(idx <= numel(time_steps));
            if isempty(idx)
                continue;
            end
            ydata = mats{sp}(m, idx);
            h = loglog(time_steps(idx), ydata, line_specs{m}{:});
            plot_handles(end + 1) = h;
            legend_entries{end + 1} = method_names{m};
            current_data = [current_data, ydata];
        end

        setLogYTicks(gca, current_data, 'minmax_every3');
        set(gca, 'XMinorGrid', 'off', 'YMinorGrid', 'off', 'FontName', tennisFont);
        ylabel(ylabs{sp}, 'Interpreter', 'latex', 'FontName', tennisFont);
        if sp == 3
            xlabel('Time step [ms]', 'FontName', tennisFont);
            legend(plot_handles, legend_entries, 'Location', 'northwest', 'FontName', tennisFont);
        end
    end
end

%% RMS norm of angular velocity error vs time step 
figure;
set(gcf, 'Name', 'RMS angular velocity norm error');
hold on;
grid on;
set(gca, 'XScale', 'log', 'YScale', 'log', ...
    'XMinorGrid', 'off', 'YMinorGrid', 'off', ...
    'FontName', tennisFont);

plot_handles = gobjects(0);
legend_entries = {};
current_data = [];

for m = 1:num_methods
    valid_points = ~isnan(w_errors_rms_norm(m, :));
    if any(valid_points)
        ydata = w_errors_rms_norm(m, valid_points);
        h = loglog(time_steps(valid_points), ydata, line_specs{m}{:});
        plot_handles(end + 1) = h;
        legend_entries{end + 1} = method_names{m};
        current_data = [current_data, ydata];
    end
end

setLogYTicks(gca, current_data, 'minmax_every3');
set(gca, 'XMinorGrid', 'off', 'YMinorGrid', 'off', 'FontName', tennisFont);
xlabel('Time step [ms]', 'FontName', tennisFont);
ylabel('RMS angular velocity norm error [rad/s]', 'FontName', tennisFont);
legend(plot_handles, legend_entries, 'Location', 'northwest', 'FontName', tennisFont);

%% Computational time vs time step
figure;
hold on;
current_data = [];
plot_handles = gobjects(num_methods, 1);
for m = 1:num_methods
    valid_points = ~isnan(computation_times(m, :));
    if any(valid_points)
        ydata = computation_times(m, valid_points);
        plot_handles(m) = loglog(time_steps(valid_points), ydata, line_specs{m}{:});
        current_data = [current_data, ydata];
    end
end
grid on;
set(gca, 'XScale', 'log', 'YScale', 'log', ...
    'XMinorGrid', 'off', 'YMinorGrid', 'off', ...
    'FontName', tennisFont);
setLogYTicks(gca, current_data, 'auto');
set(gca, 'XMinorGrid', 'off', 'YMinorGrid', 'off', 'FontName', tennisFont);
xlabel('Time step [ms]', 'FontName', tennisFont);
ylabel('Computational time [s]', 'FontName', tennisFont);
valid_methods = any(~isnan(computation_times), 2);
legend(plot_handles(valid_methods), method_names(valid_methods), ...
    'Location', 'northwest', 'FontName', tennisFont);

%% RMS angular velocity norm error vs computational time 
figure;
hold on;
grid on;
set(gca, 'XScale', 'log', 'YScale', 'log', ...
    'XMinorGrid', 'off', 'YMinorGrid', 'off', ...
    'FontName', tennisFont);
for m = 1:num_methods
    valid_idx = ~isnan(w_errors_rms_norm(m, :)) & ~isnan(computation_times(m, :));
    if any(valid_idx)
        loglog(w_errors_rms_norm(m, valid_idx), computation_times(m, valid_idx), line_specs{m}{:});
    end
end
setLogYTicks(gca, computation_times(~isnan(computation_times)), 'auto');
set(gca, 'XMinorGrid', 'off', 'YMinorGrid', 'off', 'FontName', tennisFont);
xlabel('RMS angular velocity norm error [rad/s]', 'FontName', tennisFont);
ylabel('Computational time [s]', 'FontName', tennisFont);
legend(method_names, 'Location', 'northwest', 'FontName', tennisFont);

%% Energy magnitude at final time step vs time step
figure;
hold on;
current_data = [];
plot_handles = gobjects(num_methods, 1);
for m = 1:num_methods
    valid_points = ~isnan(energy_last_step(m, :)) & energy_last_step(m, :) ~= 0;
    if any(valid_points)
        ydata = abs(energy_last_step(m, valid_points));
        plot_handles(m) = loglog(time_steps(valid_points), ydata, line_specs{m}{:});
        current_data = [current_data, ydata];
    end
end
grid on;
set(gca, 'XScale', 'log', 'YScale', 'log', ...
    'XMinorGrid', 'off', 'YMinorGrid', 'off', ...
    'FontName', tennisFont);
setLogYTicks(gca, current_data, 'every3');
set(gca, 'XMinorGrid', 'off', 'YMinorGrid', 'off', 'FontName', tennisFont);
xlabel('Time step [ms]', 'FontName', tennisFont);
ylabel('Energy magnitude at final step [J]', 'FontName', tennisFont);
valid_methods = any(~isnan(energy_last_step), 2);
legend(plot_handles(valid_methods), method_names(valid_methods), ...
    'Location', 'northwest', 'FontName', tennisFont);

%% Angular momentum norm at final time step vs time step
figure;
hold on;
current_data = [];
plot_handles = gobjects(num_methods, 1);
for m = 1:num_methods
    valid_points = ~isnan(angular_momentum_last_step(m, :));
    if any(valid_points)
        ydata = abs(angular_momentum_last_step(m, valid_points));
        plot_handles(m) = loglog(time_steps(valid_points), ydata, line_specs{m}{:});
        current_data = [current_data, ydata];
    end
end
grid on;
set(gca, 'XScale', 'log', 'YScale', 'log', ...
    'XMinorGrid', 'off', 'YMinorGrid', 'off', ...
    'FontName', tennisFont);
setLogYTicks(gca, current_data, 'every3');
set(gca, 'XMinorGrid', 'off', 'YMinorGrid', 'off', 'FontName', tennisFont);
xlabel('Time step [ms]', 'FontName', tennisFont);
ylabel('Angular momentum norm at final step [J s]', 'FontName', tennisFont);
valid_methods = any(~isnan(angular_momentum_last_step), 2);
legend(plot_handles(valid_methods), method_names(valid_methods), ...
    'Location', 'northwest', 'FontName', tennisFont);

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
    'XMinorGrid', 'off', 'YMinorGrid', 'off', ...
    'FontName', tennisFont);
setLogYTicks(gca, current_data, 'every3');
set(gca, 'XMinorGrid', 'off', 'YMinorGrid', 'off', 'FontName', tennisFont);
xlabel('Time step [ms]', 'FontName', tennisFont);
ylabel('EP unit constraint deviation [-]', 'FontName', tennisFont);
valid_methods = isgraphics(plot_handles);
legend(plot_handles(valid_methods), ep_method_names(valid_methods), ...
    'Location', 'best', 'FontName', tennisFont);
