function createComparisonPlot(error_data_x, error_data_y, error_data_z, time_steps, method_names, line_specs, plot_title, error_type)
    % createComparisonPlot - Creates a 3-subplot comparison plot for position/velocity/acceleration errors
    %
    % Inputs:
    %   error_data_x, error_data_y, error_data_z - Error matrices [methods x time_steps] for each component
    %   time_steps - Vector of time step values
    %   method_names - Cell array of method names for legend
    %   line_specs - Cell array of line specifications for each method
    %   plot_title - Main title for the figure
    %   error_type - String describing error type (e.g., 'RMS Error', 'MAX ABS Error')
    
    figure;
    titles = {'Position Error', 'Velocity Error', 'Acceleration Error'};
    ylabels = {'(m)', '(m/s)', '(m/s²)'};
    
    % Determine which error matrices to use for each subplot
    error_matrices = {error_data_x, error_data_y, error_data_z};
    
    for sp = 1:3
        subplot(3,1,sp);
        hold on;
        grid on;
        set(gca, 'XScale', 'log', 'YScale', 'log');
        
        plot_handles = gobjects(0);   % Initialize empty array for plot handles
        legend_entries = {};          % Initialize cell array for legend entries
        current_data = [];            % Collect data for tick calculation
        
        % Plot each method's data
    for m = 1:numel(method_names)
        % Find valid time-step indices for this method
        idx = find(~isnan(error_matrices{sp}(m,:)));
        % Only keep indices that actually exist in time_steps
idx = idx(idx <= numel(time_steps));
        if isempty(idx)
            continue;
        end
        % Extract and plot
        ydata = error_matrices{sp}(m, idx);
        h = loglog(time_steps(idx), ydata, line_specs{m}{:});
        plot_handles(end+1) = h;
        legend_entries{end+1} = method_names{m};
        current_data = [current_data, ydata];
    end
        
        % Apply custom tick formatting
        setLogYTicks(gca, current_data);
        
        % Set titles and labels
        title(titles{sp});
        ylabel(ylabels{sp});
        if sp == 3
            xlabel('Time Step (ms)');
            % Add legend only once after plotting all subplots
            legend(plot_handles, legend_entries, 'Location', 'northwest');
        end
    end
    
   % sgtitle([error_type ' - ' plot_title]);
end