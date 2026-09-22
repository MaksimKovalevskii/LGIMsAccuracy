function createNormComparisonPlot(error_data_norm, time_steps, method_names, line_specs, plot_title, error_type)
    % createNormComparisonPlot - Creates a 3-subplot comparison plot for norm (magnitude) errors
    %
    % Inputs:
    %   error_data_norm - Cell array with 3 error matrices {pos_norm, vel_norm, acc_norm} [methods x time_steps]
    %   time_steps - Vector of time step values
    %   method_names - Cell array of method names for legend
    %   line_specs - Cell array of line specifications for each method
    %   plot_title - Descriptor (e.g. 'The norm (magnitude)') for documentation; norm is spelled out in axis labels
    %   error_type - String describing error type (e.g. 'RMS Error', 'MAX ABS Error')

    figure;
    set(gcf, 'Name', sprintf('%s - %s', error_type, plot_title));
    qty_label = {'Position', 'Velocity', 'Acceleration'};
    units = {'m', 'm/s', 'm/s^2'};
    if strcmp(error_type, 'RMS Error')
        err_label = 'RMS error';
    elseif strcmp(error_type, 'MAX ABS Error')
        err_label = 'Maximum absolute error';
    else
        err_label = error_type;
    end

    % Shared (unified) error label for all 3 subplots
    annotation(gcf, 'textbox', [0.005, 0.5, 0.02, 0.1], ...
        'String', err_label, 'EdgeColor', 'none', 'Rotation', 90, ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
        'FontName', 'Times New Roman');

    for sp = 1:3
        subplot(3, 1, sp);
        hold on;
        grid on;
        set(gca, 'XScale', 'log', 'YScale', 'log', ...
            'XMinorGrid', 'off', 'YMinorGrid', 'off', ...
            'FontName', 'Times New Roman');

        plot_handles = gobjects(0);
        legend_entries = {};
        current_data = [];

        for m = 1:size(error_data_norm{sp}, 1)
            valid_points = ~isnan(error_data_norm{sp}(m, :));
            if any(valid_points)
                data_to_plot = error_data_norm{sp}(m, valid_points);
                h = loglog(time_steps(valid_points), data_to_plot, line_specs{m}{:});
                plot_handles(end + 1) = h;
                legend_entries{end + 1} = method_names{m};
                current_data = [current_data, data_to_plot];
            end
        end

        setLogYTicks(gca, current_data, 'minmax_every3');
        set(gca, 'XMinorGrid', 'off', 'YMinorGrid', 'off', 'FontName', 'Times New Roman');

        ylabel(sprintf('%s norm [%s]', qty_label{sp}, units{sp}), 'FontName', 'Times New Roman');
        if sp == 3
            xlabel('Time step [ms]', 'FontName', 'Times New Roman');
            legend(plot_handles, legend_entries, 'Location', 'northwest', 'FontName', 'Times New Roman');
        end
    end
end
