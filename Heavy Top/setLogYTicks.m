function setLogYTicks(ax, data_matrix)
    % setLogYTicks - Sets logarithmic y-axis ticks with at least 3 ticks
    % 
    % Inputs:
    %   ax - axis handle (e.g., gca)
    %   data_matrix - matrix of data values to determine tick range
    %                 can be a vector or matrix (will find global min/max)
    
    % Handle different input types
    if isempty(data_matrix)
        return;
    end
    
    % Flatten matrix and remove NaN/zero values
    data_vector = data_matrix(:);
    data_vector = data_vector(~isnan(data_vector) & data_vector > 0);
    
    if isempty(data_vector)
        return;
    end
    
    % Find data range
    ymin = min(data_vector);
    ymax = max(data_vector);
    
    % Create tick range
    tick_range = floor(log10(ymin)):ceil(log10(ymax));
    
    % Ensure at least 3 ticks
    if length(tick_range) < 3
        tick_range = (ceil(log10(ymax))-2):ceil(log10(ymax));
    end
    
    % Set ticks
    ytick_values = 10.^tick_range;
    set(ax, 'YTick', ytick_values);
    set(ax, 'YTickLabelMode', 'auto');
    
    % Make sure log scale is set
    set(ax, 'YScale', 'log');
    
    % Optional: reduce grid density
    set(ax, 'GridAlpha', 0.3);
    set(ax, 'YMinorGrid', 'off');  % ✅ Correct property
set(ax, 'XMinorGrid', 'off');  % ✅ Correct property
end