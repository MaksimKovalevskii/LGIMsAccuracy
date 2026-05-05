function setLogYTicks(ax, data_matrix, labelMode)
    % setLogYTicks - Sets logarithmic y-axis ticks and optional sparse labels
    %
    % Inputs:
    %   ax - axis handle (e.g., gca)
    %   data_matrix - matrix of data values to determine tick range
    %                 can be a vector or matrix (will find global min/max)
    %   labelMode (optional) - controls YTickLabel sparsity while keeping
    %                          decade ticks for grid:
    %                          'auto' (default): label all decades
    %                          'minmidmax': label only min/mid/max decade
    %                          'every3': label every 3 decades (…, -6, -3, 0, …)
    %                          'minmax_every3': label min/max and every 3 decades in between

    if nargin < 3 || isempty(labelMode)
        labelMode = 'auto';
    end

    if isempty(data_matrix)
        return;
    end

    data_vector = data_matrix(:);
    data_vector = data_vector(~isnan(data_vector) & data_vector > 0);

    if isempty(data_vector)
        return;
    end

    ymin = min(data_vector);
    ymax = max(data_vector);

    minPow = floor(log10(ymin));
    maxPow = ceil(log10(ymax));
    tick_range = minPow:maxPow;

    if length(tick_range) < 3
        tick_range = (ceil(log10(ymax)) - 2):ceil(log10(ymax));
    end

    ytick_values = 10.^tick_range;
    set(ax, 'YTick', ytick_values);
    set(ax, 'TickLabelInterpreter', 'tex');

    switch lower(string(labelMode))
        case "auto"
            set(ax, 'YTickLabelMode', 'auto');
        case "every3"
            showPow = tick_range(mod(tick_range, 3) == 0);
            labels = repmat({''}, size(tick_range));
            for i = 1:numel(tick_range)
                if any(tick_range(i) == showPow)
                    labels{i} = sprintf('10^{%d}', tick_range(i));
                end
            end
            set(ax, 'YTickLabel', labels);
        case "minmidmax"
            if isempty(tick_range)
                return;
            end
            midPow = round((minPow + maxPow) / 2);
            showPow = unique([minPow, midPow, maxPow], 'stable');
            labels = repmat({''}, size(tick_range));
            for i = 1:numel(tick_range)
                if any(tick_range(i) == showPow)
                    labels{i} = sprintf('10^{%d}', tick_range(i));
                end
            end
            set(ax, 'YTickLabel', labels);
        case "minmax_every3"
            showPow = unique([minPow:3:maxPow, maxPow], 'stable');
            labels = repmat({''}, size(tick_range));
            for i = 1:numel(tick_range)
                if any(tick_range(i) == showPow)
                    labels{i} = sprintf('10^{%d}', tick_range(i));
                end
            end
            set(ax, 'YTickLabel', labels);
        otherwise
            set(ax, 'YTickLabelMode', 'auto');
    end

    set(ax, 'YScale', 'log');
    lm = lower(char(labelMode));
    if ~strcmp(lm, 'auto')
        set(ax, 'YLim', [10^minPow, 10^maxPow]);
    end

    set(ax, 'GridAlpha', 0.3);
    set(ax, 'YMinorGrid', 'off', 'XMinorGrid', 'off');
end
