function figYAxis(fig,y_name,y_limit,varargin)
if isempty(fig)
    fig = gca;
end
ylabel(y_name);
if ~isempty(y_limit)
    ylim(y_limit);
end
if length(varargin) > 0
    y_tick = varargin{1};
    fig.YTick = y_tick;
    if length(varargin) > 1
        y_tick_label = varargin{2};
        fig.YTickLabel = y_tick_label;
    else
        fig.YTickLabel = y_tick;
    end
end
end