function figXAxis(fig,x_name,x_limit,varargin)
if isempty(fig)
    fig = gca;
end
xlabel(x_name);
if ~isempty(x_limit)
    xlim(x_limit);
end
if length(varargin) > 0
    x_tick = varargin{1};
    fig.XTick = x_tick;
    if length(varargin) > 1
        x_tick_label = varargin{2};
        fig.XTickLabel = x_tick_label;
    else
        fig.XTickLabel = x_tick;
    end
end
end