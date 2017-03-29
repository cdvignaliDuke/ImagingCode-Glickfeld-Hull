function figYAxis(fig,y_name,y_limit,varargin)

ylabel(y_name);
ylim(y_limit);
if length(varargin) > 2
    y_tick = varargin{3};
    fig.YTick = y_tick;
    if length(varargin) > 3
        y_tick_label = varargin{4};
        fig.YTickLabel = y_tick_label;
    else
        fig.YTickLabel = y_tick;
    end
end
end