function figYAxis(fig,y_name,y_limit,varargin)

ylabel(y_name);
ylim(y_limit);
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