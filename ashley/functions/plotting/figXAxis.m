function figXAxis(fig,x_name,x_limit,varargin)

xlabel(x_name);
xlim(x_limit);
if length(varargin) > 2
    x_tick = varargin{3};
    if length(varargin) > 3
        x_tick_label = varargin{4};
        fig.XTickLabel = x_tick_label;
    else
        fig.XTickLabel = x_tick;
    end
end
end