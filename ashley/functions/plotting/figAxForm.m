function figAxForm(varargin)
if isempty(varargin)
    fig = gca;
elseif isempty(varargin{1})
    fig = gca;
end
if length(varargin) <= 1 % default
    fig.TickDir = 'out';
    fig.Box = 'off';
    axis square
else
    doAxisSquare = varargin{2};
    if doAxisSquare
        fig.TickDir = 'out';
        fig.Box = 'off';
        axis square
    else  
        fig.TickDir = 'out';
        fig.Box = 'off';
    end
end
end