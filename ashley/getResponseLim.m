function responseLim = getResponseLim(varargin)

    nVars = length(varargin);
    mins = nan(1,nVars);
    maxes = nan(1,nVars);
    for i = 1:nVars
        mins(i) = min(varargin{i});
        maxes(i) = max(varargin{i});
    end
    mins = mins - (abs(mins)*.25);
    maxes = maxes + (abs(maxes)*.25);
    responseLim = [min(mins) max(maxes)];
end