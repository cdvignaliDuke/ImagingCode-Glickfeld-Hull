function [m,err] = binnedMean(data,grpVar,binEdges,varargin)

if ~isempty(varargin)
    minN = varargin{1};
else
    minN = 2;
end

if size(data,1) > 1
    if size(data,2) > 1
        error('data must be a vector')
    elseif size(data,2) == 1
        data = data';
    else
        error('?')
    end
end

n = length(binEdges)-1;
grpID = discretize(grpVar,binEdges);

m = nan(1,n);
err = nan(1,n);
for i = 1:n
    ind = grpID == i;
    if sum(ind) >= minN
        m(i) = mean(data(ind),2);
        err(i) = ste(data(ind),2);
    end
end



end