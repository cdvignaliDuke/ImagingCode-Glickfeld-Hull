function [subplotRows, subplotColumns] = optimizeSubplotDim(nPlots)

subplotColumns = ceil(sqrt(nPlots));
subplotRows = ceil(nPlots/subplotColumns);

end