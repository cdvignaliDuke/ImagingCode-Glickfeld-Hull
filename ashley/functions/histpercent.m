function [] = histpercent(data,M,color)
%    plot a histogram with M bins of vector data as percents rather than
%    counts (i.e. as hist(data,M) does

L = length(data);
[dataN databins] = hist(data,M);
datapercent = dataN/L;

bar(datapercent,color)
alpha(0.25)
xlim([1 M]);
% set(gca,'XTick',1:10:M)
% set(gca,'XTickLabel',databins(1:10:end))
end

