function [groupA, groupB] = twoRandGroups(allData)
% randomly segragate data into two groups along the first dimension
n = size(allData,1);

indA = randsample(n,round(n/2));
indB = setdiff(1:n,indA);

groupA = allData(indA,:);
groupB = allData(indB,:);
end