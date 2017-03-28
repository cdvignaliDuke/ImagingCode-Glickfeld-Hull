function [data_tc_subnnp, np_tc, np_w] = getWeightedNeuropilTimeCourse(data_reg, data_tc, mask_cell,buf,np)

% get neuropil timecourses
nCells = size(data_tc,2);

neuropil = imCellNeuropil(mask_cell,buf,np);
np_tc = zeros(size(data_tc));
for i = 1:nCells
    tempNPmask = squeeze(neuropil(:,:,i));
    if sum(sum(tempNPmask)) > 0
    np_tc(:,i) = stackGetTimeCourses(data_reg,tempNPmask);
    end
end

%get weights by maximizing skew
ii= 0.01:0.01:1;
x = zeros(length(ii), nCells);
for i = 1:100
    x(i,:) = skewness(data_tc-tcRemoveDC(np_tc*ii(i)));
end
[max_skew, ind] =  max(x,[],1);
np_w = 0.01*ind;
data_tc_subnnp = data_tc-bsxfun(@times,tcRemoveDC(np_tc),np_w);

end