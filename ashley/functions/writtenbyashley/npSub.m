function [npTC_weighted data_TCsubNP neuropil] = npSub(data_TC,mask_cell,data_reg,varargin)
%Get the timecourse of the neuropil, npTC, by finding the time-course of a disk
%around each ROI in mask_cell. This timecourse is weighted by the skewness
%of the fluorescence in the neuropil ROI to limit subtraction of real
%activity from the cell of interest. data_TCsubNP is the new neuropil subtracted 
%timecourse of all ROIs.
%Default buf = 4; np = 6;

nCells = size(data_TC,2);

if isempty(varargin)
    buf = 4;
    np = 6;
else
   buf = varargin{1};
   np = varargin{2};
end
    
neuropil = imCellNeuropil(mask_cell,buf,np);

npTC = zeros(size(data_TC));
for i = 1:size(data_TC,2)
    tempNPmask = squeeze(neuropil(:,:,i));
    if sum(sum(tempNPmask)) > 0
    npTC(:,i) = stackGetTimeCourses(data_reg,tempNPmask);
    end
end

%get weights by maximizing skew
ii= 0.01:0.01:1;
a = zeros(length(ii), nCells);
for i = 1:100
    a(i,:) = skewness(data_TC-tcRemoveDC(npTC*ii(i)));
end
[max_skew indA] =  max(a,[],1);
% skew(buf,:) = max_skew;
np_w = 0.01*indA;

npTC_weighted = bsxfun(@times,tcRemoveDC(npTC),np_w);
if nargout > 1
    data_TCsubNP = data_TC-npTC_weighted;
end
end