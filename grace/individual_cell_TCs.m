%% neuropil subtraction
% day 1
height = 24; width = 24;
sz = size(data_reg);
down = 5;
np_tc = zeros(sz(3),nCells1);
data_tc_down = zeros(sz(3),nCells1);
np_tc_down = zeros(floor(sz(3)./down), nCells1);
cell_mask_np = zeros(25,25,nCells1);
for iCell = 1:nCells1
    xCenter = round(cell_stats(iCell).Centroid(2));
    yCenter = round(cell_stats(iCell).Centroid(1));
    xCenter2 = round(cell_stats2(iCell).Centroid(2));
    yCenter2 = round(cell_stats2(iCell).Centroid(1));
    xLeft = (xCenter - width/2);
    yBottom = (yCenter - height/2);
    xLeft2 = (xCenter2 - width/2);
    yBottom2 = (yCenter2 - height/2);
%  if xLeft > 12 && xLeft < 488 && yBottom > 12 && yBottom < 772
    cell_reg = data_reg(xLeft:(xLeft+width),yBottom:(height+yBottom),:);
%     cell_reg2 = data_reg(xLeft2:(xLeft2+width),yBottom2:(height+yBottom2),:);
%     [reg shift] = shift_opt(cell_reg2,cell_reg,4);
    cell_mask_cell = mask_cell(xLeft:(xLeft+width),yBottom:(height+yBottom));
    cell_mask_np(:,:,iCell) = mask_np(xLeft:(xLeft+width),yBottom:(height+yBottom),iCell);
    data_tc = stackGetTimeCourses(cell_reg);
    data_tc_down = stackGetTimeCourses(stackGroupProject(cell_reg,5));
    data_reg_down = stackGroupProject(cell_reg,down);
    np_tc(:,iCell) = stackGetTimeCourses(cell_reg,cell_mask_np(:,:,iCell));
    np_tc_down(:,iCell) = stackGetTimeCourses(data_reg_down,cell_mask_np(:,:,iCell));
    fprintf(['Cell #' num2str(iCell) '%s/n']) 
%  end
end
ii= 0.01:0.01:1;
x = zeros(length(ii), nCells1);
for i = 1:100
    x(i,:) = skewness(data_tc_down-tcRemoveDC(np_tc_down*ii(i)));
end
[max_skew ind] =  max(x,[],1);
np_w = 0.01*ind;
npSub_tc = data_tc-bsxfun(@times,tcRemoveDC(np_tc),np_w);
