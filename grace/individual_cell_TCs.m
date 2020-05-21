%% neuropil subtraction
mouse = 'i1316';
ref_date = '200106';
day2 = '200108';
day3 = '200109';
ImgFolder = strvcat('003');
ImgFolder2 = strvcat('003');
nrun = size(ImgFolder,1);
nrun2 = size(ImgFolder2,1);
run_str = catRunName(ImgFolder, nrun);
run_str2 = catRunName(ImgFolder2, nrun2);
fnout = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\grace\Analysis\2P';

% loading data
maskD1 = load(fullfile(fnout, [ref_date '_' mouse], [ref_date '_' mouse '_' run_str], [ref_date '_' mouse '_' run_str '_mask_cell.mat']));
TCs_D1 = load(fullfile(fnout, [ref_date '_' mouse], [ref_date '_' mouse '_' ref_str], [ref_date '_' mouse '_' ref_str '_TCs.mat']));
npSub_tc1 = TCs_D1.npSub_tc;
nCells1 = size(npSub_tc1,2);
mask_cell = maskD1.mask_cell;
cell_list = intersect(1:nCells1, unique(mask_cell));
cell_stats = regionprops(mask_cell);

maskD2 = load(fullfile(fnout, [day2 '_' mouse], [day2 '_' mouse '_' run_str2], [day2 '_' mouse '_' run_str2 '_mask_cell.mat']));
mask_cell2 = maskD2.mask_cell;
cell_list2 = intersect(1:nCells1, unique(mask_cell2));
cell_stats2 = regionprops(mask_cell2);

maskD3 = load(fullfile(fnout, [day3 '_' mouse], [day3 '_' mouse '_' run_str], [day3 '_' mouse '_' run_str '_mask_cell.mat']));
mask_cell3 = maskD3.mask_cell;
cell_list3 = intersect(1:nCells1, unique(mask_cell3));
cell_stats3 = regionprops(mask_cell3);

%% day 2
height = 24; width = 24;
sz = size(data_reg);
down = 5;
np_tc = zeros(sz(3),nCells1);
data_tc = zeros(sz(3),nCells1);
data_tc_down = zeros(floor(sz(3)./down),nCells1);
np_tc_down = zeros(floor(sz(3)./down), nCells1);
cell_mask_np = zeros(25,25,nCells1);
npSub_tc = zeros(sz(3),nCells1);
np_w = zeros(1,nCells1);
ind = zeros(1,nCells1);
ii= 0.01:0.01:1;
x = zeros(length(ii), nCells1);
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
    cell_reg = double(cell_reg);
    cell_reg2 = data_reg2(xLeft2:(xLeft2+width),yBottom2:(height+yBottom2),:);
    cell_reg2 = double(cell_reg2);
    [reg shift] = shift_opt(cell_reg2,cell_reg,3);
    cell_mask_cell = mask_cell(xLeft:(xLeft+width),yBottom:(height+yBottom));
    cell_mask_np(:,:,iCell) = mask_np(xLeft:(xLeft+width),yBottom:(height+yBottom),iCell);
    data_tc(:,iCell) = stackGetTimeCourses(reg);
    data_tc_down(:,iCell) = stackGetTimeCourses(stackGroupProject(reg,5));
    data_reg_down = stackGroupProject(reg,down);
    np_tc(:,iCell) = stackGetTimeCourses(reg,cell_mask_np(:,:,iCell));
    np_tc_down(:,iCell) = stackGetTimeCourses(data_reg_down,cell_mask_np(:,:,iCell));
    fprintf(['Cell #' num2str(iCell) '%s/n']) 
%  end
    for i = 1:100
        x(i,iCell) = skewness(data_tc_down(:,iCell)-tcRemoveDC(np_tc_down(:,iCell)*ii(i)));
    end
    [max_skew ind(:,iCell)] =  max(x(:,iCell),[],1);
    np_w(:,iCell) = 0.01*ind(:,iCell);
    npSub_tc(:,iCell) = data_tc(:,iCell)-bsxfun(@times,tcRemoveDC(np_tc(:,iCell)),np_w(:,iCell));
end

save(fullfile(fnout, [day2 '_' mouse], [day2 '_' mouse '_' run_str], [day2 '_' mouse '_' run_str '_TCs.mat']), 'data_tc', 'np_tc', 'npSub_tc')
save(fullfile(fnout, [day2 '_' mouse], [day2 '_' mouse '_' run_str], [day2 '_' mouse '_' run_str '_input.mat']), 'input')

%% day 3
height = 24; width = 24;
sz = size(data_reg);
down = 5;
np_tc = zeros(sz(3),nCells1);
data_tc = zeros(sz(3),nCells1);
data_tc_down = zeros(floor(sz(3)./down),nCells1);
np_tc_down = zeros(floor(sz(3)./down), nCells1);
cell_mask_np = zeros(25,25,nCells1);
npSub_tc = zeros(sz(3),nCells1);
np_w = zeros(1,nCells1);
ind = zeros(1,nCells1);
ii= 0.01:0.01:1;
x = zeros(length(ii), nCells1);
for iCell = 1:nCells1
    xCenter = round(cell_stats(iCell).Centroid(2));
    yCenter = round(cell_stats(iCell).Centroid(1));
    xCenter3 = round(cell_stats3(iCell).Centroid(2));
    yCenter3 = round(cell_stats3(iCell).Centroid(1));
    xLeft = (xCenter - width/2);
    yBottom = (yCenter - height/2);
    xLef32 = (xCenter3 - width/2);
    yBottom3 = (yCenter3 - height/2);
%  if xLeft > 12 && xLeft < 488 && yBottom > 12 && yBottom < 772
    cell_reg = data_reg(xLeft:(xLeft+width),yBottom:(height+yBottom),:);
    cell_reg = double(cell_reg);
    cell_reg3 = data_reg3(xLeft3:(xLeft3+width),yBottom3:(height+yBottom3),:);
    cell_reg3 = double(cell_reg3);
    [reg shift] = shift_opt(cell_reg3,cell_reg,4);
    cell_mask_cell = mask_cell(xLeft:(xLeft+width),yBottom:(height+yBottom));
    cell_mask_np(:,:,iCell) = mask_np(xLeft:(xLeft+width),yBottom:(height+yBottom),iCell);
    data_tc(:,iCell) = stackGetTimeCourses(reg);
    data_tc_down(:,iCell) = stackGetTimeCourses(stackGroupProject(reg,5));
    data_reg_down = stackGroupProject(reg,down);
    np_tc(:,iCell) = stackGetTimeCourses(reg,cell_mask_np(:,:,iCell));
    np_tc_down(:,iCell) = stackGetTimeCourses(data_reg_down,cell_mask_np(:,:,iCell));
    fprintf(['Cell #' num2str(iCell) '%s/n']) 
%  end
    for i = 1:100
        x(i,iCell) = skewness(data_tc_down(:,iCell)-tcRemoveDC(np_tc_down(:,iCell)*ii(i)));
    end
    [max_skew ind(:,iCell)] =  max(x(:,iCell),[],1);
    np_w(:,iCell) = 0.01*ind(:,iCell);
    npSub_tc(:,iCell) = data_tc(:,iCell)-bsxfun(@times,tcRemoveDC(np_tc(:,iCell)),np_w(:,iCell));
end

save(fullfile(fnout, [day3 '_' mouse], [day3 '_' mouse '_' run_str], [day3 '_' mouse '_' run_str '_TCs.mat']), 'data_tc', 'np_tc', 'npSub_tc')
save(fullfile(fnout, [day3 '_' mouse], [day3 '_' mouse '_' run_str], [day3 '_' mouse '_' run_str '_input.mat']), 'input')
