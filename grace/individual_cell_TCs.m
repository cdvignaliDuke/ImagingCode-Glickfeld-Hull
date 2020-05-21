%% neuropil subtraction
mouse = 'i1316';
day2 = '200108';
day3 = '200109';
ImgFolder = strvcat('003');
ImgFolder2 = strvcat('003');
ref_date = '200106';
ref_run = strvcat('003');
nrun = size(ImgFolder,1);
nrun2 = size(ImgFolder2,1);
frame_rate = 15.5;
run_str = catRunName(ImgFolder, nrun);
run_str2 = catRunName(ImgFolder2, nrun2);
run_str3 = catRunName(ImgFolder, nrun);
ref_str = catRunName(ref_run, size(ref_run,1));
fnout = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\grace\Analysis\2P';

% loading data
maskD1 = load(fullfile(fnout, [ref_date '_' mouse], [ref_date '_' mouse '_' ref_str], [ref_date '_' mouse '_' ref_str '_mask_cell.mat']));
mask_cell = maskD1.mask_cell;
dfof = maskD1.data_dfof_max;
reg_shifts = load(fullfile(fnout, [ref_date '_' mouse], [ref_date '_' mouse '_' ref_str], [ref_date '_' mouse '_' ref_str '_reg_shifts.mat']));
reg1 = reg_shifts.data_reg_avg;
reg1(find(reg1>7000)) = 0;
reg1 = (reg1./max(max(abs(reg1))));
cell_list = intersect(1:nCells1, unique(mask_cell));
cell_stats = regionprops(mask_cell);

maskD2 = load(fullfile(fnout, [day2 '_' mouse], [day2 '_' mouse '_' run_str2], [day2 '_' mouse '_' run_str2 '_mask_cell.mat']));
mask_cell2 = maskD2.mask_cell;
reg_shifts2 = load(fullfile(fnout, [day2 '_' mouse], [day2 '_' mouse '_' run_str2], [day2 '_' mouse '_' run_str2 '_transform.mat']));
dfof2 = reg_shifts2.r2rFGTA_dfof;
reg2 = reg_shifts2.r2rFGTA;
reg2(find(reg2>7000)) = 0;
reg2 = (reg2./max(max(abs(reg2))));
cell_list2 = intersect(1:nCells1, unique(mask_cell2));
cell_stats2 = regionprops(mask_cell2);

maskD3 = load(fullfile(fnout, [day3 '_' mouse], [day3 '_' mouse '_' run_str], [day3 '_' mouse '_' run_str '_mask_cell.mat']));
mask_cell3 = maskD3.mask_cell;
reg_shifts3 = load(fullfile(fnout, [day3 '_' mouse], [day3 '_' mouse '_' run_str], [day3 '_' mouse '_' run_str '_transform.mat']));
dfof3 = reg_shifts3.r2rFGTA_dfof;
reg3 = reg_shifts3.r2rFGTA;
reg3(find(reg3>7000)) = 0;
reg3 = (reg3./max(max(abs(reg3))));
cell_list3 = intersect(1:nCells1, unique(mask_cell3));
cell_stats3 = regionprops(mask_cell3);

%% day 1
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
%     cell_reg2 = data_reg(xLeft2:(xLeft2+width),yBottom2:(height+yBottom2),:);
%     [cell_reg shift] = shift_opt(cell_reg2,cell_reg,4);
    cell_mask_cell = mask_cell(xLeft:(xLeft+width),yBottom:(height+yBottom));
    cell_mask_np(:,:,iCell) = mask_np(xLeft:(xLeft+width),yBottom:(height+yBottom),iCell);
    data_tc(:,iCell) = stackGetTimeCourses(cell_reg);
    data_tc_down(:,iCell) = stackGetTimeCourses(stackGroupProject(cell_reg,5));
    data_reg_down = stackGroupProject(cell_reg,down);
    np_tc(:,iCell) = stackGetTimeCourses(cell_reg,cell_mask_np(:,:,iCell));
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

save(fullfile(fnout, [ref_date '_' mouse], [ref_date '_' mouse '_' ref_str], [ref_date '_' mouse '_' ref_str '_TCs.mat']), 'data_tc', 'np_tc', 'npSub_tc')
save(fullfile(fnout, [ref_date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_input.mat']), 'input')
