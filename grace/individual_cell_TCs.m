%% get path names D1
date = '200118';
ImgFolder = strvcat('002');
time = strvcat('1508');
mouse = 'i1312';
nrun = size(ImgFolder,1);
frame_rate = 15.5;
run_str = catRunName(ImgFolder, nrun);
gl_fn = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lindsey\Data\2P_images';
fnout = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lindsey\Analysis\2P';
behav_fn = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\Behavior\Data';
%% load and register
data = [];
clear temp
trial_n = [];
offset = 0;
for irun = 1:nrun
    CD = fullfile(gl_fn, [mouse '\' date '\' ImgFolder(irun,:)]);
    cd(CD);
    imgMatFile = [ImgFolder(irun,:) '_000_000.mat'];
    load(imgMatFile);
    fName = fullfile(behav_fn, ['data-' mouse '-' date '-' time(irun,:) '.mat']);
    load(fName);

    nframes = info.config.frames;
    fprintf(['Reading run ' num2str(irun) '- ' num2str(nframes) ' frames \r\n'])
    data_temp = sbxread([ImgFolder(irun,:) '_000_000'],0,nframes);
    
    
    temp(irun) = input;
    if isfield(input, 'nScansOn')
        nOn = temp(irun).nScansOn;
        nOff = temp(irun).nScansOff;
        ntrials = size(temp(irun).tGratingDirectionDeg,2);

        data_temp = squeeze(data_temp);
        if nframes>ntrials*(nOn+nOff)
            data_temp = data_temp(:,:,1:ntrials*(nOn+nOff));
        elseif nframes<ntrials*(nOn+nOff)
            temp(irun) = trialChopper(temp(irun),1:ceil(nframes./(nOn+nOff)));
        end
    end
    
    offset = offset+nframes;

    data_temp = squeeze(data_temp);
    data = cat(3,data,data_temp);
    trial_n = [trial_n nframes];
end
input = concatenateDataBlocks(temp);
clear data_temp
clear temp

%% Choose register interval
t = 2000;
nep = floor(size(data,3)./t);
[n n2] = subplotn(nep);
figure; for i = 1:nep; subplot(n,n2,i); imagesc(mean(data(:,:,1+((i-1)*t):500+((i-1)*t)),3)); title([num2str(1+((i-1)*t)) '-' num2str(500+((i-1)*t))]); end

%% Register data

data_avg = mean(data(:,:,6001:6500),3);

if exist(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str]))
    load(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_reg_shifts.mat']))
    save(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_input.mat']), 'input')
    [outs, data_reg]=stackRegister_MA(data,[],[],out);
    clear out outs
else
    [out, data_reg] = stackRegister(data,data_avg);
    data_reg_avg = mean(data_reg(:,:,1:10000),3);
    reg = data_reg_avg;
    mkdir(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str]))
    save(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_reg_shifts.mat']), 'data_reg_avg', 'out', 'data_avg')
    save(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_input.mat']), 'input')
end

data_reg1 = data_reg;

clear data date ImgFolder data_reg data_avg data_reg_avg reg time

%% get path names D2
date = '200120';
ImgFolder = strvcat('003');
time = strvcat('1214');
nrun = size(ImgFolder,1);
frame_rate = 15.5;
run_str = catRunName(ImgFolder, nrun);
gl_fn = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lindsey\Data\2P_images';
fnout = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lindsey\Analysis\2P';
behav_fn = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\Behavior\Data';
%% load and register
data = [];
clear temp
trial_n = [];
offset = 0;
for irun = 1:nrun
    CD = fullfile(gl_fn, [mouse '\' date '\' ImgFolder(irun,:)]);
    cd(CD);
    imgMatFile = [ImgFolder(irun,:) '_000_000.mat'];
    load(imgMatFile);
    fName = fullfile(behav_fn, ['data-' mouse '-' date '-' time(irun,:) '.mat']);
    load(fName);

    nframes = info.config.frames;
    fprintf(['Reading run ' num2str(irun) '- ' num2str(nframes) ' frames \r\n'])
    data_temp = sbxread([ImgFolder(irun,:) '_000_000'],0,nframes);
    
    
    temp(irun) = input;
    if isfield(input, 'nScansOn')
        nOn = temp(irun).nScansOn;
        nOff = temp(irun).nScansOff;
        ntrials = size(temp(irun).tGratingDirectionDeg,2);

        data_temp = squeeze(data_temp);
        if nframes>ntrials*(nOn+nOff)
            data_temp = data_temp(:,:,1:ntrials*(nOn+nOff));
        elseif nframes<ntrials*(nOn+nOff)
            temp(irun) = trialChopper(temp(irun),1:ceil(nframes./(nOn+nOff)));
        end
    end
    
    offset = offset+nframes;

    data_temp = squeeze(data_temp);
    data = cat(3,data,data_temp);
    trial_n = [trial_n nframes];
end
input = concatenateDataBlocks(temp);
clear data_temp
clear temp

%% Choose register interval
t = 2000;
nep = floor(size(data,3)./t);
[n n2] = subplotn(nep);
figure; for i = 1:nep; subplot(n,n2,i); imagesc(mean(data(:,:,1+((i-1)*t):500+((i-1)*t)),3)); title([num2str(1+((i-1)*t)) '-' num2str(500+((i-1)*t))]); end

%% Register data

data_avg = mean(data(:,:,8001:8500),3);

if exist(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str]))
    load(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_reg_shifts.mat']))
    save(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_input.mat']), 'input')
    [outs, data_reg]=stackRegister_MA(data,[],[],out);
    clear out outs
else
    [out, data_reg] = stackRegister(data,data_avg);
    data_reg_avg2 = mean(data_reg(:,:,1:10000),3);
    reg2 = data_reg_avg2;
    mkdir(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str]))
    save(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_reg_shifts.mat']), 'data_reg_avg', 'out', 'data_avg')
    save(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_input.mat']), 'input')
end

data_reg2 = data_reg;

clear data date ImgFolder data_reg data_avg data_reg_avg2 reg time


%% Neuropil Subtraction 
ref_date = '200118';
day2 = '200120';
day3 = '200201';
mouse = 'i1312';
ImgFolder = strvcat('002');
ImgFolder2 = strvcat('003');
nrun = size(ImgFolder,1);
nrun2 = size(ImgFolder2,1);
ref_str = catRunName(ImgFolder, nrun);
run_str = catRunName(ImgFolder, nrun);
run_str2 = catRunName(ImgFolder2, nrun2);
fnout = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lindsey\Analysis\2P';

% loading data
maskD1 = load(fullfile(fnout, [ref_date '_' mouse], [ref_date '_' mouse '_' ref_str], [ref_date '_' mouse '_' ref_str '_mask_cell.mat']));
TCs_D1 = load(fullfile(fnout, [ref_date '_' mouse], [ref_date '_' mouse '_' ref_str], [ref_date '_' mouse '_' ref_str '_TCs.mat']));
pixD1 = load(fullfile(fnout, [ref_date '_' mouse], [ref_date '_' mouse '_' ref_str], [ref_date '_' mouse '_' ref_str '_pixel.mat']));
shiftsD1 = load(fullfile(fnout, [ref_date '_' mouse], [ref_date '_' mouse '_' ref_str], [ref_date '_' mouse '_' ref_str '_reg_shifts.mat']));
data_reg_avg = shiftsD1.data_reg_avg;
pixel1 = pixD1.pix;
npSub_tc1 = TCs_D1.npSub_tc;
nCells1 = size(npSub_tc1,2);
mask_cell = maskD1.mask_cell;
mask_np = maskD1.mask_np;
cell_list = intersect(1:nCells1, unique(mask_cell));
cell_stats = regionprops(mask_cell);

transD2 = load(fullfile(fnout, [day2 '_' mouse], [day2 '_' mouse '_' run_str2], [day2 '_' mouse '_' run_str2 '_transform.mat']));
fgta2 = transD2.fitGeoTAf;
% data_reg_avg with FGTA transformation
data_reg_avg2 = transD2.r2rFGTA;
% data_reg_avg without FGTA transformation
reg2 = transD2.reg;
pixD2 = load(fullfile(fnout, [day2 '_' mouse], [day2 '_' mouse '_' run_str2], [day2 '_' mouse '_' run_str2 '_pixel.mat']));
pix_fgta2 = pixD2.pix;
maskD2 = load(fullfile(fnout, [day2 '_' mouse], [day2 '_' mouse '_' run_str2], [day2 '_' mouse '_' run_str2 '_mask_cell.mat']));
mask_cell2 = maskD2.mask_cell;
mask_np2 = maskD2.mask_np;
cell_list2 = intersect(1:nCells1, unique(mask_cell2));
cell_stats2 = regionprops(mask_cell2);

% transD3 = load(fullfile(fnout, [day3 '_' mouse], [day3 '_' mouse '_' run_str], [day3 '_' mouse '_' run_str '_transform.mat']));
% fgta3 = transD3.fitGeoTAf;
% data_reg_avg3 = transD3.r2rFGTA;
% pixD3 = load(fullfile(fnout, [day3 '_' mouse], [day3 '_' mouse '_' run_str], [day3 '_' mouse '_' run_str '_pixel.mat']));
% pix_fgta3 = pixD3.pix_fgta;
% maskD3 = load(fullfile(fnout, [day3 '_' mouse], [day3 '_' mouse '_' run_str], [day3 '_' mouse '_' run_str '_mask_cell.mat']));
% mask_cell3 = maskD3.mask_cell;
% mask_np3 = maskD3.mask_np;
% cell_list3 = intersect(1:nCells1, unique(mask_cell3));
% cell_stats3 = regionprops(mask_cell3);

%% transform data_reg
sz = size(data_reg2);
data_reg22 = NaN(sz(1),sz(2),nframes);
for i = 1:nframes
    data_reg22(:,:,i) = imwarp(double(data_reg2(:,:,i)),fgta2, 'OutputView', imref2d(size(reg2)));
    if rem(i,50) == 0
        fprintf([num2str(i) '/n'])
    end
end
data_reg2_avg = mean(data_reg22,3);
figure;
subplot(1,2,1)
imagesc(data_reg2_avg)
title('avg of shifted data stack')
axis image
subplot(1,2,2)
imagesc(data_reg_avg2)
title('shifted avg of data stack')
axis image
print(fullfile(fnout, [ref_date '_' mouse], [ref_date '_' mouse '_' ref_str], ['AcrossAllDays'], ['CellMaps'], [ref_date '_' mouse '_regavg_vs_avgreg.pdf']),'-dpdf', '-bestfit')


for i = 1:nframes
    data_reg3(:,:,i) = imwarp(double(data_reg3(:,:,i)),fgta3, 'OutputView', imref2d(size(data_reg_avg3)));
    if rem(i,50) == 0
        fprintf([num2str(i) '/n'])
    end
end

%% day 2
height = 30; width = 30;
cells = NaN(nCells1,1);
for iCell = 1:nCells1
    xCenter = round(cell_stats(iCell).Centroid(2));
    yCenter = round(cell_stats(iCell).Centroid(1));
    xCenter2 = round(cell_stats2(iCell).Centroid(2));
    yCenter2 = round(cell_stats2(iCell).Centroid(1));
    xLeft = (xCenter - width/2);
    yBottom = (yCenter - height/2);
    xLeft2 = (xCenter2 - width/2);
    yBottom2 = (yCenter2 - height/2);
 if xLeft > 30 && xLeft < 482 && yBottom > 30 && yBottom < 766
 cells(iCell) = iCell;
 else 
 cells(iCell) = 0;
 end
end
cells(cells==0)=[];

%% day2 
height = 30;width = 30;
sz = size(data_reg22);
nframes = size(data_reg22,3);
down = 5;
data_tc = zeros(length(cells),sz(3));
data_tc_down = zeros(length(cells),floor(sz(3)./down));
np_tc = zeros(length(cells),sz(3));
np_tc_down = zeros(length(cells),floor(sz(3)./down));
% cell_mask_np = NaN(51,51,nCells1);
% np_w = zeros(nCells1,1);
% ind = zeros(nCells1,1);
ii= 0.01:0.01:1;
r = zeros(length(cells),1);
p = zeros(length(cells),1);
npSub_tc = zeros(length(cells),sz(3));
for iC = 1:15
    iCell = cells(iC);
    xCenter = round(cell_stats(iCell).Centroid(2));
    yCenter = round(cell_stats(iCell).Centroid(1));
    xCenter2 = round(cell_stats2(iCell).Centroid(2));
    yCenter2 = round(cell_stats2(iCell).Centroid(1));
    xLeft = (xCenter - width/2);
    yBottom = (yCenter - height/2);
    xLeft2 = (xCenter2 - width/2);
    yBottom2 = (yCenter2 - height/2);
    cell_reg_avg = data_reg_avg(xLeft:(xLeft+width),yBottom:(height+yBottom));
    pix1 = pixel1(xLeft:(xLeft+width),yBottom:(height+yBottom));
    cell_reg2 = data_reg22(xLeft2:(xLeft2+width),yBottom2:(height+yBottom2),:);
    cell_reg2 = double(cell_reg2);
    cell_reg2_avg = data_reg_avg2(xLeft2:(xLeft2+width),yBottom2:(height+yBottom2));
    pix2 = pix_fgta2(xLeft2:(xLeft2+width),yBottom2:(height+yBottom2));
    [reg shift] = shift_opt(cell_reg2_avg,cell_reg_avg,4);
    [reg2 shift2] = shift_opt(pix2,pix1,4);
    r(iCell) = triu2vec(corrcoef(cell_reg_avg(:),reg(:)));
    p(iCell) = triu2vec(corrcoef(pix1(:),reg2(:)));
    r(isnan(r))=0;
    p(isnan(p))=0;
%     fine shift to course-shifted cell squares
if r(iCell)>0.8 && p(iCell)>0.4
    [outs, reg_cell2] = stackRegister_MA(cell_reg2,[],[],repmat(shift,[nframes 1]));
else
    [outs, reg_cell2] = stackRegister_MA(cell_reg2,[],[],repmat(shift2,[nframes 1]));
end
    cell_mask = mask_cell2(xLeft2:(xLeft2+width),yBottom2:(height+yBottom2));
    cell_mask_np = mask_np2(xLeft2:(xLeft2+width),yBottom2:(height+yBottom2),iCell);
    cell_reg2_long = reshape(reg_cell2, [31.*31 nframes]);
    cell_reg2_down  = reshape(stackGroupProject(reg_cell2,5), [31.*31 nframes/5]);
    cell_mask_long = reshape(cell_mask, [31.*31 1]);
    mask_np_long = reshape(cell_mask_np, [31.*31 1]);
    data_tc(iCell,:) = mean(cell_reg2_long(find(cell_mask_long==iCell),:),1);
    data_tc_down(iCell,:) = mean(cell_reg2_down(find(cell_mask_long==iCell),:),1);
    np_tc(iCell,:) = mean(cell_reg2_long(find(mask_np_long),:),1);
    np_tc_down(iCell,:) = mean(cell_reg2_down(find(mask_np_long),:),1);
    fprintf(['Cell #' num2str(iCell) '%s/n']) 
    for i = 1:100
        x = skewness(data_tc_down(iCell,:)-tcRemoveDC(np_tc_down(iCell,:)*ii(i)));
    end
    [max_skew ind] =  max(x,[],2);
    np_w = 0.01*ind;
    npSub_tc(iCell,:) = data_tc(iCell,:)-bsxfun(@times,tcRemoveDC(np_tc(iCell,:)),np_w);
end
cells2 = find(r>0.8&p>.4);
cells3 = find(r>0.4&p>.6);
cells4 = [cells2;cells3];
cells_all = reshape(cells4,[],1);
npSub_tc = npSub_tc';
% save(fullfile(fnout, [day2 '_' mouse], [day2 '_' mouse '_' run_str2], [day2 '_' mouse '_' run_str2 '_TCs.mat']), 'data_tc', 'np_tc', 'npSub_tc')
% save(fullfile(fnout, [day2 '_' mouse], [day2 '_' mouse '_' run_str2], [day2 '_' mouse '_' run_str2 '_input.mat']), 'input')

%% mask images
start = 1;
figure;
for iC = 1:10
    iCell = cells(iC);
    xCenter = round(cell_stats(iCell).Centroid(2));
    yCenter = round(cell_stats(iCell).Centroid(1));
    xCenter2 = round(cell_stats2(iCell).Centroid(2));
    yCenter2 = round(cell_stats2(iCell).Centroid(1));
    xLeft = (xCenter - width/2);
    yBottom = (yCenter - height/2);
    xLeft2 = (xCenter2 - width/2);
    yBottom2 = (yCenter2 - height/2);
    cell_reg_avg = data_reg_avg(xLeft:(xLeft+width),yBottom:(height+yBottom));
    pix1 = pixel1(xLeft:(xLeft+width),yBottom:(height+yBottom));
    cell_reg2 = data_reg22(xLeft2:(xLeft2+width),yBottom2:(height+yBottom2),:);
    cell_reg2 = double(cell_reg2);
    cell_reg2_avg = data_reg_avg2(xLeft2:(xLeft2+width),yBottom2:(height+yBottom2));
    pix2 = pix_fgta2(xLeft2:(xLeft2+width),yBottom2:(height+yBottom2));
    [reg shift] = shift_opt(cell_reg2_avg,cell_reg_avg,4);
    [reg2 shift2] = shift_opt(pix2,pix1,4);
    r = triu2vec(corrcoef(cell_reg_avg(:),reg(:)));
    p = triu2vec(corrcoef(pix1(:),reg2(:)));
    r(isnan(r))=0;
    p(isnan(p))=0;
%     fine shift to course-shifted cell squares
if r>0.8 && p>0.4
    [outs, reg_cell2] = stackRegister_MA(cell_reg2,[],[],repmat(shift,[nframes 1]));
elseif p>0.6 && r>0.4 && r<0.8
    [outs, reg_cell2] = stackRegister_MA(cell_reg2,[],[],repmat(shift2,[nframes 1]));
else
    reg_cell2 = NaN(51,51,nframes);
end
reg_cell2_avg = mean(reg_cell2,3);
cell_mask = mask_cell2(xLeft2:(xLeft2+width),yBottom2:(height+yBottom2));
subplot(10,8,start)
imagesc(reg_cell2_avg)
pos = get(gca, 'Position');
pos(1) = 0.1;
pos(3) = 0.05;
set(gca, 'Position', pos)
axis square
axis off
hold on
bound = cell2mat(bwboundaries(cell_mask(:,:,1)));
plot(bound(:,2),bound(:,1),'.','color','r','MarkerSize',0.5);
subplot(10,8,start+1)
imagesc(cell_mask)
pos = get(gca, 'Position');
pos(1) = 0.15;
pos(3) = 0.05;
set(gca, 'Position', pos)
axis square
axis off
subplot(10,8,start+2)
imagesc(reg)
pos = get(gca, 'Position');
pos(1) = 0.2;
pos(3) = 0.05;
set(gca, 'Position', pos)
axis square
axis off
hold on
bound = cell2mat(bwboundaries(cell_mask(:,:,1)));
plot(bound(:,2),bound(:,1),'.','color','r','MarkerSize',0.5);
subplot(10,8,start+3)
imagesc(reg2)
pos = get(gca, 'Position');
pos(1) = 0.25;
pos(3) = 0.05;
set(gca, 'Position', pos)
axis square
axis off
hold on
bound = cell2mat(bwboundaries(cell_mask(:,:,1)));
plot(bound(:,2),bound(:,1),'.','color','r','MarkerSize',0.5);

    xCenter = round(cell_stats(iCell+10).Centroid(2));
    yCenter = round(cell_stats(iCell+10).Centroid(1));
    xCenter2 = round(cell_stats2(iCell+10).Centroid(2));
    yCenter2 = round(cell_stats2(iCell+10).Centroid(1));
    xLeft = (xCenter - width/2);
    yBottom = (yCenter - height/2);
    xLeft2 = (xCenter2 - width/2);
    yBottom2 = (yCenter2 - height/2);
    cell_reg_avg = data_reg_avg(xLeft:(xLeft+width),yBottom:(height+yBottom));
    pix1 = pixel1(xLeft:(xLeft+width),yBottom:(height+yBottom));
    cell_reg2 = data_reg22(xLeft2:(xLeft2+width),yBottom2:(height+yBottom2),:);
    cell_reg2 = double(cell_reg2);
    cell_reg2_avg = data_reg_avg2(xLeft2:(xLeft2+width),yBottom2:(height+yBottom2));
    pix2 = pix_fgta2(xLeft2:(xLeft2+width),yBottom2:(height+yBottom2));
    [reg shift] = shift_opt(cell_reg2_avg,cell_reg_avg,4);
    [reg2 shift2] = shift_opt(pix2,pix1,4);
    r = triu2vec(corrcoef(cell_reg_avg(:),reg(:)));
    p = triu2vec(corrcoef(pix1(:),reg2(:)));
    r(isnan(r))=0;
    p(isnan(p))=0;
%     fine shift to course-shifted cell squares
if r>0.8 && p>0.4
    [outs, reg_cell2] = stackRegister_MA(cell_reg2,[],[],repmat(shift,[nframes 1]));
elseif p>0.6 && r>0.4 && r<0.8
    [outs, reg_cell2] = stackRegister_MA(cell_reg2,[],[],repmat(shift2,[nframes 1]));
else
    reg_cell2 = NaN(51,51,nframes);
end
reg_cell2_avg = mean(reg_cell2,3);
cell_mask = mask_cell2(xLeft2:(xLeft2+width),yBottom2:(height+yBottom2));
subplot(10,8,start+4)
imagesc(reg_cell2_avg)
pos = get(gca, 'Position');
pos(1) = 0.4;
pos(3) = 0.05;
set(gca, 'Position', pos)
axis square
axis off
hold on
bound = cell2mat(bwboundaries(cell_mask(:,:,1)));
plot(bound(:,2),bound(:,1),'.','color','r','MarkerSize',0.5);
subplot(10,8,start+5)
imagesc(cell_mask)
pos = get(gca, 'Position');
pos(1) = 0.45;
pos(3) = 0.05;
set(gca, 'Position', pos)
axis square
axis off
subplot(10,8,start+6)
imagesc(reg)
pos = get(gca, 'Position');
pos(1) = 0.5;
pos(3) = 0.05;
set(gca, 'Position', pos)
axis square
axis off
hold on
bound = cell2mat(bwboundaries(cell_mask(:,:,1)));
plot(bound(:,2),bound(:,1),'.','color','r','MarkerSize',0.5);
subplot(10,8,start+7)
imagesc(reg2)
pos = get(gca, 'Position');
pos(1) = 0.55;
pos(3) = 0.05;
set(gca, 'Position', pos)
axis square
axis off
hold on
bound = cell2mat(bwboundaries(cell_mask(:,:,1)));
plot(bound(:,2),bound(:,1),'.','color','r','MarkerSize',0.5);

start = start+8;
end
% print(fullfile(fnout, [ref_date '_' mouse], [ref_date '_' mouse '_' ref_str], ['AcrossAllDays'], ['CellMaps'], [ref_date '_' mouse '_mask_outline.pdf']),'-dpdf', '-bestfit')


%% day 3
height = 50; width = 50;
sz = size(data_reg);
down = 5;
np_tc = NaN(sz(3),nCells1);
data_tc = NaN(sz(3),nCells1);
data_tc_down = NaN(floor(sz(3)./down),nCells1);
np_tc_down = NaN(floor(sz(3)./down), nCells1);
cell_mask_np = NaN(51,51,nCells1);
npSub_tc = NaN(sz(3),nCells1);
np_w = NaN(1,nCells1);
ind = NaN(1,nCells1);
ii= 0.01:0.01:1;
x = zeros(length(ii), nCells1);

for iCell = 1:nCells1
    xCenter = round(cell_stats(iCell).Centroid(2));
    yCenter = round(cell_stats(iCell).Centroid(1));
    xCenter2 = round(cell_stats3(iCell).Centroid(2));
    yCenter2 = round(cell_stats3(iCell).Centroid(1));
    xLeft = (xCenter - width/2);
    yBottom = (yCenter - height/2);
    xLeft3 = (xCenter3 - width/2);
    yBottom3 = (yCenter3 - height/2);
    cell_reg = data_reg1(xLeft:(xLeft+width),yBottom:(height+yBottom),:);
    cell_reg = double(cell_reg);
    cell_reg_avg = data_reg_avg(xLeft:(xLeft+width),yBottom:(height+yBottom));
    pix1 = pixel1(xLeft:(xLeft+width),yBottom:(height+yBottom));
    cell_reg3 = data_reg2(xLeft3:(xLeft3+width),yBottom3:(height+yBottom3),:);
    cell_reg3 = double(cell_reg3);
    cell_reg3_avg = data_reg_avg2(xLeft3:(xLeft3+width),yBottom3:(height+yBottom3));
    pix3 = pix_fgta2(xLeft3:(xLeft3+width),yBottom3:(height+yBottom3));
    [reg shift] = shift_opt(cell_reg3_avg,cell_reg_avg,4);
    [reg2 shift2] = shift_opt(pix3,pix1,4);
    r = triu2vec(corrcoef(cell_reg_avg(:),reg(:)));
    p = triu2vec(corrcoef(pix1(:),reg2(:)));
    r(isnan(r))=0;
    p(isnan(p))=0;
%     fine shift to course-shifted cell squares
if r>0.6 && p>0.4
    [outs, reg_cell3] = stackRegister_MA(cell_reg3,[],[],repmat(shift,[nframes 1]));
elseif p>0.8 && r>0.4 && r<0.6
    [outs, reg_cell3] = stackRegister_MA(cell_reg3,[],[],repmat(shift2,[nframes 1]));
end
    cell_reg3_long = reshape(reg_cell3, [51.*51 nframes]);
    cell_reg3_down  = reshape(stackGroupProject(reg_cell3,5), [51.*51 nframes/5]);
    cell_mask = mask_cell3(xLeft3:(xLeft3+width),yBottom3:(height+yBottom3));
    cell_mask_long = reshape(cell_mask, [51.*51 1]);
    cell_mask_np(:,:,iCell) = mask_np3(xLeft3:(xLeft3+width),yBottom3:(height+yBottom3),iCell);
    mask_np_long = reshape(cell_mask_np(:,:,iCell), [51.*51 1]);
    data_tc(:,iCell) = mean(cell_reg3_long(find(cell_mask_long==iCell),:),1);
    data_tc_down(:,iCell) = mean(cell_reg3_down(find(cell_mask_long==iCell),:),1);
    np_tc(:,iCell) = mean(cell_reg3_long(find(mask_np_long),:),1);
    np_tc_down(:,iCell) = mean(cell_reg3_down(find(mask_np_long),:),1);
    fprintf(['Cell #' num2str(iCell) '%s/n']) 
    for i = 1:100
        x(i,iCell) = skewness(data_tc_down(:,iCell)-tcRemoveDC(np_tc_down(:,iCell)*ii(i)));
    end
    [max_skew ind(:,iCell)] =  max(x(:,iCell),[],1);
    np_w(:,iCell) = 0.01*ind(:,iCell);
    npSub_tc(:,iCell) = data_tc(:,iCell)-bsxfun(@times,tcRemoveDC(np_tc(:,iCell)),np_w(:,iCell));
end
save(fullfile(fnout, [day3 '_' mouse], [day3 '_' mouse '_' run_str], [day3 '_' mouse '_' run_str '_TCs.mat']), 'data_tc', 'np_tc', 'npSub_tc')
save(fullfile(fnout, [day3 '_' mouse], [day3 '_' mouse '_' run_str], [day3 '_' mouse '_' run_str '_input.mat']), 'input')
