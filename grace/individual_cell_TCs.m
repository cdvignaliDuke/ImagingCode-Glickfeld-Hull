%% neuropil subtraction - load data_reg1 and data_reg2/3 first!!
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
maskD1 = load(fullfile(fnout, [ref_date '_' mouse], [ref_date '_' mouse '_' ref_str], [ref_date '_' mouse '_' ref_str '_mask_cell.mat']));
TCs_D1 = load(fullfile(fnout, [ref_date '_' mouse], [ref_date '_' mouse '_' ref_str], [ref_date '_' mouse '_' ref_str '_TCs.mat']));
pixD1 = load(fullfile(fnout, [ref_date '_' mouse], [ref_date '_' mouse '_' ref_str], [ref_date '_' mouse '_' ref_str '_pixel.mat']));
shiftsD1 = load(fullfile(fnout, [ref_date '_' mouse], [ref_date '_' mouse '_' ref_str], [ref_date '_' mouse '_' ref_str '_reg_shifts.mat']));
data_reg_avg = shiftsD1.data_reg_avg;
pixel1 = pixD1.pix_3hz;
npSub_tc1 = TCs_D1.npSub_tc;
nCells1 = size(npSub_tc1,2);
mask_cell = maskD1.mask_cell;
cell_list = intersect(1:nCells1, unique(mask_cell));
cell_stats = regionprops(mask_cell);

transD2 = load(fullfile(fnout, [day2 '_' mouse], [day2 '_' mouse '_' run_str2], [day2 '_' mouse '_' run_str2 '_transform.mat']));
fgta2 = transD2.fitGeoTAf;
% data_reg_avg with FGTA transformation
data_reg_avg2 = transD2.r2rFGTA;
% data_reg_avg without FGTA transformation
reg2 = transD2.reg;
pixD2 = load(fullfile(fnout, [day2 '_' mouse], [day2 '_' mouse '_' run_str2], [day2 '_' mouse '_' run_str2 '_pixel.mat']));
pix_fgta2 = pixD2.pix_fgta;
maskD2 = load(fullfile(fnout, [day2 '_' mouse], [day2 '_' mouse '_' run_str2], [day2 '_' mouse '_' run_str2 '_mask_cell.mat']));
mask_cell2 = maskD2.mask_cell;
mask_np2 = maskD2.mask_np;
cell_list2 = intersect(1:nCells1, unique(mask_cell2));
cell_stats2 = regionprops(mask_cell2);

transD3 = load(fullfile(fnout, [day3 '_' mouse], [day3 '_' mouse '_' run_str], [day3 '_' mouse '_' run_str '_transform.mat']));
fgta3 = transD3.fitGeoTAf;
data_reg_avg3 = transD3.r2rFGTA;
pixD3 = load(fullfile(fnout, [day3 '_' mouse], [day3 '_' mouse '_' run_str], [day3 '_' mouse '_' run_str '_pixel.mat']));
pix_fgta3 = pixD3.pix_fgta;
maskD3 = load(fullfile(fnout, [day3 '_' mouse], [day3 '_' mouse '_' run_str], [day3 '_' mouse '_' run_str '_mask_cell.mat']));
mask_cell3 = maskD3.mask_cell;
mask_np3 = maskD3.mask_np;
cell_list3 = intersect(1:nCells1, unique(mask_cell3));
cell_stats3 = regionprops(mask_cell3);

%% transform data_reg

for i = 1:nframes
    data_reg2(:,:,i) = imwarp(double(data_reg2(:,:,i)),fgta2, 'OutputView', imref2d(size(data_reg_avg2)));
    if rem(i,50) == 0
        fprintf([num2str(i) '/n'])
    end
end
data_reg2_avg = mean(data_reg2,3);
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
    xCenter2 = round(cell_stats2(iCell).Centroid(2));
    yCenter2 = round(cell_stats2(iCell).Centroid(1));
    xLeft = (xCenter - width/2);
    yBottom = (yCenter - height/2);
    xLeft2 = (xCenter2 - width/2);
    yBottom2 = (yCenter2 - height/2);
    cell_reg_avg = data_reg_avg(xLeft:(xLeft+width),yBottom:(height+yBottom));
    pix1 = pixel1(xLeft:(xLeft+width),yBottom:(height+yBottom));
    cell_reg2 = data_reg2(xLeft2:(xLeft2+width),yBottom2:(height+yBottom2),:);
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
if r>0.6 && p>0.4
    [outs, reg_cell2] = stackRegister_MA(cell_reg2,[],[],repmat(shift,[nframes 1]));
elseif p>0.8 && r>0.4 && r<0.6
    [outs, reg_cell2] = stackRegister_MA(cell_reg2,[],[],repmat(shift2,[nframes 1]));
else
    reg_cell2 = NaN(51,51,nframes);
end
    cell_reg2_long = reshape(reg_cell2, [51.*51 nframes]);
    cell_reg2_down  = reshape(stackGroupProject(reg_cell2,5), [51.*51 nframes/5]);
    cell_mask = mask_cell2(xLeft2:(xLeft2+width),yBottom2:(height+yBottom2));
    cell_mask_long = reshape(cell_mask, [51.*51 1]);
    cell_mask_np(:,:,iCell) = mask_np2(xLeft2:(xLeft2+width),yBottom2:(height+yBottom2),iCell);
    mask_np_long = reshape(cell_mask_np(:,:,iCell), [51.*51 1]);
    data_tc(:,iCell) = mean(cell_reg2_long(find(cell_mask_long==iCell),:),1);
    data_tc_down(:,iCell) = mean(cell_reg2_down(find(cell_mask_long==iCell),:),1);
    np_tc(:,iCell) = mean(cell_reg2_long(find(mask_np_long),:),1);
    np_tc_down(:,iCell) = mean(cell_reg2_down(find(mask_np_long),:),1);
    fprintf(['Cell #' num2str(iCell) '%s/n']) 
    for i = 1:100
        x(i,iCell) = skewness(data_tc_down(:,iCell)-tcRemoveDC(np_tc_down(:,iCell)*ii(i)));
    end
    [max_skew ind(:,iCell)] =  max(x(:,iCell),[],1);
    np_w(:,iCell) = 0.01*ind(:,iCell);
    npSub_tc(:,iCell) = data_tc(:,iCell)-bsxfun(@times,tcRemoveDC(np_tc(:,iCell)),np_w(:,iCell));
end
save(fullfile(fnout, [day2 '_' mouse], [day2 '_' mouse '_' run_str], [day2 '_' mouse '_' run_str '_TCs.mat']), 'data_tc', 'np_tc', 'npSub_tc')
save(fullfile(fnout, [day2 '_' mouse], [day2 '_' mouse '_' run_str], [day2 '_' mouse '_' run_str '_input.mat']), 'input')

%% mask images
start = 1;
figure;
for iCell = 1:10
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
    cell_reg2 = data_reg2(xLeft2:(xLeft2+width),yBottom2:(height+yBottom2),:);
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
if r>0.6 && p>0.4
    [outs, reg_cell2] = stackRegister_MA(cell_reg2,[],[],repmat(shift,[nframes 1]));
elseif p>0.8 && r>0.4 && r<0.6
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
    cell_reg2 = data_reg2(xLeft2:(xLeft2+width),yBottom2:(height+yBottom2),:);
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
if r>0.6 && p>0.4
    [outs, reg_cell2] = stackRegister_MA(cell_reg2,[],[],repmat(shift,[nframes 1]));
elseif p>0.8 && r>0.4 && r<0.6
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
print(fullfile(fnout, [ref_date '_' mouse], [ref_date '_' mouse '_' ref_str], ['AcrossAllDays'], ['CellMaps'], [ref_date '_' mouse '_mask_outline.pdf']),'-dpdf', '-bestfit')


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
