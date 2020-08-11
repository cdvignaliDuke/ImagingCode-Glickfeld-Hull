mouse = 'i739';
date = '170711';
ImgFolder = strvcat('002','003');
nrun = size(ImgFolder,1);
run_str = catRunName(ImgFolder, nrun);
irun = 1;
clear global
CD = ['Z:\home\lindsey\Data\2P_images\' date '_' mouse '\' ImgFolder(irun,:)];
cd(CD);
imgMatFile = [ImgFolder(irun,:) '_000_000.mat'];
load(imgMatFile);
nframes = 1000;
data_temp = sbxread([ImgFolder(irun,:) '_000_000'],0,nframes);

data = squeeze(data_temp);
load(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_reg_shifts.mat']))

[out reg] = stackRegister(data,data_avg);
figure; imagesc(mean(reg,3)); truesize; colormap gray; clim([0 3000])

load(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_mask_cell.mat']))
mask_cell_dig = mask_cell;
mask_cell_dig(find(mask_cell>1)) = 1;
rgb_img(:,:,1) = mask_cell_dig;
rgb_img(:,:,2) = mean(reg,3)./max(max(mean(reg,3),[],2),[],1);
rgb_img(:,:,3) = mask_cell_dig;
figure; imagesc(rgb_img); truesize
lim = [334 423];
vline(lim, '-c')
title([mouse ' ' date])
print(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_FOVareaPreDivide.pdf']),'-dpdf','-bestfit')

nCells = max(unique(mask_cell),[],1);
cell_area = cell(1,nCells);
for icell = 1:nCells
    [i j] = find(mask_cell == icell);
    ind1 = find(j<lim(1));
    ind2 = find(j>lim(2));
    if length(ind1) == length(j)
        cell_area{1,icell} = 'LM';
    elseif length(ind2) == length(j)
        cell_area{1,icell} = 'AL';
    else
        cell_area{1,icell} = 'NaN';
    end
end

ind_LM = find(strcmp(cell_area,'LM')==0);
if length(ind_LM)>0
    cell_img_LM = mask_cell_dig;
    cell_img_LM(find(mask_cell>=min(ind_LM,[],2))) = 0;
else
    cell_img_LM = zeros(size(mask_cell_dig));
end
ind_AL = find(strcmp(cell_area,'AL'));
if length(ind_AL)>0
    cell_img_AL = mask_cell_dig;
    cell_img_AL(find(mask_cell<min(ind_AL,[],2))) = 0;
else
    cell_img_AL = zeros(size(mask_cell_dig));
end
cell_rgb(:,:,1) = cell_img_AL;
cell_rgb(:,:,2) = mean(reg,3)./max(max(mean(reg,3),[],2),[],1);
cell_rgb(:,:,3) = cell_img_LM;
figure; imagesc(cell_rgb); truesize; axis off
vline(lim, '-c')
title([mouse ' ' date '; LM- blue; AL- red'])
print(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_FOVareaPostDivide.pdf']),'-dpdf','-bestfit')
save(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_areaDivide.mat']), 'lim','cell_rgb','cell_area')



