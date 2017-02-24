% edit Load_SBXdataset_fast.m
%%

awFSAVdatasets_V1
iexp = 17;

for irun = 1:expt(iexp).nrun;

SubNum = expt(iexp).SubNum;
date = expt(iexp).date;
time = expt(iexp).time_mat(irun,:);
ImgFolder = expt(iexp).runs(irun,:);
mouse = expt(iexp).mouse;
fName = [ImgFolder '_000_000'];

Load_SBXdataPlusMWorksData    
    
data_sub = data-min(min(min(data,[],1),[],2),[],3);
data = data_sub;
clear data_sub

% data_avg = mean(data(:,:,2000:2100),3);
% figure; imagesq(data_avg); colormap(gray)

%get direction tuning registration image and get cells
dirFolder = expt(iexp).dirtuning;
fileDirMasks = fullfile('Z:\analysis\',mouse,'two-photon imaging', date, dirFolder);
cd(fileDirMasks);
load('regImg.mat');
load('maskROIcombo.mat');

mask_cell = mask_cell_combo;
neuropil = neuropil_combo;

[out data_reg] = stackRegister(data, data_avg);
clear data


% get timecourses
cellsUsed = unique(mask_cell(:));
cellsUsed = cellsUsed(2:end);
dataTC = zeros(size(data_reg,3),length(cellsUsed));
for iCell = 1:length(cellsUsed)
    
    maskLin = mask_cell(:);
    mask_temp = zeros(size(maskLin));
    mask_temp(maskLin == cellsUsed(iCell)) = 1;
    mask_temp = reshape(mask_temp,size(mask_cell));
    dataTC(:,iCell) = stackGetTimeCourses(data_reg, mask_temp);
end

% get neuropil timecourses
buf = 4;
np = 6;
nCells = size(dataTC,2);


npTC = zeros(size(dataTC));
for i = 1:nCells
    tempNPmask = squeeze(neuropil(:,:,i));
    if sum(sum(tempNPmask)) > 0
    npTC(:,i) = stackGetTimeCourses(data_reg,tempNPmask);
    end
end

% npTC = stackGetTimeCourses(data_reg,neuropil);


dataTC_mavg = tsmovavg(dataTC,'s',10,1);
npTC_mavg = tsmovavg(npTC,'s',10,1);

% down = 10;
% dataTC_down = reshape(mean(reshape(dataTC',size(dataTC',1),down,size(dataTC',2)/down),2),size(dataTC',1),size(dataTC',2)/down)';
% npTC_down = reshape(mean(reshape(npTC',size(npTC',1),down,size(npTC',2)/down),2),size(npTC',1),size(npTC',2)/down)';

ii= 0.01:0.01:1;
x = zeros(length(ii), nCells);
for i = 1:100
    x(i,:) = skewness(dataTC_mavg-tcRemoveDC(npTC_mavg*ii(i)));
end
[max_skew ind] =  max(x,[],1);
% skew(buf,:) = max_skew;
np_w = 0.01*ind;
dataTCsub = dataTC-bsxfun(@times,tcRemoveDC(npTC),np_w);


%save in
try
    fileSave = fullfile('Z:\analysis\',mouse,'two-photon imaging', date, ImgFolder);
    cd(fileSave);
catch
    fileSave = fullfile('Z:\analysis\',mouse,'two-photon imaging', date);
    cd(fileSave);
    mkdir(ImgFolder);
    fileSave = fullfile('Z:\analysis\',mouse,'two-photon imaging', date, ImgFolder);
    cd(fileSave);
end
dataTimecourse.dataTC = dataTC;
dataTimecourse.dataTCsub = dataTCsub;
dataTimecourse.npilTC = npTC;
dataTimecourse.mouse = mouse;
dataTimecourse.date = date;
dataTimecourse.dataset = ImgFolder;
save('Timecourses.mat','dataTimecourse')
clear data_reg dataTimecourse
end
