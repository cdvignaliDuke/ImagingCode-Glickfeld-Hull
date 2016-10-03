% iexp = 1;
awFSAVdatasets_V1;
datasetStr = '_V1';
rc = behavConstsAV;

for iexp = 10:size(expt,2)
fn = fullfile(rc.ashleyAnalysis,expt(iexp).mouse,expt(iexp).folder,expt(iexp).date,expt(iexp).dirtuning,'cells only');

%% Direction Tuning
% load data 
mworks = ['data-' 'i' expt(iexp).SubNum '-' expt(iexp).date '-' expt(iexp).dirtuning_time]; 
load (fullfile(rc.behavData,mworks));

load(fullfile(fn,'Timecourses.mat'))

% set variables
down = 10;
nON = (input.nScansOn)./down;
nOFF = (input.nScansOff)./down;
nStim = input.gratingDirectionStepN;
nRep = size(data_TC,1)./((nON+nOFF)*nStim);
nTrials = (nStim.*nRep);
DirectionDeg = cell2mat(input.tGratingDirectionDeg);
if length(DirectionDeg > nTrials)
    DirectionDeg = DirectionDeg(1:nTrials);
end
Dirs = unique(DirectionDeg);

% dF/F, F per trial
stimOFF_ind = 1:nOFF+nON:size(data_TC,1);

dF_data = zeros(size(data_TC));
dFoverF_data = zeros(size(data_TC));
for i = 1:nTrials
    indAll = stimOFF_ind(i):stimOFF_ind(i)+(nOFF+nON-1);
    indF = stimOFF_ind(i)+5:stimOFF_ind(i)+(nOFF-1);
    dF_data(indAll,:) = bsxfun(@minus,data_TC(indAll,:),mean(data_TC(indF,:),1));
    dFoverF_data(indAll,:) = bsxfun(@rdivide,dF_data(indAll,:),mean(data_TC(indF,:),1));
end

% find on indices for the first frame of each stimulus start period and iti (Off) period

stimON_ind = nOFF+1:nOFF+nON:size(dFoverF_data,1);

% sort data_TC into 20 frame (10 pre, 10 post) traces around stimON 

dFoverFCellsTrials = zeros(10+nON,size(dFoverF_data,2),nTrials);
for i = 1:nTrials
    dFoverFCellsTrials(:,:,i) = dFoverF_data(stimON_ind(i)-10:stimON_ind(i)+(nON-1),:);
end

dFoverF_meanDirResp = zeros(size(dFoverFCellsTrials,1),size(dFoverFCellsTrials,2),nStim);
for i = 1:nStim
    trials = find(DirectionDeg(:,1:nTrials) == Dirs(i));
    dFoverF_meanDirResp(:,:,i) = mean(dFoverFCellsTrials(:,:,trials),3);
end

figure;
for i = 1:nStim
    plot(dFoverF_meanDirResp(:,1,i));
    hold on
end
title([expt(iexp).mouse expt(iexp).date])

% find magnitude of response to stim
dFoverF_meanOFFDirResp = (squeeze(mean(dFoverF_meanDirResp(1:10,:,:),1)));
DirRespPerCell = (squeeze(mean(dFoverF_meanDirResp(11:end,:,:),1)));

figure;
plot(DirRespPerCell(1,:))

% find direction preference
DirRespPerCell_sub = DirRespPerCell-min(min(min(DirRespPerCell,[],1),[],2));
[dirPref_val,dirPref_ind] = max(DirRespPerCell_sub,[],2);

pref2orth_dir = [nStim/2+1:nStim 1:nStim/2];
dirOrth_ind = zeros(size(dirPref_ind));
for i = 1:nStim
    cells = find(dirPref_ind == i);
    dirOrth_ind(cells) = pref2orth_dir(i);
end

dirOrth_val = zeros(size(dirPref_ind)); 
for i = 1:size(data_TC,2)
    ind = dirOrth_ind(i);
    cell = DirRespPerCell_sub(i,:);
    dirOrth_val(i) = cell(ind);
end

cellDSI = zeros(size(dirPref_ind));
cellDSI = (dirPref_val - dirOrth_val)./(dirPref_val + dirOrth_val);
dirSlctvCells = find(cellDSI > 0.3);
figure;hist(cellDSI,10);title('DSI')

% find orientation preference
dFoverF_meanOriResp = zeros(size(dFoverFCellsTrials,1),size(dFoverFCellsTrials,2),nStim/2);
for i = 1:nStim/2
    trials = find(DirectionDeg(:,1:nTrials) == Dirs(i) | DirectionDeg(:,1:nTrials) == Dirs(i+nStim/2));
    dFoverF_meanOriResp(:,:,i) = mean(dFoverFCellsTrials(:,:,trials),3);
end
dFoverF_meanOFFOriResp = (squeeze(mean(dFoverF_meanOriResp(1:10,:,:),1)))+1;
dFoverF_meanONOriResp = (squeeze(mean(dFoverF_meanOriResp(11:end,:,:),1)))+1;
OriRespPerCell = dFoverF_meanONOriResp;

figure;
plot(OriRespPerCell(1,:))

OriRespPerCell_sub = OriRespPerCell-min(min(min(OriRespPerCell,[],1),[],2));
[oriPref_val,oriPref_ind] = max(OriRespPerCell_sub,[],2);

pref2orth_ori = [(nStim/2)/2+1:nStim/2 1:(nStim/2)/2];
oriOrth_ind = zeros(size(oriPref_ind));
for i = 1:nStim/2
    cells = find(oriPref_ind == i);
    oriOrth_ind(cells) = pref2orth_ori(i);
end

oriOrth_val = zeros(size(oriPref_ind)); 
for i = 1:size(data_TC,2)
    ind = oriOrth_ind(i);
    cell = OriRespPerCell_sub(i,:);
    oriOrth_val(i) = cell(ind);
end

cellOSI = zeros(size(oriPref_ind));
cellOSI = (oriPref_val - oriOrth_val)./(oriPref_val + oriOrth_val);
oriSlctvCells = find(cellOSI > 0.3);
figure;hist(cellOSI,10);title('OSI')

% get tuning curves

dFoverFDirResp = zeros(nStim,size(data_TC,2));
errbar = zeros(nStim,size(data_TC,2));
for i = 1:nStim 
    trials = find(DirectionDeg == Dirs(i));
    dFoverFDirResp(i,:) = squeeze(mean(mean(dFoverFCellsTrials(11:16,:,trials),1),3));
    errbar(i,:) = std(mean(dFoverFCellsTrials(11:16,:,trials),1),[],3)/sqrt(size(dFoverFCellsTrials(11:16,:,trials),3));
end

% save tuning info
save(fullfile(fn,'TuningPreferences.mat'),'oriPref_ind','dirPref_ind','dirSlctvCells','oriSlctvCells','dFoverFDirResp','dFoverF_meanDirResp')

%% register behavior data
for irun = 1:expt(iexp).nrun;
fnout = fullfile(rc.ashleyAnalysis,expt(iexp).mouse,expt(iexp).folder,expt(iexp).date, expt(iexp).runs(irun,:),'cells only');

% load data 
mworks = ['data-' 'i' expt(iexp).SubNum '-' expt(iexp).date '-' expt(iexp).time_mat(irun,:)]; 
load (fullfile(rc.behavData,mworks));

% Set current directory to imaging data location
CD = fullfile('Z:\data\', expt(iexp).mouse, expt(iexp).folder,expt(iexp).date,expt(iexp).runs(irun,:));
cd(CD);
fName = [expt(iexp).runs(irun,:) '_000_000'];
imgMatFile = [fName '.mat'];
load(imgMatFile);

nframes = info.config.frames;
% nframes = 13577
tic
data = sbxread(fName,0,nframes);
toc
data = squeeze(data);   
    
data_sub = data-min(min(min(data,[],1),[],2),[],3);
data = data_sub;
clear data_sub

% data_avg = mean(data(:,:,2000:2100),3);
% figure; imagesq(data_avg); colormap(gray)

%get direction tuning registration image and get cells
dirFolder = fullfile(rc.ashleyAnalysis,expt(iexp).mouse,expt(iexp).folder,expt(iexp).date,expt(iexp).dirtuning);
load(fullfile(dirFolder,'regImg.mat'));
load(fullfile(fn,'mask&TCDir.mat'));
load(fullfile(fn,'neuropil.mat'));
clear npTC npSubTC data_TC

[out data_reg] = stackRegister(data, data_avg);
clear data


% get timecourses
dataTC = stackGetTimeCourses(data_reg, mask_cell);

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

dataTimecourse.dataTC = dataTC;
dataTimecourse.dataTCsub = dataTCsub;
dataTimecourse.npilTC = npTC;
mkdir(fnout)
save(fullfile(fnout,'Timecourses.mat'),'dataTimecourse')
clear data_reg dataTimecourse
end

end
