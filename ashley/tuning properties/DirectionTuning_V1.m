%load data

SubNum = '516';
date = '141003';
time = '1210';
ImgFolder = '005';
mouse = '516';
fName = '005_000_000';

% load MWorks file
% load MWorks file
CD = ['Z:\data\' mouse '\mworks\' date];
cd(CD);
mworks = ['data-' 'i' SubNum '-' date '-' time]; 
load (mworks);


% analysis folder
try
    filedir = fullfile('Z:\analysis\',mouse,'two-photon imaging', date, ImgFolder);
    cd(filedir);
catch
    filedir = fullfile('Z:\analysis\',mouse,'two-photon imaging');
    cd(filedir)
    mkdir(date,ImgFolder)
    filedir = fullfile('Z:\analysis\',mouse,'two-photon imaging', date, ImgFolder);
    cd(filedir);
end

data_reg = readtiff('DirectionTuning_V1.tif');
%used load_SBXdataset_fast.m to load data
%% Parameters
down = 10;
nON = (input.nScansOn)./down;
nOFF = (input.nScansOff)./down;
nStim = input.gratingDirectionStepN;

%% downsample and register data

%average signals in time
data_down = stackGroupProject(data,down);
clear data

%remove negative data (by addition)
data_sub = data_down-min(min(min(data_down,[],1),[],2),[],3);
clear data_down

% register
data_avg = mean(data_sub(:,:,90:100),3);
figure; imagesq(data_avg); colormap(gray)

[out data_reg] = stackRegister(data_sub, data_avg);
clear data_sub

%save data_reg
writetiff(data_reg, 'DirectionTuning_V1');

%%
nRep = size(data_reg,3)./((nON+nOFF)*nStim);
nTrials = (nStim.*nRep);
%% create dF/F stack

nOFF_ind = zeros(1,(nOFF*nStim*nRep));
start = 1;
for iStim = 1:(nRep*nStim)
    nOFF_ind(1, start:start+nOFF-1) = 1+((iStim-1)*(nOFF+nON)):nOFF + ((iStim-1)*(nOFF+nON));
    start = start+nOFF;
end

nON_ind = setdiff(1:size(data_reg,3),nOFF_ind);
nON_avg = mean(data_reg(:,:,nON_ind),3);
nOFF_avg = mean(data_reg(:,:,nOFF_ind),3);

% dF average F
dF_data = bsxfun(@minus,data_reg, nOFF_avg);

max_dF = max(dF_data,[],3);
figure; imagesq(max_dF); colormap(gray)

bwout = imCellEditInteractive(max_dF);
mask_cell = bwlabel(bwout);

data_TC = stackGetTimeCourses(data_reg,mask_cell);
figure; tcOffsetPlot(data_TC)

negImg_maxDF = mask_cell<1;

%% use correlation dF/F to find ROIS

% %use max dF if too many cells
% b = 5;
% siz = size(data_reg);
% corr_map = zeros(siz(1),siz(2));
% for ix = b:siz(2)-b
%     for iy = b:siz(1)-b
%         TC = data_reg(iy,ix,:);
%         surround = (data_reg(iy-1,ix-1,:)+data_reg(iy-1,ix,:)+data_reg(iy-1,ix+1,:)+data_reg(iy,ix-1,:)+data_reg(iy,ix+1,:)+data_reg(iy+1,ix-1,:)+data_reg(iy+1,ix,:)+data_reg(iy+1,ix+1,:))/8;
%         R = corrcoef(TC,surround);
%         corr_map(iy,ix) = R(1,2);
%     end
% end
% figure; imagesq(corr_map); colormap(gray)
% 
% bwout = imCellEditInteractive(corr_map);
% mask_cell = bwlabel(bwout);

%timecourses
data_TC = stackGetTimeCourses(data_reg,mask_cell);
figure; tcOffsetPlot(data_TC)
%%
fileSave = fullfile('Z:\analysis\',mouse,'two-photon imaging', date, ImgFolder);
cd(fileSave);
save('mask&TCDir.mat','mask_cell','data_TC');
% load('mask&TCDir.mat')


%%
orig_rate = 30;
final_rate = 3;
down = orig_rate./final_rate;
nON = (input.nScansOn)./down;
nOFF = (input.nScansOff)./down;
nStim = input.gratingDirectionStepN;
nRep = size(data_TC,1)./((nON+nOFF)*nStim);
nTrials = (nStim.*nRep);
DirectionDeg = cell2mat(input.tGratingDirectionDeg);
Dirs = unique(DirectionDeg);

%% dF/F by trial
% F per trial
stimOFF_ind = 1:nOFF+nON:size(data_TC,1);

dF_data = zeros(size(data_TC));
dFoverF_data = zeros(size(data_TC));
for i = 1:nTrials
    indAll = stimOFF_ind(i):stimOFF_ind(i)+(nOFF+nON-1);
    indF = stimOFF_ind(i)+5:stimOFF_ind(i)+(nOFF-1);
    dF_data(indAll,:) = bsxfun(@minus,data_TC(indAll,:),mean(data_TC(indF,:),1));
    dFoverF_data(indAll,:) = bsxfun(@rdivide,dF_data(indAll,:),mean(data_TC(indF,:),1));
end
%% dF/F (by cell) for each stimulus type

% find on indices for the first frame of each stimulus start period and iti (Off) period

stimON_ind = nOFF+1:nOFF+nON:size(dFoverF_data,1);

% sort data_TC into 20 frame (10 pre, 10 post) traces around stimON 

dFoverFCellsTrials = zeros(10+nON,size(dFoverF_data,2),nTrials);
for i = 1:nTrials
    dFoverFCellsTrials(:,:,i) = dFoverF_data(stimON_ind(i)-10:stimON_ind(i)+(nON-1),:);
end

dFoverF_meanDirResp = zeros(size(dFoverFCellsTrials,1),size(dFoverFCellsTrials,2),nStim);
for i = 1:nStim
    trials = find(DirectionDeg == Dirs(i));
    dFoverF_meanDirResp(:,:,i) = mean(dFoverFCellsTrials(:,:,trials),3);
end

figure;
for i = 1:nStim
    plot(dFoverF_meanDirResp(:,10,i));
    hold on
end

%% find magnitude of response to stim

dFoverF_meanOFFDirResp = (squeeze(mean(dFoverF_meanDirResp(1:10,:,:),1)));

dFoverF_meanONDirResp = (squeeze(mean(dFoverF_meanDirResp(11:end,:,:),1)));
DirRespPerCell = dFoverF_meanONDirResp;
% DirRespPerCell = dFoverF_meanONDirResp./dFoverF_meanOFFDirResp;

figure;
plot(DirRespPerCell(10,:))

%% find direction preference

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

% DirRespPerCell_abvZero = DirRespPerCell;
% DirRespPerCell_abvZero(DirRespPerCell_abvZero < 0) = 0;
% 
% [dirPref_val,dirPref_ind] = max(DirRespPerCell_abvZero,[],2);
% 
% pref2orth_dir = [nStim/2+1:nStim 1:nStim/2];
% dirOrth_ind = zeros(size(dirPref_ind));
% for i = 1:nStim
%     cells = find(dirPref_ind == i);
%     dirOrth_ind(cells) = pref2orth_dir(i);
% end
% 
% dirOrth_val = zeros(size(dirPref_ind)); 
% for i = 1:size(data_TC,2)
%     ind = dirOrth_ind(i);
%     cell = DirRespPerCell_abvZero(i,:);
%     dirOrth_val(i) = cell(ind);
% end
% 
% cellDSI = zeros(size(dirPref_ind));
% cellDSI = (dirPref_val - dirOrth_val)./(dirPref_val + dirOrth_val);
% dirSlctvCells = find(cellDSI > 0.3);

%% find orientation preference
dFoverF_meanOriResp = zeros(size(dFoverFCellsTrials,1),size(dFoverFCellsTrials,2),nStim/2);
for i = 1:nStim/2
    trials = find(DirectionDeg == Dirs(i) | DirectionDeg == Dirs(i+nStim/2));
    dFoverF_meanOriResp(:,:,i) = mean(dFoverFCellsTrials(:,:,trials),3);
end

dFoverF_meanOFFOriResp = (squeeze(mean(dFoverF_meanOriResp(1:10,:,:),1)))+1;

dFoverF_meanONOriResp = (squeeze(mean(dFoverF_meanOriResp(11:end,:,:),1)))+1;

OriRespPerCell = dFoverF_meanONOriResp;

figure;
plot(OriRespPerCell(24,:))

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
% [oriPref_val,oriPref_ind] = max(OriRespPerCell,[],2);
% 
% pref2orth_ori = [(nStim/2)/2+1:nStim/2 1:(nStim/2)/2];
% oriOrth_ind = zeros(size(oriPref_ind));
% for i = 1:nStim/2
%     cells = find(oriPref_ind == i);
%     oriOrth_ind(cells) = pref2orth_ori(i);
% end
% 
% oriOrth_val = zeros(size(oriPref_ind)); 
% for i = 1:size(data_TC,2)
%     ind = oriOrth_ind(i);
%     cell = OriRespPerCell(i,:);
%     oriOrth_val(i) = cell(ind);
% end
% 
% cellOSI = zeros(size(oriPref_ind));
% cellOSI = (oriPref_val - oriOrth_val)./(oriPref_val + oriOrth_val);

%save orientation and direction preferences
save('oriTuningPreferences.mat','oriPref_ind','dirPref_ind','dirSlctvCells','oriSlctvCells')
