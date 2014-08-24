%% Parameters
SubNum = '001';
date = '140815';
time = '006';
ImgFolder = '006';
mouse = 'AW01';
fName = '006_000_000';
experiment = 'Flashing Stim';

%%

%load retinotopy info
retfile = ['Z:\2P imaging\Analysis\' mouse '\' experiment '\' date '\Retinotopy\analysis']; 
load (retfile);
clear data_reg

%load dataset and mworks file with 'Load_SBXdataset_fast.m'
edit Load_SBXdataset_fast.m

%remove negative data (by addition)
data_sub = data-min(min(min(data,[],1),[],2),[],3);
data = data_sub
clear data_sub

[out data_reg] = stackRegister(data, data_avg);
clear data

save 'data_reg.mat' data_reg
%% call some mworks variables

nTrials = input.trialSinceReset;
cLeverDown = cell2mat(input.cLeverDown);
cTargetOn = cell2mat(input.cTargetOn);
cLeverUp = cell2mat(input.cLeverUp);

%% Find off and on indices for each trial

start = 1;
start1 = 1;
for itrial = 1:nTrials
    start2 = cLeverDown(itrial)-start +1;
    nOFF_ind (1,start1:start1+start2-1) =  start:cLeverDown(itrial);
    start1 = start1+start2;
    start = cLeverUp(itrial);
    
end
siz = size(nOFF_ind,2);
nOFF_ind = nOFF_ind(:,1:siz-1);
% start1 = 1;
% for itrial = 1:nTrials
%     start = cLeverDown(itrial);
%     start2 = cLeverUp(itrial)-start +1;
%     nON_ind (1,start1:start1+start2) = start:cLeverUp(itrial);
%     start1 = start1+start2;
% end

%F=average fluorescence for all iti
F1 = mean(data_reg(:,:,nOFF_ind),3);

%% Find dFoverF for each trial for F1 iti

dF1_data = bsxfun(@minus,data_reg, F1);
dF1overF1_data = bsxfun(@rdivide, dF1_data, F1);
max_dF1 = max(dF1overF1_data,[],3);
figure; imagesq(max_dF1); colormap(gray)   

save 'dF1_data.mat' dF1_data
%%
%using Retinotopy_V1 ROI mask, find ROIs in flashingstim dataset
data_TC_FS = stackGetTimeCourses(dF1overF1_data,mask_cell);
figure; tcOffsetPlot(data_TC_FS)

%save data_reg and ROI timecourses
save 'timecourses.mat' data_TC_FS;

%%
%plot data_TC_FS for just cells that have the preferred retinotopy
CellPrefTimecourses = figure; tcOffsetPlot(data_TC_FS(:,FSpos_ind))
saveas(CellPrefTimecourses,'CellPrefTimecourses.fig')    
    

%%
%Plot data_TC_FS in ms instead of frames

