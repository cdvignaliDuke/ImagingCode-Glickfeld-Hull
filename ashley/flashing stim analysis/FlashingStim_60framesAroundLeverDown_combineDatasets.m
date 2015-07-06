%combine two datasets to have 3 trial types - vis only, aud only, and
%vis+aud
SubNum = '608';
mouse = 'AW08';
date = '150209';

%% first dataset (_1) - vis only (V) and aud only (A)
time = '1636';
ImgFolder = '003';

% load MWorks file
CD = ['Z:\data\' mouse '\mworks\' date];
cd(CD);
mworks = ['data-' 'i' SubNum '-' date '-' time]; 
load (mworks);

%load timecourse
fileSave = fullfile('Z:\analysis\',mouse,'two-photon imaging', date, ImgFolder);
cd(fileSave);
load('dataTC.mat');
dataTC_1 = dataTimecourse.dataTC;

cLeverDown_1 = cell2mat_padded(input.cLeverDown);
cTargetOn_1 = cell2mat_padded(input.cTargetOn);
tCyclesOn_1 = cell2mat_padded(input.tCyclesOn);
block2_1 = cell2mat_padded(input.tBlock2TrialNumber);
cycles_1 = unique(tCyclesOn_1);
cycTime_1 = input.nFramesOn + input.nFramesOff;
nTrials_1 = input.trialSinceReset;

V_ind_1 = find(block2_1 == 0);
A_ind = find(block2_1 == 1);

%dF/F
Data_1 = zeros(60,size(dataTC_1,2),nTrials_1);
DataDF_1 = zeros(60,size(dataTC_1,2),nTrials_1);
DFoverF_1 = zeros(60,size(dataTC_1,2),nTrials_1);
for itrial = 1:nTrials_1
    Data_1(:,:,itrial) = dataTC_1(cLeverDown_1(itrial)-30:cLeverDown_1(itrial)+29,:);
    DataDF_1(:,:,itrial) = bsxfun(@minus,Data_1(:,:,itrial),mean(Data_1(1:30,:,itrial),1));
    DFoverF_1(:,:,itrial) = bsxfun(@rdivide,DataDF_1(:,:,itrial),mean(Data_1(1:30,:,itrial),1));
end

% aud only trials
DFoverFmeanA = mean(mean(DFoverF_1(:,:,A_ind),3),2);

%% second dataset (_2) - vis only (V) and vis + aud only (AV)
time = '1653';
ImgFolder = '004';

% load MWorks file
CD = ['Z:\data\' mouse '\mworks\' date];
cd(CD);
mworks = ['data-' 'i' SubNum '-' date '-' time]; 
load (mworks);

%load timecourse
fileSave = fullfile('Z:\analysis\',mouse,'two-photon imaging', date, ImgFolder);
cd(fileSave);
load('dataTC.mat');
dataTC_2 = dataTimecourse.dataTC;

cLeverDown_2 = cell2mat_padded(input.cLeverDown);
cTargetOn_2 = cell2mat_padded(input.cTargetOn);
tCyclesOn_2 = cell2mat_padded(input.tCyclesOn);
block2_2 = cell2mat_padded(input.tBlock2TrialNumber);
cycles_2 = unique(tCyclesOn_2);
cycTime_2 = input.nFramesOn + input.nFramesOff;
nTrials_2 = input.trialSinceReset;

V_ind_2 = find(block2_2 == 0);
AV_ind = find(block2_2 ==1);

%dF/F
Data_2 = zeros(60,size(dataTC_2,2),nTrials_2);
DataDF_2 = zeros(60,size(dataTC_2,2),nTrials_2);
DFoverF_2 = zeros(60,size(dataTC_2,2),nTrials_2);
for itrial = 1:nTrials_2
    Data_2(:,:,itrial) = dataTC_2(cLeverDown_2(itrial)-30:cLeverDown_2(itrial)+29,:);
    DataDF_2(:,:,itrial) = bsxfun(@minus,Data_2(:,:,itrial),mean(Data_2(1:30,:,itrial),1));
    DFoverF_2(:,:,itrial) = bsxfun(@rdivide,DataDF_2(:,:,itrial),mean(Data_2(1:30,:,itrial),1));
end

% vis+aud trials
DFoverFmeanAV = mean(mean(DFoverF_2(:,:,AV_ind),3),2);

%vis only trials
DFoverF_V = cat(3,DFoverF_1(:,:,V_ind_1),DFoverF_2(:,:,V_ind_2));
DFoverFmeanV = mean(mean(DFoverF_V,3),2);

%% plot average all trials all cells for each condition

figure;
plot(DFoverFmeanV,'g');
hold on
plot(DFoverFmeanA,'r');
hold on
plot(DFoverFmeanAV, 'm');
hold on
vline(30,'k');
hold on
numCyc = ceil(30/cycTime_2);
for i = 1:numCyc
    line = (i-1)*cycTime_2+30;
    vline(line,'k:')
    hold on
end
errbar_V = std(mean(DFoverF_V,2),[],3)/sqrt(size(DFoverF_V,3));
errbar_A = std(mean(DFoverF_1(:,:,A_ind),2),[],3)/sqrt(size(DFoverF_1(:,:,A_ind),3));
errbar_AV = std(mean(DFoverF_2(:,:,AV_ind),2),[],3)/sqrt(size(DFoverF_2(:,:,AV_ind),3));
shadedErrorBar([],DFoverFmeanV,errbar_V, 'g', 1);
shadedErrorBar([],DFoverFmeanA,errbar_A, 'r', 1);
shadedErrorBar([],DFoverFmeanAV,errbar_AV, 'm', 1);
hold on
title('Avg all trials, all cells; 60 frames around trial start');
ylabel('dF/F');
xlabel('frames');





