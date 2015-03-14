%combine two datasets to have 3 trial types - vis only, aud only, and
%vis+aud
SubNum = '607';
mouse = 'AW07';
date = '150109';

%% first dataset (_1) - vis only (V) and aud only (A)
time = '1734';
ImgFolder = '004';

% load MWorks file
CD = ['Z:\data\' mouse '\mworks\' date];
cd(CD);
mworks = ['data-' 'i' SubNum '-' date '-' time]; 
load (mworks);
input1 = input;

%load timecourse
fileSave = fullfile('Z:\analysis\',mouse,'two-photon imaging', date, ImgFolder);
cd(fileSave);
load('dataTC.mat');
dataTC_1 = dataTimecourse.dataTC;

clear input dataTimecourse
%% second dataset (_2) - vis only (V) and vis + aud only (AV)
time = '1717';
ImgFolder = '003';

% load MWorks file
CD = ['Z:\data\' mouse '\mworks\' date];
cd(CD);
mworks = ['data-' 'i' SubNum '-' date '-' time]; 
load (mworks);
input2 = input;

%load timecourse
fileSave = fullfile('Z:\analysis\',mouse,'two-photon imaging', date, ImgFolder);
cd(fileSave);
load('dataTC.mat');
dataTC_2 = dataTimecourse.dataTC;

clear input dataTimecourse
%% combine dataTCs
dataTC = cat(1,dataTC_1,dataTC_2);
%% mworks variables
addInd = size(dataTC_1,1);
cLeverDown = cat(1,cell2mat_padded(input1.cLeverDown),(cell2mat_padded(input2.cLeverDown)+addInd));
cTargetOn = cat(1,cell2mat_padded(input1.cTargetOn),(cell2mat_padded(input2.cTargetOn)+addInd));
tCyclesOn = cat(1,cell2mat_padded(input1.tCyclesOn),cell2mat_padded(input2.tCyclesOn));
block2 = cat(1,cell2mat_padded(input1.tBlock2TrialNumber),cell2mat_padded(input2.tBlock2TrialNumber));
cycles = unique(tCyclesOn);
cycTime = input1.nFramesOn + input1.nFramesOff;
nTrials = input1.trialSinceReset+input2.trialSinceReset;

V_ind = find(block2 == 0);
A_ind = find(block2(1:input1.trialSinceReset,:));
AV_ind = find(block2(input1.trialSinceReset:end,:));

F = zeros(30,size(dataTC,2),nTrials);
Data = zeros(60,size(dataTC,2),nTrials);
DataDF = zeros(60,size(dataTC,2),nTrials);
DFoverF = zeros(60,size(dataTC,2),nTrials);
for itrial = 1:nTrials
    F(:,:,itrial) = dataTC(cLeverDown(itrial)-30:cLeverDown(itrial)-1,:);
    Data(:,:,itrial) = dataTC(cTargetOn(itrial)-30:cTargetOn(itrial)+29,:);
    DataDF(:,:,itrial) = bsxfun(@minus,Data(:,:,itrial),mean(F(:,:,itrial),1));
    DFoverF(:,:,itrial) = bsxfun(@rdivide,DataDF(:,:,itrial),mean(F(:,:,itrial),1));
end

% trial-type averages
DFoverFmeanAV = mean(mean(DFoverF(:,:,AV_ind),3),2);
DFoverFmeanV = mean(mean(DFoverF(:,:,V_ind),3),2);
DFoverFmeanA = mean(mean(DFoverF(:,:,A_ind),3),2);

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
numCyc = ceil(30/cycTime);
for i = 1:numCyc
    line = (i-1)*cycTime+30;
    vline(line,'k:')
    hold on
end
errbar_V = std(mean(DFoverF(:,:,V_ind),2),[],3)/sqrt(size(DFoverF(:,:,V_ind),3));
errbar_A = std(mean(DFoverF(:,:,A_ind),2),[],3)/sqrt(size(DFoverF(:,:,A_ind),3));
errbar_AV = std(mean(DFoverF(:,:,AV_ind),2),[],3)/sqrt(size(DFoverF(:,:,AV_ind),3));
shadedErrorBar([],DFoverFmeanV,errbar_V, 'g', 1);
shadedErrorBar([],DFoverFmeanA,errbar_A, 'r', 1);
shadedErrorBar([],DFoverFmeanAV,errbar_AV, 'm', 1);
hold on
title('Avg all trials, all cells; 60 frames around trial start');
ylabel('dF/F');
xlabel('frames');
