%combine two datasets to have 3 trial types - vis only, aud only, and
%vis+aud
SubNum = '613';
mouse = 'AW13';
date = '150511';

%% first dataset (_1) - vis only (V) and aud only (A)
time = '1420';
ImgFolder = '001';

% load MWorks file
CD = ['\\CRASH.dhe.duke.edu\data\home\ashley\data\' mouse '\mworks\' date];
cd(CD);
mworks = ['data-' 'i' SubNum '-' date '-' time]; 
load (mworks);
input1 = input;

%load timecourse
fileSave = fullfile('\\CRASH.dhe.duke.edu\data\home\ashley\analysis\',mouse,'two-photon imaging', date, ImgFolder);
cd(fileSave);
% load('dataTC.mat');
load('Timecourses.mat')
dataTC_1 = dataTimecourse.dataTCsub;
% dataTC_1 = dataTimecourse.dataTC;

clear input dataTimecourse
%% second dataset (_2) - vis only (V) and vis + aud only (AV)
time = '1437';
ImgFolder = '002';

% load MWorks file
CD = ['\\CRASH.dhe.duke.edu\data\home\ashley\data\' mouse '\mworks\' date];
cd(CD);
mworks = ['data-' 'i' SubNum '-' date '-' time]; 
load (mworks);
input2 = input;

%load timecourse
fileSave = fullfile('\\CRASH.dhe.duke.edu\data\home\ashley\analysis\',mouse,'two-photon imaging', date, ImgFolder);
cd(fileSave);
% load('dataTC.mat');
load('Timecourses.mat')
dataTC_2 = dataTimecourse.dataTCsub;
% dataTC_2 = dataTimecourse.dataTC;

clear input dataTimecourse
%% third dataset (_3) 
time = '1453';
ImgFolder = '003';

% load MWorks file
CD = ['\\CRASH.dhe.duke.edu\data\home\ashley\data\' mouse '\mworks\' date];
cd(CD);
mworks = ['data-' 'i' SubNum '-' date '-' time]; 
load (mworks);
input3 = input;

%load timecourse
fileSave = fullfile('\\CRASH.dhe.duke.edu\data\home\ashley\analysis\',mouse,'two-photon imaging', date, ImgFolder);
cd(fileSave);
% load('dataTC.mat');
load('Timecourses.mat')
dataTC_3 = dataTimecourse.dataTCsub;

clear input dataTimecourse
%% combine dataTCs
dataTC = cat(1,dataTC_1,dataTC_2,dataTC_3);

%% mworks variables - 3 datasets
addInd_2 = size(dataTC_1,1);
addInd_3 = size(dataTC_1,1)+addInd_2;
cLeverDown = cat(1,cell2mat_padded(input1.cLeverDown),(cell2mat_padded(input2.cLeverDown)+addInd_2),(cell2mat_padded(input3.cLeverDown)+addInd_3));
cLeverUp = cat(1,cell2mat_padded(input1.cLeverUp),(cell2mat_padded(input2.cLeverUp)+addInd_2),(cell2mat_padded(input3.cLeverUp)+addInd_3));
cTargetOn = cat(1,cell2mat_padded(input1.cTargetOn),(cell2mat_padded(input2.cTargetOn)+addInd_2),(cell2mat_padded(input3.cTargetOn)+addInd_3));
cCatchOn = cat(1,cell2mat_padded(input1.cCatchOn),(cell2mat_padded(input2.cCatchOn)+addInd_2),(cell2mat_padded(input3.cCatchOn)+addInd_3));
cCatchOn(cCatchOn == addInd_2 | cCatchOn == addInd_3) = 0;
tCyclesOn = cat(1,cell2mat_padded(input1.tCyclesOn),cell2mat_padded(input2.tCyclesOn),cell2mat_padded(input3.tCyclesOn));
nCyclesOn = cat(1,cell2mat_padded(input1.nCyclesOn),cell2mat_padded(input2.nCyclesOn),cell2mat_padded(input3.nCyclesOn));
isFA = cat(1,cell2mat_padded(input1.tFalseAlarm),cell2mat_padded(input2.tFalseAlarm),cell2mat_padded(input3.tFalseAlarm));
catchCycle = cat(1,cell2mat_padded(input1.catchCyclesOn),cell2mat_padded(input2.catchCyclesOn),cell2mat_padded(input3.catchCyclesOn));
cycles = unique(tCyclesOn);
cycTime = input1.nFramesOn + input1.nFramesOff;
frameRateS = 30; %hard-coded for now, but should be available in scanbox-yeti datasets' info file
RateFRperMS = frameRateS/1000;
% cycTime = ceil((input1.stimOnTimeMs+input1.stimOffTimeMs)*RateFRperMS);
nTrials = input1.trialSinceReset+input2.trialSinceReset+input3.trialSinceReset;
trialOutcome = cat(2,input1.trialOutcomeCell,input2.trialOutcomeCell,input3.trialOutcomeCell);
DirectionDeg = cell2mat_padded(cat(2,input1.gratingDirectionDeg,input2.gratingDirectionDeg,input3.gratingDirectionDeg));
catchDirectionDeg = cell2mat_padded(cat(2,input1.tCatchGratingDirectionDeg,input2.tCatchGratingDirectionDeg,input3.tCatchGratingDirectionDeg));
Dirs = unique(DirectionDeg);
catchDirs = unique(catchDirectionDeg);
isCatchTrial = catchDirectionDeg > 0;
tooFastTime = input1.nFramesTooFast;
maxReactTime = input1.nFramesReact;
minCyclesOn = input1.minCyclesOn;
maxCyclesOn = input1.maxCyclesOn;



catchTrialOutcome = num2cell(NaN(length(nCyclesOn),1));
catchIndex = find(isCatchTrial == 1);
for i = 1:sum(isCatchTrial)
    if isFA(catchIndex(i)) == 1
        catchTrialOutcome{catchIndex(i),1} = 'FA';
    elseif cCatchOn(catchIndex(i)) == 0
        catchTrialOutcome{catchIndex(i),1} = 'failure';
    elseif (cLeverUp(catchIndex(i)) - cCatchOn(catchIndex(i))) < tooFastTime
        catchTrialOutcome{catchIndex(i),1} = 'failure';
    elseif (cLeverUp(catchIndex(i)) - cCatchOn(catchIndex(i))) > maxReactTime
        catchTrialOutcome{catchIndex(i),1} = 'CR';
    end
end

block2 = cat(1,cell2mat_padded(input1.tBlock2TrialNumber),cell2mat_padded(input2.tBlock2TrialNumber),cell2mat_padded(input3.tBlock2TrialNumber));
V_ind = find(block2 == 0);
% AorAV_ind = find(block2 == 1);
% A_dataset_ind = AorAV_ind <= input1.trialSinceReset;
% A_ind = AorAV_ind(A_dataset_ind == 1);
% AV_ind = AorAV_ind(A_dataset_ind == 0);
AV_ind = find(block2 == 1);
%% new folder to save in
% 
% ImgFolder = '004+005'
% try
%     fileSave = fullfile('Z:\analysis\',mouse,'two-photon imaging', date, ImgFolder);
%     cd(fileSave);
% catch
%     fileSave = fullfile('Z:\analysis\',mouse,'two-photon imaging', date);
%     cd(fileSave);
%     mkdir(ImgFolder)
%     fileSave = fullfile('Z:\analysis\',mouse,'two-photon imaging', date, ImgFolder);
%     cd(fileSave);
% end
% 
% % save('dataTC.mat','dataTC');

%% divide up data by cycle- align to lever down

for icyc = 1:length(cycles)
    ind = find(tCyclesOn == cycles(icyc));
    Data = zeros(cycTime.*(cycles(icyc)+1)+60,size(dataTC,2),length(ind));
    DataDF = zeros(cycTime.*(cycles(icyc)+1)+60,size(dataTC,2),length(ind));
    DataDFoverF = zeros(cycTime.*(cycles(icyc)+1)+60,size(dataTC,2),length(ind));
    if cLeverDown(end,1)+30 > size(dataTC,1)
        for itrial = 1:length(ind)-1
            Data(:,:,itrial) = dataTC(cLeverDown(ind(itrial))-30:cLeverDown(ind(itrial))+29+(cycTime.*(cycles(icyc)+1)),:);
            DataDF(:,:,itrial) = bsxfun(@minus, Data(:,:,itrial), mean(Data(1:30,:,itrial),1));
            DataDFoverF(:,:,itrial) = bsxfun(@rdivide, DataDF(:,:,itrial), mean(Data(1:30,:,itrial),1));
        end
    else
        for itrial = 1:length(ind)
            Data(:,:,itrial) = dataTC(cLeverDown(ind(itrial))-30:cLeverDown(ind(itrial))+29+(cycTime.*(cycles(icyc)+1)),:);
            DataDF(:,:,itrial) = bsxfun(@minus, Data(:,:,itrial), mean(Data(1:30,:,itrial),1));
            DataDFoverF(:,:,itrial) = bsxfun(@rdivide, DataDF(:,:,itrial), mean(Data(1:30,:,itrial),1));
        end
    end
    cycData{icyc} = Data;
    cycDataDF{icyc} = DataDF;
    cycDataDFoverF{icyc} = DataDFoverF;
end

%% Align data to lever up
Data = zeros(105,size(dataTC,2),length(trialOutcome));
DataDF = zeros(105,size(dataTC,2),length(trialOutcome));
DataDFoverF = zeros(105,size(dataTC,2),length(trialOutcome));
if cLeverUp(end,1)+30 > size(dataTC,1)
    for itrial = 1:length(trialOutcome)-1
        Data(:,:,itrial) = dataTC(cLeverUp(itrial)-30:cLeverUp(itrial)+74,:);
        DataDF(:,:,itrial) = bsxfun(@minus, Data(:,:,itrial), mean(Data(1:30,:,itrial),1));
        DataDFoverF(:,:,itrial) = bsxfun(@rdivide, DataDF(:,:,itrial), mean(Data(1:30,:,itrial),1));
    end
else
    for itrial = 1:length(trialOutcome)
        Data(:,:,itrial) = dataTC(cLeverUp(itrial)-30:cLeverUp(itrial)+74,:);
        DataDF(:,:,itrial) = bsxfun(@minus, Data(:,:,itrial), mean(Data(1:30,:,itrial),1));
        DataDFoverF(:,:,itrial) = bsxfun(@rdivide, DataDF(:,:,itrial), mean(Data(1:30,:,itrial),1));
    end
end

DataDFoverFavg = squeeze(mean(DataDFoverF,2));
FIx = find(strcmp(trialOutcome, 'failure'));
SIx = find(strcmp(trialOutcome, 'success'));
FIxlong = intersect(find(tCyclesOn>3), FIx);
SIxlong = intersect(find(tCyclesOn>3), SIx);
Fb1Ix = intersect(V_ind, FIxlong);
Fb2Ix = intersect(AV_ind, FIxlong);
Sb1Ix = intersect(V_ind, SIxlong);
Sb2Ix = intersect(AV_ind, SIxlong);

figure;
tt = [-30:74].*(1000/30);
subplot(2,2,1)
shadedErrorBar(tt, nanmean(DataDFoverFavg(:,Fb1Ix),2), nanstd(DataDFoverFavg(:,Fb1Ix),[],2)/sqrt(length(Fb1Ix)), 'r');
hold on
shadedErrorBar(tt, nanmean(DataDFoverFavg(:,Sb1Ix),2), nanstd(DataDFoverFavg(:,Sb1Ix),[],2)/sqrt(length(Sb1Ix)), 'k'); 
title(['Visual trials: ' num2str(length(Sb1Ix)) ' Successes; ' num2str(length(Fb1Ix)) ' Earlies']) 
xlim([-500 1500])
subplot(2,2,2)
shadedErrorBar(tt, nanmean(DataDFoverFavg(:,Fb2Ix),2), nanstd(DataDFoverFavg(:,Fb2Ix),[],2)/sqrt(length(Fb2Ix)), 'r');
hold on
shadedErrorBar(tt, nanmean(DataDFoverFavg(:,Sb2Ix),2), nanstd(DataDFoverFavg(:,Sb2Ix),[],2)/sqrt(length(Sb2Ix)), 'k'); 
title(['Auditory trials: ' num2str(length(Sb2Ix)) ' Successes; ' num2str(length(Fb2Ix)) ' Earlies']) 
xlim([-500 1500])
subplot(2,2,3)
shadedErrorBar(tt, nanmean(DataDFoverFavg(:,Fb1Ix),2), nanstd(DataDFoverFavg(:,Fb1Ix),[],2)/sqrt(length(Fb1Ix)), 'k');
hold on
shadedErrorBar(tt, nanmean(DataDFoverFavg(:,Fb2Ix),2), nanstd(DataDFoverFavg(:,Fb2Ix),[],2)/sqrt(length(Fb2Ix)), 'g'); 
title(['Early trials: ' num2str(length(Fb1Ix)) ' Visual; ' num2str(length(Fb2Ix)) ' Auditory']) 
xlim([-500 1500])
subplot(2,2,4)
shadedErrorBar(tt, nanmean(DataDFoverFavg(:,Sb1Ix),2), nanstd(DataDFoverFavg(:,Sb1Ix),[],2)/sqrt(length(Sb1Ix)), 'k');
hold on
shadedErrorBar(tt, nanmean(DataDFoverFavg(:,Sb2Ix),2), nanstd(DataDFoverFavg(:,Sb2Ix),[],2)/sqrt(length(Sb2Ix)), 'g'); 
alignYaxes
title(['Success trials: ' num2str(length(Sb1Ix)) ' Visual; ' num2str(length(Sb2Ix)) ' Auditory']) 
xlim([-500 1500])



