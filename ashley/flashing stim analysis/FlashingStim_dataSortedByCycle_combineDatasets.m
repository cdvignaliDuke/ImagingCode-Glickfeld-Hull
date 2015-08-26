%combine two datasets to have 3 trial types - vis only, aud only, and
%vis+aud
SubNum = '614';
mouse = 'AW14';
date = '150623';

%% first dataset (_1) - vis only (V) and aud only (A)
time = '1142';
ImgFolder = '002';

% load MWorks file
CD = ['Z:\data\' mouse '\mworks\' date];
cd(CD);
mworks = ['data-' 'i' SubNum '-' date '-' time]; 
load (mworks);
input1 = input;

%load timecourse
fileSave = fullfile('Z:\analysis\',mouse,'two-photon imaging', date, ImgFolder);
cd(fileSave);
% load('dataTC.mat');
load('Timecourses.mat')
dataTC_1 = dataTimecourse.dataTCsub;
% dataTC_1 = dataTimecourse.dataTC;

clear input dataTimecourse
%% second dataset (_2) - vis only (V) and vis + aud only (AV)
time = '1158';
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
% load('dataTC.mat');
load('Timecourses.mat')
dataTC_2 = dataTimecourse.dataTCsub;
% dataTC_2 = dataTimecourse.dataTC;

clear input dataTimecourse
%% third dataset (_3) 
time = '1217';
ImgFolder = '004';

% load MWorks file
CD = ['Z:\data\' mouse '\mworks\' date];
cd(CD);
mworks = ['data-' 'i' SubNum '-' date '-' time]; 
load (mworks);
input3 = input;

%load timecourse
fileSave = fullfile('Z:\analysis\',mouse,'two-photon imaging', date, ImgFolder);
cd(fileSave);
% load('dataTC.mat');
load('Timecourses.mat')
dataTC_3 = dataTimecourse.dataTCsub;

clear input dataTimecourse
%% combine dataTCs
dataTC = cat(1,dataTC_1,dataTC_2,dataTC_3);
% dataTC = cat(1,dataTC_1,dataTC_2);
%% mworks variables - 2 datasets
% addInd_2 = size(dataTC_1,1);
% cLeverDown = cat(1,cell2mat_padded(input1.cLeverDown),(cell2mat_padded(input2.cLeverDown)+addInd_2));
% cTargetOn = cat(1,cell2mat_padded(input1.cTargetOn),(cell2mat_padded(input2.cTargetOn)+addInd_2));
% tCyclesOn = cat(1,cell2mat_padded(input1.tCyclesOn),cell2mat_padded(input2.tCyclesOn));
% cycles = unique(tCyclesOn);
% cycTime = input1.nFramesOn + input1.nFramesOff;
% frameRateS = 30; %hard-coded for now, but should be available in scanbox-yeti datasets' info file
% RateFRperMS = frameRateS/1000;
% % cycTime = ceil((input1.stimOnTimeMs+input1.stimOffTimeMs)*RateFRperMS);
% nTrials = input1.trialSinceReset+input2.trialSinceReset;
% trialOutcome = cat(2,input1.trialOutcomeCell,input2.trialOutcomeCell);
% DirectionDeg = cell2mat_padded(cat(2,input1.gratingDirectionDeg,input2.gratingDirectionDeg));
% tooFastTime = input1.tooFastTimeMs*RateFRperMS;
% maxReactTime = input1.reactTimeMs*RateFRperMS;
% reactTimesMs = cell2mat_padded(cat(2,input1.reactTimesMs,input2.reactTimesMs));
% successReactTime = reactTimesMs(strcmp(trialOutcome,'success'))*RateFRperMS;
% meanSuccessReactTime = mean(successReactTime,1);
% 
% 
% block2 = cat(1,cell2mat_padded(input1.tBlock2TrialNumber),cell2mat_padded(input2.tBlock2TrialNumber));
% V_ind = find(block2 == 0);
% % AorAV_ind = find(block2 == 1);
% % A_dataset_ind = AorAV_ind <= input1.trialSinceReset;
% % A_ind = AorAV_ind(A_dataset_ind == 1);
% % AV_ind = AorAV_ind(A_dataset_ind == 0);
% AV_ind = find(block2 == 1);

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

%% analysis

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

figure;
for icyc = 1:length(cycles)
    DataDFoverF = cycDataDFoverF{icyc};
    trials = find(tCyclesOn == cycles(icyc));
    V_cycInd = find(ismember(trials,V_ind));
%     A_cycInd = find(ismember(trials,A_ind));
    AV_cycInd = find(ismember(trials,AV_ind));
    V_avg = mean(mean(DataDFoverF(:,:,V_cycInd),3),2);
%     A_avg = mean(mean(DataDFoverF(:,:,A_cycInd),3),2);
    AV_avg = mean(mean(DataDFoverF(:,:,AV_cycInd),3),2);
    subplot(3,4,icyc);
    plot(V_avg,'g');
    hold on
%     plot(A_avg,'r');
    hold on
    plot(AV_avg,'m');
    hold on
%     title([num2str(cycles(icyc)) ' cycles; ' num2str(length(V_ind)) ' visual trials; ' num2str(length(A_ind)) ' auditory trials']);
    hold on
    vline(30,'k');
    hold on
    for i = 1:cycles(icyc)-1
        L = (i*cycTime)+30;
        vline(L,'k:');
        hold on
    end
    vline((cycles(icyc)*cycTime+30),'c');
    hold on
end

L = zeros(size(cycles));
for icyc = 1:length(cycles)
    dataDFoverF_cmlvNoTarget = [];
    V_indAll = [];
%     A_indAll = [];
    AV_indAll = [];
    running_ind = 0;
    L = ceil(30+ (cycles(icyc))*cycTime);
    cyc_ind = icyc:length(cycles);
    for i = cyc_ind
        dataDFoverF = cycDataDFoverF{i};
        dataDFoverF_NoTarget = zeros(L,size(dataDFoverF,2),size(dataDFoverF,3));
        dataDFoverF_NoTarget = dataDFoverF(1:L,:,:);
        dataDFoverF_cmlvNoTarget = cat(3,dataDFoverF_cmlvNoTarget,dataDFoverF_NoTarget);
        trials = find(tCyclesOn == cycles(i));
        V_cycInd = find(ismember(trials,V_ind));
%         A_cycInd = find(ismember(trials,A_ind));
        AV_cycInd = find(ismember(trials,AV_ind));
        V_indAll = cat(1,V_indAll, V_cycInd+running_ind);
%         A_indAll = cat(1,A_indAll, A_cycInd+running_ind);
        AV_indAll = cat(1,AV_indAll,AV_cycInd+running_ind);
        running_ind = length(trials)+running_ind;
    end        
    cycDataDFoverF_cmlvNoTarget{icyc} = dataDFoverF_cmlvNoTarget; 
    cycV_ind{icyc} = V_indAll;
%     cycA_ind{icyc} = A_indAll;
    cycAV_ind{icyc} = AV_indAll;
end

for icyc = 1:length(cycles)
    ind = sort(cat(1,cycV_ind{icyc},cycAV_ind{icyc}));   
    cycTrialOutcome{icyc} = trialOutcome(ind);
    cycDirectionDeg{icyc} = DirectionDeg(ind);
    if sum(catchIndex) > 0
    cycCatchDirectionDeg{icyc} = catchDirectionDeg(ind);
    cycCatchTrialOutcome{icyc} = catchTrialOutcome(ind);
    cycCatchCycle{icyc} = catchCycle(ind);
    end
end

% legendinfo = {'vis only','aud only','vis+aud'};
figure;
for icyc = 1:length(cycles)
    dataDFoverF = cycDataDFoverF_cmlvNoTarget{icyc};
    V_cycInd = cycV_ind{icyc};
%     A_cycInd = cycA_ind{icyc};
    AV_cycInd = cycAV_ind{icyc};
    V_avg = mean(mean(dataDFoverF(:,:,V_cycInd),3),2);
%     A_avg = mean(mean(dataDFoverF(:,:,A_cycInd),3),2);
    AV_avg = mean(mean(dataDFoverF(:,:,AV_cycInd),3),2);
    errbar_V = (std(mean(dataDFoverF(:,:,V_cycInd),2),[],3))/(sqrt(size(dataDFoverF(:,:,V_cycInd),3)));
%     errbar_A = (std(mean(dataDFoverF(:,:,A_cycInd),2),[],3))/(sqrt(size(dataDFoverF(:,:,A_cycInd),3)));
    errbar_AV = (std(mean(dataDFoverF(:,:,AV_cycInd),2),[],3))/(sqrt(size(dataDFoverF(:,:,AV_cycInd),3)));
    subplot(3,4,icyc);
%     plot(V_avg(20:end,:),'g');
%     hold on
%     plot(A_avg(20:end,:),'r');
%     hold on
%     plot(AV_avg(20:end,:),'m');
%     hold on
    errorbar(V_avg(20:end,:),errbar_V(20:end,:),'g')
    hold on
%     errorbar(A_avg(20:end,:),errbar_A(20:end,:),'b')
    hold on
    errorbar(AV_avg(20:end,:),errbar_AV(20:end,:),'k')
    hold on
    vline(10,'k')
    hold on
    for i = 1:cycles(icyc)-1
        L = (i*cycTime)+11;
        vline(L,'k:');
        hold on
    end
    vline((cycles(icyc)*cycTime+11),'c');
    hold on
    if icyc == 1
        title([num2str(size(dataDFoverF,2)) ' cells'])
    else
%     title([num2str(length(V_cycInd)) ' visual trials; ' num2str(length(A_cycInd)) ' auditory trials; ' num2str(length(AV_cycInd)) ' vis+aud trials'])
    title([num2str(length(V_cycInd)) ' visual trials; ' num2str(length(AV_cycInd)) ' vis+aud trials'])
    end
    hold on
    xlim([0 length(V_avg(20:end,:))+5])
    ylim([-0.05 0.05])
%     if icyc == 1
%         legend(legendinfo,'Location','SouthEast')
%     end
end


% fileSave = fullfile('Z:\analysis\',mouse,'two-photon imaging', date, ImgFolder);
% cd(fileSave);
% save('cycDataDFoverF_cmlvNoTarget.mat', 'cycDataDFoverF_cmlvNoTarget', 'cycA_ind', 'cycV_ind');

