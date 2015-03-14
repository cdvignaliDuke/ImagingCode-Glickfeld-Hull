%combine two datasets to have 3 trial types - vis only, aud only, and
%vis+aud
SubNum = '516';
mouse = '516';
date = '141003';

%% first dataset (_1) 
time = '1156';
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
load('dataTC.mat');
dataTC_1 = dataTimecourse.dataTC;

clear input dataTimecourse
%% second dataset (_2)
time = '1159';
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
%% third dataset (_3)
time = '1159';
ImgFolder = '003';

% load MWorks file
CD = ['Z:\data\' mouse '\mworks\' date];
cd(CD);
mworks = ['data-' 'i' SubNum '-' date '-' time]; 
load (mworks);
input3 = input;

%load timecourse
fileSave = fullfile('Z:\analysis\',mouse,'two-photon imaging', date, ImgFolder);
cd(fileSave);
load('dataTC.mat');
dataTC_3 = dataTimecourse.dataTC;

clear input dataTimecourse

%% combine dataTCs
dataTC = cat(1,dataTC_1,dataTC_2,dataTC_3);
%% mworks variables
addInd_2 = size(dataTC_1,1);
addInd_3 = size(dataTC_2,2)+addInd_2;
cLeverDown = cat(1,cell2mat_padded(input1.cLeverDown),(cell2mat_padded(input2.cLeverDown)+addInd_2),(cell2mat_padded(input3.cLeverDown))+addInd_3);
cTargetOn = cat(1,cell2mat_padded(input1.cTargetOn),cell2mat_padded(input2.cTargetOn)+addInd_2,cell2mat_padded(input2.cTargetOn)+addInd_3);
tCyclesOn = cat(1,cell2mat_padded(input1.tCyclesOn),cell2mat_padded(input2.tCyclesOn),cell2mat_padded(input3.tCyclesOn));
cycles = unique(tCyclesOn);
% cycTime = input1.nFramesOn + input1.nFramesOff;
frameRateS = 30; %hard-coded for now, but should be available in scanbox-yeti datasets' info file
RateFRperMS = frameRateS/1000;
cycTime = ceil((input1.stimOnTimeMs+input1.stimOffTimeMs)*RateFRperMS);
nTrials = input1.trialSinceReset+input2.trialSinceReset+input3.trialSinceReset;
trialOutcome = cat(2,input1.trialOutcomeCell,input2.trialOutcomeCell,input3.trialOutcomeCell);
DirectionDeg = cell2mat_padded(cat(2,input1.gratingDirectionDeg,input2.gratingDirectionDeg,input3.gratingDirectionDeg));
tooFastTime = input1.tooFastTimeMs*RateFRperMS;
maxReactTime = input1.reactTimeMs*RateFRperMS;
reactTimesMs = cell2mat_padded(cat(2,input1.reactTimesMs,input2.reactTimesMs,input3.reactTimesMs));
successReactTime = reactTimesMs(strcmp(trialOutcome,'success'))*RateFRperMS;
meanSuccessReactTime = mean(successReactTime,1);

block2 = cat(1,cell2mat_padded(input1.tBlock2TrialNumber),cell2mat_padded(input2.tBlock2TrialNumber),cell2mat_padded(input3.tBlock2TrialNumber));
V_ind = find(block2 == 0);
AV_ind = find(block2 == 1);

%% new folder to save in
ImgFolder = '002+003+004'
try
    fileSave = fullfile('Z:\analysis\',mouse,'two-photon imaging', date, ImgFolder);
    cd(fileSave);
catch
    fileSave = fullfile('Z:\analysis\',mouse,'two-photon imaging', date);
    cd(fileSave);
    mkdir(ImgFolder)
    fileSave = fullfile('Z:\analysis\',mouse,'two-photon imaging', date, ImgFolder);
    cd(fileSave);
end

% save('dataTC.mat','dataTC');

%% analysis

for icyc = 1:length(cycles)
    ind = find(tCyclesOn == cycles(icyc));
    Data = zeros(cycTime.*(cycles(icyc)+1)+60,size(dataTC,2),length(ind));
    DataDF = zeros(cycTime.*(cycles(icyc)+1)+60,size(dataTC,2),length(ind));
    DataDFoverF = zeros(cycTime.*(cycles(icyc)+1)+60,size(dataTC,2),length(ind));
    for itrial = 1:length(ind)
        Data(:,:,itrial) = dataTC(cLeverDown(ind(itrial))-30:cLeverDown(ind(itrial))+29+(cycTime.*(cycles(icyc)+1)),:);
        DataDF(:,:,itrial) = bsxfun(@minus, Data(:,:,itrial), mean(Data(1:30,:,itrial),1));
        DataDFoverF(:,:,itrial) = bsxfun(@rdivide, DataDF(:,:,itrial), mean(Data(1:30,:,itrial),1));
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
    AV_cycInd = find(ismember(trials,AV_ind));
    V_avg = mean(mean(DataDFoverF(:,:,V_cycInd),3),2);
    AV_avg = mean(mean(DataDFoverF(:,:,AV_cycInd),3),2);
    subplot(3,3,icyc);
    plot(V_avg,'g');
    hold on
    plot(AV_avg,'m');
    hold on
    title([num2str(cycles(icyc)) ' cycles; ' num2str(length(V_cycInd)) ' visual trials; ' num2str(length(AV_cycInd)) ' auditory trials']);
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
        AV_cycInd = find(ismember(trials,AV_ind));
        V_indAll = cat(1,V_indAll, V_cycInd+running_ind);
        AV_indAll = cat(1,AV_indAll,AV_cycInd+running_ind);
        running_ind = length(trials)+running_ind;
    end        
    cycDataDFoverF_cmlvNoTarget{icyc} = dataDFoverF_cmlvNoTarget; 
    cycV_ind{icyc} = V_indAll;
    cycAV_ind{icyc} = AV_indAll;
end

figure;
for icyc = 1:length(cycles)
    dataDFoverF = cycDataDFoverF_cmlvNoTarget{icyc};
    V_cycInd = cycV_ind{icyc};
    AV_cycInd = cycAV_ind{icyc};
    V_avg = mean(mean(dataDFoverF(:,:,V_cycInd),3),2);
    AV_avg = mean(mean(dataDFoverF(:,:,AV_cycInd),3),2);
    subplot(3,3,icyc);
    plot(V_avg,'g');
    hold on
    plot(AV_avg,'m');
    hold on
    vline(30,'k')
    hold on
    for i = 1:cycles(icyc)-1
        L = (i*cycTime)+30;
        vline(L,'k:');
        hold on
    end
    vline((cycles(icyc)*cycTime+30),'c');
    hold on
    title([num2str(cycles(icyc)) ' cycles; ' num2str(length(V_cycInd)) ' visual trials; ' num2str(length(AV_cycInd)) ' auditory trials'])
    hold on
end


fileSave = fullfile('Z:\analysis\',mouse,'two-photon imaging', date, ImgFolder);
cd(fileSave);
% save('cycDataDFoverF_cmlvNoTarget.mat', 'cycDataDFoverF_cmlvNoTarget', 'cycAV_ind', 'cycV_ind');

