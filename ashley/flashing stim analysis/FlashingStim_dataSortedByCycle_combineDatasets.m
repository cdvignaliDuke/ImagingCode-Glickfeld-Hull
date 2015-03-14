%combine two datasets to have 3 trial types - vis only, aud only, and
%vis+aud
SubNum = '607';
mouse = 'AW07';
date = '141215';

%% first dataset (_1) - vis only (V) and aud only (A)
time = '1655';
ImgFolder = '005';

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
time = '1635';
ImgFolder = '004';

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
addInd_2 = size(dataTC_1,1);
cLeverDown = cat(1,cell2mat_padded(input1.cLeverDown),(cell2mat_padded(input2.cLeverDown)+addInd_2));
cTargetOn = cat(1,cell2mat_padded(input1.cTargetOn),(cell2mat_padded(input2.cTargetOn)+addInd_2));
tCyclesOn = cat(1,cell2mat_padded(input1.tCyclesOn),cell2mat_padded(input2.tCyclesOn));
cycles = unique(tCyclesOn);
cycTime = input1.nFramesOn + input1.nFramesOff;
% frameRateS = 30; %hard-coded for now, but should be available in scanbox-yeti datasets' info file
% RateFRperMS = frameRateS/1000;
% cycTime = ceil((input1.stimOnTimeMs+input1.stimOffTimeMs)*RateFRperMS);
nTrials = input1.trialSinceReset+input2.trialSinceReset;


block2 = cat(1,cell2mat_padded(input1.tBlock2TrialNumber),cell2mat_padded(input2.tBlock2TrialNumber));
V_ind = find(block2 == 0);
AorAV_ind = find(block2 == 1);
A_dataset_ind = AorAV_ind <= input1.trialSinceReset;
A_ind = AorAV_ind(A_dataset_ind == 1);
AV_ind = AorAV_ind(A_dataset_ind == 0);

%% new folder to save in
ImgFolder = '004+005'
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
    A_cycInd = find(ismember(trials,A_ind));
    AV_cycInd = find(ismember(trials,AV_ind));
    V_avg = mean(mean(DataDFoverF(:,:,V_cycInd),3),2);
    A_avg = mean(mean(DataDFoverF(:,:,A_cycInd),3),2);
    AV_avg = mean(mean(DataDFoverF(:,:,AV_cycInd),3),2);
    subplot(3,3,icyc);
    plot(V_avg,'g');
    hold on
    plot(A_avg,'r');
    hold on
    plot(AV_avg,'m');
    hold on
    title([num2str(cycles(icyc)) ' cycles; ' num2str(length(V_ind)) ' visual trials; ' num2str(length(A_ind)) ' auditory trials']);
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
    A_indAll = [];
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
        A_cycInd = find(ismember(trials,A_ind));
        AV_cycInd = find(ismember(trials,AV_ind));
        V_indAll = cat(1,V_indAll, V_cycInd+running_ind);
        A_indAll = cat(1,A_indAll, A_cycInd+running_ind);
        AV_indAll = cat(1,AV_indAll,AV_cycInd+running_ind);
        running_ind = length(trials)+running_ind;
    end        
    cycDataDFoverF_cmlvNoTarget{icyc} = dataDFoverF_cmlvNoTarget; 
    cycV_ind{icyc} = V_indAll;
    cycA_ind{icyc} = A_indAll;
    cycAV_ind{icyc} = AV_indAll;
end

figure;
for icyc = 1:length(cycles)
    dataDFoverF = cycDataDFoverF_cmlvNoTarget{icyc};
    V_cycInd = cycV_ind{icyc};
    A_cycInd = cycA_ind{icyc};
    AV_cycInd = cycAV_ind{icyc};
    V_avg = mean(mean(dataDFoverF(:,:,V_cycInd),3),2);
    A_avg = mean(mean(dataDFoverF(:,:,A_cycInd),3),2);
    AV_avg = mean(mean(dataDFoverF(:,:,AV_cycInd),3),2);
    subplot(3,3,icyc);
    plot(V_avg,'g');
    hold on
    plot(A_avg,'r');
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
    title([num2str(cycles(icyc)) ' cycles; ' num2str(length(V_ind)) ' visual trials; ' num2str(length(A_ind)) ' auditory trials'])
    hold on
end


fileSave = fullfile('Z:\analysis\',mouse,'two-photon imaging', date, ImgFolder);
cd(fileSave);
% save('cycDataDFoverF_cmlvNoTarget.mat', 'cycDataDFoverF_cmlvNoTarget', 'cycA_ind', 'cycV_ind');

