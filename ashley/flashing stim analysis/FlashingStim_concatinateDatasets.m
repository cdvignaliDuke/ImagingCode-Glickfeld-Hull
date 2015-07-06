%concatinate datasets
SubNum = '004';
mouse = 'AW04';
date = '140923';

%% first dataset (_1) 
time = '1509';
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
%% second dataset (_2) 
time = '1527';
ImgFolder = '007';

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
block2 = cat(1,cell2mat_padded(input1.tBlock2TrialNumber),cell2mat_padded(input2.tBlock2TrialNumber));
cycles = unique(tCyclesOn);
% cycTime = input.nFramesOn + input.nFramesOff;
frameRateS = 30; %hard-coded for now, but should be available in scanbox-yeti datasets' info file
RateFRperMS = frameRateS/1000;
cycTime = ceil((input1.stimOnTimeMs+input1.stimOffTimeMs)*RateFRperMS);
nTrials = input1.trialSinceReset+input2.trialSinceReset;

V_ind = find(block2 ==0);
AV_ind = find(block2 ==1);

%% dF/F
Data = zeros(60,size(dataTC,2),nTrials);
DataDF = zeros(60,size(dataTC,2),nTrials);
DFoverF = zeros(60,size(dataTC,2),nTrials);
for itrial = 1:nTrials
    Data(:,:,itrial) = dataTC(cLeverDown(itrial)-30:cLeverDown(itrial)+29,:);
    DataDF(:,:,itrial) = bsxfun(@minus,Data(:,:,itrial),mean(Data(1:30,:,itrial),1));
    DFoverF(:,:,itrial) = bsxfun(@rdivide,DataDF(:,:,itrial),mean(Data(1:30,:,itrial),1));
end

DFoverFmeanV = mean(mean(DFoverF(:,:,V_ind),3),2);
DFoverFmeanAV = mean(mean(DFoverF(:,:,AV_ind),3),2);

figure;
plot(DFoverFmeanV,'g');
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
errbar_AV = std(mean(DFoverF(:,:,AV_ind),2),[],3)/sqrt(size(DFoverF(:,:,AV_ind),3));
shadedErrorBar([],DFoverFmeanV,errbar_V, 'g', 1);
shadedErrorBar([],DFoverFmeanAV,errbar_AV, 'm', 1);
hold on
title('Avg all trials, all cells; 60 frames around trial start');
ylabel('dF/F');
xlabel('frames');

%% cyc sorting
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
    V_ind = find(block2(trials) == 0);
    AV_ind = find(block2(trials) == 1);
    V_avg = mean(mean(DataDFoverF(:,:,V_ind),3),2);
    AV_avg = mean(mean(DataDFoverF(:,:,AV_ind),3),2);
    subplot(3,3,icyc);
    plot(V_avg,'g');
    hold on
    plot(AV_avg,'m');
    hold on
    title([num2str(cycles(icyc)) ' cycles; ' num2str(length(V_ind)) ' visual trials; ' num2str(length(AV_ind)) ' auditory trials']);
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

% plot cumulative average for each trial length (# of cycles), average
% up to, but not including the target

L = zeros(size(cycles));
for icyc = 1:length(cycles)
    dataDFoverF_cmlvNoTarget = [];
    V_indAll = [];
    A_indAll = [];
    running_ind = 0;
    L = ceil(30+ (cycles(icyc))*cycTime);
    cyc_ind = icyc:length(cycles);
    for i = cyc_ind
        dataDFoverF = cycDataDFoverF{i};
        dataDFoverF_NoTarget = zeros(L,size(dataDFoverF,2),size(dataDFoverF,3));
        dataDFoverF_NoTarget = dataDFoverF(1:L,:,:);
        dataDFoverF_cmlvNoTarget = cat(3,dataDFoverF_cmlvNoTarget,dataDFoverF_NoTarget);
        trials = find(tCyclesOn == cycles(i));
        V_ind = find(block2(trials) == 0);
        AV_ind = find(block2(trials) == 1);
        V_indAll = cat(1,V_indAll, V_ind+running_ind);
        A_indAll = cat(1,A_indAll, AV_ind+running_ind);
        running_ind = length(trials)+running_ind;
    end        
    cycDataDFoverF_cmlvNoTarget{icyc} = dataDFoverF_cmlvNoTarget; 
    cycV_ind{icyc} = V_indAll;
    cycA_ind{icyc} = A_indAll;
end

figure;
for icyc = 1:length(cycles)
    dataDFoverF = cycDataDFoverF_cmlvNoTarget{icyc};
    V_ind = cycV_ind{icyc};
    AV_ind = cycA_ind{icyc};
    V_avg = mean(mean(dataDFoverF(:,:,V_ind),3),2);
    AV_avg = mean(mean(dataDFoverF(:,:,AV_ind),3),2);
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
    title([num2str(cycles(icyc)) ' cycles; ' num2str(length(V_ind)) ' visual trials; ' num2str(length(AV_ind)) ' auditory trials'])
    hold on
end

%% plot cells
preStimResp_V = zeros(size(V_ind,2),size(DFoverF,2));
for itrial =1:size(V_ind,1);
    for icell = 1:size(DFoverF,2)
        preStimResp_V(itrial,icell) = mean(DFoverF(1:30,icell,V_ind(itrial)),1);
    end
end

baselineStimResp_V = zeros(size(V_ind,2),size(DFoverF,2));
for itrial = 1:size(V_ind,1);
    for icell = 1:size(DFoverF,2)
        baselineStimResp_V(itrial,icell) = mean(DFoverF(36:40,icell,V_ind(itrial)),1);
    end
end

baselineStimRespTtest_V= ttest(preStimResp_V,baselineStimResp_V,'alpha', 0.01);
baselineStimRespIndex_V = find(baselineStimRespTtest_V == 1);

preStimResp_A = zeros(size(AV_ind,2),size(DFoverF,2));
for itrial =1:size(AV_ind,1);
    for icell = 1:size(DFoverF,2)
        preStimResp_A(itrial,icell) = mean(DFoverF(1:30,icell,AV_ind(itrial)),1);
    end
end

baselineStimResp_A = zeros(size(AV_ind,2),size(DFoverF,2));
for itrial = 1:size(AV_ind,1);
    for icell = 1:size(DFoverF,2)
        baselineStimResp_A(itrial,icell) = mean(DFoverF(36:40,icell,AV_ind(itrial)),1);
    end
end

baselineStimRespTtest_A= ttest(preStimResp_A,baselineStimResp_A,'alpha', 0.01);
baselineStimRespIndex_A = find(baselineStimRespTtest_A == 1);


%plot visually responsive cells
start = 1;
figure;
for iplot = 1:20
    cell = baselineStimRespIndex_V(start);
    subplot(4,5,iplot)
    plot(mean(DFoverF(:,cell,V_ind),3),'g');
    hold on
    plot(mean(DFoverF(:,cell,AV_ind),3),'r');
    hold on
    vline(30,'k');
    hold on
    for i = 1:cycles(1)-1
        L = (i*cycTime)+30;
        vline(L,'k:');
        hold on
    end
    vline((cycles(1)*cycTime+30),'c');
    hold on
    title(['Cell ' num2str(cell)]);
    hold on
    start = start+1;
    xlim([0 length(mean(DFoverF(:,cell,AV_ind),3))]);
end

