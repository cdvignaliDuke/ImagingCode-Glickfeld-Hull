
SubNum = '608';
mouse = 'AW08';
date = '150217';
time = '1620';
ImgFolder = '005';

% load MWorks file
CD = ['Z:\data\' mouse '\mworks\' date];
cd(CD);
mworks = ['data-' 'i' SubNum '-' date '-' time]; 
load (mworks);

%load timecourse
fileSave = fullfile('Z:\analysis\',mouse,'two-photon imaging', date, ImgFolder);
cd(fileSave);
load('dataTC.mat');
dataTC = dataTimecourse.dataTC;

%variables from mworks
cLeverDown = cell2mat_padded(input.cLeverDown);
cTargetOn = cell2mat_padded(input.cTargetOn);
tCyclesOn = cell2mat_padded(input.tCyclesOn);
block2 = cell2mat_padded(input.tBlock2TrialNumber);

cycles = unique(tCyclesOn);
cycTime = input.nFramesOn + input.nFramesOff;
% frameRateS = 30; %hard-coded for now, but should be available in scanbox-yeti datasets' info file
% RateFRperMS = frameRateS/1000;
% cycTime = ceil((input.stimOnTimeMs+input.stimOffTimeMs)*RateFRperMS);
nTrials = input.trialSinceReset;

% % special case where last trial take place within last 2 frames collected
% % by mworks, but not scanbox
% cLeverDown = cLeverDown(1:end-2,:);
% cTargetOn = cTargetOn(1:end-2,:);
% tCyclesOn = tCyclesOn(1:end-2,:);
% block2 = block2(1:end-2,:);
% nTrials = nTrials-2;

%dFoverF
Data = zeros(60,size(dataTC,2),nTrials);
DataDF = zeros(60,size(dataTC,2),nTrials);
DFoverF = zeros(60,size(dataTC,2),nTrials);
for itrial = 1:nTrials
    Data(:,:,itrial) = dataTC(cLeverDown(itrial)-30:cLeverDown(itrial)+29,:);
    DataDF(:,:,itrial) = bsxfun(@minus,Data(:,:,itrial),mean(Data(1:10,:,itrial),1));
    DFoverF(:,:,itrial) = bsxfun(@rdivide,DataDF(:,:,itrial),mean(Data(1:10,:,itrial),1));
end

%plot average all trials all cells
visualIndex = find(block2==0);
auditoryIndex = find(block2==1);
DFoverFmeanV = mean(mean(DFoverF(:,:,visualIndex),3),2);
DFoverFmeanA = mean(mean(DFoverF(:,:,auditoryIndex),3),2);
figure;
plot(DFoverFmeanV,'g');
hold on
plot(DFoverFmeanA,'r');
hold on
vline(30,'k');

%plot cells
cell = 1;
figure;
for iplot = 1:16
    subplot(4,4,iplot)
    plot(mean(DFoverF(:,cell,visualIndex),3),'g');
    hold on
    plot(mean(DFoverF(:,cell,auditoryIndex),3),'r');
    hold on
    vline(30,'k');
    cell = cell+1;
end

%ttest responsive cells
preStimResp_V = zeros(size(visualIndex,2),size(DFoverF,2));
for itrial =1:size(visualIndex,2);
    for icell = 1:size(DFoverF,2)
        preStimResp_V(itrial,icell) = mean(DFoverF(1:30,icell,visualIndex(itrial)),1);
    end
end

baselineStimResp_V = zeros(size(visualIndex,2),size(DFoverF,2));
for itrial = 1:size(visualIndex,2);
    for icell = 1:size(DFoverF,2)
        baselineStimResp_V(itrial,icell) = mean(DFoverF(36:40,icell,visualIndex(itrial)),1);
    end
end

baselineStimRespTtest_V= ttest(preStimResp_V,baselineStimResp_V,'alpha', 0.01);
baselineStimRespIndex_V = find(baselineStimRespTtest_V == 1);

preStimResp_A = zeros(size(auditoryIndex,2),size(DFoverF,2));
for itrial =1:size(auditoryIndex,2);
    for icell = 1:size(DFoverF,2)
        preStimResp_A(itrial,icell) = mean(DFoverF(1:30,icell,auditoryIndex(itrial)),1);
    end
end

baselineStimResp_A = zeros(size(auditoryIndex,2),size(DFoverF,2));
for itrial = 1:size(auditoryIndex,2);
    for icell = 1:size(DFoverF,2)
        baselineStimResp_A(itrial,icell) = mean(DFoverF(31:end,icell,auditoryIndex(itrial)),1);
    end
end

baselineStimRespTtest_A = ttest(preStimResp_A,baselineStimResp_A);
baselineStimRespIndex_A = find(baselineStimRespTtest_A == 1);

baselineStimRespIndex_both = intersect(baselineStimRespIndex_A,baselineStimRespIndex_V);
baselineStimRespIndex_Aonly = setdiff(baselineStimRespIndex_A,baselineStimRespIndex_both);
baselineStimRespIndex_Vonly = setdiff(baselineStimRespIndex_V,baselineStimRespIndex_both);

%plot responsive cells
cell = 1;
figure;
for iplot = 1:16
    plotcell = baselineStimRespIndex_V(cell);
    subplot(4,4,iplot)
    plot(mean(DFoverF(:,plotcell,visualIndex),3),'g');
    hold on
    plot(mean(DFoverF(:,plotcell,auditoryIndex),3),'r');
    hold on
    vline(30,'k');
    cell = cell+1
end    

% use mean instead of std and ste
figure;
for iplot = 1:16
    subplot(4,4,iplot)
    errbar_V = std(std(DFoverF(:,baselineStimRespIndex_V,visualIndex),[], 2),[],3);
    shadedErrorBar([],mean(mean(DFoverF(:,baselineStimRespIndex_V,visualIndex),2),3),errbar_V, 'g', 1);
    hold on
    errbar_A = std(std(DFoverF(:,baselineStimRespIndex_V,auditoryIndex),[], 2),[],3);
    shadedErrorBar([],mean(mean(DFoverF(:,baselineStimRespIndex_V,auditoryIndex),2),3),errbar_A,'r', 1);
    
    hold on
    vline(30,'k');
end 


save('DFoverF_60framesTrialStart.mat', 'DFoverF');

