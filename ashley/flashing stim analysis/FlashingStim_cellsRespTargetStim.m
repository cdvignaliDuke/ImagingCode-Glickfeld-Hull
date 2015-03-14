SubNum = '604';
date = '141022';
time = '1600';
ImgFolder = '002';
mouse = 'AW04';

% load MWorks file
CD = ['Z:\data\' mouse '\mworks\' date];
cd(CD);
mworks = ['data-' 'i' SubNum '-' date '-' time]; 
load (mworks);

% load dataTC
fileSave = fullfile('Z:\analysis\',mouse,'two-photon imaging', date, ImgFolder);
cd(fileSave);
load('DFoverF_60framesTargetOn.mat')

%variables from mworks
cLeverDown = cell2mat_padded(input.cLeverDown);
cTargetOn = cell2mat_padded(input.cTargetOn);
tCyclesOn = cell2mat_padded(input.tCyclesOn);
block2 = cell2mat_padded(input.tBlock2TrialNumber);
cycles = unique(tCyclesOn);
% cycTime = input.nFramesOn + input.nFramesOff;
frameRateS = 30; %hard-coded for now, but should be available in scanbox-yeti datasets' info file
RateFRperMS = frameRateS/1000;
cycTime = ceil((input.stimOnTimeMs+input.stimOffTimeMs)*RateFRperMS);
nTrials = input.trialSinceReset;

% special case where last trial take place within last 2 frames collected
% by mworks, but not scanbox
cLeverDown = cLeverDown(1:end-10,:);
cTargetOn = cTargetOn(1:end-10,:);
tCyclesOn = tCyclesOn(1:end-10,:);
block2 = block2(1:end-10,:);
nTrials = nTrials-10;

%plot cells
v_ind = find(block2==0);
a_ind = find(block2==1);

cell = 1;
figure;
for iplot = 1:16
    subplot(4,4,iplot)
    plot(mean(DFoverF(:,cell,v_ind),3),'g');
    hold on
    plot(mean(DFoverF(:,cell,a_ind),3),'r');
    hold on
    vline(30,'c');
    hold on
    vline(30-cycTime,'k:');
    hold on
    title(['Cell ' num2str(cell)]);
    hold on
    cell = cell+1;
end

% find cells that respond to target stimulus, visual and auditory
% conditions, signal compare to trace just before target comes on
preStimResp_V = zeros(size(v_ind,2),size(DFoverF,2));
for itrial =1:size(v_ind,1);
    for icell = 1:size(DFoverF,2)
        preStimResp_V(itrial,icell) = mean(DFoverF(28:32,icell,v_ind(itrial)),1);
    end
end

baselineStimResp_V = zeros(size(v_ind,2),size(DFoverF,2));
for itrial = 1:size(v_ind,1);
    for icell = 1:size(DFoverF,2)
        baselineStimResp_V(itrial,icell) = mean(DFoverF(36:40,icell,v_ind(itrial)),1);
    end
end

baselineStimRespTtest_V= ttest(preStimResp_V,baselineStimResp_V,'alpha', 0.01);
baselineStimRespIndex_V = find(baselineStimRespTtest_V == 1);

preStimResp_A = zeros(size(a_ind,2),size(DFoverF,2));
for itrial =1:size(a_ind,1);
    for icell = 1:size(DFoverF,2)
        preStimResp_A(itrial,icell) = mean(DFoverF(28:32,icell,a_ind(itrial)),1);
    end
end

baselineStimResp_A = zeros(size(v_ind,2),size(DFoverF,2));
for itrial = 1:size(a_ind,1);
    for icell = 1:size(DFoverF,2)
        baselineStimResp_A(itrial,icell) = mean(DFoverF(36:40,icell,a_ind(itrial)),1);
    end
end

baselineStimRespTtest_A= ttest(preStimResp_A,baselineStimResp_A,'alpha', 0.01);
baselineStimRespIndex_A = find(baselineStimRespTtest_A == 1);


%plot visually responsive cells
start = 1;
figure;
for iplot = 1:16
    cell = baselineStimRespIndex_V(start);
    subplot(4,4,iplot)
    plot(mean(DFoverF(:,cell,v_ind),3),'g');
    hold on
    plot(mean(DFoverF(:,cell,a_ind),3),'r');
    hold on
    vline(30,'k');
    hold on
    vline(30-cycTime,'k:');
    hold on
    title(['Cell ' num2str(cell)]);
    hold on
    start = start+1;
end
