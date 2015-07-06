SubNum = '608';
date = '150209';
time = '1744';
ImgFolder = '009';
mouse = 'AW08';

% load MWorks file
CD = ['Z:\data\' mouse '\mworks\' date];
cd(CD);
mworks = ['data-' 'i' SubNum '-' date '-' time]; 
load (mworks);

% load dataTC
fileSave = fullfile('Z:\analysis\',mouse,'two-photon imaging', date, ImgFolder);
cd(fileSave);
load('cycDataDFoverF_cmlvNoTarget.mat')

tCyclesOn = cell2mat_padded(input.tCyclesOn);
cycles = unique(tCyclesOn);
cycTime = input.nFramesOn + input.nFramesOff;
% frameRateS = 30; %hard-coded for now, but should be available in scanbox-yeti datasets' info file
% RateFRperMS = frameRateS/1000;
% cycTime = ceil((input.stimOnTimeMs+input.stimOffTimeMs)*RateFRperMS);
nTrials = input.trialSinceReset;

% special case where last trial take place within last 2 frames collected
% by mworks, but not scanbox

% tCyclesOn = tCyclesOn(1:end-10,:);
% nTrials = nTrials-10;

dataTrialStart = cycDataDFoverF_cmlvNoTarget{1};
v_ind = cycV_ind{1};
a_ind = cycA_ind{1};

%plot cells
cell = 1;
figure;
for iplot = 1:16
    subplot(4,4,iplot)
    plot(mean(dataTrialStart(:,cell,v_ind),3),'g');
    hold on
    plot(mean(dataTrialStart(:,cell,a_ind),3),'r');
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
    cell = cell+1;
end

% find cells that respond to first stimulus, visual and auditory conditions
preStimResp_V = zeros(size(v_ind,2),size(dataTrialStart,2));
for itrial =1:size(v_ind,1);
    for icell = 1:size(dataTrialStart,2)
        preStimResp_V(itrial,icell) = mean(dataTrialStart(1:30,icell,v_ind(itrial)),1);
    end
end

baselineStimResp_V = zeros(size(v_ind,2),size(dataTrialStart,2));
for itrial = 1:size(v_ind,1);
    for icell = 1:size(dataTrialStart,2)
        baselineStimResp_V(itrial,icell) = mean(dataTrialStart(36:40,icell,v_ind(itrial)),1);
    end
end

baselineStimRespTtest_V= ttest(preStimResp_V,baselineStimResp_V,'alpha', 0.01);
baselineStimRespIndex_V = find(baselineStimRespTtest_V == 1);

preStimResp_A = zeros(size(a_ind,2),size(dataTrialStart,2));
for itrial =1:size(a_ind,1);
    for icell = 1:size(dataTrialStart,2)
        preStimResp_A(itrial,icell) = mean(dataTrialStart(1:30,icell,a_ind(itrial)),1);
    end
end

baselineStimResp_A = zeros(size(v_ind,2),size(dataTrialStart,2));
for itrial = 1:size(a_ind,1);
    for icell = 1:size(dataTrialStart,2)
        baselineStimResp_A(itrial,icell) = mean(dataTrialStart(36:40,icell,a_ind(itrial)),1);
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
    plot(mean(dataTrialStart(:,cell,v_ind),3),'g');
    hold on
    plot(mean(dataTrialStart(:,cell,a_ind),3),'r');
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
    xlim([0 length(mean(dataTrialStart(:,cell,v_ind),3))]);
end

%plot cells with zero deg ori pref
dirFolder = '007';
CD = ['Z:\analysis\' mouse '\two-photon imaging\' date '\' dirFolder];
cd(CD);
load('oriTuningPreferences.mat')


fileSave = fullfile('Z:\analysis\',mouse,'two-photon imaging', date, ImgFolder);
cd(fileSave);

cellsPrefZero = find(oriPref_ind == 1);
cellsPrefNinety = find(oriPref_ind == 3);

start = 1;
figure;
for iplot = 1:20
    cell = cellsPrefZero(start);
    subplot(4,5,iplot)
    plot(mean(dataTrialStart(:,cell,v_ind),3),'g');
    hold on
    plot(mean(dataTrialStart(:,cell,a_ind),3),'r');
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
    xlim([0 length(mean(dataTrialStart(:,cell,v_ind),3))]);
end


