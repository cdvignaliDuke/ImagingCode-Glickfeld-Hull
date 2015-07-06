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

%load timecourse
fileSave = fullfile('Z:\analysis\',mouse,'two-photon imaging', date, ImgFolder);
cd(fileSave);
load('dataTC.mat');
dataTC = dataTimecourse.dataTC;

%% mworks variables

cLeverDown = cell2mat_padded(input.cLeverDown);
cTargetOn = cell2mat_padded(input.cTargetOn);
tCyclesOn = cell2mat_padded(input.tCyclesOn);
cycles = unique(tCyclesOn);
% cycTime = input1.nFramesOn + input1.nFramesOff;
frameRateS = 30; %hard-coded for now, but should be available in scanbox-yeti datasets' info file
RateFRperMS = frameRateS/1000;
cycTime = ceil((input.stimOnTimeMs+input.stimOffTimeMs)*RateFRperMS);
nTrials = input.trialSinceReset;
trialOutcome = input.trialOutcomeCell;
DirectionDeg = cell2mat_padded(input.gratingDirectionDeg);
tooFastTime = input.tooFastTimeMs*RateFRperMS;
maxReactTime = input.reactTimeMs*RateFRperMS;
reactTimesMs = cell2mat_padded(input.reactTimesMs);
successReactTime = reactTimesMs(strcmp(trialOutcome,'success'))*RateFRperMS;
meanSuccessReactTime = mean(successReactTime,1);

block2 = cell2mat_padded(input.tBlock2TrialNumber);
V_ind = find(block2 == 0);
AV_ind = find(block2 == 1);

Success_ind = find(strcmp('success',trialOutcome));
Miss_ind = find(strcmp('ignore',trialOutcome));
Early_ind = find(strcmp('failure',trialOutcome));

Dirs = unique(DirectionDeg);

L = 60;
data_aroundTarget = zeros(L,size(dataTC,2),length(Success_ind));
dF_aroundTarget = zeros(L,size(dataTC,2),length(Success_ind)); 
dFoverF_aroundTarget = zeros(L,size(dataTC,2),length(Success_ind));
for i = 1:length(Success_ind)
    trial = Success_ind(i);
    data_aroundTarget(:,:,i) = dataTC(cTargetOn(trial)-(L/2):cTargetOn(trial)+((L/2)-1),:);
    dF_aroundTarget(:,:,i) = bsxfun(@minus, data_aroundTarget(:,:,i), mean(data_aroundTarget((L/2)-5:(L/2),:,i),1));
    dFoverF_aroundTarget(:,:,i) = bsxfun(@rdivide, dF_aroundTarget(:,:,i),mean(data_aroundTarget((L/2)-5:(L/2),:,i),1));
end

cellsAvg = squeeze(mean(dFoverF_aroundTarget,2));

V_Success_ind = find(ismember(Success_ind,intersect(Success_ind,V_ind)));
AV_Success_ind = find(ismember(Success_ind,intersect(Success_ind,AV_ind)));

V_successAvg = mean(mean(dFoverF_aroundTarget(:,:,V_Success_ind),3),2);
AV_successAvg = mean(mean(dFoverF_aroundTarget(:,:,AV_Success_ind),3),2);

start = 1;
figure;
for i = 1:18
    subplot(4,5,i)
    if ismember(start,V_Success_ind) == 1
        plot(cellsAvg(:,start),'g')
    elseif ismember(start,AV_Success_ind) == 1
        plot(cellsAvg(:,start), 'r')
    end
    hold on
    plot(V_successAvg,'g:')
    hold on
    plot(AV_successAvg,'r:')
    hold on
    vline((L/2),'c')
    hold on
    for i = 1:2
        vline((L/2)-(cycTime*i),'k:');
    end
    hold on
    vline((L/2)+tooFastTime,'k')
    vline((L/2) + maxReactTime,'k')
    vline((L/2) + meanSuccessReactTime, 'b')
%     if strcmp(trialOutcome(1,start),'success') == 1
%         set(subplot(5,5,i),'color',[0.9 0.9 0.9])
%     end
    xlim([0 length(cellsAvg(:,start))])
    title(['Success Trial ' num2str(start)])
    start = start+1;
end
