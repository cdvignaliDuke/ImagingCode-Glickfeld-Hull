%load dataset
SubNum = '510';
date = '140903';
time = '1257';
ImgFolder = '003';
mouse = '510';
fName = '003_000_000';
experiment = 'Flashing Stim';

% load MWorks file
CD = ['Z:\2P imaging\MWorks\' mouse '\' date];
cd(CD);
mworks = ['data-' 'i' SubNum '-' date '-' time]; 
load (mworks);


%set current directory for saving
Save = ['Z:\2P imaging\Analysis\' mouse '\' experiment '\' date '\FlashingStimAnalysis\Rsp2VisStim'];
cd(Save)
%%
%name and convert some mworks variables
cLeverDown = double(cell2mat(input.cLeverDown));
nTrials = (input.trialSinceReset)-1;
RateFRperMS = 30./1000;
Block2ON = double(cell2mat(input.tBlock2TrialNumber));

%%
%find timecourse for iti right before start of trial
meanFrameData = squeeze(mean(mean(data_reg,1),2));

ITItimecourse = zeros(nTrials,15);
for itrial = 1:nTrials
    ind = cLeverDown(itrial)-15:cLeverDown(itrial)-1;
    ITItimecourse(itrial,:) = meanFrameData(ind);
end
ITItimecourse_mean = mean(ITItimecourse,1);

%find timecourse for first 15 frames of trial
TrStarttimecourse = zeros(nTrials,15);
for itrial = 1:nTrials
    ind = cLeverDown(itrial):cLeverDown(itrial)+14;
    TrStarttimecourse(itrial,:) = meanFrameData(ind);
end
TrStarttimecourse_mean = mean(TrStarttimecourse,1);

%align frames for just those around the start of each trial
framesaroundleverdown = zeros(size(data_reg,1),size(data_reg,2),nTrials*30);
start = 1;
L = 30;
for itrial = 1:nTrials
    ind = cLeverDown(itrial)-15:cLeverDown(itrial)+14;
    framesaroundleverdown(:,:,start:start+L-1) = data_reg(:,:,ind);
    start = start+L;
end
writetiff(framesaroundleverdown,'PrePostTrialStart.tif');

PreANDPostTrStartTimecourses = figure;
title('Pre (b) and Post (r) Trial Start');
xlabel('raw fluorescence')
ylabel('frame')
plot(ITItimecourse_mean,'b');
hold on
plot(TrStarttimecourse_mean,'r');