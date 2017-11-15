%quick behavior imaging check
clear all
close all
awFSAVdatasets_naive100ms
iexp = 1;
irun = 1;
quickBehaviorSummary(expt,iexp)

frameRateHz = expt(iexp).frame_rate;
nTrialFrames = frameRateHz*2;
nBaselineFrames = frameRateHz;

subnum = expt(iexp).SubNum;
mouse = expt(iexp).mouse;
expDate = expt(iexp).date;
expTime = expt(iexp).time_mat(irun,:);
expFolder = expt(iexp).runs(irun,:);
fname = [expFolder '_000_000'];

mworks = loadMworksFile(subnum,expDate,expTime);

data = loadsbx_choosepmt(1,mouse,expDate,expFolder,fname);

trialStart = double(cell2mat(mworks.cFirstStim));
trialFramesInd = linspaceNDim(trialStart - nBaselineFrames+1,...
    trialStart+nTrialFrames,nBaselineFrames+nTrialFrames)';
ntrials = length(trialStart);

[ypix,xpix,~] = size(data);
trialFrames = data(:,:,trialFramesInd(:));
d = double(reshape(trialFrames,ypix,xpix,nBaselineFrames+nTrialFrames,ntrials));

F = mean(d(:,:,1:nBaselineFrames,:),3);
dFF = bsxfun(@rdivide,bsxfun(@minus,d,F),F);

trialMeanFrames = squeeze(mean(dFF(...
    :,:,nBaselineFrames+1:nBaselineFrames+nTrialFrames,:),3));

maxDFF = max(trialMeanFrames,[],3);

bwout = imCellEditInteractive(maxDFF);
mask = bwlabel(bwout);

tc = stackGetTimeCourses(data,mask);

%%
if iscell(mworks.nFramesOn)
    cycLengthFrames = unique(cell2mat(mworks.nFramesOn)) + ...
        unique(cell2mat(mworks.nFramesOff));
else
    cycLengthFrames = mworks.nFramesOn + mworks.nFramesOff;
end
nTrialFrames = cycLengthFrames*4;
outcome = mworks.trialOutcomeCell;
hitsAndMisses = strcmp(outcome,'success') | strcmp(outcome,'ignore');
nCycPerTrial = cell2mat(mworks.tCyclesOn);
trialInd = nCycPerTrial >=4 & hitsAndMisses;
trialStartInd = trialStart(trialInd);
trialFramesInd = linspaceNDim(trialStartInd - nBaselineFrames+1,...
    trialStartInd+nTrialFrames,nBaselineFrames+nTrialFrames)';
ntrials = length(trialStartInd);

ncells = size(tc,2);

trialTC = nan(nBaselineFrames+(cycLengthFrames*4),ncells,ntrials);
for itrial = 1:ntrials
    tcThisTrial = tc(trialFramesInd(:,itrial),:);
    trialTC(:,:,itrial) = tcThisTrial;
end

trialF = mean(trialTC(1:nBaselineFrames,:,:));
trialDFF = bsxfun(@rdivide,bsxfun(@minus,trialTC,trialF),trialF);

meanDFFAllCells = mean(trialDFF,3);

tt = (-nBaselineFrames+1:cycLengthFrames*4)./frameRateHz;
figure;
plot(tt,meanDFFAllCells)
title({'hit and miss trials, aud and vis, selected cells';[num2str(ntrials) ' trials']})
figXAxis([],'time (s)',[-0.25 tt(end)])
figYAxis([],'dF/F',[])


