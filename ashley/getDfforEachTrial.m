function dfff = getDFFforEachTrial(data,trType,tTrialStart, tTargetOn)
frRateHz = 30;
nBaselineFrames = frRateHz;
nTargetFrames = frRateHz;
nStartFrames = frRateHz*8;
nTrials = length(trType);
nCells = size(data,2);

dff_target = nan(nBaselineFrames+nTargetFrames,nCells,nTrials);
for itrial = 1:nTrials
    f0 = mean(data(tTargetOn(itrial)-nBaselineFrames+1:tTargetOn(itrial),:),1);
    f1 = data(tTargetOn(itrial)-nBaselineFrames+1:tTargetOn(itrial)+nTargetFrames,:);
    dff_target(:,:,itrial) = bsxfun(@rdivide, bsxfun(@minus, f1, f0), f0);
end

dff_start = nan(nBaselineFrames+nStartFrames,nCells,nTrials);
for itrial = 1:nTrials
    f0 = mean(data(tTrialStart(itrial)-nBaselineFrames+1:tTrialStart(itrial),:),1);
    f1 = data(tTrialStart(itrial)-nBaselineFrames+1:tTrialStart(itrial)+nStartFrames,:);
    dff_start(:,:,itrial) = bsxfun(@rdivide, bsxfun(@minus, f1, f0), f0);
end

dff_targetWithStartBL = nan(nBaselineFrames+nTargetFrames,nCells,nTrials);
for itrial = 1:nTrials
    f0 = mean(data(tTrialStart(itrial)-nBaselineFrames+1:tTrialStart(itrial),:),1);
    f1 = data(tTargetOn(itrial)-nBaselineFrames+1:tTargetOn(itrial)+nTargetFrames,:);
    dff_targetWithStartBL(:,:,itrial) = bsxfun(@rdivide, bsxfun(@minus, f1, f0), f0);
end
