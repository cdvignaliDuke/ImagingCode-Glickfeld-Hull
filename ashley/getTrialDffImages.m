function dff = getTrialDffImages(imageStack,trialStartInd,...
    nBaselineFrames,nTrialFrames)

    [ypix,xpix,nframes] = size(imageStack);
    ntrials = length(trialStartInd);
    if trialStartInd(end)+nTrialFrames > nframes
        extraTrials = trialStartInd+nTrialFrames > nframes;
        ntrials = ntrials-sum(extraTrials);
        trialStartInd = trialStartInd(1:ntrials);
    end

    dff = nan(ypix,xpix,nBaselineFrames+nTrialFrames,ntrials);
    for itrial = 1:ntrials
        baselineFramesThisTrialInd = trialStartInd(itrial)-nBaselineFrames+1:...
            trialStartInd(itrial);
        framesThisTrialInd = trialStartInd(itrial)-nBaselineFrames+1:...
            trialStartInd(itrial)+nTrialFrames;
        f0 = mean(imageStack(:,:,baselineFramesThisTrialInd),3);
        f1 = imageStack(:,:,framesThisTrialInd);
        dff(:,:,:,itrial) = bsxfun(@rdivide, bsxfun(@minus, f1, f0), f0);
    end
end