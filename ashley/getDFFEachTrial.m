function dff = getDFFEachTrial(data,trialBaselineInd,trialStimInd,...
    nStimFrames,frRateHz)
% data is imaging data, frames x neurons, trial indicies are frame indicies
% for each trial, nStimFrames is length of time in frames you want the
% time-course to be after the stimulus start
% length of baseline is always 1s
nBaselineFrames = frRateHz;
nTrials = length(trialBaselineInd);
[nFrames, nCells] = size(data);

dff = nan(nBaselineFrames+nStimFrames,nCells,nTrials);
for itrial = 1:nTrials
    f0 = mean(data(trialBaselineInd(itrial)-nBaselineFrames+1:trialBaselineInd(itrial),:),1);
    frame_ind = trialStimInd(itrial)-nBaselineFrames+1:trialStimInd(itrial)+nStimFrames;
    if frame_ind(end) > nFrames
        extraframes = sum(frame_ind > nFrames);
        frame_ind_chop = frame_ind(1):frame_ind(frame_ind == nFrames);
        padtrialnan = nan(extraframes,nCells);
        f1 = cat(1, data(frame_ind_chop,:), padtrialnan);
    else
        f1 = data(frame_ind,:);
    end
    dff(:,:,itrial) = bsxfun(@rdivide, bsxfun(@minus, f1, f0), f0);
end
