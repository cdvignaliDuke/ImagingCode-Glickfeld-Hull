function tcs = scriptProcessTimeCourses(raw, fc, artifacts,Fs)
%SCRIPTPROCESSTIMECOURSES
% TCS = SCRIPTPROCESSTIMECOURSES(RAW, FC, ARTIFACTS)
%
% DEPRECATED

%% process
tcs.raw = raw;
[tcs.nsamples,tcs.ncells]=size(tcs.raw);
tcs.fc = fc;
tcs.baseline = repmat(mean(tcs.raw,1),tcs.nsamples,1);
tcs.df = tcs.raw-tcs.baseline;

% high pass filter
period = 1/fc(1)*Fs; % in samples
tcs.lowcut = tcLowCut (tcs.df, period, 'gaussian', 1);

tcs.dfof = tcs.lowcut./tcs.baseline*100;

%% low pass filter
wn = min(fc(2)/(Fs/2),0.99);
[b,a]=butter(5,wn);
tcs.smooth = filtfilt(b,a,tcs.dfof);

if nargin > 2 | ~isempty(artifacts)
    tcs.correct = tcRemoveArtifacts(tcs.smooth,artifacts);
end

% tcs.trial = tcTrialAverage(tcs.lowcut,trialdur);
% tcs.shuffleCorrect = tcShuffleCorrect(tcs.dfof,trialdur);

return;