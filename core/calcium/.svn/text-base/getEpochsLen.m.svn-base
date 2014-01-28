function len =  getEpochsLen(epochs)
%GETEPOCHLEN
% len =  getEpochLen(epochs)

[nTrials,nStim]=size(epochs);

len = inf;

for iStim = 1:nStim
    for iTrial = 1:nTrials
        len=min(len,length(epochs{iTrial,iStim}));
    end
end

return;
