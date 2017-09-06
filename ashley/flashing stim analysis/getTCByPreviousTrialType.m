function tc = getTCByPreviousTrialType(tcInfo,nCycles,previousType)
hitTrials = 1;
missTrials = 2;
%hits & misses
trialLengthInd = cat(2,tcInfo.outcome(hitTrials).tcyc >= nCycles,...
    tcInfo.outcome(missTrials).tcyc >= nCycles);
previousTrialInd = cat(2,tcInfo.outcome(hitTrials).prevTrType,...
    tcInfo.outcome(missTrials).prevTrType) == previousType;


trialComboInd = trialLengthInd & previousTrialInd;

tcAllTrials = cat(3,tcInfo.outcome(hitTrials).resp, ...
    tcInfo.outcome(missTrials).resp);
tc = tcAllTrials(:,:,trialComboInd);
end