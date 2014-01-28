function trials = tcSortTrials(timeCourses,epochs)
%TCSORTTRIALS 
% TRIALS = TCSORTTRIALS(TIMECOURSES,EPOCHS) where EPOCHS is a cell array of
% frame indices and 
% TRIALS is amultidimensional array of size [EPOCHLEN,NCELLS,NTRIALS,NSTIM]

eLen = getEpochsLen(epochs);

[nSamples,nCells] = size(timeCourses);

[nTrials,nStim]=size(epochs);

trials = zeros(eLen,nCells,nTrials,nStim);

for iTrial = 1:nTrials
    for iStim = 1:nStim
        
        ind = epochs{iTrial,iStim}(1:eLen);        
        
        trials(:,:,iTrial,iStim) = timeCourses(ind,:);
        
    end
end

return;
