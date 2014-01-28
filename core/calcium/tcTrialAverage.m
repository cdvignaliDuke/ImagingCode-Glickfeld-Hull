function out = tcTrialAverage(timeCourses,epochs,blanks,win)
%TCTRIALAVERAGE
% AV = TCTRIALAVERAGE(TIMECOURSES,EPOCHS,BLANKS,WIN)

eLen = getEpochsLen(epochs);
bLen = getEpochsLen(blanks);

if win < 0
    win = 1:win;
else
    win = bLen-win+1:bLen;
end

[nSamples,nCells] = size(timeCourses);

[nTrials,nStim]=size(epochs);

out.trials.F = zeros(eLen,nCells,nTrials,nStim);
out.trials.dF = zeros(eLen,nCells,nTrials,nStim);
out.trials.ratio = zeros(eLen,nCells,nTrials,nStim);

result = cell(nCells,1);

for iTrial = 1:nTrials
    for iStim = 1:nStim
        
        ind = epochs{iTrial,iStim}(1:eLen);        
        
        out.trials.F(:,:,iTrial,iStim) = timeCourses(ind,:);
        
        if bLen > 0
            baseind = blanks{iTrial,iStim}(win);
        else
            baseind = epochs{iTrial,iStim}(1:eLen);
        end
        
        baseline = mean(timeCourses(baseind,:),1);
        
        out.trials.dF(:,:,iTrial,iStim) = bsxfun(@minus,out.trials.F(:,:,iTrial,iStim),baseline);
        
        out.trials.ratio(:,:,iTrial,iStim) = bsxfun(@rdivide,out.trials.dF(:,:,iTrial,iStim),baseline)*100;
        
    end
end

out.av.F = squeeze(mean(out.trials.F,3));
out.av.dF = squeeze(mean(out.trials.dF,3));
out.av.ratio = squeeze(mean(out.trials.ratio,3));

return;

traces = squeeze(mean(results.trials.dF,1));