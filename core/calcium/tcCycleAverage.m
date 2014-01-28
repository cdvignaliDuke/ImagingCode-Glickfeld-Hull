function av = tcCycleAverage(timeCourses,trialdur)
%TCCYCLEAVERAGE Average time course over cycle
% AV = TCCYCLEAVERAGE(TIMECOURSES,TRIALDUR)
%
% See also tcEpochAverage, tcTrialAverage.
%

[nsamples,ncells]=size(timeCourses);

ntrials = floor(nsamples/trialdur);

shuffled=reshape(timeCourses(:),trialdur,ntrials*ncells);

av = zeros(trialdur,ncells);

for icell =1:ncells
    sel = (icell-1)*ntrials+1:icell*ntrials;
    av(:,icell) = mean(shuffled(:,sel),2);
end

return