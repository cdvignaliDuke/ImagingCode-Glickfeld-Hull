function [relabeledTrialOri, meanOriEaBin] = binAndRelabelTrialOrientation(...
    trialOrientation,oriBins)
trialOrientation = round(trialOrientation,2,'significant');
nt = length(trialOrientation);
nTargetBins = length(oriBins)-2;
oriBinID = discretize(trialOrientation,oriBins,'IncludedEdge','right');

relabeledTrialOri = nan(1,nt);
relabeledTrialOri(oriBinID == 1) = 0;

meanOriEaBin = nan(1,nTargetBins);
for i = 1:nTargetBins
    ind = oriBinID == (i+1);
    relabeledTrialOri(ind) = oriBins(i+2);
    if nargout > 1
        meanOriEaBin(i) = mean(trialOrientation(ind));
    end
end
end