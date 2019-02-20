function matchTrialInd = getMatchedOutcomeTrialIndex(...
    trOut, trStimID, minTrN)

h = strcmp(trOut,'h');
m = strcmp(trOut,'m');

nHit = histcounts(trStimID(h & trStimID > 1),2:4);
nMiss = histcounts(trStimID(m & trStimID > 1),2:4);

matchTrialInd = cell(1,2);
if any(nHit+nMiss < minTrN)
    ind = find(any(nHit+nMiss < minTrN))+1;
    matchTrialInd = (h|m) & trStimID~=ind;
else
    for i = 1:2
        if nHit(i) == 0 || nMiss(i) == 0
            continue
        elseif nHit(i) > nMiss(i)
            ind = trStimID == i+1 & m;
            ind2 = randsample(find(trStimID == i+1 & h),sum(ind));
            ind(ind2) = true;
        elseif nHit(1) < nMiss(i)
            ind = trStimID == i+1 & h;
            ind2 = randsample(find(trStimID == i+1 & m),sum(ind));
            ind(ind2) = true;
        else
            ind = h | m;
        end
        matchTrialInd{i} = ind;
    end
    matchTrialInd = sum(cell2mat(matchTrialInd(~isempty(matchTrialInd))'));
end
end

