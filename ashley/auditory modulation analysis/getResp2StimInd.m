function resp2StimIndEaCond = getResp2StimInd(neuronTuning,stimInd)

    nconditions = size(neuronTuning,2);
    ncells = length(stimInd);
    resp2StimIndEaCond = cell(1,nconditions);
    for icond = 1:nconditions
        resp2StimInd = nan(1,ncells);
        for icell = 1:ncells
            thisCellInd = stimInd(icell);
            thisCellResp = neuronTuning(icond).rawTuningCurve(thisCellInd,icell);
            resp2StimInd(icell) = thisCellResp;
        end
        resp2StimIndEaCond{icond} = resp2StimInd;
    end

end