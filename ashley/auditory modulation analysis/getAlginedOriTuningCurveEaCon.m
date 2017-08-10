function alignedTuningCurveEaCond = getAlginedOriTuningCurveEaCond(...
    neuronTuning,dataField,stimPref,stimuli)
    % stimuli should be symmetrical
    
    ndelay = size(neuronTuning,2);
    ncells = length(stimPref);
    centerInd = find(stimuli == 90);
    
    for idelay = 1:ndelay
        delayData = eval(['neuronTuning(idelay).' dataField]);
        for icell = 1:ncells
            thisCellPref = stimPref(icell);
            nshift = thisCellPref - centerInd;
            thisCellTuningCurve = delayData(:,icell);
            thisCellCenteredCurve = circshift(thisCellTuningCurve,nshift);
        end
    end    

end