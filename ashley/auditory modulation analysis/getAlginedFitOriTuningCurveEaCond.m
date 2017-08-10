function alignedTuningCurveEaCond = getAlginedRawOriTuningCurveEaCond(...
    neuronTuning,dataField,stimPref,stimuli)

    stimuli = [stimuli 180];
    
    ndelay = size(neuronTuning,2);
    ncells = length(stimPref);
    nstim = length(stimuli);
    centerInd = find(stimuli == 90);
    alignedTuningCurveEaCond = cell(1,ndelay);
    for idelay = 1:ndelay
        delayData = eval(['neuronTuning(idelay).' dataField]);
        delayDataCircular = cat(1,delayData,delayData(1,:));
        centeredCurves = nan(nstim,ncells);
        for icell = 1:ncells
            thisCellPref = stimPref(icell);
            if thisCellPref < centerInd
                nshift = centerInd - thisCellPref;
            else
                nshift = thisCellPref - centerInd;
            end
            thisCellTuningCurve = delayDataCircular(:,icell);
            thisCellCenteredCurve = circshift(thisCellTuningCurve,nshift);
            centeredCurves(:,icell) = thisCellCenteredCurve;
        end
        alignedTuningCurveEaCond{idelay} = centeredCurves;
    end    

end