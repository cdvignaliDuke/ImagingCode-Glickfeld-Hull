function centeredTuningCurves = getOriTuningCurvesCenteredOnPref(neuronTuningCurves, neuronPref, oris)
% neuronTuningCurves should be a matrix of size number of stim x number of
% neurons, orientations should match responses of the rows of this matrix
if oris(end) ~= 180
    orisAll = [oris 180];
    tuningCurves = cat(1, neuronTuningCurves, neuronTuningCurves(1,:));
else
    orisAll = oris;
    tuningCurves = neuronTuningCurves;
end

[nStim,nCell] = size(tuningCurves);

center_ind = find(orisAll == 90);

centeredTuningCurves = nan(nStim,nCell);
for icell = 1:nCell
    pref_ind = neuronPref(icell);
    if ~isnan(pref_ind)
        nOris2shift = center_ind-pref_ind;
        thisCellTuningCurve = tuningCurves(:,icell);
        centeredTuningCurves(:,icell) = circshift(thisCellTuningCurve,nOris2shift);
    end
end

end
