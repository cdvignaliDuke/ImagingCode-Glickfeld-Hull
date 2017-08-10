function centeredTuningCurves = getOriTuningCurvesCenteredOnPref(neuronTuningCurves, neuronPref, oris)
% neuronTuningCurves should be a matrix of size number of stim x number of
% neurons, orientations should match responses of the rows of this matrix
if oris(1) ~= oris(end)+180
    orisAll = [oris 180];
    tuningCurves = cat(1, neuronTuningCurves, neuronTuningCurves(1,:);
else
    orisAll = oris;
    tuningCurves = neuronTuningCurves;
end

[nStim,nCell] = size(neuronTuningCurves);

for icell = 1:nCell
    pref_ind = neuronPref(icell);
end

end
