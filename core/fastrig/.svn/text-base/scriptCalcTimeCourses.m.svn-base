function tc = scriptCalcTimeCourses(array,expt)

tc.neurons = calcTimeCourses(array,expt.maskNeurons);
tc.glia = calcTimeCourses(array,expt.maskGlia);
tc.neuropil = calcTimeCourses(array,expt.maskNeuropil);

return;

function tc = calcTimeCourses(array,mask)

tc.raw = stackGetTimeCourses(array, mask);
[tc.dur,tc.ncells]=size(tc.raw);
tc.lowcut = tcLowCut (tc.neurons, expt.trialdur*5, 'gaussian', 1);
tc.dc = repmat(mean(filtered,1),dur,1);
tc.contrast = (filtered-tc.dc)./tc.dc*100;

tc.
trialAv = tcTrialAverage(contrast,expt.trialdur);
shufflecorrected = tcShuffleCorrect(contrast,expt.trialdur);

return;
