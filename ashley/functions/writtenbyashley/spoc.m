function steps = spoc(maxVal, spo, nSteps)
% calculate vector of values for a given max value, MAXVAL, number of
% steps, NSTEPS, and steps per octave, SPO.

steps = maxVal./(2.^(([1:nSteps]-1)./spo));
end