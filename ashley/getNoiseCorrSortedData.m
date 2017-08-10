function noiseCorr = getNoiseCorrSortedData(data)
% data should be timepoints x cells x trials
meanResp = squeeze(mean(data,1))';
noiseCorr = corrcoef(meanResp);
end