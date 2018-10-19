function tuningBootPctSuccess = getTuningReliability(tOri,stimResp,tuning,nboot)
% stimResp should be cells x trials, tOri should be n trials long

% nboot = 10;
    orientations = unique(tOri);
    nOri = length(orientations);
    nCells = size(stimResp,1);

    tuning_boot = nan(nOri,nCells,nboot);
    for iori = 1:nOri
        ind = tOri == orientations(iori);
        for iboot = 1:nboot
            sampleInd = randsample(find(ind),sum(ind),true);
            tuning_boot(iori,:,iboot) = mean(stimResp(:,sampleInd),2);
        end    
    end

    tuningPeakID_boot = nan(nCells,nboot);
    for icell = 1:nCells
        [~,tuningPeakID_boot(icell,:)] = max(squeeze(tuning_boot(:,icell,:)),[],1);
    end

    [~,tuningPeakID] = max(tuning,[],1);
    tuningDiff = tuningPeakID_boot - tuningPeakID';
    tuningBootPctSuccess = sum(tuningDiff == 0,2)./nboot;
end