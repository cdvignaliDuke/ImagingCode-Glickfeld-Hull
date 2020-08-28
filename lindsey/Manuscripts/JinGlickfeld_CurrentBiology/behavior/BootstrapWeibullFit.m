function [bootStats]=BootstrapWeibullFit(trialVec, pctCorrect,NBootstrapReps, intensV,DoClampAtZero, use50Thresh) 
nTotB = squeeze(trialVec);
pctCorrB = pctCorrect;
pctCorrB(pctCorrB>=1) = 1-10*eps;
pctCorrB(pctCorrB<=0) = 0+10*eps; 
threshB = repmat(NaN, [NBootstrapReps, 1]);
coefEstsB = repmat(NaN, [NBootstrapReps, 4]);
weightsB = repmat(NaN, [NBootstrapReps length(nTotB)]);
 for iR = 1:NBootstrapReps
        tNCorr = binornd(nTotB, pctCorrB);
    
        tPctCorr = tNCorr./nTotB;
        tPctCorr(tPctCorr >= 1) = 1-10*eps;
        tPctCorr(tPctCorr == 0) = 0+10*eps;
        
        % do the fit
        bootFitS = weibullFitLG(intensV, tPctCorr, DoClampAtZero, use50Thresh, ...
            {'nTrials', nTotB });  
       

        threshB(iR) = bootFitS.thresh;
        coefEstsB(iR,:) = bootFitS.coefEsts;
        weightsB(iR,:) = bootFitS.fitWeights';
        
        if mod(iR, 25) == 0
            statusbar(0, '%s done %4g%%: %d of %d bootstrap reps', ...
                mfilename, iR/NBootstrapReps*100, iR, NBootstrapReps);
        end
 end
 
 nNaN = mean(sum(isnan(coefEstsB)));
    if nNaN > 0
        disp(sprintf('** Errors/Reps: %d/%d ', nNaN, NBootstrapReps));
    end
    
    slopeB = coefEstsB(:,2);

    bootStats.threshMean = nanmean(threshB);  % should be near the true value
    bootStats.threshStd = nanstd(threshB);  % should be near the true value
    bootStats.ci95 = prctile(threshB, [0 100] + 5*[1 -1]/2);
    bootStats.ci99 = prctile(threshB, [0 100] + 1*[1 -1]/2);
    
  
    bootStats.slopeMean = nanmean(slopeB);
    bootStats.slopeStd = std(slopeB);
    bootStats.slopeCi95 = prctile(slopeB, [0 100] + 5*[1 -1]/2);
    bootStats.slopeCi99 = prctile(slopeB, [0 100] + 1*[1 -1]/2);
 
 
 
 
end
