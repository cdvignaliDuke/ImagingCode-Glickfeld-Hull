function pctCorr = cmbPctCorr(nTrials1,rate1,nTrials2,rate2)
    pctCorr = ((nTrials1*rate1)+(nTrials2*rate2))./(nTrials1+nTrials2);
end