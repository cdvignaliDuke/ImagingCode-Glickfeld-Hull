% testing curve fitting
% here trying different fitting methods

i=28; % cell

figure(1);clf;
iCon = 1; % con 1 (0.2)
for iCon = 1:nCon
    fprintf(['\nCell# ' num2str(i) ', bIndex ' num2str(bIndex(i))])

    fprintf(['\nContrast: ' num2str(cons(iCon))])

    dumdum = zeros(1,8);%[0]; % define zero point for dF/F
    szs0 = zeros(1,8);%[0.1]; % and size
    nPts = 8;
    for iSz = 1:nSize
        nPts = nPts + length(sizeTune{iSz,iCon,i});
        dumdum = [dumdum sizeTune{iSz,iCon,i}'];
        szs0 = [szs0 szs(iSz)*ones(1,length(sizeTune{iSz,iCon,i}))];
    end

    SStot = sum((dumdum-mean(dumdum)).^2);
    fprintf([' (SStot: ' num2str(SStot) ')'])

    % make initial guesses
    maxMean = max(sizeMean(:,iCon,i));
    
    % single sigmoid model
    fprintf('\nEvaluating fit 1...\n')
    guess1 = [0.5*maxMean 0.35 10]; % guesses for 3 params (Ae, ke, xe)
    lb = [0 0 0];
    ub = [inf 1 80]; %[inf inf inf]; % here limit steepness to 1 (above this is way too steep)
    %[fitx,resnorm1] = lsqcurvefit(logfit1,guess,szs0,dumdum,lb,ub,opts);
    %guess = lsqcurvefit(logfit1,guess,szs0,dumdum,lb,ub,opts);
    [fitx,OF1] = fminsearchbnd(@(c)sigevalpen1(c,szs0,dumdum),guess1,lb,ub,opts);
    resnorm1 = sigeval1(fitx,szs0,dumdum);
    
    fprintf('Fit 1 readout:\n')
    fprintf(['Ex amp: ' num2str(fitx(1)) ', steepness: ' num2str(fitx(2)) ', center: ' num2str(fitx(3))])
    r2 = 1 - resnorm1/SStot; % resnorm is SSE
    fprintf(['\nFit 1 R-sq:' num2str(r2) '\n'])
    fprintf(['Fit 1 OF:' num2str(OF1) ', SSE:' num2str(resnorm1) ', Pen:' num2str(OF1-resnorm1) '\n'])
    sizeFits.Rsq1(i, iCon) = r2;
    sizeFits.coefs1(i,iCon,:) = fitx;

    fitout1 = logfit1(fitx,szRng);
    maxResp1 = max(fitout1);
    prefSize1 = szRng(find(fitout1>(0.9*maxResp1),1));
    suppInd1 = 0;
    fprintf(['Pref size: ' num2str(prefSize1) ' (dF/F: ' num2str(0.9*maxResp1) ')\n'])
    fprintf(['Suppression Index: ' num2str(suppInd1) ' (no suppression in single sigmoid model)\n'])

    % double sigmoid model
    fprintf('\nEvaluating fit 2...\n')
    guess2 = [maxMean 0.4 10 maxMean-mean(dumdum) 0.3 20]; % guesses for (Ae, k1, x1, Ai, k2 x2)
    lb = [0 0 0 0 0 0];
    ub = [inf 1 80 inf 1 80]; %[inf inf inf inf inf inf]; % here limit steepness to 1 (above this is way too steep)
    %[fitx,resnorm2] = lsqcurvefit(logfit2,guess,szs0,dumdum,lb,ub,opts);
    %guess = lsqcurvefit(logfit2,guess,szs0,dumdum,lb,ub,opts);
    [fitx,OF2] = fminsearchbnd(@(c)sigevalpen2(c,szs0,dumdum),guess2,lb,ub,opts);
    resnorm2 = sigeval2(fitx,szs0,dumdum);

    fprintf('Fit 2 readout:\n')
    fprintf(['Ex amp: ' num2str(fitx(1)) ', steepness: ' num2str(fitx(2)+fitx(5)) ', center: ' num2str(fitx(3))])
    fprintf(['\nInh amp: ' num2str(fitx(4)) ', steepness: ' num2str(fitx(5)) ', center: ' num2str(fitx(3)+fitx(6))])
    r2 = 1 - resnorm2/SStot; % resnorm is SSE
    fprintf(['\nFit 2 R-sq: ' num2str(r2) '\n'])
    fprintf(['Fit 2 OF:' num2str(OF2) ', SSE:' num2str(resnorm2) ', Pen:' num2str(OF2-resnorm2) '\n'])
    sizeFits.Rsq2(i,iCon) = r2;
    sizeFits.coefs2(i,iCon,:) = fitx;

    fitout2 = logfit2(fitx,szRng);
    maxResp2 = max(fitout2);
    prefSize2 = szRng(find(fitout2==maxResp2,1));
    suppInd2 = 1 - fitout2(end)/maxResp2;
    fprintf(['Pref size: ' num2str(prefSize2) ' (dF/F: ' num2str(maxResp2) ')\n'])
    fprintf(['Suppression Index: ' num2str(suppInd2) '\n'])

    % plot all contrasts
    subplot(2,nCon,iCon)
    errorbar([0 szs],[0 sizeMean(:,iCon,i)'],[0 sizeSEM(:,iCon,i)'])
    hold on
    plot(szs0,dumdum,'.b')
    plot(szRng,fitout1,'-')
    plot(szRng,logfit1(guess1,szRng),'g--')
    hold off
    ylim([min([-0.5*maxResp1 min(sizeMean(:,iCon,i))]) 1.4*max([maxResp2 max(sizeMean(:,iCon,i))])])
    title(['Cell #' num2str(i) ', Con ' num2str(cons(iCon)) ', R^2=' num2str(sizeFits.Rsq1(i,iCon))])
    xlabel('Stimulus Size')
    ylabel('dF/F')
    legend('mean','data','fit','guess')

    subplot(2,nCon,iCon+nCon)
    errorbar([0 szs],[0 sizeMean(:,iCon,i)'],[0 sizeSEM(:,iCon,i)'])
    hold on
    plot(szs0,dumdum,'.b')
    plot(szRng,fitout2,'-')
    plot(szRng,logfit2(guess2,szRng),'g--')
    hold off
    ylim([min([-0.5*maxResp1 min(sizeMean(:,iCon,i))]) 1.4*max([maxResp2 max(sizeMean(:,iCon,i))])])
    title(['Cell #' num2str(i) ', Con ' num2str(cons(iCon)) ', R^2=' num2str(sizeFits.Rsq2(i,iCon))])
    xlabel('Stimulus Size')
    ylabel('dF/F')
    legend('mean','data','fit','guess')
    
    % now F-test
    fprintf([num2str(nPts) ' points\n'])
    dfS = nPts - 3; % single sigmoid model df = #sizes - 3 model parameters
    dfD = nPts - 6; % double sigmoid model df = #sizes - 6 model parameters
    Fcrit = finv(0.95,dfD,dfS); % measure critical F at a=0.05
    sizeFits.Fscore(i,iCon) = (resnorm2 - resnorm1)/(dfD-dfS) ./ (resnorm2/dfD);
    Ftest = sizeFits.Fscore(i,iCon)>Fcrit;
    sizeFits.Ftest(i,iCon) = Ftest;
end

fprintf(['\nF-tests for cell i=' num2str(i) ', distIndex ' num2str(bIndex(i)) ':\nF-scores: ' num2str(sizeFits.Fscore(i,:)) '\nChoose model 2? ' num2str(sizeFits.Ftest(i,:)) '\n'])

