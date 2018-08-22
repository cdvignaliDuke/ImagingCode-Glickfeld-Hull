close all
clear all

load('Z:\Analysis\FSAV Choice\FSAV_decodeData.mat')
fnout = 'Z:\Analysis\FSAV Choice';
nexp = length(dcExpt);

doExptPlots = 0;
%%
minTrN = 20;
theta90Threshold = 11.25;
decisionVariable = 0.5;
detectThresholdHR = 0.8;
detectThresholdOri = 45;
oriBinsFine = [0 1 16 22.5 32 45 64 90];
oriBins = [0 1 32 90];
ampBins = [0 0.0001 0.1 1];
weightsDiffZeroAlpha = 0.05/2;
%%
targetWeightsAll_vis = [];
detectWeightsAll_vis = [];

targetCorrsAll_vis = [];
detectCorrsAll_vis = [];

targetWeightsAll_shuf_vis = [];
detectWeightsAll_shuf_vis = [];

targetCorrsAll_shuf_vis = [];
detectCorrsAll_shuf_vis = [];

targetWeightsAll_aud = [];
detectWeightsAll_aud = [];

targetCorrsAll_aud = [];
detectCorrsAll_aud = [];

targetWeightsAll_shuf_aud = [];
detectWeightsAll_shuf_aud = [];

targetCorrsAll_shuf_aud = [];
detectCorrsAll_shuf_aud = [];

nCells = zeros(1,nexp);
oriPrefAll = [];
% pcv = cell(1,nexp);
exptName = cell(1,nexp);
for iexp=1:nexp
    disp([dcExpt(iexp).mouse '-' dcExpt(iexp).date])
   trOut=dcExpt(iexp).trialOutcome;
    [Y_yeses, Y_targets] = getStimAndBehaviorYs(trOut);
    
    
    tOri = dcExpt(iexp).trialOrientation;
    tOri = binAndRelabelTrialOrientation(tOri,oriBins);
    nOriBin = histcounts(discretize(tOri,oriBins,'IncludedEdge','right'));
    nOriBin = nOriBin(nOriBin > 0);
    orientations = unique(tOri);
    nOri = length(orientations);
    minBinN = min(nOriBin);
    if minBinN < minTrN
        minBinN = min(nOriBin(nOriBin >= minTrN));
    end
    matchTrialsInd = [];
    for iori = 1:nOri
        ind = find(tOri == orientations(iori));
        if length(ind) >= minTrN
            if iori == 1
                if minBinN == nOriBin(iori)
                    error('not enough FA/CR trials')
                else
                    n = (nOri-1).*minBinN;
                    if n > length(ind)
                        error('not enough FA/CR trials')
                    else
                        indSample = randsample(ind,n);
                        matchTrialsInd = cat(2,matchTrialsInd,indSample);
                    end
                end
            else
                indSample = randsample(ind,minBinN);
                matchTrialsInd = cat(2,matchTrialsInd,indSample);
            end
        end
    end
    

   X_targets=dcExpt(iexp).stimResp';
   
   cellIdx = dcExpt(iexp).signifResponsiveCells' & ...
       dcExpt(iexp).oriFitTheta90 < theta90Threshold;
   nCells(iexp) = sum(cellIdx);

   X_targets=X_targets(matchTrialsInd,cellIdx);
   
   Y_yeses=Y_yeses(matchTrialsInd);   
   Y_targets=Y_targets(matchTrialsInd); 
   
%    [c,score,latent] = pca(X_targets);
   
%    pcv{iexp} = cumsum(laten)./sum(latent);
%    exptName{iexp} = [dcExpt(iexp).mouse '-' dcExpt(iexp).date];
   
    Y_yeses_shuf = Y_yeses(randperm(length(Y_yeses)));
    Y_targets_shuf = Y_targets(randperm(length(Y_targets)));
    
%    X_targets = zscore(X_targets);
 
   detectCorr = corr(Y_yeses,X_targets)';
   targetCorr = corr(Y_targets,X_targets)';
   
   detectCorr_shuf = corr(Y_yeses_shuf,X_targets)';
   targetCorr_shuf = corr(Y_targets_shuf,X_targets)';
   
   detectCorrsAll_vis = cat(1,detectCorrsAll_vis,detectCorr);
   targetCorrsAll_vis = cat(1,targetCorrsAll_vis,targetCorr);
   
   detectCorrsAll_shuf_vis = cat(1,detectCorrsAll_shuf_vis,detectCorr_shuf);
   targetCorrsAll_shuf_vis = cat(1,targetCorrsAll_shuf_vis,targetCorr_shuf);
   
   X_targets=bsxfun(@plus,X_targets,-mean(X_targets));
   X_targets=bsxfun(@times,X_targets,1./std(X_targets));
   C=eye(size(X_targets,2));%C=inv(sqrtm(cov(X)));
   X_targets=X_targets*C;
   p=1;
   [temp1,dev1,targetGLM]=glmfit(X_targets,Y_targets,'binomial');
   [temp2,dev2,detectGLM]=glmfit(X_targets,Y_yeses,'binomial');
   idx=find(targetGLM.p>p|detectGLM.p>p);
   temp1(idx)=0;
   temp2(idx)=0;
   B_targets = C*temp1(2:end,1);
   B_yeses = C*temp2(2:end,1);
   
   [temp1,dev1,targetGLM]=glmfit(X_targets,Y_targets_shuf,'binomial');
   [temp2,dev2,detectGLM]=glmfit(X_targets,Y_yeses_shuf,'binomial');
   idx=find(targetGLM.p>p|detectGLM.p>p);
   temp1(idx)=0;
   temp2(idx)=0;
   B_targets_shuf = C*temp1(2:end,1);
   B_yeses_shuf = C*temp2(2:end,1);
   
   targetWeightsAll_vis = cat(1,targetWeightsAll_vis,B_targets);
   detectWeightsAll_vis = cat(1,detectWeightsAll_vis,B_yeses);
  
   targetWeightsAll_shuf_vis = cat(1,targetWeightsAll_shuf_vis,B_targets_shuf);
   detectWeightsAll_shuf_vis = cat(1,detectWeightsAll_shuf_vis,B_yeses_shuf);
  
%    SI=dcExpt(iexp).baseTargRatio;
   
   oriPref = dcExpt(iexp).oriPref(cellIdx);
   oriPrefAll = cat(2,oriPrefAll,oriPref);
   
   %***AUDITORY TRIALS
   trOut=dcExpt(iexp).audTrialOutcome;
    [Y_yeses, Y_targets] = getStimAndBehaviorYs(trOut);
    
   tAmp = dcExpt(iexp).audTrialAmp;
    tAmp = binAndRelabelTrialOrientation(tAmp,ampBins);
    nAmpBin = histcounts(discretize(tAmp,ampBins,'IncludedEdge','right'));
    nAmpBin = nAmpBin(nAmpBin > 0);
    amplitudes = unique(tAmp);
    nAmp = length(amplitudes);
    minBinN = min(nAmpBin);
    if minBinN < minTrN
        minBinN = min(nAmpBin(nAmpBin >= minTrN));
    end
    matchTrialsInd = [];
    for iamp = 1:nAmp
        ind = find(tAmp == amplitudes(iamp));
        if length(ind) >= minTrN
            if iamp == 1
                if minBinN == nAmpBin(iamp)
                    error('not enough FA/CR trials')
                else
                    n = (nAmp-1).*minBinN;
                    if n > length(ind)
                        error('not enough FA/CR trials')
                    else
                        indSample = randsample(ind,n);
                        matchTrialsInd = cat(2,matchTrialsInd,indSample);
                    end
                end
            else
                indSample = randsample(ind,minBinN);
                matchTrialsInd = cat(2,matchTrialsInd,indSample);
            end
        end
    end 
    

   X_targets=dcExpt(iexp).audStimResp';
   
   cellIdx = dcExpt(iexp).signifResponsiveCells' & ...
       dcExpt(iexp).oriFitTheta90 < theta90Threshold;

   X_targets=X_targets(matchTrialsInd,cellIdx);
   
   Y_yeses=Y_yeses(matchTrialsInd);   
   Y_targets=Y_targets(matchTrialsInd);   
   
    Y_yeses_shuf = Y_yeses(randperm(length(Y_yeses)));
    Y_targets_shuf = Y_targets(randperm(length(Y_targets)));
   
   detectCorr = corr(Y_yeses,X_targets)';
   targetCorr = corr(Y_targets,X_targets)';
   
   detectCorr_shuf = corr(Y_yeses_shuf,X_targets)';
   targetCorr_shuf = corr(Y_targets_shuf,X_targets)';
   
   detectCorrsAll_aud = cat(1,detectCorrsAll_aud,detectCorr);
   targetCorrsAll_aud = cat(1,targetCorrsAll_aud,targetCorr);
   
   detectCorrsAll_shuf_aud = cat(1,detectCorrsAll_shuf_aud,detectCorr_shuf);
   targetCorrsAll_shuf_aud = cat(1,targetCorrsAll_shuf_aud,targetCorr_shuf);
   
   X_targets=bsxfun(@plus,X_targets,-mean(X_targets));
   X_targets=bsxfun(@times,X_targets,1./std(X_targets));
   C=eye(size(X_targets,2));%C=inv(sqrtm(cov(X)));
   X_targets=X_targets*C;
   p=1;
   [temp1,dev1,targetGLM]=glmfit(X_targets,Y_targets,'binomial');
   [temp2,dev2,detectGLM]=glmfit(X_targets,Y_yeses,'binomial');
   idx=find(targetGLM.p>p|detectGLM.p>p);
   temp1(idx)=0;
   temp2(idx)=0;
   B_targets = C*temp1(2:end,1);
   B_yeses = C*temp2(2:end,1);
   
   [temp1,dev1,targetGLM]=glmfit(X_targets,Y_targets_shuf,'binomial');
   [temp2,dev2,detectGLM]=glmfit(X_targets,Y_yeses_shuf,'binomial');
   idx=find(targetGLM.p>p|detectGLM.p>p);
   temp1(idx)=0;
   temp2(idx)=0;
   B_targets_shuf = C*temp1(2:end,1);
   B_yeses_shuf = C*temp2(2:end,1);
   
   targetWeightsAll_aud = cat(1,targetWeightsAll_aud,B_targets);
   detectWeightsAll_aud = cat(1,detectWeightsAll_aud,B_yeses);
  
   targetWeightsAll_shuf_aud = cat(1,targetWeightsAll_shuf_aud,B_targets_shuf);
   detectWeightsAll_shuf_aud = cat(1,detectWeightsAll_shuf_aud,B_yeses_shuf);
end

%%
ind = targetWeightsAll_vis < 0;
targetWeightsBinned_vis{1} = targetWeightsAll_vis(ind);
detectWeightsBinned_vis{1} = detectWeightsAll_vis(ind);
targetWeightsBinned_vis{2} = targetWeightsAll_vis(~ind);
detectWeightsBinned_vis{2} = detectWeightsAll_vis(~ind);

ind = targetWeightsAll_aud < 0;
targetWeightsBinned_aud{1} = targetWeightsAll_aud(ind);
detectWeightsBinned_aud{1} = detectWeightsAll_aud(ind);
targetWeightsBinned_aud{2} = targetWeightsAll_aud(~ind);
detectWeightsBinned_aud{2} = detectWeightsAll_aud(~ind);
%%
% figure
% for iexp = 1:nexp
%     h = plot(pcv{iexp});
%     hold on
%     
% end
% hline(0.95,'k--')
% figXAxis([],'N Components',[])
% figYAxis([],'% Var Expl.',[])
% figAxForm
% print(fullfile(fnout,'pca_pctVarExplained'),'-dpdf','-fillpage')
%%
setFigParams4Print('portrait')
set(0,'defaultAxesFontSize',12)
weightLim = [-2 5];
corrLim = [-1 1];
figure
suptitle('Visual Trials')
subplot 321
scatter(detectWeightsAll_vis,detectCorrsAll_vis,'bo')
hold on
hline(0,'k:')
vline(0,'k:')
figXAxis([],'Detect LR Weight',weightLim)
figYAxis([],'Detect Correlation',corrLim)
figAxForm
title('All Trials')
subplot 323
scatter(targetWeightsAll_vis,targetCorrsAll_vis,'bo')
hold on
figXAxis([],'Target LR Weight',weightLim)
figYAxis([],'Target Correlation',corrLim)
hline(0,'k:')
vline(0,'k:')
figAxForm
subplot 325   
s = scatter(targetWeightsAll_vis,detectWeightsAll_vis,'bo');
s.MarkerFaceColor = [1 1 1];
hold on
figXAxis([],'Target LR Weight',weightLim)
figYAxis([],'Detect LR Weight',weightLim)
hline(0,'k:')
vline(0,'k:')
figAxForm
subplot 326   
s = scatter(targetWeightsAll_shuf_vis,detectWeightsAll_shuf_vis,'bo');
s.MarkerFaceColor = [1 1 1];
hold on
figXAxis([],'Target LR Weight',weightLim)
figYAxis([],'Detect LR Weight',weightLim)
hline(0,'k:')
vline(0,'k:')
figAxForm

subplot 322
p = nan(1,2);
for i = 1:2
    x = mean(targetWeightsBinned_vis{i});
    xerr = ste(targetWeightsBinned_vis{i},1);
    y = detectWeightsBinned_vis{i};
    yerr = ste(y,1);
    [ttest_pass,p(i)] = ttest(y,[],'alpha',0.025);
    plot(x,y,'k.');
    hold on
    h = errorbar(x,mean(y),yerr,yerr,xerr,xerr,'ko');
    if ttest_pass
        h.Color = 'r';
    end
end
figXAxis([],'Target LR Weight',weightLim)
figYAxis([],'Detect LR Weight',weightLim)
hline(0,'k:')
vline(0,'k:')
figAxForm
title(sprintf('hi/lo target, p = %s/%s',round(p,2,'significant')))

print(fullfile(fnout,'corrXWeight_shuffle_vis'),'-dpdf','-fillpage')

figure
suptitle('Auditory Trials')
subplot 321
scatter(detectWeightsAll_aud,detectCorrsAll_aud,'bo')
hold on
figXAxis([],'Detect LR Weight',weightLim)
figYAxis([],'Detect Correlation',corrLim)
figAxForm
hline(0,'k:')
vline(0,'k:')
title('All Trials')
subplot 323
scatter(targetWeightsAll_aud,targetCorrsAll_aud,'bo')
hold on
figXAxis([],'Target LR Weight',weightLim)
figYAxis([],'Target Correlation',corrLim)
hline(0,'k:')
vline(0,'k:')
figAxForm
subplot 325   
s = scatter(targetWeightsAll_aud,detectWeightsAll_aud,'bo');
s.MarkerFaceColor = [1 1 1];
hold on
figXAxis([],'Target LR Weight',weightLim)
figYAxis([],'Detect LR Weight',weightLim)
hline(0,'k:')
vline(0,'k:')
figAxForm
subplot 326   
s = scatter(targetWeightsAll_shuf_aud,detectWeightsAll_shuf_aud,'bo');
s.MarkerFaceColor = [1 1 1];
hold on
figXAxis([],'Target LR Weight',weightLim)
figYAxis([],'Detect LR Weight',weightLim)
hline(0,'k:')
vline(0,'k:')
figAxForm


subplot 322
p = nan(1,2);
for i = 1:2
    x = mean(targetWeightsBinned_aud{i});
    xerr = ste(targetWeightsBinned_aud{i},1);
    y = detectWeightsBinned_aud{i};
    yerr = ste(y,1);
    [ttest_pass,p(i)] = ttest(y,[],'alpha',0.025);
    plot(x,y,'k.');
    hold on
    h = errorbar(x,mean(y),yerr,yerr,xerr,xerr,'ko');
    if ttest_pass
        h.Color = 'r';
    end
end
figXAxis([],'Target LR Weight',weightLim)
figYAxis([],'Detect LR Weight',weightLim)
hline(0,'k:')
vline(0,'k:')
figAxForm
title(sprintf('hi/lo target, p = %s/%s',round(p,2,'significant')))

print(fullfile(fnout,'corrXWeight_shuffle_aud'),'-dpdf','-fillpage')