close all
clear all

load('Z:\Analysis\FSAV Choice\FSAV_decodeData.mat')
fnout = 'Z:\Analysis\FSAV Choice';
nexp = length(dcExpt);

doExptPlots = 0;
%%
theta90Threshold = 11.25;
decisionVariable = 0.5;
detectThresholdHR = 0.8;
detectThresholdOri = 45;
oriBinsFine = [0 1 16 22.5 32 45 64 90];
oriBins = [0 1 45 90];
weightsDiffZeroAlpha = 0.05/2;
%%
targetWeightsAll = [];
detectWeightsAll = [];
targetWeightsVar = [];
detectWeightsVar = [];
nCells = [];
targetCorrsAll = [];
detectCorrsAll = [];
oriPrefAll = [];
for iexp=1:nexp
   trOut=dcExpt(iexp).trialOutcome;
    [Y_yeses, Y_targets] = getStimAndBehaviorYs(trOut);

   X_targets=dcExpt(iexp).stimResp';
   
   cellIdx = dcExpt(iexp).signifResponsiveCells' & ...
       dcExpt(iexp).oriFitTheta90 < theta90Threshold;

   X_targets=X_targets(:,cellIdx);
   
   idx=find(~isnan(sum(X_targets,2)));
   Y_yeses=Y_yeses(idx,1);   
   Y_targets=Y_targets(idx,1);   
   X_targets=X_targets(idx,:);
   
%    X_targets = zscore(X_targets);
 
   detectCorr = corr(Y_yeses,X_targets)';
   targetCorr = corr(Y_targets,X_targets)';
   
   detectCorrsAll = cat(1,detectCorrsAll,detectCorr);
   targetCorrsAll = cat(1,targetCorrsAll,targetCorr);
   
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
   
   targetWeightsAll = cat(1,targetWeightsAll,B_targets);
   detectWeightsAll = cat(1,detectWeightsAll,B_yeses);
   
   targetWeightsVar = cat(1,targetWeightsVar,var(B_targets));
   detectWeightsVar = cat(1,detectWeightsVar,var(B_yeses));
   nCells = cat(1,nCells,length(B_targets));
      
%    SI=dcExpt(iexp).baseTargRatio;
   
   oriPref = dcExpt(iexp).oriPref(cellIdx);
   oriPrefAll = cat(2,oriPrefAll,oriPref);

   if doExptPlots
    figure
    suptitle([dcExpt(iexp).mouse '-' dcExpt(iexp).date])
    subplot 221
    scatter(B_yeses,detectCorr,'bo')
    hold on
    hline(0,'k:')
    vline(0,'k:')
    figXAxis([],'Detect LR Weight',[])
    figYAxis([],'Detect Correlation',[])
    figAxForm
    subplot 222
    scatter(B_targets,targetCorr,'bo')
    hold on
    hline(0,'k:')
    vline(0,'k:')
    figXAxis([],'Target LR Weight',[])
    figYAxis([],'Target Correlation',[])
    figAxForm
    subplot 223   
    s = scatter(B_targets,B_yeses,'bo');
    s.MarkerFaceColor = [1 1 1];
    hold on
    hline(0,'k:')
    vline(0,'k:')
    figXAxis([],'Target LR Weight',[])
    figYAxis([],'Detect LR Weight',[])
    figAxForm
    subplot 224
    axis off
    text(0.5,1,sprintf('N Cells = %s',num2str(sum(cellIdx))))
    hold on
    text(0.5, 0.5,sprintf('N Trials = %s',num2str(length(Y_yeses))))
    figXAxis([],'',[0 2])
    figYAxis([],'',[0 2])
    
    print(fullfile(fnout,...
        [dcExpt(iexp).mouse '_' dcExpt(iexp).date '_corrXweight']),...
        '-dpdf','-fillpage')

    figure
    suptitle([dcExpt(iexp).mouse '-' dcExpt(iexp).date])
    subplot 221
    scatter(oriPref,B_targets,'bo')
    hold on
    hline(0,'k:')
    vline(0,'k:')
    figXAxis([],'Orientation Preference',[-10 190],[0:45:180],[0:45:180])
    figYAxis([],'Target LR Weight',[])
    figAxForm
    subplot 222
    scatter(oriPref,targetCorr,'bo')
    hold on
    hline(0,'k:')
    vline(0,'k:')
    figXAxis([],'Orientation Preference',[-10 190],[0:45:180],[0:45:180])
    figYAxis([],'Target Correlation',[])
    figAxForm
    subplot 223
    scatter(oriPref,B_yeses,'bo')
    hold on
    hline(0,'k:')
    vline(0,'k:')
    figXAxis([],'Orientation Preference',[-10 190],[0:45:180],[0:45:180])
    figYAxis([],'Detect LR Weight',[])
    figAxForm
    subplot 224
    scatter(oriPref,detectCorr,'bo')
    hold on
    hline(0,'k:')
    vline(0,'k:')
    figXAxis([],'Orientation Preference',[-10 190],[0:45:180],[0:45:180])
    figYAxis([],'Detect Correlation',[])
    figAxForm
    
    print(fullfile(fnout,...
        [dcExpt(iexp).mouse '_' dcExpt(iexp).date '_oriPrefXcorrXweight']),...
        '-dpdf','-fillpage')
   end
end

figure
subplot 121
scatter(nCells,targetWeightsVar,'ko')
figXAxis([],'N Cells',[])
figYAxis([],'Target Weight Variance',[0 0.7])
figAxForm
subplot 122
scatter(nCells,detectWeightsVar,'ko')
figXAxis([],'N Cells',[])
figYAxis([],'Detect Weight Variance',[0 0.7])
figAxForm
print(fullfile(fnout,'weightVarianceXn'),'-dpdf','-fillpage')

figure
subplot 221
scatter(detectWeightsAll,detectCorrsAll,'bo')
hold on
hline(0,'k:')
vline(0,'k:')
figXAxis([],'Detect LR Weight',[])
figYAxis([],'Detect Correlation',[])
figAxForm
subplot 222
scatter(targetWeightsAll,targetCorrsAll,'bo')
hold on
hline(0,'k:')
vline(0,'k:')
figXAxis([],'Target LR Weight',[])
figYAxis([],'Target Correlation',[])
figAxForm
subplot 223   
s = scatter(targetWeightsAll,detectWeightsAll,'bo');
s.MarkerFaceColor = [1 1 1];
hold on
hline(0,'k:')
vline(0,'k:')
figXAxis([],'Target LR Weight',[])
figYAxis([],'Detect LR Weight',[])
figAxForm

print(fullfile(fnout,'corrXweight'),'-dpdf','-fillpage')

figure
suptitle(sprintf('Tuning Cutoff = %s', num2str(theta90Threshold)))
subplot 221
scatter(oriPrefAll,targetWeightsAll,'bo')
hold on
hline(0,'k:')
vline(0,'k:')
figXAxis([],'Orientation Preference',[-10 190],[0:45:180],[0:45:180])
figYAxis([],'Target LR Weight',[])
figAxForm
subplot 222
scatter(oriPrefAll,targetCorrsAll,'bo')
hold on
hline(0,'k:')
vline(0,'k:')
figXAxis([],'Orientation Preference',[-10 190],[0:45:180],[0:45:180])
figYAxis([],'Target Correlation',[])
figAxForm
subplot 223
scatter(oriPrefAll,detectWeightsAll,'bo')
hold on
hline(0,'k:')
vline(0,'k:')
figXAxis([],'Orientation Preference',[-10 190],[0:45:180],[0:45:180])
figYAxis([],'Detect LR Weight',[])
figAxForm
subplot 224
scatter(oriPrefAll,detectCorrsAll,'bo')
hold on
hline(0,'k:')
vline(0,'k:')
figXAxis([],'Orientation Preference',[-10 190],[0:45:180],[0:45:180])
figYAxis([],'Detect Correlation',[])
figAxForm

print(fullfile(fnout,'oriPrefXcorrXweight'),'-dpdf','-fillpage')

%% test significance of weights
% DW = detect weight
% TW = target weight
% DC = detect corr
% TC = target corr

DWvsTWdata = cell(1,2);
DWvsTCdata = cell(1,2);
TWvsTCdata = cell(1,2);

ind = targetWeightsAll > 0;
DWvsTWdata{1} = detectWeightsAll(~ind);
DWvsTWdata{2} = detectWeightsAll(ind);
TWmean(1) = mean(targetWeightsAll(~ind));
TWmean(2) = mean(targetWeightsAll(ind));

ind = targetCorrsAll > 0;
DWvsTCdata{1} = detectWeightsAll(~ind);
DWvsTCdata{2} = detectWeightsAll(ind);
TCmean(1) = mean(targetCorrsAll(~ind));
TCmean(2) = mean(targetCorrsAll(ind));

TWvsTCdata{1} = targetWeightsAll(~ind);
TWvsTCdata{2} = targetWeightsAll(ind);

figure
suptitle(sprintf('Red is signif from zero, alpha = %s',num2str(weightsDiffZeroAlpha)))
subplot 221
p = nan(1,2);
for i = 1:2
    x = TWmean(i);
    y = DWvsTWdata{i};
    yerr = ste(y,1);
    plot(x,y,'k.')
    hold on
    [ttest_pass,p(i)] = ttest(y,[],'alpha',weightsDiffZeroAlpha);
    if ttest_pass
        h = errorbar(x,mean(y),yerr,'ro');
        h.MarkerFaceColor = [1 1 1];
    else
        h = errorbar(x,mean(y),yerr,'ko');
        h.MarkerFaceColor = [1 1 1];
    end
end
figXAxis([],'Target Weight',[-1 1],TWmean,{'< 0';'> 0'})
figYAxis([],'Detect Weight',[])
figAxForm
hline(0,'k--')
vline(0,'k--')
subplot 222
p = nan(1,2);
for i = 1:2
    x = TCmean(i);
    y = DWvsTCdata{i};
    yerr = ste(y,1);
    plot(x,y,'k.')
    hold on
    [ttest_pass,p(i)] = ttest(y,[],'alpha',weightsDiffZeroAlpha);
    if ttest_pass
        h = errorbar(x,mean(y),yerr,'ro');
        h.MarkerFaceColor = [1 1 1];
    else
        h = errorbar(x,mean(y),yerr,'ko');
        h.MarkerFaceColor = [1 1 1];
    end
end
figXAxis([],'Target Correlation',[-1 1],TWmean,{'< 0';'> 0'})
figYAxis([],'Detect Weight',[])
figAxForm
hline(0,'k--')
vline(0,'k--')
subplot 223
p = nan(1,2);
for i = 1:2
    x = TCmean(i);
    y = TWvsTCdata{i};
    yerr = ste(y,1);
    plot(x,y,'k.')
    hold on
    [ttest_pass,p(i)] = ttest(y,[],'alpha',weightsDiffZeroAlpha);
    if ttest_pass
        h = errorbar(x,mean(y),yerr,'ro');
        h.MarkerFaceColor = [1 1 1];
    else
        h = errorbar(x,mean(y),yerr,'ko');
        h.MarkerFaceColor = [1 1 1];
    end
end
figXAxis([],'Target Correlation',[-1 1],TWmean,{'< 0';'> 0'})
figYAxis([],'Target Weight',[])
figAxForm
hline(0,'k--')
vline(0,'k--')

print(fullfile(fnout,'testPosAndNegWeights'),'-dpdf','-fillpage')
%% fit with hard trials

targetWeightsAll_hard = [];
detectWeightsAll_hard = [];
targetCorrsAll_hard = [];
detectCorrsAll_hard = [];
oriPrefAll = [];
for iexp = 1:nexp
   trOut=dcExpt(iexp).trialOutcome;
    tOri = dcExpt(iexp).trialOrientation;
    tOri = binAndRelabelTrialOrientation(tOri,oriBinsFine);
    
    hardTrialIndex = tOri <= detectThresholdOri;
    
    [Y_yeses, Y_targets] = getStimAndBehaviorYs(trOut(hardTrialIndex));

   X_targets=dcExpt(iexp).stimResp(:,hardTrialIndex)';
   
   cellIdx = dcExpt(iexp).signifResponsiveCells' & ...
       dcExpt(iexp).oriFitTheta90 < theta90Threshold;

   X_targets=X_targets(:,cellIdx);
   
   idx=find(~isnan(sum(X_targets,2)));
   Y_yeses=Y_yeses(idx,1);   
   Y_targets=Y_targets(idx,1);   
   X_targets=X_targets(idx,:);
   
%    X_targets = zscore(X_targets);
 
   detectCorr = corr(Y_yeses,X_targets)';
   targetCorr = corr(Y_targets,X_targets)';
   
   detectCorrsAll_hard = cat(1,detectCorrsAll_hard,detectCorr);
   targetCorrsAll_hard = cat(1,targetCorrsAll_hard,targetCorr);
   
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
   
   targetWeightsAll_hard = cat(1,targetWeightsAll_hard,B_targets);
   detectWeightsAll_hard = cat(1,detectWeightsAll_hard,B_yeses);
   
   oriPref = dcExpt(iexp).oriPref(cellIdx);
   oriPrefAll = cat(2,oriPrefAll,oriPref);
end
figure
suptitle(sprintf('Targets < %s deg',num2str(detectThresholdOri)))
subplot 221
scatter(oriPrefAll,targetWeightsAll_hard,'bo')
hold on
hline(0,'k:')
vline(0,'k:')
figXAxis([],'Orientation Preference',[-10 190],[0:45:180],[0:45:180])
figYAxis([],'Target LR Weight',[])
figAxForm
subplot 222
scatter(oriPrefAll,targetCorrsAll_hard,'bo')
hold on
hline(0,'k:')
vline(0,'k:')
figXAxis([],'Orientation Preference',[-10 190],[0:45:180],[0:45:180])
figYAxis([],'Target Correlation',[])
figAxForm
subplot 223
scatter(oriPrefAll,detectWeightsAll_hard,'bo')
hold on
hline(0,'k:')
vline(0,'k:')
figXAxis([],'Orientation Preference',[-10 190],[0:45:180],[0:45:180])
figYAxis([],'Detect LR Weight',[])
figAxForm
subplot 224
scatter(oriPrefAll,detectCorrsAll_hard,'bo')
hold on
hline(0,'k:')
vline(0,'k:')
figXAxis([],'Orientation Preference',[-10 190],[0:45:180],[0:45:180])
figYAxis([],'Detect Correlation',[])
figAxForm

print(fullfile(fnout,'oriPrefXcorrXweight_hardTrialsOri'),'-dpdf','-fillpage')

figure
subplot 221
scatter(detectWeightsAll_hard,detectCorrsAll_hard,'bo')
hold on
hline(0,'k:')
vline(0,'k:')
figXAxis([],'Detect LR Weight',[])
figYAxis([],'Detect Correlation',[])
figAxForm
subplot 222
scatter(targetWeightsAll_hard,targetCorrsAll_hard,'bo')
hold on
hline(0,'k:')
vline(0,'k:')
figXAxis([],'Target LR Weight',[])
figYAxis([],'Target Correlation',[])
figAxForm
subplot 223   
s = scatter(targetWeightsAll_hard,detectWeightsAll_hard,'bo');
s.MarkerFaceColor = [1 1 1];
hold on
hline(0,'k:')
vline(0,'k:')
figXAxis([],'Target LR Weight',[])
figYAxis([],'Detect LR Weight',[])
figAxForm
suptitle(sprintf('Targets < %s deg',num2str(detectThresholdOri)))

print(fullfile(fnout,'corrXweight_hardTrialsOri'),'-dpdf','-fillpage')

weightLim = [-1 2.5];
figure
suptitle(sprintf('Targets < %s deg',num2str(detectThresholdOri)))
subplot 121
scatter(detectWeightsAll,detectWeightsAll_hard)
hold on
hline(0,'k:')
vline(0,'k:')
plot(weightLim,weightLim,'k--')
title('Detect Weights')
figXAxis([],'All Trials',weightLim)
figYAxis([],'Hard Trials',weightLim)
figAxForm
subplot 122
scatter(targetWeightsAll,targetWeightsAll_hard)
hold on
hline(0,'k:')
vline(0,'k:')
plot(weightLim,weightLim,'k--')
title('Target Weights')
figXAxis([],'All Trials',weightLim)
figYAxis([],'Hard Trials',weightLim)
figAxForm

print(fullfile(fnout,'hardXallTrialsWeights'),'-dpdf','-fillpage')

%% predictions and cross-validation, binned orientation

pctCorrOriBinned_detect_allData = cell(1,nexp);
pctCorrOriBinned_detect_holdout = cell(1,nexp);
pctCorrOriBinned_target_allData = cell(1,nexp);
pctCorrOriBinned_target_holdout = cell(1,nexp);
exptOris = cell(1,nexp);
for iexp = 1:nexp

   trOut=dcExpt(iexp).trialOutcome;

   [Y_yeses, Y_targets] = getStimAndBehaviorYs(trOut);
   
   X_targets=dcExpt(iexp).stimResp';
   
   cellIdx = dcExpt(iexp).signifResponsiveCells' & ...
       dcExpt(iexp).oriFitTheta90 < theta90Threshold;

   X_targets=X_targets(:,cellIdx);
   
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
   
   yhat_target = glmval(targetGLM.beta,X_targets,'logit');
   yhat_detect = glmval(detectGLM.beta,X_targets,'logit');
   nTrials_all = size(X_targets,1);
   
   pctCorr_target = sum(yhat_target > decisionVariable)./nTrials_all;
   [~,pctCorr_ci] = binofit(sum(yhat_target > decisionVariable),nTrials_all);
   pctCorrTarget_err(:,1) = pctCorr_target' - pctCorr_ci(:,1);
   pctCorrTarget_err(:,2) = pctCorr_ci(:,2) - pctCorr_target';
   pctCorr_detect = sum(yhat_detect > decisionVariable)./nTrials_all;
   [~,pctCorr_ci] = binofit(sum(yhat_detect > decisionVariable),nTrials_all);
   pctCorrDetect_err(:,1) = pctCorr_detect' - pctCorr_ci(:,1);
   pctCorrDetect_err(:,2) = pctCorr_ci(:,2) - pctCorr_detect';
   
   
    yhat_targetHoldouts = nan(1,nTrials_all);
    yhat_detectHoldouts = nan(1,nTrials_all);
    for i = 1:nTrials_all
        X_holdout = X_targets(i,:);
        othersInd = [1:(i-1),(i+1):nTrials_all];
        X_others = X_targets(othersInd,:);
        
        [~,~,targetOthersGLM] = glmfit(X_others,Y_targets(othersInd),'binomial');
        yhat_targetHoldouts(i) = glmval(targetOthersGLM.beta,X_holdout,'logit');
        
        [~,~,detectOthersGLM] = glmfit(X_others,Y_yeses(othersInd),'binomial');
        yhat_detectHoldouts(i) = glmval(detectOthersGLM.beta,X_holdout,'logit');        
    end
    pctCorr_targetXVal = sum(yhat_targetHoldouts > decisionVariable)./nTrials_all;
    pctCorr_detectXVal = sum(yhat_detectHoldouts > decisionVariable)./nTrials_all;
    
    tOri = dcExpt(iexp).trialOrientation;
    [tOri,meanBinnedOri] = binAndRelabelTrialOrientation(tOri,oriBins);
    orientations = unique(tOri);
    nOri = length(orientations);
    
    pctCorrOri_detect = nan(1,nOri);
    pctCorrXValOri_detect = nan(1,nOri);
    pctCorrOri_target = nan(1,nOri);
    pctCorrXValOri_target = nan(1,nOri);    
    nOriHits_detect = nan(1,nOri);
    nOriHits_target = nan(1,nOri);
    nOriTrials = nan(1,nOri);
    for iori = 1:nOri
        trInd = find(tOri == orientations(iori));
        X = X_targets(trInd,:);
        
        yhat_target = glmval(targetGLM.beta,X,'logit');
        pctCorrOri_target(iori) = mean(yhat_target > decisionVariable);
        nOriHits_target(iori) = sum(yhat_target> decisionVariable);
        
        yhat_detect = glmval(detectGLM.beta,X,'logit');
        pctCorrOri_detect(iori) = mean(yhat_detect > decisionVariable);
        nOriHits_detect(iori) = sum(yhat_detect> decisionVariable);
        
        nTrials = length(trInd);
        nOriTrials(iori) = nTrials;
        yhat_detectHoldouts = nan(1,nTrials);
        yhat_targetHoldouts = nan(1,nTrials);
        for i = 1:nTrials
            X_holdout = X(i,:);
            oriTrInd = trInd(i);
            othersInd = [1:(oriTrInd-1),(oriTrInd+1):nTrials_all];
            X_others = X_targets(othersInd,:);
            
            [~,~,detectOthersGLM] = glmfit(...
                X_others,Y_yeses(othersInd),'binomial');
            yhat_detectHoldouts(i) = glmval(...
                detectOthersGLM.beta,X_holdout,'logit');
            
            [~,~,targetOthersGLM] = glmfit(...
                X_others,Y_targets(othersInd),'binomial');
            yhat_targetHoldouts(i) = glmval(...
                targetOthersGLM.beta,X_holdout,'logit');
            
        end
        pctCorrXValOri_detect(iori) = mean(yhat_detectHoldouts > decisionVariable);
        pctCorrXValOri_target(iori) = mean(yhat_targetHoldouts > decisionVariable);
        
    end
    [~,pctCorrOri_ci] = binofit(nOriHits_detect,nOriTrials);
    pctCorrOri_detectErr = nan(length(nOriTrials),2);
    pctCorrOri_detectErr(:,1) = pctCorrOri_detect' - pctCorrOri_ci(:,1);
    pctCorrOri_detectErr(:,2) = pctCorrOri_ci(:,2) - pctCorrOri_detect';
    pctCorr_detectErr = cat(1,pctCorrOri_detectErr,pctCorrDetect_err);
    [~,pctCorrOri_ci] = binofit(nOriHits_target,nOriTrials);
    pctCorrOri_targetErr = nan(length(nOriTrials),2);
    pctCorrOri_targetErr(:,1) = pctCorrOri_target' - pctCorrOri_ci(:,1);
    pctCorrOri_targetErr(:,2) = pctCorrOri_ci(:,2) - pctCorrOri_target';
    pctCorr_targetErr = cat(1,pctCorrOri_targetErr,pctCorrTarget_err);
    
    figure
    subplot 121
    x = [orientations 100];
    y = [pctCorrOri_detect.*100 pctCorr_detect.*100];
    yuerr = pctCorr_detectErr(:,1)'.*100;
    ylerr = pctCorr_detectErr(:,2)'.*100;
    h = errorbar(x,y,yuerr,ylerr,'ko-');
    h.MarkerFaceColor = [0 0 0];
    hold on
    y = [pctCorrXValOri_detect.*100 pctCorr_detectXVal.*100];
    h = plot(x,y,'ko');
    h.MarkerFaceColor = [1 1 1];
    xLabel = [cellfun(@(y)...
        num2str(round(y,2,'significant')),num2cell(x(1:end-1)),'unif',0),...
        {'All Trials'}];
    figXAxis([],'Orientation (deg)',[-10 110],x,xLabel)
    figYAxis([],'% Yes',[0 110])
    figAxForm
    legend({'All Trials';'Holdout X Val'},'location','northeastoutside')
    title({'Detect GLM';[dcExpt(iexp).mouse '-' dcExpt(iexp).date]})
   
    subplot 122
    x = [orientations 100];
    y = [pctCorrOri_target.*100 pctCorr_target.*100];
    yuerr = pctCorr_targetErr(:,1)'.*100;
    ylerr = pctCorr_targetErr(:,2)'.*100;
    h = errorbar(x,y,yuerr,ylerr,'ko-');
    h.MarkerFaceColor = [0 0 0];
    hold on
    y = [pctCorrXValOri_target.*100 pctCorr_targetXVal.*100];
    h = plot(x,y,'ko');
    h.MarkerFaceColor = [1 1 1];
    xLabel = [cellfun(@(y)...
        num2str(round(y,2,'significant')),num2cell(x(1:end-1)),'unif',0),...
        {'All Trials'}];
    figXAxis([],'Orientation (deg)',[-10 110],x,xLabel)
    figYAxis([],'% Target Present',[0 110])
    figAxForm
    legend({'All Trials';'Holdout X Val'},'location','northeastoutside')
    title({'Target GLM';[dcExpt(iexp).mouse '-' dcExpt(iexp).date]})

    pctCorrOriBinned_detect_allData{iexp} = [pctCorrOri_detect.*100 pctCorr_detect.*100];
    pctCorrOriBinned_detect_holdout{iexp} = [pctCorrXValOri_detect.*100 pctCorr_detectXVal.*100];
    
    pctCorrOriBinned_target_allData{iexp} = [pctCorrOri_target.*100 pctCorr_target.*100];
    pctCorrOriBinned_target_holdout{iexp} = [pctCorrXValOri_target.*100 pctCorr_targetXVal.*100];
    
    exptOris{iexp} = x;
    
    print(fullfile(fnout,[dcExpt(iexp).mouse '-' dcExpt(iexp).date '_xValidation_binnedOri']),'-dpdf','-fillpage')
end
    
figure
suptitle('Summary Across Experiments, Binned Target Orientation')
subplot 231
allOris = [];
for i = 1:nexp
    plot(exptOris{i},pctCorrOriBinned_detect_allData{i},'o-')
    hold on
    allOris = unique([allOris,exptOris{i}]);
end
xLabel = [cellfun(@(y)...
    num2str(round(y,2,'significant')),num2cell(allOris(1:end-1)),'unif',0),...
    {'All Trials'}];
figXAxis([],'Orientation (deg)',[-10 110],allOris,xLabel)
figYAxis([],'% Yes',[0 110])
figAxForm
title('Detect GLM, All Trials')
pctCorrOriBinned_catAllData = cell(1,length(allOris));
for iexp = 1:nexp
    for i = 1:length(allOris)
        ind = exptOris{iexp} == allOris(i);
        pctCorrOriBinned_catAllData{i} = ...
            [pctCorrOriBinned_catAllData{i} pctCorrOriBinned_detect_allData{iexp}(ind)];
    end
end
subplot 232
x = 1:length(allOris);
y1data = cell2mat(pctCorrOriBinned_detect_allData');
y1 = mean(y1data,1);
y1err = ste(y1data,1);
y2data = cell2mat(pctCorrOriBinned_detect_holdout');
y2 = mean(y2data,1);
y2err = ste(y2data,1);
h = bar(x,[y1' y2'],'grouped');
h(1).EdgeColor = 'k';
h(1).FaceColor = 'k';
h(2).EdgeColor = 'k';
h(2).FaceColor = [1 1 1];
hold on
errorbar(x-0.15,y1,y1err,'k.');
errorbar(x+0.15,y2,y2err,'k.');
legend(h,{'Predicted';'X-Valid'},'location','northwest')
figXAxis([],'Orientation (deg)',[0 length(x)+1],x,xLabel)
figYAxis([],'% Yes',[0 110])
figAxForm
subplot 233
x = 1:length(allOris);
y1data = cell2mat(pctCorrOriBinned_detect_allData');
y2data = cell2mat(pctCorrOriBinned_detect_holdout');
ydata = y1data-y2data;
y = mean(ydata,1);
yerr = ste(ydata,1);
h = errorbar(x,y,yerr,'ko');
h.MarkerFaceColor = 'k';
hold on
hline(0,'k:')
for i = 1:nexp
    y = ydata(:,i);
    plot(x(i),y,'k.')
end
legend(h,{'Predicted';'X-Valid'},'location','northwest')
figXAxis([],'Orientation (deg)',[0 length(x)+1],x,xLabel)
figYAxis([],'% Yes',[0 110])
figAxForm


subplot 234
allOris = [];
for i = 1:nexp
    plot(exptOris{i},pctCorrOriBinned_target_allData{i},'o-')
    hold on
    allOris = unique([allOris,exptOris{i}]);
end
xLabel = [cellfun(@(y)...
    num2str(round(y,2,'significant')),num2cell(allOris(1:end-1)),'unif',0),...
    {'All Trials'}];
figXAxis([],'Orientation (deg)',[-10 110],allOris,xLabel)
figYAxis([],'% Target Present',[0 110])
figAxForm
title('Target GLM, All Trials')
pctCorrOriBinned_catAllData = cell(1,length(allOris));
for iexp = 1:nexp
    for i = 1:length(allOris)
        ind = exptOris{iexp} == allOris(i);
        pctCorrOriBinned_catAllData{i} = ...
            [pctCorrOriBinned_catAllData{i} pctCorrOriBinned_target_allData{iexp}(ind)];
    end
end
subplot 235
x = 1:length(allOris);
y1data = cell2mat(pctCorrOriBinned_target_allData');
y1 = mean(y1data,1);
y1err = ste(y1data,1);
y2data = cell2mat(pctCorrOriBinned_target_holdout');
y2 = mean(y2data,1);
y2err = ste(y2data,1);
h = bar(x,[y1' y2'],'grouped');
h(1).EdgeColor = 'k';
h(1).FaceColor = 'k';
h(2).EdgeColor = 'k';
h(2).FaceColor = [1 1 1];
hold on
errorbar(x-0.15,y1,y1err,'k.');
errorbar(x+0.15,y2,y2err,'k.');
legend(h,{'Predicted';'X-Valid'},'location','northwest')
figXAxis([],'Orientation (deg)',[0 length(x)+1],x,xLabel)
figYAxis([],'% Target Present',[0 110])
figAxForm

print(fullfile(fnout,'pctCorr_byOriBinned'),'-dpdf','-fillpage')

%% Auditory trials

targetWeightsAll_aud = [];
detectWeightsAll_aud = [];
targetWeightsVar = [];
detectWeightsVar = [];
nCells = [];
targetCorrsAll_aud = [];
detectCorrsAll_aud = [];
for iexp=1:nexp
   trOut=dcExpt(iexp).audTrialOutcome;
    [Y_yeses, Y_targets] = getStimAndBehaviorYs(trOut);

   X_targets=dcExpt(iexp).audStimResp';
   
   cellIdx = dcExpt(iexp).signifResponsiveCells' & ...
       dcExpt(iexp).oriFitTheta90 < theta90Threshold;

   X_targets=X_targets(:,cellIdx);
   
   idx=find(~isnan(sum(X_targets,2)));
   Y_yeses=Y_yeses(idx,1);   
   Y_targets=Y_targets(idx,1);   
   X_targets=X_targets(idx,:);
   
%    X_targets = zscore(X_targets);
 
   detectCorr = corr(Y_yeses,X_targets)';
   targetCorr = corr(Y_targets,X_targets)';
   
   detectCorrsAll_aud = cat(1,detectCorrsAll_aud,detectCorr);
   targetCorrsAll_aud = cat(1,targetCorrsAll_aud,targetCorr);
   
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
   
   targetWeightsAll_aud = cat(1,targetWeightsAll_aud,B_targets);
   detectWeightsAll_aud = cat(1,detectWeightsAll_aud,B_yeses);
   
   targetWeightsVar = cat(1,targetWeightsVar,var(B_targets));
   detectWeightsVar = cat(1,detectWeightsVar,var(B_yeses));
   nCells = cat(1,nCells,length(B_targets));
      
%    SI=dcExpt(iexp).baseTargRatio;
%    
   oriPref = dcExpt(iexp).oriPref(cellIdx);
%    oriPrefAll = cat(2,oriPrefAll,oriPref);

   if doExptPlots
    figure
    suptitle([dcExpt(iexp).mouse '-' dcExpt(iexp).date])
    subplot 221
    scatter(B_yeses,detectCorr,'bo')
    hold on
    hline(0,'k:')
    vline(0,'k:')
    figXAxis([],'Detect LR Weight',[])
    figYAxis([],'Detect Correlation',[])
    figAxForm
    subplot 222
    scatter(B_targets,targetCorr,'bo')
    hold on
    hline(0,'k:')
    vline(0,'k:')
    figXAxis([],'Target LR Weight',[])
    figYAxis([],'Target Correlation',[])
    figAxForm
    subplot 223   
    s = scatter(B_targets,B_yeses,'bo');
    s.MarkerFaceColor = [1 1 1];
    hold on
    hline(0,'k:')
    vline(0,'k:')
    figXAxis([],'Target LR Weight',[])
    figYAxis([],'Detect LR Weight',[])
    figAxForm
    subplot 224
    axis off
    text(0.5,1,sprintf('N Cells = %s',num2str(sum(cellIdx))))
    hold on
    text(0.5, 0.5,sprintf('N Trials = %s',num2str(length(Y_yeses))))
    figXAxis([],'',[0 2])
    figYAxis([],'',[0 2])
    
    print(fullfile(fnout,...
        [dcExpt(iexp).mouse '_' dcExpt(iexp).date '_audTrials_corrXweight']),...
        '-dpdf','-fillpage')

    figure
    suptitle([dcExpt(iexp).mouse '-' dcExpt(iexp).date])
    subplot 221
    scatter(oriPref,B_targets,'bo')
    hold on
    hline(0,'k:')
    vline(0,'k:')
    figXAxis([],'Orientation Preference',[-10 190],[0:45:180],[0:45:180])
    figYAxis([],'Target LR Weight',[])
    figAxForm
    subplot 222
    scatter(oriPref,targetCorr,'bo')
    hold on
    hline(0,'k:')
    vline(0,'k:')
    figXAxis([],'Orientation Preference',[-10 190],[0:45:180],[0:45:180])
    figYAxis([],'Target Correlation',[])
    figAxForm
    subplot 223
    scatter(oriPref,B_yeses,'bo')
    hold on
    hline(0,'k:')
    vline(0,'k:')
    figXAxis([],'Orientation Preference',[-10 190],[0:45:180],[0:45:180])
    figYAxis([],'Detect LR Weight',[])
    figAxForm
    subplot 224
    scatter(oriPref,detectCorr,'bo')
    hold on
    hline(0,'k:')
    vline(0,'k:')
    figXAxis([],'Orientation Preference',[-10 190],[0:45:180],[0:45:180])
    figYAxis([],'Detect Correlation',[])
    figAxForm
    
    print(fullfile(fnout,...
        [dcExpt(iexp).mouse '_' dcExpt(iexp).date '_audTrials_oriPrefXcorrXweight']),...
        '-dpdf','-fillpage')
   end
end

figure
subplot 221
scatter(detectWeightsAll_aud,detectCorrsAll_aud,'bo')
hold on
hline(0,'k:')
vline(0,'k:')
figXAxis([],'Detect LR Weight',[])
figYAxis([],'Detect Correlation',[])
figAxForm
subplot 222
scatter(targetWeightsAll_aud,targetCorrsAll_aud,'bo')
hold on
hline(0,'k:')
vline(0,'k:')
figXAxis([],'Target LR Weight',[])
figYAxis([],'Target Correlation',[])
figAxForm
subplot 223   
s = scatter(targetWeightsAll_aud,detectWeightsAll_aud,'bo');
s.MarkerFaceColor = [1 1 1];
hold on
hline(0,'k:')
vline(0,'k:')
figXAxis([],'Target LR Weight',[])
figYAxis([],'Detect LR Weight',[])
figAxForm

print(fullfile(fnout,'audTrials_corrXweight'),'-dpdf','-fillpage')

figure
suptitle(sprintf('Tuning Cutoff = %s', num2str(theta90Threshold)))
subplot 221
scatter(oriPrefAll,targetWeightsAll_aud,'bo')
hold on
hline(0,'k:')
vline(0,'k:')
figXAxis([],'Orientation Preference',[-10 190],[0:45:180],[0:45:180])
figYAxis([],'Target LR Weight',[])
figAxForm
subplot 222
scatter(oriPrefAll,targetCorrsAll_aud,'bo')
hold on
hline(0,'k:')
vline(0,'k:')
figXAxis([],'Orientation Preference',[-10 190],[0:45:180],[0:45:180])
figYAxis([],'Target Correlation',[])
figAxForm
subplot 223
scatter(oriPrefAll,detectWeightsAll_aud,'bo')
hold on
hline(0,'k:')
vline(0,'k:')
figXAxis([],'Orientation Preference',[-10 190],[0:45:180],[0:45:180])
figYAxis([],'Detect LR Weight',[])
figAxForm
subplot 224
scatter(oriPrefAll,detectCorrsAll_aud,'bo')
hold on
hline(0,'k:')
vline(0,'k:')
figXAxis([],'Orientation Preference',[-10 190],[0:45:180],[0:45:180])
figYAxis([],'Detect Correlation',[])
figAxForm

print(fullfile(fnout,'_audTrials_oriPrefXcorrXweight'),'-dpdf','-fillpage')


DWvsTWdata = cell(1,2);
DWvsTCdata = cell(1,2);
TWvsTCdata = cell(1,2);

ind = targetWeightsAll_aud > 0;
DWvsTWdata{1} = detectWeightsAll_aud(~ind);
DWvsTWdata{2} = detectWeightsAll_aud(ind);
TWmean(1) = mean(targetWeightsAll_aud(~ind));
TWmean(2) = mean(targetWeightsAll_aud(ind));

ind = targetCorrsAll_aud > 0;
DWvsTCdata{1} = detectWeightsAll_aud(~ind);
DWvsTCdata{2} = detectWeightsAll_aud(ind);
TCmean(1) = mean(targetCorrsAll_aud(~ind));
TCmean(2) = mean(targetCorrsAll_aud(ind));

TWvsTCdata{1} = targetWeightsAll_aud(~ind);
TWvsTCdata{2} = targetWeightsAll_aud(ind);

figure
suptitle(sprintf('Red is signif from zero, alpha = %s',num2str(weightsDiffZeroAlpha)))
subplot 221
p = nan(1,2);
for i = 1:2
    x = TWmean(i);
    y = DWvsTWdata{i};
    yerr = ste(y,1);
    plot(x,y,'k.')
    hold on
    [ttest_pass,p(i)] = ttest(y,[],'alpha',weightsDiffZeroAlpha);
    if ttest_pass
        h = errorbar(x,mean(y),yerr,'ro');
        h.MarkerFaceColor = [1 1 1];
    else
        h = errorbar(x,mean(y),yerr,'ko');
        h.MarkerFaceColor = [1 1 1];
    end
end
figXAxis([],'Target Weight',[-1 1],TWmean,{'< 0';'> 0'})
figYAxis([],'Detect Weight',[])
figAxForm
hline(0,'k--')
vline(0,'k--')
subplot 222
p = nan(1,2);
for i = 1:2
    x = TCmean(i);
    y = DWvsTCdata{i};
    yerr = ste(y,1);
    plot(x,y,'k.')
    hold on
    [ttest_pass,p(i)] = ttest(y,[],'alpha',weightsDiffZeroAlpha);
    if ttest_pass
        h = errorbar(x,mean(y),yerr,'ro');
        h.MarkerFaceColor = [1 1 1];
    else
        h = errorbar(x,mean(y),yerr,'ko');
        h.MarkerFaceColor = [1 1 1];
    end
end
figXAxis([],'Target Correlation',[-1 1],TWmean,{'< 0';'> 0'})
figYAxis([],'Detect Weight',[])
figAxForm
hline(0,'k--')
vline(0,'k--')
subplot 223
p = nan(1,2);
for i = 1:2
    x = TCmean(i);
    y = TWvsTCdata{i};
    yerr = ste(y,1);
    plot(x,y,'k.')
    hold on
    [ttest_pass,p(i)] = ttest(y,[],'alpha',weightsDiffZeroAlpha);
    if ttest_pass
        h = errorbar(x,mean(y),yerr,'ro');
        h.MarkerFaceColor = [1 1 1];
    else
        h = errorbar(x,mean(y),yerr,'ko');
        h.MarkerFaceColor = [1 1 1];
    end
end
figXAxis([],'Target Correlation',[-1 1],TWmean,{'< 0';'> 0'})
figYAxis([],'Target Weight',[])
figAxForm
hline(0,'k--')
vline(0,'k--')

print(fullfile(fnout,'audTrials_testPosAndNegWeights'),'-dpdf','-fillpage')

corrLim = [-0.4 1];
weightLim = [-1.5 3];
figure
subplot 221
x = targetCorrsAll;
y = targetCorrsAll_aud;
scatter(x,y,'ko');
hold on
plot(corrLim,corrLim,'k--')
figXAxis([],'Visual Target Correlation',corrLim);
figYAxis([],'Auditory Target Correlation',corrLim);
figAxForm
hline(0,'k--')
vline(0,'k--')
subplot 222
x = detectCorrsAll;
y = detectCorrsAll_aud;
scatter(x,y,'ko');
hold on
plot(corrLim,corrLim,'k--')
figXAxis([],'Visual Detect Correlation',corrLim);
figYAxis([],'Auditory Detect Correlation',corrLim);
figAxForm
hline(0,'k--')
vline(0,'k--')

subplot 223
x = targetWeightsAll;
y = targetWeightsAll_aud;
scatter(x,y,'ko');
hold on
plot(weightLim,weightLim,'k--')
figXAxis([],'Visual Target Weight',weightLim);
figYAxis([],'Auditory Target Weight',weightLim);
figAxForm
hline(0,'k--')
vline(0,'k--')
subplot 224
x = detectWeightsAll;
y = detectWeightsAll_aud;
scatter(x,y,'ko');
hold on
plot(weightLim,weightLim,'k--')
figXAxis([],'Visual Detect Weight',weightLim);
figYAxis([],'Auditory Detect Weight',weightLim);
figAxForm
hline(0,'k--')
vline(0,'k--')

print(fullfile(fnout,'compareVisAudCorrAndWeight'),'-dpdf','-fillpage')