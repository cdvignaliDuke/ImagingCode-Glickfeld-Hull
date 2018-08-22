close all
clear all

load('Z:\Analysis\FSAV Choice\FSAV_decodeData.mat')
dataLabel = 'allExpt_';
fnout = 'Z:\Analysis\FSAV Choice\';
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
hardTrialCutoff = 40;
%%
targetWeightsAll_cells_vis = [];
detectWeightsAll_cells_vis = [];

targetCorrsAll_cells_vis = [];
detectCorrsAll_cells_vis = [];

targetWeightsAll_pca_vis = [];
detectWeightsAll_pca_vis = [];

targetCorrsAll_pca_vis = [];
detectCorrsAll_pca_vis = [];

nCells = zeros(1,nexp);
pcNumber = [];
pve = cell(1,nexp);
exptName = cell(1,nexp);
for iexp=1:nexp
   trOut=dcExpt(iexp).trialOutcome;
    [Y_yeses, Y_targets] = getStimAndBehaviorYs(trOut);
    

   X_targets=dcExpt(iexp).stimResp';
   
   cellIdx = dcExpt(iexp).signifResponsiveCells' & ...
       dcExpt(iexp).oriFitTheta90 < theta90Threshold;
   nCells(iexp) = sum(cellIdx);

   X_targets=X_targets(:,cellIdx);
   
   idx=find(~isnan(sum(X_targets,2)));
   Y_yeses=Y_yeses(idx,1);   
   Y_targets=Y_targets(idx,1);   
   X_targets=X_targets(idx,:);
   
   [c,score,latent] = pca(X_targets);
   
   pve_expt = cumsum(latent)./sum(latent);
   pve{iexp} = pve_expt;
%    dim75 = find(pve_expt > 0.75,1);
    dim75 = length(pve_expt);
   X_prime = score(:,1:dim75);
   X_prime = zscore(X_prime);
   pcNumber = cat(2,pcNumber,1:dim75);
   c = c(1:dim75,1:dim75);
   
   
   exptName{iexp} = [dcExpt(iexp).mouse '-' dcExpt(iexp).date];
   
   detectCorr = corr(Y_yeses,X_targets)';
   targetCorr = corr(Y_targets,X_targets)';
   
   detectCorrsAll_cells_vis = cat(1,detectCorrsAll_cells_vis,detectCorr);
   targetCorrsAll_cells_vis = cat(1,targetCorrsAll_cells_vis,targetCorr);
   
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
   
   targetWeightsAll_cells_vis = cat(1,targetWeightsAll_cells_vis,B_targets);
   detectWeightsAll_cells_vis = cat(1,detectWeightsAll_cells_vis,B_yeses);
  
   
   detectCorr = corr(Y_yeses,X_prime)';
   targetCorr = corr(Y_targets,X_prime)';
   
   detectCorrsAll_pca_vis = cat(1,detectCorrsAll_pca_vis,detectCorr);
   targetCorrsAll_pca_vis = cat(1,targetCorrsAll_pca_vis,targetCorr);
   
   C=eye(size(X_prime,2));%C=inv(sqrtm(cov(X)));
   X_prime=X_prime*C;
   p=1;
   [temp1,dev1,targetGLM]=glmfit(X_prime,Y_targets,'binomial');
   [temp2,dev2,detectGLM]=glmfit(X_prime,Y_yeses,'binomial');
   idx=find(targetGLM.p>p|detectGLM.p>p);
   temp1(idx)=0;
   temp2(idx)=0;
   B_targets = C*temp1(2:end,1);
   B_yeses = C*temp2(2:end,1);
   
   targetWeightsAll_pca_vis = cat(1,targetWeightsAll_pca_vis,B_targets);
   detectWeightsAll_pca_vis = cat(1,detectWeightsAll_pca_vis,B_yeses);
end

%%
figure
for iexp = 1:nexp
    h = plot(pve{iexp});
    hold on    
end
hline(0.95,'k--')
figXAxis([],'N Components',[])
figYAxis([],'% Var Expl.',[])
figAxForm
print(fullfile(fnout,'pca_pctVarExplained'),'-dpdf')

figure
suptitle('Visual Trials, weights and corrrlations are abs val')
x = pcNumber;
subplot 221
y = abs(targetCorrsAll_pca_vis);
scatter(x,y,'ko');
figXAxis([],'PC Number',[]);
figYAxis([],'Correlation',[])
figAxForm
title('Target')

subplot 222
y = abs(detectCorrsAll_pca_vis);
scatter(x,y,'ko');
figXAxis([],'PC Number',[]);
figYAxis([],'Correlation',[])
figAxForm
title('Detect')

subplot 223
y = abs(targetWeightsAll_pca_vis);
scatter(x,y,'ko');
figXAxis([],'PC Number',[]);
figYAxis([],'LR Weight',[])
figAxForm
title('Target')

subplot 224
y = abs(detectWeightsAll_pca_vis);
scatter(x,y,'ko');
figXAxis([],'PC Number',[]);
figYAxis([],'LR Weight',[])
figAxForm
title('Detect')

print(fullfile(fnout,'weightOrCorrPerPC'),'-dpdf')
%% 

setFigParams4Print('portrait')
set(0,'defaultAxesFontSize',12)
weightLim = [-2.2 4.4];
corrLim = [-1 1];
figure
suptitle('Visual Trials')
subplot 321
scatter(detectWeightsAll_cells_vis,detectCorrsAll_cells_vis,'ko')
hold on
hline(0,'k:')
vline(0,'k:')
figXAxis([],'Detect LR Weight',weightLim)
figYAxis([],'Detect Correlation',corrLim)
figAxForm
title('Cells as Predictors')
subplot 323
scatter(targetWeightsAll_cells_vis,targetCorrsAll_cells_vis,'ko')
hold on
hline(0,'k:')
vline(0,'k:')
figXAxis([],'Target LR Weight',weightLim)
figYAxis([],'Target Correlation',corrLim)
figAxForm
subplot 325   
s = scatter(targetWeightsAll_cells_vis,detectWeightsAll_cells_vis,'ko');
s.MarkerFaceColor = [1 1 1];
hold on
hline(0,'k:')
vline(0,'k:')
figXAxis([],'Target LR Weight',weightLim)
figYAxis([],'Detect LR Weight',weightLim)
figAxForm

ind1 = pcNumber == 1;
ind2 = pcNumber == 2;
subplot 322
scatter(detectWeightsAll_pca_vis,detectCorrsAll_pca_vis,'ko')
hold on
hline(0,'k:')
vline(0,'k:')
scatter(detectWeightsAll_pca_vis(ind1),detectCorrsAll_pca_vis(ind1),'ro')
scatter(detectWeightsAll_pca_vis(ind2),detectCorrsAll_pca_vis(ind2),'bo')
legend({'';'PC1';'PC2'},'location','southeast')
figXAxis([],'Detect LR Weight',weightLim)
figYAxis([],'Detect Correlation',corrLim)
figAxForm
title('PCs as Predictors')
subplot 324
scatter(targetWeightsAll_pca_vis,targetCorrsAll_pca_vis,'ko')
hold on
hline(0,'k:')
vline(0,'k:')
scatter(targetWeightsAll_pca_vis(ind1),targetCorrsAll_pca_vis(ind1),'ro')
scatter(targetWeightsAll_pca_vis(ind2),targetCorrsAll_pca_vis(ind2),'bo')
figXAxis([],'Target LR Weight',weightLim)
figYAxis([],'Target Correlation',corrLim)
figAxForm
subplot 326   
s = scatter(targetWeightsAll_pca_vis,detectWeightsAll_pca_vis,'ko');
s.MarkerFaceColor = [1 1 1];
hold on
s = scatter(targetWeightsAll_pca_vis(ind1),detectWeightsAll_pca_vis(ind1),'ro');
s = scatter(targetWeightsAll_pca_vis(ind2),detectWeightsAll_pca_vis(ind2),'bo');
hline(0,'k:')
vline(0,'k:')
figXAxis([],'Target LR Weight',weightLim)
figYAxis([],'Detect LR Weight',weightLim)
figAxForm

print(fullfile(fnout,'corrAndWeight_PCvsCells'),'-dpdf','-fillpage')

%% cross-validation - cells
pctCorrOri_detect_all = cell(1,nexp);
pctCorrOri_detect_holdout = cell(1,nexp);
pctCorrOri_target_all = cell(1,nexp);
pctCorrOri_target_holdout = cell(1,nexp);

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
   nTrials_all = length(Y_yeses);
   
    [~,~,targetGLM] = glmfit(X_targets,Y_targets,'binomial');
    yhat_target = glmval(targetGLM.beta,X_targets,'logit') > decisionVariable;
    pctCorr_target = mean((Y_targets-yhat_target) == 0);
    [~,pctCorrCI_target] = binofit(sum((Y_targets-yhat_target) == 0), ...
        length(Y_targets));
    pctCorrErr_target = nan(1,2);
    pctCorrErr_target(1) = pctCorr_target - pctCorrCI_target(1);
    pctCorrErr_target(2) = pctCorrCI_target(2) - pctCorr_target;
    
    [~,~,detectGLM] = glmfit(X_targets,Y_yeses,'binomial');
    yhat_detect = glmval(detectGLM.beta,X_targets,'logit') > decisionVariable;
    pctCorr_detect = mean((Y_yeses-yhat_detect) == 0);
    [~,pctCorrCI_detect] = binofit(sum((Y_yeses-yhat_detect) == 0), ...
        length(Y_yeses));
    pctCorrErr_detect = nan(1,2);
    pctCorrErr_detect(1) = pctCorr_detect - pctCorrCI_detect(1);
    pctCorrErr_detect(2) = pctCorrCI_detect(2) - pctCorr_detect;

    ysub_targetHoldouts = nan(1,nTrials_all);
    ysub_detectHoldouts = nan(1,nTrials_all);
    for i = 1:nTrials_all
        X_holdout = X_targets(i,:);
        othersInd = [1:(i-1),(i+1):nTrials_all];
        X_others = X_targets(othersInd,:);
        
        [~,~,targetOthersGLM] = glmfit(X_others,Y_targets(othersInd),'binomial');        
        Y_holdout = Y_targets(i);
        yhat = glmval(targetOthersGLM.beta,X_holdout,'logit');
        ysub_targetHoldouts(i) = Y_holdout - (yhat > decisionVariable);
        
        [~,~,detectOthersGLM] = glmfit(X_others,Y_yeses(othersInd),'binomial');
        Y_holdout = Y_yeses(i);
        yhat = glmval(detectOthersGLM.beta,X_holdout,'logit');
        ysub_detectHoldouts(i) = Y_holdout - (yhat > decisionVariable);       
    end
    pctCorr_targetXVal = mean(ysub_targetHoldouts == 0);
    pctCorr_detectXVal = mean(ysub_detectHoldouts == 0);
    
    tOri = dcExpt(iexp).trialOrientation;
    tOri = binAndRelabelTrialOrientation(tOri,oriBinsFine);
    orientations = unique(tOri);
    nOri = length(orientations);
    
    pctCorrXValOri_detect = nan(1,nOri);
    pctCorrXValOri_target = nan(1,nOri);    
    nOriHits_detect = nan(1,nOri);
    nOriHits_target = nan(1,nOri);
    nOriTrials = nan(1,nOri);
    pctCorrOri_detect = nan(1,nOri);
    pctCorrOri_target = nan(1,nOri);
    pctTargetOriErr = nan(nOri,2);
    pctDetectOriErr = nan(nOri,2);
    for iori = 1:nOri
        trInd = find(tOri == orientations(iori));
        X = X_targets(trInd,:);
        
        yhat = glmval(detectGLM.beta,X,'logit') > decisionVariable;
        pctCorrOri_detect(iori) = mean((Y_yeses(trInd) - yhat) == 0);
        [~,pctDetectOriCI] = binofit(sum((Y_yeses(trInd) - yhat) == 0),...
            length(Y_yeses(trInd)));
        pctDetectOriErr(iori,1) = pctCorrOri_detect(iori) - pctDetectOriCI(1);
        pctDetectOriErr(iori,2) = pctDetectOriCI(2) - pctDetectOriErr(iori);
        
        yhat = glmval(targetGLM.beta,X,'logit') > decisionVariable;
        pctCorrOri_target(iori) = mean((Y_targets(trInd) - yhat) == 0);        
        [~,pctTargetOriCI] = binofit(sum((Y_targets(trInd) - yhat) == 0),...
            length(Y_targets(trInd)));
        pctTargetOriErr(iori,1) = pctCorrOri_target(iori) - pctTargetOriCI(1);
        pctTargetOriErr(iori,2) = pctTargetOriCI(2) - pctTargetOriErr(iori);
        
        nTrials = length(trInd);
        nOriTrials(iori) = nTrials;
        ysub_detectHoldouts = nan(1,nTrials);
        ysub_targetHoldouts = nan(1,nTrials);
        for i = 1:nTrials
            X_holdout = X(i,:);
            oriTrInd = trInd(i);
            othersInd = [1:(oriTrInd-1),(oriTrInd+1):nTrials_all];
            X_others = X_targets(othersInd,:);
            
            [~,~,detectOthersGLM] = glmfit(...
                X_others,Y_yeses(othersInd),'binomial');
            yhat_detectHoldout = glmval(...
                detectOthersGLM.beta,X_holdout,'logit') > decisionVariable;
            Y_holdout = Y_yeses(oriTrInd);
            ysub_detectHoldouts(i) = Y_holdout - yhat_detectHoldout;
            
            [~,~,targetOthersGLM] = glmfit(...
                X_others,Y_targets(othersInd),'binomial');
            yhat_targetHoldouts = glmval(...
                targetOthersGLM.beta,X_holdout,'logit') > decisionVariable;
            Y_holdout = Y_targets(oriTrInd);
            ysub_targetHoldouts(i) = Y_holdout - yhat_targetHoldouts;
        end
        pctCorrXValOri_detect(iori) = mean(ysub_detectHoldouts == 0);
        pctCorrXValOri_target(iori) = mean(ysub_targetHoldouts == 0);
    end

    pctCorrOri_detect_all{iexp} = [pctCorrOri_detect.*100 pctCorr_detect.*100];
    pctCorrOri_detect_holdout{iexp} = [pctCorrXValOri_detect.*100 pctCorr_detectXVal.*100];

    pctCorrOri_target_all{iexp} = [pctCorrOri_target.*100 pctCorr_target.*100];
    pctCorrOri_target_holdout{iexp} = [pctCorrXValOri_target.*100 pctCorr_targetXVal.*100];

    exptOris{iexp} = [orientations 100];
end

pctCorrOri_detect_sub = cellfun(@(x,y) x-y,pctCorrOri_detect_all,...
    pctCorrOri_detect_holdout,'unif',0);
pctCorrOri_target_sub = cellfun(@(x,y) x-y,pctCorrOri_target_all,...
    pctCorrOri_target_holdout,'unif',0);

set(0,'defaultAxesFontSize',12)
figure
suptitle('Summary Across Experiments')
subplot 231
allOris = [];
for i = 1:nexp
    x = exptOris{i};
    allOris = unique([allOris,x]);
end
allExpt = nan(length(allOris),nexp);
for i = 1:nexp
    x = exptOris{i};
    y = pctCorrOri_detect_all{i};
    h = plot(x,y,'-');
    h.Color = [0.5 0.5 0.5];
    hold on
    for iori = 1:length(x)
        ind = allOris == x(iori);
        allExpt(ind,i) = y(iori);
    end
end
h = errorbar(allOris,nanmean(allExpt,2),ste(allExpt,2),'ko-');
h.MarkerFaceColor = 'k';
xLabel = [cellfun(@(y)...
    num2str(round(y,2,'significant')),num2cell(allOris(1:end-1)),'unif',0),...
    {'All Trials'}];
figXAxis([],'Orientation (deg)',[-10 110],allOris,xLabel)
figYAxis([],'% Correct',[0 110])
figAxForm
title('Detect GLM - All Trials')

subplot 232
allExpt = nan(length(allOris),nexp);
for i = 1:nexp
    x = exptOris{i};
    y = pctCorrOri_detect_holdout{i};
    h = plot(x,y,'-');
    h.Color = [0.5 0.5 0.5];
    hold on
    for iori = 1:length(x)
        ind = allOris == x(iori);
        allExpt(ind,i) = y(iori);
    end
end
h = errorbar(allOris,nanmean(allExpt,2),ste(allExpt,2),'ko-');
h.MarkerFaceColor = 'k';
xLabel = [cellfun(@(y)...
    num2str(round(y,2,'significant')),num2cell(allOris(1:end-1)),'unif',0),...
    {'All Trials'}];
figXAxis([],'Orientation (deg)',[-10 110],allOris,xLabel)
figYAxis([],'% Correct',[0 110])
figAxForm
title('Detect GLM - Holdout Trials')

pctCorrOri_sub_catAllData = cell(1,length(allOris));
for iexp = 1:nexp
    for i = 1:length(allOris)
        ind = exptOris{iexp} == allOris(i);
        pctCorrOri_sub_catAllData{i} = ...
            [pctCorrOri_sub_catAllData{i} pctCorrOri_detect_sub{iexp}(ind)];
    end
end
subplot 233
for i = 1:length(allOris)
    h = plot(allOris(i),pctCorrOri_sub_catAllData{i},'k.');
%     h.Color = [0.5 0.5 0.5]
    hold on
end
x = allOris;
y = cellfun(@mean,pctCorrOri_sub_catAllData);
yerr = cellfun(@(z) ste(z,2),pctCorrOri_sub_catAllData);
h = errorbar(x,y,yerr,'ko');
h.MarkerFaceColor = 'k';
figXAxis([],'Orientation (deg)',[-10 110],allOris,xLabel)
figYAxis([],'% Correct (All - Holdout)',[-100 100])
figAxForm
hline(0,'k--')

subplot 234
allExpt = nan(length(allOris),nexp);
for i = 1:nexp
    x = exptOris{i};
    y = pctCorrOri_target_all{i};
    h = plot(x,y,'-');
    h.Color = [0.5 0.5 0.5];
    hold on
    for iori = 1:length(x)
        ind = allOris == x(iori);
        allExpt(ind,i) = y(iori);
    end
end
h = errorbar(allOris,nanmean(allExpt,2),ste(allExpt,2),'ko-');
h.MarkerFaceColor = 'k';
xLabel = [cellfun(@(y)...
    num2str(round(y,2,'significant')),num2cell(allOris(1:end-1)),'unif',0),...
    {'All Trials'}];
figXAxis([],'Orientation (deg)',[-10 110],allOris,xLabel)
figYAxis([],'% Correct',[0 110])
figAxForm
title('Target GLM - All Trials')

subplot 235
allExpt = nan(length(allOris),nexp);
for i = 1:nexp
    x = exptOris{i};
    y = pctCorrOri_target_holdout{i};
    h = plot(x,y,'-');
    h.Color = [0.5 0.5 0.5];
    hold on
    for iori = 1:length(x)
        ind = allOris == x(iori);
        allExpt(ind,i) = y(iori);
    end
end
h = errorbar(allOris,nanmean(allExpt,2),ste(allExpt,2),'ko-');
h.MarkerFaceColor = 'k';
xLabel = [cellfun(@(y)...
    num2str(round(y,2,'significant')),num2cell(allOris(1:end-1)),'unif',0),...
    {'All Trials'}];
figXAxis([],'Orientation (deg)',[-10 110],allOris,xLabel)
figYAxis([],'% Correct',[0 110])
figAxForm
title('Target GLM - Holdout Trials')

pctCorrOri_sub_catAllData = cell(1,length(allOris));
for iexp = 1:nexp
    for i = 1:length(allOris)
        ind = exptOris{iexp} == allOris(i);
        pctCorrOri_sub_catAllData{i} = ...
            [pctCorrOri_sub_catAllData{i} pctCorrOri_target_sub{iexp}(ind)];
    end
end
subplot 236
for i = 1:length(allOris)
    h = plot(allOris(i),pctCorrOri_sub_catAllData{i},'k.');
    hold on
end
x = allOris;
y = cellfun(@mean,pctCorrOri_sub_catAllData);
yerr = cellfun(@(z) ste(z,2),pctCorrOri_sub_catAllData);
h = errorbar(x,y,yerr,'ko');
h.MarkerFaceColor = 'k';
figXAxis([],'Orientation (deg)',[-10 110],allOris,xLabel)
figYAxis([],'% Correct (All - Holdout)',[-100 100])
figAxForm
hline(0,'k--')

print(fullfile(fnout,[dataLabel 'pctCorr_all_cells']),'-dpdf','-fillpage')

%% cross-validation - pca
pctCorrOri_detect_all = cell(1,nexp);
pctCorrOri_detect_holdout = cell(1,nexp);
pctCorrOri_target_all = cell(1,nexp);
pctCorrOri_target_holdout = cell(1,nexp);

exptOris = cell(1,nexp);
for iexp = 1:nexp

   trOut=dcExpt(iexp).trialOutcome;

   [Y_yeses, Y_targets] = getStimAndBehaviorYs(trOut);
   
   X_targets=dcExpt(iexp).stimResp';
   
   cellIdx = dcExpt(iexp).signifResponsiveCells' & ...
       dcExpt(iexp).oriFitTheta90 < theta90Threshold;

   X_targets=X_targets(:,cellIdx);
   
   [c,score,latent] = pca(X_targets);
   pve_expt = cumsum(latent)./sum(latent);
   dim75 = find(pve_expt > 0.95,1);
   X_targets = score(:,1:dim75);
   
   X_targets=bsxfun(@plus,X_targets,-mean(X_targets));
   X_targets=bsxfun(@times,X_targets,1./std(X_targets));
   C=eye(size(X_targets,2));%C=inv(sqrtm(cov(X)));
   X_targets=X_targets*C;
   p=1;
   nTrials_all = length(Y_yeses);
   
    [~,~,targetGLM] = glmfit(X_targets,Y_targets,'binomial');
    yhat_target = glmval(targetGLM.beta,X_targets,'logit') > decisionVariable;
    pctCorr_target = mean((Y_targets-yhat_target) == 0);
    [~,pctCorrCI_target] = binofit(sum((Y_targets-yhat_target) == 0), ...
        length(Y_targets));
    pctCorrErr_target = nan(1,2);
    pctCorrErr_target(1) = pctCorr_target - pctCorrCI_target(1);
    pctCorrErr_target(2) = pctCorrCI_target(2) - pctCorr_target;
    
    [~,~,detectGLM] = glmfit(X_targets,Y_yeses,'binomial');
    yhat_detect = glmval(detectGLM.beta,X_targets,'logit') > decisionVariable;
    pctCorr_detect = mean((Y_yeses-yhat_detect) == 0);
    [~,pctCorrCI_detect] = binofit(sum((Y_yeses-yhat_detect) == 0), ...
        length(Y_yeses));
    pctCorrErr_detect = nan(1,2);
    pctCorrErr_detect(1) = pctCorr_detect - pctCorrCI_detect(1);
    pctCorrErr_detect(2) = pctCorrCI_detect(2) - pctCorr_detect;

    ysub_targetHoldouts = nan(1,nTrials_all);
    ysub_detectHoldouts = nan(1,nTrials_all);
    for i = 1:nTrials_all
        X_holdout = X_targets(i,:);
        othersInd = [1:(i-1),(i+1):nTrials_all];
        X_others = X_targets(othersInd,:);
        
        [~,~,targetOthersGLM] = glmfit(X_others,Y_targets(othersInd),'binomial');        
        Y_holdout = Y_targets(i);
        yhat = glmval(targetOthersGLM.beta,X_holdout,'logit');
        ysub_targetHoldouts(i) = Y_holdout - (yhat > decisionVariable);
        
        [~,~,detectOthersGLM] = glmfit(X_others,Y_yeses(othersInd),'binomial');
        Y_holdout = Y_yeses(i);
        yhat = glmval(detectOthersGLM.beta,X_holdout,'logit');
        ysub_detectHoldouts(i) = Y_holdout - (yhat > decisionVariable);       
    end
    pctCorr_targetXVal = mean(ysub_targetHoldouts == 0);
    pctCorr_detectXVal = mean(ysub_detectHoldouts == 0);
    
    tOri = dcExpt(iexp).trialOrientation;
    tOri = binAndRelabelTrialOrientation(tOri,oriBinsFine);
    orientations = unique(tOri);
    nOri = length(orientations);
    
    pctCorrXValOri_detect = nan(1,nOri);
    pctCorrXValOri_target = nan(1,nOri);    
    nOriHits_detect = nan(1,nOri);
    nOriHits_target = nan(1,nOri);
    nOriTrials = nan(1,nOri);
    pctCorrOri_detect = nan(1,nOri);
    pctCorrOri_target = nan(1,nOri);
    pctTargetOriErr = nan(nOri,2);
    pctDetectOriErr = nan(nOri,2);
    for iori = 1:nOri
        trInd = find(tOri == orientations(iori));
        X = X_targets(trInd,:);
        
        yhat = glmval(detectGLM.beta,X,'logit') > decisionVariable;
        pctCorrOri_detect(iori) = mean((Y_yeses(trInd) - yhat) == 0);
        [~,pctDetectOriCI] = binofit(sum((Y_yeses(trInd) - yhat) == 0),...
            length(Y_yeses(trInd)));
        pctDetectOriErr(iori,1) = pctCorrOri_detect(iori) - pctDetectOriCI(1);
        pctDetectOriErr(iori,2) = pctDetectOriCI(2) - pctDetectOriErr(iori);
        
        yhat = glmval(targetGLM.beta,X,'logit') > decisionVariable;
        pctCorrOri_target(iori) = mean((Y_targets(trInd) - yhat) == 0);        
        [~,pctTargetOriCI] = binofit(sum((Y_targets(trInd) - yhat) == 0),...
            length(Y_targets(trInd)));
        pctTargetOriErr(iori,1) = pctCorrOri_target(iori) - pctTargetOriCI(1);
        pctTargetOriErr(iori,2) = pctTargetOriCI(2) - pctTargetOriErr(iori);
        
        nTrials = length(trInd);
        nOriTrials(iori) = nTrials;
        ysub_detectHoldouts = nan(1,nTrials);
        ysub_targetHoldouts = nan(1,nTrials);
        for i = 1:nTrials
            X_holdout = X(i,:);
            oriTrInd = trInd(i);
            othersInd = [1:(oriTrInd-1),(oriTrInd+1):nTrials_all];
            X_others = X_targets(othersInd,:);
            
            [~,~,detectOthersGLM] = glmfit(...
                X_others,Y_yeses(othersInd),'binomial');
            yhat_detectHoldout = glmval(...
                detectOthersGLM.beta,X_holdout,'logit') > decisionVariable;
            Y_holdout = Y_yeses(oriTrInd);
            ysub_detectHoldouts(i) = Y_holdout - yhat_detectHoldout;
            
            [~,~,targetOthersGLM] = glmfit(...
                X_others,Y_targets(othersInd),'binomial');
            yhat_targetHoldouts = glmval(...
                targetOthersGLM.beta,X_holdout,'logit') > decisionVariable;
            Y_holdout = Y_targets(oriTrInd);
            ysub_targetHoldouts(i) = Y_holdout - yhat_targetHoldouts;
        end
        pctCorrXValOri_detect(iori) = mean(ysub_detectHoldouts == 0);
        pctCorrXValOri_target(iori) = mean(ysub_targetHoldouts == 0);
    end

    pctCorrOri_detect_all{iexp} = [pctCorrOri_detect.*100 pctCorr_detect.*100];
    pctCorrOri_detect_holdout{iexp} = [pctCorrXValOri_detect.*100 pctCorr_detectXVal.*100];

    pctCorrOri_target_all{iexp} = [pctCorrOri_target.*100 pctCorr_target.*100];
    pctCorrOri_target_holdout{iexp} = [pctCorrXValOri_target.*100 pctCorr_targetXVal.*100];

    exptOris{iexp} = [orientations 100];
end

pctCorrOri_detect_sub = cellfun(@(x,y) x-y,pctCorrOri_detect_all,...
    pctCorrOri_detect_holdout,'unif',0);
pctCorrOri_target_sub = cellfun(@(x,y) x-y,pctCorrOri_target_all,...
    pctCorrOri_target_holdout,'unif',0);

setFigParams4Print('landscape')
set(0,'defaultAxesFontSize',12)
figure
suptitle('Predictions - PCs')
subplot 231
allOris = [];
for i = 1:nexp
    x = exptOris{i};
    allOris = unique([allOris,x]);
end
allExpt = nan(length(allOris),nexp);
for i = 1:nexp
    x = exptOris{i};
    y = pctCorrOri_detect_all{i};
    h = plot(x,y,'-');
    h.Color = [0.5 0.5 0.5];
    hold on
    for iori = 1:length(x)
        ind = allOris == x(iori);
        allExpt(ind,i) = y(iori);
    end
end
h = errorbar(allOris,nanmean(allExpt,2),ste(allExpt,2),'ko-');
h.MarkerFaceColor = 'k';
xLabel = [cellfun(@(y)...
    num2str(round(y,2,'significant')),num2cell(allOris(1:end-1)),'unif',0),...
    {'All Trials'}];
figXAxis([],'Orientation (deg)',[-10 110],allOris,xLabel)
figYAxis([],'% Correct',[0 110])
figAxForm
title('Detect GLM - All Trials')

subplot 232
allExpt = nan(length(allOris),nexp);
for i = 1:nexp
    x = exptOris{i};
    y = pctCorrOri_detect_holdout{i};
    h = plot(x,y,'-');
    h.Color = [0.5 0.5 0.5];
    hold on
    for iori = 1:length(x)
        ind = allOris == x(iori);
        allExpt(ind,i) = y(iori);
    end
end
h = errorbar(allOris,nanmean(allExpt,2),ste(allExpt,2),'ko-');
h.MarkerFaceColor = 'k';
xLabel = [cellfun(@(y)...
    num2str(round(y,2,'significant')),num2cell(allOris(1:end-1)),'unif',0),...
    {'All Trials'}];
figXAxis([],'Orientation (deg)',[-10 110],allOris,xLabel)
figYAxis([],'% Correct',[0 110])
figAxForm
title('Detect GLM - Holdout Trials')

pctCorrOri_sub_catAllData = cell(1,length(allOris));
for iexp = 1:nexp
    for i = 1:length(allOris)
        ind = exptOris{iexp} == allOris(i);
        pctCorrOri_sub_catAllData{i} = ...
            [pctCorrOri_sub_catAllData{i} pctCorrOri_detect_sub{iexp}(ind)];
    end
end
subplot 233
for i = 1:length(allOris)
    h = plot(allOris(i),pctCorrOri_sub_catAllData{i},'k.');
%     h.Color = [0.5 0.5 0.5]
    hold on
end
x = allOris;
y = cellfun(@mean,pctCorrOri_sub_catAllData);
yerr = cellfun(@(z) ste(z,2),pctCorrOri_sub_catAllData);
h = errorbar(x,y,yerr,'ko');
h.MarkerFaceColor = 'k';
figXAxis([],'Orientation (deg)',[-10 110],allOris,xLabel)
figYAxis([],'% Correct (All - Holdout)',[-100 100])
figAxForm
hline(0,'k--')

subplot 234
allExpt = nan(length(allOris),nexp);
for i = 1:nexp
    x = exptOris{i};
    y = pctCorrOri_target_all{i};
    h = plot(x,y,'-');
    h.Color = [0.5 0.5 0.5];
    hold on
    for iori = 1:length(x)
        ind = allOris == x(iori);
        allExpt(ind,i) = y(iori);
    end
end
h = errorbar(allOris,nanmean(allExpt,2),ste(allExpt,2),'ko-');
h.MarkerFaceColor = 'k';
xLabel = [cellfun(@(y)...
    num2str(round(y,2,'significant')),num2cell(allOris(1:end-1)),'unif',0),...
    {'All Trials'}];
figXAxis([],'Orientation (deg)',[-10 110],allOris,xLabel)
figYAxis([],'% Correct',[0 110])
figAxForm
title('Target GLM - All Trials')

subplot 235
allExpt = nan(length(allOris),nexp);
for i = 1:nexp
    x = exptOris{i};
    y = pctCorrOri_target_holdout{i};
    h = plot(x,y,'-');
    h.Color = [0.5 0.5 0.5];
    hold on
    for iori = 1:length(x)
        ind = allOris == x(iori);
        allExpt(ind,i) = y(iori);
    end
end
h = errorbar(allOris,nanmean(allExpt,2),ste(allExpt,2),'ko-');
h.MarkerFaceColor = 'k';
xLabel = [cellfun(@(y)...
    num2str(round(y,2,'significant')),num2cell(allOris(1:end-1)),'unif',0),...
    {'All Trials'}];
figXAxis([],'Orientation (deg)',[-10 110],allOris,xLabel)
figYAxis([],'% Correct',[0 110])
figAxForm
title('Target GLM - Holdout Trials')

pctCorrOri_sub_catAllData = cell(1,length(allOris));
for iexp = 1:nexp
    for i = 1:length(allOris)
        ind = exptOris{iexp} == allOris(i);
        pctCorrOri_sub_catAllData{i} = ...
            [pctCorrOri_sub_catAllData{i} pctCorrOri_target_sub{iexp}(ind)];
    end
end
subplot 236
for i = 1:length(allOris)
    h = plot(allOris(i),pctCorrOri_sub_catAllData{i},'k.');
    hold on
end
x = allOris;
y = cellfun(@mean,pctCorrOri_sub_catAllData);
yerr = cellfun(@(z) ste(z,2),pctCorrOri_sub_catAllData);
h = errorbar(x,y,yerr,'ko');
h.MarkerFaceColor = 'k';
figXAxis([],'Orientation (deg)',[-10 110],allOris,xLabel)
figYAxis([],'% Correct (All - Holdout)',[-100 100])
figAxForm
hline(0,'k--')

print(fullfile(fnout,[dataLabel 'pctCorr_all_pc']),'-dpdf','-fillpage')
