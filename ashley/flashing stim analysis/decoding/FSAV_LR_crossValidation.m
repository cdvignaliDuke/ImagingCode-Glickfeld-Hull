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
%% cross-validation
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
    
    if doExptPlots
        setFigParams4Print('landscape')
        figure
        subplot 121
        x = [orientations 100];
        y = [pctCorrOri_detect.*100 pctCorr_detect.*100];
        ylerr = cat(1,pctDetectOriErr(:,1),pctCorrErr_detect(1))'.*100;
        yuerr = cat(1,pctDetectOriErr(:,2),pctCorrErr_detect(2))'.*100;
        h = errorbar(x,y,ylerr,yuerr,'ko-');
        h.MarkerFaceColor = [0 0 0];
        hold on
        y = [pctCorrXValOri_detect.*100 pctCorr_detectXVal.*100];
        h = plot(x,y,'ko-');
        h.MarkerFaceColor = [1 1 1];
        xLabel = [cellfun(@(y)...
            num2str(round(y,2,'significant')),num2cell(x(1:end-1)),'unif',0),...
            {'All Trials'}];
        figXAxis([],'Orientation (deg)',[-10 110],x,xLabel)
        figYAxis([],'% Correct',[0 110])
        figAxForm
        legend({'All Trials';'Holdout X Val'},'location','northeastoutside')
        title({'Detect GLM';[dcExpt(iexp).mouse '-' dcExpt(iexp).date]})

        subplot 122
        x = [orientations 100];
        y = [pctCorrOri_target.*100 pctCorr_target.*100];
        ylerr = cat(1,pctTargetOriErr(:,1),pctCorrErr_target(1))'.*100;
        yuerr = cat(1,pctTargetOriErr(:,2),pctCorrErr_target(2))'.*100;
        h = errorbar(x,y,ylerr,yuerr,'ko-');
        h.MarkerFaceColor = [0 0 0];
        hold on
        y = [pctCorrXValOri_target.*100 pctCorr_targetXVal.*100];
        h = plot(x,y,'ko-');
        h.MarkerFaceColor = [1 1 1];
        xLabel = [cellfun(@(y)...
            num2str(round(y,2,'significant')),num2cell(x(1:end-1)),'unif',0),...
            {'All Trials'}];
        figXAxis([],'Orientation (deg)',[-10 110],x,xLabel)
        figYAxis([],'% Correct',[0 110])
        figAxForm
        legend({'All Trials';'Holdout X Val'},'location','northeastoutside')
        title({'Target GLM';[dcExpt(iexp).mouse '-' dcExpt(iexp).date]})

        print(fullfile(fnout,[dcExpt(iexp).mouse '-' dcExpt(iexp).date '_xValidation']),'-dpdf','-fillpage')
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

print(fullfile(fnout,[dataLabel 'pctCorr_all']),'-dpdf','-fillpage')

%% hard trials only
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
   
   
   hardTrialsInd = dcExpt(iexp).trialOrientation < hardTrialCutoff;

   X_targets=X_targets(:,cellIdx);
   
   X_targets=bsxfun(@plus,X_targets,-mean(X_targets));
   X_targets=bsxfun(@times,X_targets,1./std(X_targets));
   C=eye(size(X_targets,2));%C=inv(sqrtm(cov(X)));
   X_targets=X_targets*C;
   p=1;
   nTrials_all = length(Y_yeses);
   
   X_targets_hard = X_targets(hardTrialsInd,:);
   Y_targets_hard = Y_targets(hardTrialsInd);
   Y_yeses_hard = Y_yeses(hardTrialsInd);
   
   nTrials_hard = length(Y_targets_hard);
   
    [~,~,targetGLM_hard] = glmfit(X_targets_hard,Y_targets_hard,'binomial');
    yhat_target = glmval(targetGLM_hard.beta,X_targets,'logit') > decisionVariable;
    pctCorr_target = mean((Y_targets-yhat_target) == 0);
    [~,pctCorrCI_target] = binofit(sum((Y_targets-yhat_target) == 0), ...
        length(Y_targets));
    pctCorrErr_target = nan(1,2);
    pctCorrErr_target(1) = pctCorr_target - pctCorrCI_target(1);
    pctCorrErr_target(2) = pctCorrCI_target(2) - pctCorr_target;
    
    [~,~,detectGLM_hard] = glmfit(X_targets_hard,Y_yeses_hard,'binomial');
    yhat_detect = glmval(detectGLM_hard.beta,X_targets,'logit') > decisionVariable;
    pctCorr_detect = mean((Y_yeses-yhat_detect) == 0);
    [~,pctCorrCI_detect] = binofit(sum((Y_yeses-yhat_detect) == 0), ...
        length(Y_yeses));
    pctCorrErr_detect = nan(1,2);
    pctCorrErr_detect(1) = pctCorr_detect - pctCorrCI_detect(1);
    pctCorrErr_detect(2) = pctCorrCI_detect(2) - pctCorr_detect;

    ysub_targetHoldouts = nan(1,nTrials_all);
    ysub_detectHoldouts = nan(1,nTrials_all);
    for i = 1:nTrials_hard
        X_holdout = X_targets_hard(i,:);
        othersInd = [1:(i-1),(i+1):nTrials_hard];
        X_others = X_targets_hard(othersInd,:);
        
        [~,~,targetOthersGLM] = glmfit(X_others,Y_targets_hard(othersInd),'binomial');        
        Y_holdout = Y_targets_hard(i);
        yhat = glmval(targetOthersGLM.beta,X_holdout,'logit');
        ysub_targetHoldouts(i) = Y_holdout - (yhat > decisionVariable);
        
        [~,~,detectOthersGLM] = glmfit(X_others,Y_yeses_hard(othersInd),'binomial');
        Y_holdout = Y_yeses_hard(i);
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
    nTrialsOri{iexp} = nan(nOri,1);
    for iori = 1:nOri
        trInd = find(tOri == orientations(iori));
        nTrialsOri{iexp}(iori) = length(trInd);
        X = X_targets(trInd,:);
        
        yhat = glmval(detectGLM_hard.beta,X,'logit') > decisionVariable;
        pctCorrOri_detect(iori) = mean((Y_yeses(trInd) - yhat) == 0);
        [~,pctDetectOriCI] = binofit(sum((Y_yeses(trInd) - yhat) == 0),...
            length(Y_yeses(trInd)));
        pctDetectOriErr(iori,1) = pctCorrOri_detect(iori) - pctDetectOriCI(1);
        pctDetectOriErr(iori,2) = pctDetectOriCI(2) - pctDetectOriErr(iori);
        
        yhat = glmval(targetGLM_hard.beta,X,'logit') > decisionVariable;
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
            if ismember(oriTrInd,find(hardTrialsInd))
                oriTrInd_hard = find(find(hardTrialsInd) == oriTrInd);
                othersInd = [1:(oriTrInd_hard-1),(oriTrInd_hard+1):nTrials_hard];
            else
                othersInd = 1:nTrials_hard;
            end
            X_others = X_targets_hard(othersInd,:);
            
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
    
    if doExptPlots
        setFigParams4Print('landscape')
        figure
        subplot 121
        x = [orientations 100];
        y = [pctCorrOri_detect.*100 pctCorr_detect.*100];
        ylerr = cat(1,pctDetectOriErr(:,1),pctCorrErr_detect(1))'.*100;
        yuerr = cat(1,pctDetectOriErr(:,2),pctCorrErr_detect(2))'.*100;
        h = errorbar(x,y,ylerr,yuerr,'ko-');
        h.MarkerFaceColor = [0 0 0];
        hold on
        y = [pctCorrXValOri_detect.*100 pctCorr_detectXVal.*100];
        h = plot(x,y,'ko-');
        h.MarkerFaceColor = [1 1 1];
        xLabel = [cellfun(@(y)...
            num2str(round(y,2,'significant')),num2cell(x(1:end-1)),'unif',0),...
            {'All Trials'}];
        figXAxis([],'Orientation (deg)',[-10 110],x,xLabel)
        figYAxis([],'% Correct',[0 110])
        figAxForm
        legend({'All Trials';'Holdout X Val'},'location','northeastoutside')
        title({'Detect GLM';[dcExpt(iexp).mouse '-' dcExpt(iexp).date]})

        subplot 122
        x = [orientations 100];
        y = [pctCorrOri_target.*100 pctCorr_target.*100];
        ylerr = cat(1,pctTargetOriErr(:,1),pctCorrErr_target(1))'.*100;
        yuerr = cat(1,pctTargetOriErr(:,2),pctCorrErr_target(2))'.*100;
        h = errorbar(x,y,ylerr,yuerr,'ko-');
        h.MarkerFaceColor = [0 0 0];
        hold on
        y = [pctCorrXValOri_target.*100 pctCorr_targetXVal.*100];
        h = plot(x,y,'ko-');
        h.MarkerFaceColor = [1 1 1];
        xLabel = [cellfun(@(y)...
            num2str(round(y,2,'significant')),num2cell(x(1:end-1)),'unif',0),...
            {'All Trials'}];
        figXAxis([],'Orientation (deg)',[-10 110],x,xLabel)
        figYAxis([],'% Correct',[0 110])
        figAxForm
        legend({'All Trials';'Holdout X Val'},'location','northeastoutside')
        title({'Target GLM';[dcExpt(iexp).mouse '-' dcExpt(iexp).date]})

        print(fullfile(fnout,[dcExpt(iexp).mouse '-' dcExpt(iexp).date '_xValidation_hardTrials']),'-dpdf','-fillpage')
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

set(0,'defaultAxesFontSize',16)
figure
suptitle('Summary Across Experiments')
subplot 221
allOris = [];
for i = 1:nexp
    plot(exptOris{i},pctCorrOri_detect_all{i},'o-')
    hold on
    allOris = unique([allOris,exptOris{i}]);
end
xLabel = [cellfun(@(y)...
    num2str(round(y,2,'significant')),num2cell(allOris(1:end-1)),'unif',0),...
    {'All Trials'}];
figXAxis([],'Orientation (deg)',[-10 110],allOris,xLabel)
figYAxis([],'% Correct',[0 110])
figAxForm
title('Detect GLM')

pctCorrOri_sub_catAllData = cell(1,length(allOris));
for iexp = 1:nexp
    for i = 1:length(allOris)
        ind = exptOris{iexp} == allOris(i);
        pctCorrOri_sub_catAllData{i} = ...
            [pctCorrOri_sub_catAllData{i} pctCorrOri_detect_sub{iexp}(ind)];
    end
end
subplot 222
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

subplot 223
allOris = [];
for i = 1:nexp
    plot(exptOris{i},pctCorrOri_target_all{i},'o-')
    hold on
    allOris = unique([allOris,exptOris{i}]);
end
xLabel = [cellfun(@(y)...
    num2str(round(y,2,'significant')),num2cell(allOris(1:end-1)),'unif',0),...
    {'All Trials'}];
figXAxis([],'Orientation (deg)',[-10 110],allOris,xLabel)
figYAxis([],'% Correct',[0 110])
figAxForm
title('Target GLM')

pctCorrOri_sub_catAllData = cell(1,length(allOris));
for iexp = 1:nexp
    for i = 1:length(allOris)
        ind = exptOris{iexp} == allOris(i);
        pctCorrOri_sub_catAllData{i} = ...
            [pctCorrOri_sub_catAllData{i} pctCorrOri_target_sub{iexp}(ind)];
    end
end
subplot 224
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

print(fullfile(fnout,[dataLabel 'pctCorr_hardTrials']),'-dpdf','-fillpage')

%% predictions and cross-validation, binned orientation

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
    tOri = binAndRelabelTrialOrientation(tOri,oriBins);
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
    nTrialsOri{iexp} = nan(nOri,1);
    for iori = 1:nOri
        trInd = find(tOri == orientations(iori));
        nTrialsOri{iexp}(iori) = length(trInd);
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
    
    if doExptPlots
        setFigParams4Print('landscape')
        figure
        subplot 121
        x = [orientations 100];
        y = [pctCorrOri_detect.*100 pctCorr_detect.*100];
        ylerr = cat(1,pctDetectOriErr(:,1),pctCorrErr_detect(1))'.*100;
        yuerr = cat(1,pctDetectOriErr(:,2),pctCorrErr_detect(2))'.*100;
        h = errorbar(x,y,ylerr,yuerr,'ko-');
        h.MarkerFaceColor = [0 0 0];
        hold on
        y = [pctCorrXValOri_detect.*100 pctCorr_detectXVal.*100];
        h = plot(x,y,'ko-');
        h.MarkerFaceColor = [1 1 1];
        xLabel = [cellfun(@(y)...
            num2str(round(y,2,'significant')),num2cell(x(1:end-1)),'unif',0),...
            {'All Trials'}];
        figXAxis([],'Orientation (deg)',[-10 110],x,xLabel)
        figYAxis([],'% Correct',[0 110])
        figAxForm
        legend({'All Trials';'Holdout X Val'},'location','northeastoutside')
        title({'Detect GLM';[dcExpt(iexp).mouse '-' dcExpt(iexp).date]})

        subplot 122
        x = [orientations 100];
        y = [pctCorrOri_target.*100 pctCorr_target.*100];
        ylerr = cat(1,pctTargetOriErr(:,1),pctCorrErr_target(1))'.*100;
        yuerr = cat(1,pctTargetOriErr(:,2),pctCorrErr_target(2))'.*100;
        h = errorbar(x,y,ylerr,yuerr,'ko-');
        h.MarkerFaceColor = [0 0 0];
        hold on
        y = [pctCorrXValOri_target.*100 pctCorr_targetXVal.*100];
        h = plot(x,y,'ko-');
        h.MarkerFaceColor = [1 1 1];
        xLabel = [cellfun(@(y)...
            num2str(round(y,2,'significant')),num2cell(x(1:end-1)),'unif',0),...
            {'All Trials'}];
        figXAxis([],'Orientation (deg)',[-10 110],x,xLabel)
        figYAxis([],'% Correct',[0 110])
        figAxForm
        legend({'All Trials';'Holdout X Val'},'location','northeastoutside')
        title({'Target GLM';[dcExpt(iexp).mouse '-' dcExpt(iexp).date]})


        print(fullfile(fnout,[dcExpt(iexp).mouse '-' dcExpt(iexp).date '_xValidation_binnedOri']),'-dpdf','-fillpage')
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

set(0,'defaultAxesFontSize',16)
figure
suptitle('Summary Across Experiments')
subplot 221
allOris = [];
for i = 1:nexp
    plot(exptOris{i},pctCorrOri_detect_all{i},'o-')
    hold on
    allOris = unique([allOris,exptOris{i}]);
end
xLabel = [cellfun(@(y)...
    num2str(round(y,2,'significant')),num2cell(allOris(1:end-1)),'unif',0),...
    {'All Trials'}];
figXAxis([],'Orientation (deg)',[-10 110],allOris,xLabel)
figYAxis([],'% Correct',[0 110])
figAxForm
title('Detect GLM')

pctCorrOri_sub_catAllData = cell(1,length(allOris));
for iexp = 1:nexp
    for i = 1:length(allOris)
        ind = exptOris{iexp} == allOris(i);
        pctCorrOri_sub_catAllData{i} = ...
            [pctCorrOri_sub_catAllData{i} pctCorrOri_detect_sub{iexp}(ind)];
    end
end
subplot 222
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

subplot 223
allOris = [];
for i = 1:nexp
    plot(exptOris{i},pctCorrOri_target_all{i},'o-')
    hold on
    allOris = unique([allOris,exptOris{i}]);
end
xLabel = [cellfun(@(y)...
    num2str(round(y,2,'significant')),num2cell(allOris(1:end-1)),'unif',0),...
    {'All Trials'}];
figXAxis([],'Orientation (deg)',[-10 110],allOris,xLabel)
figYAxis([],'% Correct',[0 110])
figAxForm
title('Target GLM')

pctCorrOri_sub_catAllData = cell(1,length(allOris));
for iexp = 1:nexp
    for i = 1:length(allOris)
        ind = exptOris{iexp} == allOris(i);
        pctCorrOri_sub_catAllData{i} = ...
            [pctCorrOri_sub_catAllData{i} pctCorrOri_target_sub{iexp}(ind)];
    end
end
subplot 224
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

print(fullfile(fnout,[dataLabel 'pctCorr_byOriBinned']),'-dpdf','-fillpage')

%%
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

    
    pctCorr_targetOthers = nan(1,nTrials_all);
    pctCorr_detectOthers = nan(1,nTrials_all);
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
        
        yhat_others = glmval(targetOthersGLM.beta,X_others,'logit') > decisionVariable;
        pctCorr_targetOthers(i) = mean((yhat_others - Y_targets(othersInd)) == 0);
        
        [~,~,detectOthersGLM] = glmfit(X_others,Y_yeses(othersInd),'binomial');
        Y_holdout = Y_yeses(i);
        yhat = glmval(detectOthersGLM.beta,X_holdout,'logit');
        ysub_detectHoldouts(i) = Y_holdout - (yhat > decisionVariable);    
        
        yhat_others = glmval(detectOthersGLM.beta,X_others,'logit') > decisionVariable;
        pctCorr_detectOthers(i) = mean((yhat_others - Y_yeses(othersInd)) == 0);
    end
    pctCorr_targetXVal = mean(ysub_targetHoldouts == 0);
    pctCorr_detectXVal = mean(ysub_detectHoldouts == 0);
    
    tOri = dcExpt(iexp).trialOrientation;
    tOri = binAndRelabelTrialOrientation(tOri,oriBins);
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
    nTrialsOri{iexp} = nan(nOri,1);
    pctCorrOri_detectOthers = cell(1,nOri);
    pctCorrOri_targetOthers = cell(1,nOri);
    for iori = 1:nOri
        trInd = find(tOri == orientations(iori));
        nTrialsOri{iexp}(iori) = length(trInd);
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
        pctCorrOri_detectOthers{iori} = nan(1,nTrials);
        pctCorrOri_targetOthers{iori} = nan(1,nTrials);
        ysub_detectHoldouts = nan(1,nTrials);
        ysub_targetHoldouts = nan(1,nTrials);
        for i = 1:nTrials
            X_holdout = X(i,:);
            oriTrInd = trInd(i);
            othersInd = [1:(oriTrInd-1),(oriTrInd+1):nTrials_all];
            othersIndOri = ~ismember(trInd,oriTrInd);
            X_others = X_targets(othersInd,:);
            
            [~,~,detectOthersGLM] = glmfit(...
                X_others,Y_yeses(othersInd),'binomial');
            yhat_detectHoldout = glmval(...
                detectOthersGLM.beta,X_holdout,'logit') > decisionVariable;
            Y_holdout = Y_yeses(oriTrInd);
            ysub_detectHoldouts(i) = Y_holdout - yhat_detectHoldout;
            
            yhat_others = glmval(detectOthersGLM.beta,X(othersIndOri,:),'logit')...
                > decisionVariable;
            pctCorrOri_detectOthers{iori}(i) = mean(...
                (yhat_others - Y_yeses(trInd(othersIndOri))) == 0);
            
            [~,~,targetOthersGLM] = glmfit(...
                X_others,Y_targets(othersInd),'binomial');
            yhat_targetHoldouts = glmval(...
                targetOthersGLM.beta,X_holdout,'logit') > decisionVariable;
            Y_holdout = Y_targets(oriTrInd);
            ysub_targetHoldouts(i) = Y_holdout - yhat_targetHoldouts;
            
            yhat_others = glmval(targetOthersGLM.beta,X(othersIndOri,:),'logit')...
                > decisionVariable;
            pctCorrOri_targetOthers{iori}(i) = mean(...
                (yhat_others - Y_targets(trInd(othersIndOri))) == 0);
        end
        pctCorrXValOri_detect(iori) = mean(ysub_detectHoldouts == 0);
        pctCorrXValOri_target(iori) = mean(ysub_targetHoldouts == 0);
    end
    
    if doExptPlots
        setFigParams4Print('landscape')
        figure
        subplot 121
        x = [orientations 100];
        y = [cellfun(@mean,pctCorrOri_targetOthers).*100 ...
            mean(pctCorr_targetOthers).*100];
        yerr = [cellfun(@(a) ste(a,2),pctCorrOri_targetOthers).*100 ...
            ste(pctCorr_targetOthers,2).*100];
        h = errorbar(x,y,yerr,'ko-');
        h.MarkerFaceColor = [0 0 0];
        hold on
        y = [pctCorrXValOri_target.*100 pctCorr_targetXVal.*100];
        h = plot(x,y,'ko-');
        h.MarkerFaceColor = [1 1 1];
        xLabel = [cellfun(@(y)...
            num2str(round(y,2,'significant')),num2cell(x(1:end-1)),'unif',0),...
            {'All Trials'}];
        figXAxis([],'Orientation (deg)',[-10 110],x,xLabel)
        figYAxis([],'% Correct',[0 110])
        figAxForm
        legend({'Other Trials';'Holdout X Val'},'location','northeastoutside')
        title({'Target GLM';[dcExpt(iexp).mouse '-' dcExpt(iexp).date]})

        subplot 122
        x = [orientations 100];
        y = [cellfun(@mean,pctCorrOri_detectOthers).*100 ...
            mean(pctCorr_detectOthers).*100];
        yerr = [cellfun(@(a) ste(a,2),pctCorrOri_detectOthers).*100 ...
            ste(pctCorr_detectOthers,2).*100];
        h = errorbar(x,y,yerr,'ko-');
        h.MarkerFaceColor = [0 0 0];
        hold on
        y = [pctCorrXValOri_detect.*100 pctCorr_detectXVal.*100];
        h = plot(x,y,'ko-');
        h.MarkerFaceColor = [1 1 1];
        xLabel = [cellfun(@(y)...
            num2str(round(y,2,'significant')),num2cell(x(1:end-1)),'unif',0),...
            {'All Trials'}];
        figXAxis([],'Orientation (deg)',[-10 110],x,xLabel)
        figYAxis([],'% Correct',[0 110])
        figAxForm
        legend({'Other Trials';'Holdout X Val'},'location','northeastoutside')
        title({'Detect GLM';[dcExpt(iexp).mouse '-' dcExpt(iexp).date]})


        print(fullfile(fnout,[dcExpt(iexp).mouse '-' dcExpt(iexp).date '_xValidation_binnedOri_matched']),'-dpdf','-fillpage')
    end
        pctCorrOri_detect_all{iexp} = [cellfun(@mean,pctCorrOri_detectOthers).*100 ...
            mean(pctCorr_detectOthers).*100];
        pctCorrOri_detect_holdout{iexp} = [pctCorrXValOri_detect.*100 ...
            pctCorr_detectXVal.*100];

        pctCorrOri_target_all{iexp} = [cellfun(@mean,pctCorrOri_targetOthers).*100 ...
            mean(pctCorr_targetOthers).*100];
        pctCorrOri_target_holdout{iexp} = [pctCorrXValOri_target.*100 ...
            pctCorr_targetXVal.*100];

        exptOris{iexp} = [orientations 100];
end


pctCorrOri_detect_sub = cellfun(@(x,y) x-y,pctCorrOri_detect_all,...
    pctCorrOri_detect_holdout,'unif',0);
pctCorrOri_target_sub = cellfun(@(x,y) x-y,pctCorrOri_target_all,...
    pctCorrOri_target_holdout,'unif',0);

set(0,'defaultAxesFontSize',16)
figure
suptitle('Summary Across Experiments')
subplot 221
allOris = [];
for i = 1:nexp
    plot(exptOris{i},pctCorrOri_detect_all{i},'o-')
    hold on
    allOris = unique([allOris,exptOris{i}]);
end
xLabel = [cellfun(@(y)...
    num2str(round(y,2,'significant')),num2cell(allOris(1:end-1)),'unif',0),...
    {'All Trials'}];
figXAxis([],'Orientation (deg)',[-10 110],allOris,xLabel)
figYAxis([],'% Correct',[0 110])
figAxForm
title('Detect GLM')

pctCorrOri_sub_catAllData = cell(1,length(allOris));
for iexp = 1:nexp
    for i = 1:length(allOris)
        ind = exptOris{iexp} == allOris(i);
        pctCorrOri_sub_catAllData{i} = ...
            [pctCorrOri_sub_catAllData{i} pctCorrOri_detect_sub{iexp}(ind)];
    end
end
subplot 222
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

subplot 223
allOris = [];
for i = 1:nexp
    plot(exptOris{i},pctCorrOri_target_all{i},'o-')
    hold on
    allOris = unique([allOris,exptOris{i}]);
end
xLabel = [cellfun(@(y)...
    num2str(round(y,2,'significant')),num2cell(allOris(1:end-1)),'unif',0),...
    {'All Trials'}];
figXAxis([],'Orientation (deg)',[-10 110],allOris,xLabel)
figYAxis([],'% Correct',[0 110])
figAxForm
title('Target GLM')

pctCorrOri_sub_catAllData = cell(1,length(allOris));
for iexp = 1:nexp
    for i = 1:length(allOris)
        ind = exptOris{iexp} == allOris(i);
        pctCorrOri_sub_catAllData{i} = ...
            [pctCorrOri_sub_catAllData{i} pctCorrOri_target_sub{iexp}(ind)];
    end
end
subplot 224
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

print(fullfile(fnout,[dataLabel 'pctCorr_byOriBinned_matched']),'-dpdf','-fillpage')

%%

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

    
    pctCorr_targetOthers = nan(1,nTrials_all);
    pctCorr_detectOthers = nan(1,nTrials_all);
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
        
        yhat_others = glmval(targetOthersGLM.beta,X_others,'logit') > decisionVariable;
        pctCorr_targetOthers(i) = mean((yhat_others - Y_targets(othersInd)) == 0);
        
        [~,~,detectOthersGLM] = glmfit(X_others,Y_yeses(othersInd),'binomial');
        Y_holdout = Y_yeses(i);
        yhat = glmval(detectOthersGLM.beta,X_holdout,'logit');
        ysub_detectHoldouts(i) = Y_holdout - (yhat > decisionVariable);    
        
        yhat_others = glmval(detectOthersGLM.beta,X_others,'logit') > decisionVariable;
        pctCorr_detectOthers(i) = mean((yhat_others - Y_yeses(othersInd)) == 0);
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
    nTrialsOri{iexp} = nan(nOri,1);
    pctCorrOri_detectOthers = cell(1,nOri);
    pctCorrOri_targetOthers = cell(1,nOri);
    for iori = 1:nOri
        trInd = find(tOri == orientations(iori));
        nTrialsOri{iexp}(iori) = length(trInd);
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
        pctCorrOri_detectOthers{iori} = nan(1,nTrials);
        pctCorrOri_targetOthers{iori} = nan(1,nTrials);
        ysub_detectHoldouts = nan(1,nTrials);
        ysub_targetHoldouts = nan(1,nTrials);
        for i = 1:nTrials
            X_holdout = X(i,:);
            oriTrInd = trInd(i);
            othersInd = [1:(oriTrInd-1),(oriTrInd+1):nTrials_all];
            othersIndOri = ~ismember(trInd,oriTrInd);
            X_others = X_targets(othersInd,:);
            
            [~,~,detectOthersGLM] = glmfit(...
                X_others,Y_yeses(othersInd),'binomial');
            yhat_detectHoldout = glmval(...
                detectOthersGLM.beta,X_holdout,'logit') > decisionVariable;
            Y_holdout = Y_yeses(oriTrInd);
            ysub_detectHoldouts(i) = Y_holdout - yhat_detectHoldout;
            
            yhat_others = glmval(detectOthersGLM.beta,X(othersIndOri,:),'logit')...
                > decisionVariable;
            pctCorrOri_detectOthers{iori}(i) = mean(...
                (yhat_others - Y_yeses(trInd(othersIndOri))) == 0);
            
            [~,~,targetOthersGLM] = glmfit(...
                X_others,Y_targets(othersInd),'binomial');
            yhat_targetHoldouts = glmval(...
                targetOthersGLM.beta,X_holdout,'logit') > decisionVariable;
            Y_holdout = Y_targets(oriTrInd);
            ysub_targetHoldouts(i) = Y_holdout - yhat_targetHoldouts;
            
            yhat_others = glmval(targetOthersGLM.beta,X(othersIndOri,:),'logit')...
                > decisionVariable;
            pctCorrOri_targetOthers{iori}(i) = mean(...
                (yhat_others - Y_targets(trInd(othersIndOri))) == 0);
        end
        pctCorrXValOri_detect(iori) = mean(ysub_detectHoldouts == 0);
        pctCorrXValOri_target(iori) = mean(ysub_targetHoldouts == 0);
    end
    
    if doExptPlots
        setFigParams4Print('landscape')
        figure
        subplot 121
        x = [orientations 100];
        y = [cellfun(@mean,pctCorrOri_targetOthers).*100 ...
            mean(pctCorr_targetOthers).*100];
        yerr = [cellfun(@(a) ste(a,2),pctCorrOri_targetOthers).*100 ...
            ste(pctCorr_targetOthers,2).*100];
        h = errorbar(x,y,yerr,'ko-');
        h.MarkerFaceColor = [0 0 0];
        hold on
        y = [pctCorrXValOri_target.*100 pctCorr_targetXVal.*100];
        h = plot(x,y,'ko-');
        h.MarkerFaceColor = [1 1 1];
        xLabel = [cellfun(@(y)...
            num2str(round(y,2,'significant')),num2cell(x(1:end-1)),'unif',0),...
            {'All Trials'}];
        figXAxis([],'Orientation (deg)',[-10 110],x,xLabel)
        figYAxis([],'% Correct',[0 110])
        figAxForm
        legend({'Other Trials';'Holdout X Val'},'location','northeastoutside')
        title({'Target GLM';[dcExpt(iexp).mouse '-' dcExpt(iexp).date]})

        subplot 122
        x = [orientations 100];
        y = [cellfun(@mean,pctCorrOri_detectOthers).*100 ...
            mean(pctCorr_detectOthers).*100];
        yerr = [cellfun(@(a) ste(a,2),pctCorrOri_detectOthers).*100 ...
            ste(pctCorr_detectOthers,2).*100];
        h = errorbar(x,y,yerr,'ko-');
        h.MarkerFaceColor = [0 0 0];
        hold on
        y = [pctCorrXValOri_detect.*100 pctCorr_detectXVal.*100];
        h = plot(x,y,'ko-');
        h.MarkerFaceColor = [1 1 1];
        xLabel = [cellfun(@(y)...
            num2str(round(y,2,'significant')),num2cell(x(1:end-1)),'unif',0),...
            {'All Trials'}];
        figXAxis([],'Orientation (deg)',[-10 110],x,xLabel)
        figYAxis([],'% Correct',[0 110])
        figAxForm
        legend({'Other Trials';'Holdout X Val'},'location','northeastoutside')
        title({'Detect GLM';[dcExpt(iexp).mouse '-' dcExpt(iexp).date]})

        print(fullfile(fnout,[dcExpt(iexp).mouse '-' dcExpt(iexp).date '_xValidation_matched']),'-dpdf','-fillpage')
    end

        pctCorrOri_detect_all{iexp} = [cellfun(@mean,pctCorrOri_detectOthers).*100 ...
            mean(pctCorr_detectOthers).*100];
        pctCorrOri_detect_holdout{iexp} = [pctCorrXValOri_detect.*100 ...
            pctCorr_detectXVal.*100];

        pctCorrOri_target_all{iexp} = [cellfun(@mean,pctCorrOri_targetOthers).*100 ...
            mean(pctCorr_targetOthers).*100];
        pctCorrOri_target_holdout{iexp} = [pctCorrXValOri_target.*100 ...
            pctCorr_targetXVal.*100];
        exptOris{iexp} = [orientations 100];

end


pctCorrOri_detect_sub = cellfun(@(x,y) x-y,pctCorrOri_detect_all,...
    pctCorrOri_detect_holdout,'unif',0);
pctCorrOri_target_sub = cellfun(@(x,y) x-y,pctCorrOri_target_all,...
    pctCorrOri_target_holdout,'unif',0);

set(0,'defaultAxesFontSize',16)
figure
suptitle('Summary Across Experiments')
subplot 221
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
xLabel = [cellfun(@(a)...
    num2str(round(a,2,'significant')),num2cell(allOris(1:end-1)),'unif',0),...
    {'All Trials'}];
figXAxis([],'Orientation (deg)',[-10 110],allOris,xLabel)
figYAxis([],'% Correct',[0 110])
figAxForm
title('Detect GLM')

pctCorrOri_sub_catAllData = cell(1,length(allOris));
for iexp = 1:nexp
    for i = 1:length(allOris)
        ind = exptOris{iexp} == allOris(i);
        pctCorrOri_sub_catAllData{i} = ...
            [pctCorrOri_sub_catAllData{i} pctCorrOri_detect_sub{iexp}(ind)];
    end
end
subplot 222
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

subplot 223
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
title('Target GLM')

pctCorrOri_sub_catAllData = cell(1,length(allOris));
for iexp = 1:nexp
    for i = 1:length(allOris)
        ind = exptOris{iexp} == allOris(i);
        pctCorrOri_sub_catAllData{i} = ...
            [pctCorrOri_sub_catAllData{i} pctCorrOri_target_sub{iexp}(ind)];
    end
end
subplot 224
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

print(fullfile(fnout,[dataLabel 'pctCorr_byOri_matched']),'-dpdf','-fillpage')
