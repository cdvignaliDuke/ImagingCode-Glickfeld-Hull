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
ampBinsFine = [0 exp(linspace(log(0.0001),log(1),7))];
% ampBins = [0 exp(linspace(log(0.0001),log(1),3))];
ampBins = [0 0.0001 0.1 1];
weightsDiffZeroAlpha = 0.05/2;
nPC = 10;
maxCells = 8;
weightLim = [-2.2 5];
corrLim = [-1 1];
%%
amplitudes = [];
orientations = [];
for iexp = 1:nexp
    tOri = dcExpt(iexp).trialOrientation;
    tOri = binAndRelabelTrialOrientation(tOri,oriBins);
    orientations = cat(2,orientations,unique(tOri));
    tAmp = dcExpt(iexp).audTrialAmp;
    tAmp = binAndRelabelTrialOrientation(tAmp,ampBins);
    amplitudes = cat(2,amplitudes,unique(tAmp));
end
allOris = unique(orientations);
allAmps = unique(amplitudes);

%% weights, predictions, cross validation; Visual & Auditory Trials
targetWeightsAll_vis = [];
detectWeightsAll_vis = [];

targetCorrsAll_vis = [];
detectCorrsAll_vis = [];

targetWeightsAll_aud = [];
detectWeightsAll_aud = [];

targetCorrsAll_aud = [];
detectCorrsAll_aud = [];

nCells = zeros(1,nexp);
nTrialsExpt = zeros(1,nexp);
nTrialsExpt_aud = zeros(1,nexp);
pcNumber = [];
exptName = cell(1,nexp);

pctCorrTarget_visXori = nan(nexp,length(allOris)+1);
pctCorrDetect_visXori = nan(nexp,length(allOris)+1);
pctCorrTarget_visXori_XVal = nan(nexp,length(allOris)+1);
pctCorrDetect_visXori_XVal = nan(nexp,length(allOris)+1);

pctCorrTarget_audXamp = nan(nexp,length(allAmps)+1);
pctCorrDetect_audXamp = nan(nexp,length(allAmps)+1);
pctCorrTarget_audXamp_XVal = nan(nexp,length(allAmps)+1);
pctCorrDetect_audXamp_XVal = nan(nexp,length(allAmps)+1);

for iexp=1:nexp
    
   disp([dcExpt(iexp).mouse '-' dcExpt(iexp).date])
    %****VISUAL TRIALS
   trOut=dcExpt(iexp).trialOutcome;
    [Y_yeses_vis, Y_targets_vis] = getStimAndBehaviorYs(trOut);
    
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
   nTrialsExpt(iexp) = length(matchTrialsInd);

   X_targets=X_targets(matchTrialsInd,cellIdx);
   tOri = tOri(matchTrialsInd);
   
%    idx=find(~isnan(sum(X_targets,2)));
   Y_yeses_vis=Y_yeses_vis(matchTrialsInd);   
   Y_targets_vis=Y_targets_vis(matchTrialsInd);   
%    X_targets=X_targets(idx,:);
   
   [c,score,latent] = pca(X_targets);
   if sum(cellIdx) < nPC
       npc = sum(cellIdx);
   else
       npc = nPC;
   end
   X_prime_vis = score(:,1:npc);
   X_prime_vis = zscore(X_prime_vis);
   pcNumber = cat(2,pcNumber,1:npc);
   
   %**Cells or PCs?
   X_targets = zscore(X_targets);
   X_training = X_targets;
    
   detectCorr = corr(Y_yeses_vis,X_training)';
   targetCorr = corr(Y_targets_vis,X_training)';
   
   detectCorrsAll_vis = cat(1,detectCorrsAll_vis,detectCorr);
   targetCorrsAll_vis = cat(1,targetCorrsAll_vis,targetCorr);
   
   C=eye(size(X_training,2));
   p=1;
   [temp1,dev1,targetGLM_vis]=glmfit(X_training,Y_targets_vis,'binomial');
   [temp2,dev2,detectGLM_vis]=glmfit(X_training,Y_yeses_vis,'binomial');
   idx=find(targetGLM_vis.p>p|detectGLM_vis.p>p);
   temp1(idx)=0;
   temp2(idx)=0;
   B_targets = C*temp1(2:end,1);
   B_yeses = C*temp2(2:end,1);
   
   targetWeightsAll_vis = cat(1,targetWeightsAll_vis,B_targets);
   detectWeightsAll_vis = cat(1,detectWeightsAll_vis,B_yeses);
   
   dv_target = mean(Y_targets_vis);
   dv_detect = mean(Y_yeses_vis);

    yhat_target = glmval(targetGLM_vis.beta,X_training,'logit') > dv_target;
    pctCorr_target_vis = mean((Y_targets_vis-yhat_target) == 0);
   
    yhat_detect = glmval(detectGLM_vis.beta,X_training,'logit') > dv_detect;
    pctCorr_detect_vis = mean((Y_yeses_vis-yhat_detect) == 0);
    
    pctCorrTarget_visXori(iexp,end) = pctCorr_target_vis;
    pctCorrDetect_visXori(iexp,end) = pctCorr_detect_vis;
    
    nTrials_all = size(X_training,1);
   
    ysub_targetHoldouts = nan(1,nTrials_all);
    ysub_detectHoldouts = nan(1,nTrials_all);
    for i = 1:nTrials_all
        X_holdout = X_training(i,:);
        othersInd = [1:(i-1),(i+1):nTrials_all];
        X_others = X_training(othersInd,:);
        
        [~,~,targetOthersGLM] = glmfit(X_others,Y_targets_vis(othersInd),'binomial');        
        Y_holdout = Y_targets_vis(i);
        yhat = glmval(targetOthersGLM.beta,X_holdout,'logit');
        ysub_targetHoldouts(i) = Y_holdout - (yhat > dv_target);
        
        [~,~,detectOthersGLM] = glmfit(X_others,Y_yeses_vis(othersInd),'binomial');
        Y_holdout = Y_yeses_vis(i);
        yhat = glmval(detectOthersGLM.beta,X_holdout,'logit');
        ysub_detectHoldouts(i) = Y_holdout - (yhat > dv_detect);       
    end
    pctCorr_targetXVal_vis = mean(ysub_targetHoldouts == 0);
    pctCorr_detectXVal_vis = mean(ysub_detectHoldouts == 0);
        
    pctCorrTarget_visXori_XVal(iexp,end) = pctCorr_targetXVal_vis;
    pctCorrDetect_visXori_XVal(iexp,end) = pctCorr_detectXVal_vis;
   
    
    pctCorrXValOri_detect = nan(1,nOri);
    pctCorrXValOri_target = nan(1,nOri);  
    pctCorrOri_detect = nan(1,nOri);
    pctCorrOri_target = nan(1,nOri);
   for iori = 1:nOri
        trInd = find(tOri == orientations(iori));
        X = X_training(trInd,:);
        
        yhat = glmval(detectGLM_vis.beta,X,'logit') > dv_detect;
        pctCorrOri_detect(iori) = mean((Y_yeses_vis(trInd) - yhat) == 0);
        
        yhat = glmval(targetGLM_vis.beta,X,'logit') > dv_target;
        pctCorrOri_target(iori) = mean((Y_targets_vis(trInd) - yhat) == 0);        
        
        nTrials = length(trInd);
        ysub_detectHoldouts = nan(1,nTrials);
        ysub_targetHoldouts = nan(1,nTrials);
        for i = 1:nTrials
            X_holdout = X(i,:);
            oriTrInd = trInd(i);
            othersInd = [1:(oriTrInd-1),(oriTrInd+1):nTrials_all];
            X_others = X_training(othersInd,:);
            
            [~,~,detectOthersGLM] = glmfit(...
                X_others,Y_yeses_vis(othersInd),'binomial');
            yhat_detectHoldout = glmval(...
                detectOthersGLM.beta,X_holdout,'logit') > dv_detect;
            Y_holdout = Y_yeses_vis(oriTrInd);
            ysub_detectHoldouts(i) = Y_holdout - yhat_detectHoldout;
            
            [~,~,targetOthersGLM] = glmfit(...
                X_others,Y_targets_vis(othersInd),'binomial');
            yhat_targetHoldouts = glmval(...
                targetOthersGLM.beta,X_holdout,'logit') > dv_target;
            Y_holdout = Y_targets_vis(oriTrInd);
            ysub_targetHoldouts(i) = Y_holdout - yhat_targetHoldouts;
        end
        pctCorrXValOri_detect(iori) = mean(ysub_detectHoldouts == 0);
        pctCorrXValOri_target(iori) = mean(ysub_targetHoldouts == 0);
   end
    
   ind = ismember(allOris,orientations);
    pctCorrTarget_visXori(iexp,ind) = pctCorrOri_target;
    pctCorrDetect_visXori(iexp,ind) = pctCorrOri_detect;
   pctCorrTarget_visXori_XVal(iexp,ind) = pctCorrXValOri_target;
    pctCorrDetect_visXori_XVal(iexp,ind) = pctCorrXValOri_detect;

%****AUDITORY TRIALS   
   trOut=dcExpt(iexp).audTrialOutcome;
    [Y_yeses_aud, Y_targets_aud] = getStimAndBehaviorYs(trOut);
    
    
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
    
   nTrialsExpt_aud(iexp) = length(matchTrialsInd);

   X_targets=dcExpt(iexp).audStimResp';
   
   cellIdx = dcExpt(iexp).signifResponsiveCells' & ...
       dcExpt(iexp).oriFitTheta90 < theta90Threshold;

   X_targets=X_targets(matchTrialsInd,cellIdx);
   tAmp = tAmp(matchTrialsInd);
   
   [c,score,latent] = pca(X_targets);
   
   pve_expt = cumsum(latent)./sum(latent);
   if sum(cellIdx) < nPC
       npc = sum(cellIdx);
   else
       npc = nPC;
   end
   X_prime_aud = score(:,1:npc);
   X_prime_aud = zscore(X_prime_aud);
   pcNumber = cat(2,pcNumber,1:npc);
    
%    idx=find(~isnan(sum(X_targets,2)));
   Y_yeses_aud=Y_yeses_aud(matchTrialsInd);   
   Y_targets_aud=Y_targets_aud(matchTrialsInd);   
%    X_targets=X_targets(matchTrialsInd,:);

   %**Cells or PCs?
   X_training = X_prime_aud;
   
   detectCorr = corr(Y_yeses_aud,X_training)';
   targetCorr = corr(Y_targets_aud,X_training)';
   
   detectCorrsAll_aud = cat(1,detectCorrsAll_aud,detectCorr);
   targetCorrsAll_aud = cat(1,targetCorrsAll_aud,targetCorr);
   
   C=eye(size(X_training,2));
   p=1;
   [temp1,dev1,targetGLM_aud]=glmfit(X_training,Y_targets_aud,'binomial');
   [temp2,dev2,detectGLM_aud]=glmfit(X_training,Y_yeses_aud,'binomial');
   idx=find(targetGLM_aud.p>p|detectGLM_aud.p>p);
   temp1(idx)=0;
   temp2(idx)=0;
   B_targets = C*temp1(2:end,1);
   B_yeses = C*temp2(2:end,1);
   
   targetWeightsAll_aud = cat(1,targetWeightsAll_aud,B_targets);
   detectWeightsAll_aud = cat(1,detectWeightsAll_aud,B_yeses);
   
   dv_target = mean(Y_targets_aud);
   dv_detect = mean(Y_yeses_aud);
   
    yhat_target = glmval(targetGLM_aud.beta,X_training,'logit') > dv_target;
    pctCorr_target_aud = mean((Y_targets_aud-yhat_target) == 0);
   
    yhat_detect = glmval(detectGLM_aud.beta,X_training,'logit') > dv_detect;
    pctCorr_detect_aud = mean((Y_yeses_aud-yhat_detect) == 0);
    
    pctCorrTarget_audXamp(iexp,end) = pctCorr_target_aud;
    pctCorrDetect_audXamp(iexp,end) = pctCorr_detect_aud;
    
    nTrials_all = size(X_training,1);
   
    ysub_targetHoldouts = nan(1,nTrials_all);
    ysub_detectHoldouts = nan(1,nTrials_all);
    for i = 1:nTrials_all
        X_holdout = X_training(i,:);
        othersInd = [1:(i-1),(i+1):nTrials_all];
        X_others = X_training(othersInd,:);
        
        [~,~,targetOthersGLM] = glmfit(X_others,Y_targets_aud(othersInd),'binomial');        
        Y_holdout = Y_targets_aud(i);
        yhat = glmval(targetOthersGLM.beta,X_holdout,'logit');
        ysub_targetHoldouts(i) = Y_holdout - (yhat > dv_target);
        
        [~,~,detectOthersGLM] = glmfit(X_others,Y_yeses_aud(othersInd),'binomial');
        Y_holdout = Y_yeses_aud(i);
        yhat = glmval(detectOthersGLM.beta,X_holdout,'logit');
        ysub_detectHoldouts(i) = Y_holdout - (yhat > dv_detect);       
    end
    pctCorr_targetXVal_aud = mean(ysub_targetHoldouts == 0);
    pctCorr_detectXVal_aud = mean(ysub_detectHoldouts == 0);
    
    pctCorrTarget_audXamp_XVal(iexp,end) = pctCorr_targetXVal_aud;
    pctCorrDetect_audXamp_XVal(iexp,end) = pctCorr_detectXVal_aud;
   
%     tAmp = dcExpt(iexp).audTrialAmp;
%     tAmp = binAndRelabelTrialOrientation(tAmp,ampBinsFine);
%     amplitudes = unique(tAmp);
%     nAmp = length(amplitudes);
    
    pctCorrXValAmp_detect = nan(1,nAmp);
    pctCorrXValAmp_target = nan(1,nAmp);  
    pctCorrAmp_detect = nan(1,nAmp);
    pctCorrAmp_target = nan(1,nAmp);
   for iamp = 1:nAmp
        trInd = find(tAmp == amplitudes(iamp));
        X = X_training(trInd,:);
        
        yhat = glmval(detectGLM_aud.beta,X,'logit') > dv_detect;
        pctCorrAmp_detect(iamp) = mean((Y_yeses_aud(trInd) - yhat) == 0);
        
        yhat = glmval(targetGLM_aud.beta,X,'logit') > dv_target;
        pctCorrAmp_target(iamp) = mean((Y_targets_aud(trInd) - yhat) == 0);        
        
        nTrials = length(trInd);
        ysub_detectHoldouts = nan(1,nTrials);
        ysub_targetHoldouts = nan(1,nTrials);
        for i = 1:nTrials
            X_holdout = X(i,:);
            oriTrInd = trInd(i);
            othersInd = [1:(oriTrInd-1),(oriTrInd+1):nTrials_all];
            X_others = X_training(othersInd,:);
            
            [~,~,detectOthersGLM] = glmfit(...
                X_others,Y_yeses_aud(othersInd),'binomial');
            yhat_detectHoldout = glmval(...
                detectOthersGLM.beta,X_holdout,'logit') > dv_detect;
            Y_holdout = Y_yeses_aud(oriTrInd);
            ysub_detectHoldouts(i) = Y_holdout - yhat_detectHoldout;
            
            [~,~,targetOthersGLM] = glmfit(...
                X_others,Y_targets_aud(othersInd),'binomial');
            yhat_targetHoldouts = glmval(...
                targetOthersGLM.beta,X_holdout,'logit') > dv_target;
            Y_holdout = Y_targets_aud(oriTrInd);
            ysub_targetHoldouts(i) = Y_holdout - yhat_targetHoldouts;
        end
        pctCorrXValAmp_detect(iamp) = mean(ysub_detectHoldouts == 0);
        pctCorrXValAmp_target(iamp) = mean(ysub_targetHoldouts == 0);
   end
    
   ind = ismember(allAmps,amplitudes);
    pctCorrTarget_audXamp(iexp,ind) = pctCorrAmp_target;
    pctCorrDetect_audXamp(iexp,ind) = pctCorrAmp_detect;
   pctCorrTarget_audXamp_XVal(iexp,ind) = pctCorrXValAmp_target;
    pctCorrDetect_audXamp_XVal(iexp,ind) = pctCorrXValAmp_detect;
       
end

%% percent correct and cross valiation by orienation or amplitude

%**Visual Trials
pctCorrDetect_visXori_sub = pctCorrDetect_visXori-pctCorrDetect_visXori_XVal;
pctCorrTarget_visXori_sub = pctCorrTarget_visXori-pctCorrTarget_visXori_XVal;

setFigParams4Print('landscape')
set(0,'defaultAxesFontSize',12)
figure
suptitle('Summary Across Experiments - Visual')
subplot 231
% allExpt = nan(length(allOris),nexp);
for i = 1:nexp
    x = [allOris 100];
    y = pctCorrDetect_visXori(i,:).*100;
    ind = ~isnan(y);
    h = plot(x(ind),y(ind),'.');
    h.Color = [0.5 0.5 0.5];
    hold on
end
h = errorbar(x,nanmean(pctCorrDetect_visXori.*100,1),...
    ste(pctCorrDetect_visXori.*100,1),'ko');
h.MarkerFaceColor = [1 1 1];
xLabel = [cellfun(@(a)...
    num2str(round(a,2,'significant')),num2cell(allOris),'unif',0),...
    {'All Trials'}];
figXAxis([],'Orientation (deg)',[-10 110],x,xLabel)
figYAxis([],'% Correct',[0 110])
figAxForm
title('Detect - All Trials')

subplot 232
% allExpt = nan(length(allOris),nexp);
for i = 1:nexp
    x = [allOris 100];
    y = pctCorrDetect_visXori_XVal(i,:).*100;
    ind = ~isnan(y);
    h = plot(x(ind),y(ind),'.');
    h.Color = [0.5 0.5 0.5];
    hold on
end
h = errorbar(x,nanmean(pctCorrDetect_visXori_XVal.*100,1),...
    ste(pctCorrDetect_visXori_XVal.*100,1),'ko');
h.MarkerFaceColor = [1 1 1];
xLabel = [cellfun(@(a)...
    num2str(round(a,2,'significant')),num2cell(allOris),'unif',0),...
    {'All Trials'}];
figXAxis([],'Orientation (deg)',[-10 110],x,xLabel)
figYAxis([],'% Correct',[0 110])
figAxForm
title('Detect - Holdout Trials')

subplot 233
% allExpt = nan(length(allOris),nexp);
for i = 1:nexp
    x = [allOris 100];
    y = pctCorrDetect_visXori_sub(i,:).*100;
    ind = ~isnan(y);
    h = plot(x(ind),y(ind),'.');
    h.Color = [0.5 0.5 0.5];
    hold on
end
h = errorbar(x,nanmean(pctCorrDetect_visXori_sub.*100,1),...
    ste(pctCorrDetect_visXori_sub.*100,1),'ko');
h.MarkerFaceColor = [1 1 1];
xLabel = [cellfun(@(a)...
    num2str(round(a,2,'significant')),num2cell(allOris),'unif',0),...
    {'All Trials'}];
figXAxis([],'Orientation (deg)',[-10 110],x,xLabel)
figYAxis([],'% Correct (All-HO)',[-30 30])
figAxForm
title('Detect - Subtraction')
hline(0,'k--')

subplot 234
% allExpt = nan(length(allOris),nexp);
for i = 1:nexp
    x = [allOris 100];
    y = pctCorrTarget_visXori(i,:).*100;
    ind = ~isnan(y);
    h = plot(x(ind),y(ind),'.');
    h.Color = [0.5 0.5 0.5];
    hold on
end
h = errorbar(x,nanmean(pctCorrTarget_visXori.*100,1),...
    ste(pctCorrTarget_visXori.*100,1),'ko');
h.MarkerFaceColor = [1 1 1];
xLabel = [cellfun(@(a)...
    num2str(round(a,2,'significant')),num2cell(allOris),'unif',0),...
    {'All Trials'}];
figXAxis([],'Orientation (deg)',[-10 110],x,xLabel)
figYAxis([],'% Correct',[0 110])
figAxForm
title('Target - All Trials')

subplot 235
% allExpt = nan(length(allOris),nexp);
for i = 1:nexp
    x = [allOris 100];
    y = pctCorrTarget_visXori_XVal(i,:).*100;
    ind = ~isnan(y);
    h = plot(x(ind),y(ind),'.');
    h.Color = [0.5 0.5 0.5];
    hold on
end
h = errorbar(x,nanmean(pctCorrTarget_visXori_XVal.*100,1),...
    ste(pctCorrTarget_visXori_XVal.*100,1),'ko');
h.MarkerFaceColor = [1 1 1];
xLabel = [cellfun(@(a)...
    num2str(round(a,2,'significant')),num2cell(allOris),'unif',0),...
    {'All Trials'}];
figXAxis([],'Orientation (deg)',[-10 110],x,xLabel)
figYAxis([],'% Correct',[0 110])
figAxForm
title('Target - Holdout Trials')

subplot 236
% allExpt = nan(length(allOris),nexp);
for i = 1:nexp
    x = [allOris 100];
    y = pctCorrTarget_visXori_sub(i,:).*100;
    ind = ~isnan(y);
    h = plot(x(ind),y(ind),'.');
    h.Color = [0.5 0.5 0.5];
    hold on
end
h = errorbar(x,nanmean(pctCorrTarget_visXori_sub.*100,1),...
    ste(pctCorrTarget_visXori_sub.*100,1),'ko');
h.MarkerFaceColor = [1 1 1];
xLabel = [cellfun(@(a)...
    num2str(round(a,2,'significant')),num2cell(allOris),'unif',0),...
    {'All Trials'}];
figXAxis([],'Orientation (deg)',[-10 110],x,xLabel)
figYAxis([],'% Correct (All-HO)',[-30 30])
figAxForm
title('Target - Subtraction')
hline(0,'k--')

print(fullfile(fnout,'pctCorrectWithXVal_matched_vis'),'-dpdf','-fillpage')

% **Auditory Trials
pctCorrDetect_audXamp_sub = pctCorrDetect_audXamp-pctCorrDetect_audXamp_XVal;
pctCorrTarget_audXamp_sub = pctCorrTarget_audXamp-pctCorrTarget_audXamp_XVal;
allAmps_plot = [0.001 allAmps(2:end) 5].*100;

set(0,'defaultAxesFontSize',12)
figure
suptitle('Summary Across Experiments - Auditory')
subplot 231
% allExpt = nan(length(allOris),nexp);
for i = 1:nexp
    x = allAmps_plot;
    y = pctCorrDetect_audXamp(i,:).*100;
    ind = ~isnan(y);
    h = plot(x(ind),y(ind),'.');
    h.Color = [0.5 0.5 0.5];
    hold on
end
h = errorbar(x,nanmean(pctCorrDetect_audXamp.*100,1),...
    ste(pctCorrDetect_audXamp.*100,1),'ko');
h.MarkerFaceColor = [1 1 1];
figXAxis([],'Amplitude (% Max)',[0 1100])
ax = gca;
set(ax,'xscale','log')
figYAxis([],'% Correct',[0 110])
figAxForm
title('Detect - All Trials')

subplot 232
% allExpt = nan(length(allAmps),nexp);
for i = 1:nexp
    x = allAmps_plot;
    y = pctCorrDetect_audXamp_XVal(i,:).*100;
    ind = ~isnan(y);
    h = plot(x(ind),y(ind),'.');
    h.Color = [0.5 0.5 0.5];
    hold on
end
h = errorbar(x,nanmean(pctCorrDetect_audXamp_XVal.*100,1),...
    ste(pctCorrDetect_audXamp_XVal.*100,1),'ko');
h.MarkerFaceColor = [1 1 1];
figXAxis([],'Amplitude (% Max)',[0 1100])
ax = gca;
set(ax,'xscale','log')
figYAxis([],'% Correct',[0 110])
figAxForm
title('Detect - Holdout Trials')

subplot 233
% allExpt = nan(length(allAmps),nexp);
for i = 1:nexp
    x = allAmps_plot;
    y = pctCorrDetect_audXamp_sub(i,:).*100;
    ind = ~isnan(y);
    h = plot(x(ind),y(ind),'.');
    h.Color = [0.5 0.5 0.5];
    hold on
end
h = errorbar(x,nanmean(pctCorrDetect_audXamp_sub.*100,1),...
    ste(pctCorrDetect_audXamp_sub.*100,1),'ko');
h.MarkerFaceColor = [1 1 1];
figXAxis([],'Amplitude (% Max)',[0 1100])
ax = gca;
set(ax,'xscale','log')
figYAxis([],'% Correct',[-30 30])
figAxForm
title('Detect - Subtraction')
hline(0,'k--')

subplot 234
% allExpt = nan(length(allAmps),nexp);
for i = 1:nexp
    x = allAmps_plot;
    y = pctCorrTarget_audXamp(i,:).*100;
    ind = ~isnan(y);
    h = plot(x(ind),y(ind),'.');
    h.Color = [0.5 0.5 0.5];
    hold on
end
h = errorbar(x,nanmean(pctCorrTarget_audXamp.*100,1),...
    ste(pctCorrTarget_audXamp.*100,1),'ko');
h.MarkerFaceColor = [1 1 1];
figXAxis([],'Amplitude (% Max)',[0 1100])
ax = gca;
set(ax,'xscale','log')
figYAxis([],'% Correct',[0 110])
figAxForm
title('Target - All Trials')

subplot 235
% allExpt = nan(length(allAmps),nexp);
for i = 1:nexp
    x = allAmps_plot;
    y = pctCorrTarget_audXamp_XVal(i,:).*100;
    ind = ~isnan(y);
    h = plot(x(ind),y(ind),'.');
    h.Color = [0.5 0.5 0.5];
    hold on
end
h = errorbar(x,nanmean(pctCorrTarget_audXamp_XVal.*100,1),...
    ste(pctCorrTarget_audXamp_XVal.*100,1),'ko');
h.MarkerFaceColor = [1 1 1];
figXAxis([],'Amplitude (% Max)',[0 1100])
ax = gca;
set(ax,'xscale','log')
figYAxis([],'% Correct',[0 110])
figAxForm
title('Target - Holdout Trials')

subplot 236
% allExpt = nan(length(allAmps),nexp);
for i = 1:nexp
    x = allAmps_plot;
    y = pctCorrTarget_audXamp_sub(i,:).*100;
    ind = ~isnan(y);
    h = plot(x(ind),y(ind),'.');
    h.Color = [0.5 0.5 0.5];
    hold on
end
h = errorbar(x,nanmean(pctCorrTarget_audXamp_sub.*100,1),...
    ste(pctCorrTarget_audXamp_sub.*100,1),'ko');
h.MarkerFaceColor = [1 1 1];
figXAxis([],'Amplitude (% Max)',[0 1100])
ax = gca;
set(ax,'xscale','log')
figYAxis([],'% Correct',[-30 30])
figAxForm
title('Target - Subtraction')
hline(0,'k--')

print(fullfile(fnout,'pctCorrectWithXVal_matched_aud'),'-dpdf','-fillpage')

%%

