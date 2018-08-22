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
ampBinsFine = [0 exp(linspace(log(0.0001),log(1),7))];
ampBins = [0 exp(linspace(log(0.0001),log(1),3))];
weightsDiffZeroAlpha = 0.05/2;
nPC = 15;
weightLim = [-2.2 4];
corrLim = [-1 1];
%%
amplitudes = [];
orientations = [];
for iexp = 1:nexp
    tOri = dcExpt(iexp).trialOrientation;
    tOri = binAndRelabelTrialOrientation(tOri,oriBinsFine);
    orientations = cat(2,orientations,unique(tOri));
    tAmp = dcExpt(iexp).audTrialAmp;
    tAmp = binAndRelabelTrialOrientation(tAmp,ampBinsFine);
    amplitudes = cat(2,amplitudes,unique(tAmp));
end
allOris = unique(orientations);
allAmps = unique(amplitudes);

%% weights, predictions, cross validation
targetWeightsAll_vis = [];

nCells = zeros(1,nexp);
pcNumber = [];
exptName = cell(1,nexp);

pctCorrTarget_visXori = nan(nexp,length(allOris)+1);
pctCorrTarget_visXori_XVal = nan(nexp,length(allOris)+1);

for iexp=1:nexp
    
   disp([dcExpt(iexp).mouse '-' dcExpt(iexp).date])
    %****VISUAL TRIALS
   trOut=dcExpt(iexp).trialOutcome;
    [~, Y_targets_vis] = getStimAndBehaviorYs(trOut);
    
   X_targets=dcExpt(iexp).stimResp';
   
   cellIdx = dcExpt(iexp).signifResponsiveCells' & ...
       dcExpt(iexp).oriFitTheta90 < theta90Threshold;
   nCells(iexp) = sum(cellIdx);

   X_targets=X_targets(:,cellIdx);
   
   idx=find(~isnan(sum(X_targets,2)));
   Y_targets_vis=Y_targets_vis(idx,1);   
   X_targets=X_targets(idx,:);
   
   X_targets = zscore(X_targets);
%    
%    [c,score,latent] = pca(X_targets);
%    if sum(cellIdx) < nPC
%        npc = sum(cellIdx);
%        nCells(iexp) = sum(cellIdx);
%    else
%        npc = nPC;
%        nCells(iexp) =nPC;
%    end
%    X_prime_vis = score(:,1:npc);
%    X_prime_vis = zscore(X_prime_vis);
%    pcNumber = cat(2,pcNumber,1:npc);
   
    nTrials_all = size(X_targets,1);
   
    tOri = dcExpt(iexp).trialOrientation;
    tOri = binAndRelabelTrialOrientation(tOri,oriBinsFine);
    orientations = unique(tOri);
    nOri = length(orientations);
    
    stimPredictors = nan(nTrials_all,nOri);
    for iori = 1:nOri
        trInd = tOri == orientations(iori);
        if iori == 1
            stimPredictors(:,iori) = zeros(nTrials_all,1);
        else
            stimPredictors(:,iori) = trInd;
        end
    end
    
    X_targets_visStim = cat(2,X_targets,stimPredictors);
    
   [B_vis,stats_vis]=lassoglm(X_targets_visStim,Y_targets_vis,'binomial');
    [~,loc]=min(stats_vis.Deviance);

    B_best = B_vis(1:nCells(iexp)+nOri,loc);
    B_intc = stats_vis.Intercept(loc);
    
   p = 1;
   [temp1,dev1,targetGLM_vis]=glmfit(X_targets_visStim,Y_targets_vis,'binomial');
   idx=find(targetGLM_vis.p>p);
   temp1(idx)=0;
   C=eye(size(X_targets_visStim,2));
   B_targets = C*temp1(2:end,1);
   
    
%    [B_vis,stats_vis]=lassoglm(X_targets_visStim,Y_targets_vis,'binomial');
  
   
   targetWeightsAll_vis = cat(1,targetWeightsAll_vis,B_best);
%    detectWeightsAll_vis = cat(1,detectWeightsAll_vis,B_yeses);

    yhat_target = glmval([B_intc; B_best],X_targets,'logit') > decisionVariable;
    pctCorr_target_vis = mean((Y_targets_vis-yhat_target) == 0);
   
    pctCorrTarget_visXori(iexp,end) = pctCorr_target_vis;
   
    ysub_targetHoldouts = nan(1,nTrials_all);
    for i = 1:nTrials_all
        X_holdout = X_targets(i,:);
        othersInd = [1:(i-1),(i+1):nTrials_all];
        X_others = X_targets_visStim(othersInd,:);
       
       [B_vis_xval,stats_vis_xval]=lassoglm(X_others,Y_targets_vis(othersInd),...
           'binomial');
       [m,loc]=min(stats_vis.Deviance);

        B_best_xval = B_vis_xval(1:nCells(iexp),loc);
        B_intc_xval = stats_vis_xval.Intercept(loc);
        Y_holdout = Y_targets_vis(i);
        yhat = glmval([B_intc_xval; B_best_xval],X_holdout,'logit');
        ysub_targetHoldouts(i) = Y_holdout - (yhat > decisionVariable);      
    end
    pctCorr_targetXVal_vis = mean(ysub_targetHoldouts == 0);
        
    pctCorrTarget_visXori_XVal(iexp,end) = pctCorr_targetXVal_vis;
    
    pctCorrXValOri_target = nan(1,nOri);  
    pctCorrOri_target = nan(1,nOri);
   for iori = 1:nOri
        trInd = find(tOri == orientations(iori));
        X = X_targets(trInd,:);
                
        yhat = glmval([B_intc; B_best],X,'logit') > decisionVariable;
        pctCorrOri_target(iori) = mean((Y_targets_vis(trInd) - yhat) == 0);        
        
        nTrials = length(trInd);
        ysub_targetHoldouts = nan(1,nTrials);
        for i = 1:nTrials
            X_holdout = X(i,:);
            oriTrInd = trInd(i);
            othersInd = [1:(oriTrInd-1),(oriTrInd+1):nTrials_all];
            X_others = X_targets_visStim(othersInd,:);
            
            [B_vis_xval,stats_vis_xval]=lassoglm(...
                X_others,Y_targets_vis(othersInd),'binomial');
            [~,loc]=min(stats_vis.Deviance);

            B_best_xval = B_vis_xval(1:nCells(iexp),loc);
            B_intc_xval = stats_vis_xval.Intercept(loc);
            yhat_targetHoldouts = glmval(...
                [B_intc_xval; B_best_xval],X_holdout,'logit') > decisionVariable;
            Y_holdout = Y_targets_vis(oriTrInd);
            ysub_targetHoldouts(i) = Y_holdout - yhat_targetHoldouts;
        end
        pctCorrXValOri_target(iori) = mean(ysub_targetHoldouts == 0);
   end
    
   ind = ismember(allOris,orientations);
    pctCorrTarget_visXori(iexp,ind) = pctCorrOri_target;
   pctCorrTarget_visXori_XVal(iexp,ind) = pctCorrXValOri_target;
end

%%
targetWeightsAll_noReg_vis = [];

pctCorrTarget_noReg_visXori = nan(nexp,length(allOris)+1);
pctCorrTarget_noReg_visXori_XVal = nan(nexp,length(allOris)+1);

for iexp=1:nexp
    
   disp([dcExpt(iexp).mouse '-' dcExpt(iexp).date])
    %****VISUAL TRIALS
   trOut=dcExpt(iexp).trialOutcome;
    [~, Y_targets_vis] = getStimAndBehaviorYs(trOut);
    

   X_targets=dcExpt(iexp).stimResp';
   
   cellIdx = dcExpt(iexp).signifResponsiveCells' & ...
       dcExpt(iexp).oriFitTheta90 < theta90Threshold;
   nCells(iexp) = sum(cellIdx);

   X_targets=X_targets(:,cellIdx);
   
   idx=find(~isnan(sum(X_targets,2)));
   Y_targets_vis=Y_targets_vis(idx,1);   
   X_targets=X_targets(idx,:);
   
   [c,score,latent] = pca(X_targets);
   if sum(cellIdx) < nPC
       npc = sum(cellIdx);
   else
       npc = nPC;
   end
   X_targets = score(:,1:npc);
   X_targets = zscore(X_targets);
   pcNumber = cat(2,pcNumber,1:npc);
   
   C=eye(size(X_targets,2));
   p=1;
   [temp1,dev1,targetGLM_vis]=glmfit(X_targets,Y_targets_vis,'binomial');
   idx=find(targetGLM_vis.p>p);
   temp1(idx)=0;
   B_targets = C*temp1(2:end,1);
   
   targetWeightsAll_noReg_vis = cat(1,targetWeightsAll_noReg_vis,B_targets);
   
    yhat_target = glmval(targetGLM_vis.beta,X_targets,'logit') > decisionVariable;
    pctCorr_target_vis = mean((Y_targets_vis-yhat_target) == 0);
   
    pctCorrTarget_visXori(iexp,end) = pctCorr_target_vis;
    
    nTrials_all = size(X_targets,1);
   
    ysub_targetHoldouts = nan(1,nTrials_all);
    for i = 1:nTrials_all
        X_holdout = X_targets(i,:);
        othersInd = [1:(i-1),(i+1):nTrials_all];
        X_others = X_targets(othersInd,:);
        
        [~,~,targetOthersGLM] = glmfit(X_others,Y_targets_vis(othersInd),'binomial');        
        Y_holdout = Y_targets_vis(i);
        yhat = glmval(targetOthersGLM.beta,X_holdout,'logit');
        ysub_targetHoldouts(i) = Y_holdout - (yhat > decisionVariable);
               
    end
    pctCorr_targetXVal_vis = mean(ysub_targetHoldouts == 0);
        
    pctCorrTarget_visXori_XVal(iexp,end) = pctCorr_targetXVal_vis;
   
    tOri = dcExpt(iexp).trialOrientation;
    tOri = binAndRelabelTrialOrientation(tOri,oriBinsFine);
    orientations = unique(tOri);
    nOri = length(orientations);
    
    pctCorrXValOri_target = nan(1,nOri);  
    pctCorrOri_target = nan(1,nOri);
   for iori = 1:nOri
        trInd = find(tOri == orientations(iori));
        X = X_targets(trInd,:);
        
        yhat = glmval(targetGLM_vis.beta,X,'logit') > decisionVariable;
        pctCorrOri_target(iori) = mean((Y_targets_vis(trInd) - yhat) == 0);        
        
        nTrials = length(trInd);
        ysub_targetHoldouts = nan(1,nTrials);
        for i = 1:nTrials
            X_holdout = X(i,:);
            oriTrInd = trInd(i);
            othersInd = [1:(oriTrInd-1),(oriTrInd+1):nTrials_all];
            X_others = X_targets(othersInd,:);
            
            [~,~,targetOthersGLM] = glmfit(...
                X_others,Y_targets_vis(othersInd),'binomial');
            yhat_targetHoldouts = glmval(...
                targetOthersGLM.beta,X_holdout,'logit') > decisionVariable;
            Y_holdout = Y_targets_vis(oriTrInd);
            ysub_targetHoldouts(i) = Y_holdout - yhat_targetHoldouts;
        end
        pctCorrXValOri_target(iori) = mean(ysub_targetHoldouts == 0);
   end
    
   ind = ismember(allOris,orientations);
    pctCorrTarget_noReg_visXori(iexp,ind) = pctCorrOri_target;
   pctCorrTarget_noReg_visXori_XVal(iexp,ind) = pctCorrXValOri_target;   
end

%%
figure
suptitle('Visual Trials')
subplot 131
scatter(targetWeightsAll_vis,targetWeightsAll_noReg_vis,'ko')
hold on
plot(weightLim,weightLim,'k--')
figXAxis([],'Regularization',weightLim)
figYAxis([],'No Reg',weightLim)
figAxForm
title('Target LR Weights')

subplot 132
x = [allOris 100];
y = pctCorrTarget_visXori.*100;
h = errorbar(x,nanmean(y,1),ste(y,1),'ko');
h.MarkerFaceColor = 'k';
hold on
y = pctCorrTarget_noReg_visXori.*100;
h = errorbar(x,nanmean(y,1),ste(y,1),'ko');
h.MarkerFaceColor = 'b';
legend({'Reg';'No Reg'},'location','southeast')
figXAxis([],'Orientation',[-10 110])
figYAxis([],'% Correct',[0 110])
figAxForm

subplot 133
x = (pctCorrTarget_visXori - pctCorrTarget_visXori_XVal).*100;
y = (pctCorrTarget_noReg_visXori - pctCorrTarget_noReg_visXori_XVal).*100;
scatter(x(:),y(:),'ko');
hold on
plot([-2 10],[-2 10],'k--')
figXAxis([],'Reg',[-2 10])
figYAxis([],'No Reg',[-2 10])
figAxForm
title('% Correct (All - XVal)')