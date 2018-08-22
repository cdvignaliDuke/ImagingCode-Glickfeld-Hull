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
targetWeightsAll = [];
targetWeightsAll_lasso = [];
targetWeightsAll_lassoLocZero = [];

nCells = zeros(1,nexp);
pcNumber = [];
exptName = cell(1,nexp);

pctCorrTarget_visXori = nan(nexp,length(allOris)+1);
pctCorrTarget_visXori_XVal = nan(nexp,length(allOris)+1);

pctCorrTarget_visXori_locZero = nan(nexp,length(allOris)+1);
pctCorrTarget_visXori_XVal_locZero = nan(nexp,length(allOris)+1);

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
   
    nTrials_all = size(X_targets,1);
    
    if isempty(X_targets)
        continue
    end
   
    tOri = dcExpt(iexp).trialOrientation;
    tOri = binAndRelabelTrialOrientation(tOri,oriBinsFine);
    orientations = unique(tOri);
    nOri = length(orientations);
% %     
%     stimPredictors = nan(nTrials_all,nOri);
%     for iori = 1:nOri
%         trInd = tOri == orientations(iori);
%         if iori == 1
%             stimPredictors(:,iori) = zeros(nTrials_all,1);
%         else
%             stimPredictors(:,iori) = trInd;
%         end
%     end
%     
%     X_targets_visStim = cat(2,X_targets,stimPredictors);

   [B_vis,stats_vis]=lassoglm(X_targets,Y_targets_vis,'binomial');    
%    [B_vis,stats_vis]=lassoglm(X_targets_visStim,Y_targets_vis,'binomial',...
%        'Offset',Y_targets_vis);
    [~,loc]=min(stats_vis.Deviance);
    [~,locZero] = find(any(B_vis == 0),1);

    B_best = B_vis(:,loc);
    B_locZero = B_vis(:,locZero-1);
    B_intc = stats_vis.Intercept(loc);
    B_intc_locZero = stats_vis.Intercept(locZero-1);
    
   p = 1;
   [temp1,dev1,targetGLM_vis]=glmfit(X_targets,Y_targets_vis,'binomial');
   idx=find(targetGLM_vis.p>p);
   temp1(idx)=0;
   C=eye(size(X_targets,2));
   B_targets = C*temp1(2:end,1);
   
   targetWeightsAll = cat(1,targetWeightsAll,B_targets);
   targetWeightsAll_lasso = cat(1,targetWeightsAll_lasso,B_best);
   targetWeightsAll_lassoLocZero = cat(1,targetWeightsAll_lassoLocZero,B_locZero);

   decisionVariable = mean(Y_targets_vis);
   
    yhat_target = glmval([B_intc; B_best],X_targets,'logit') > decisionVariable;
    pctCorr_target_vis = mean((Y_targets_vis-yhat_target) == 0);
   
    pctCorrTarget_visXori(iexp,end) = pctCorr_target_vis;
    
    yhat_target = glmval([B_intc_locZero; B_locZero],X_targets,'logit')...
        > decisionVariable;
    pctCorr_target_vis = mean((Y_targets_vis-yhat_target) == 0);
   
    pctCorrTarget_visXori_locZero(iexp,end) = pctCorr_target_vis;
   
    ysub_targetHoldouts = nan(1,nTrials_all);
    ysub_targetHoldouts_locZero = nan(1,nTrials_all);
    for i = 1:nTrials_all
        X_holdout = X_targets(i,:);
        othersInd = [1:(i-1),(i+1):nTrials_all];
        X_others = X_targets(othersInd,:);
       
       [B_vis_xval,stats_vis_xval]=lassoglm(X_others,Y_targets_vis(othersInd),...
           'binomial');
       [~,loc]=min(stats_vis.Deviance);
        [~,locZero] = find(any(B_vis == 0),1);

        B_best_xval = B_vis_xval(:,loc);
        B_intc_xval = stats_vis_xval.Intercept(loc);
                
        B_locZero_xval = B_vis_xval(:,loc);
        B_intcLocZero_xval = stats_vis_xval.Intercept(loc);
        
        Y_holdout = Y_targets_vis(i);
        yhat = glmval([B_intc_xval; B_best_xval],X_holdout,'logit');
        ysub_targetHoldouts(i) = Y_holdout - (yhat > decisionVariable); 
        
        yhat = glmval([B_intcLocZero_xval; B_locZero_xval],X_holdout,'logit');
        ysub_targetHoldouts_locZero(i) = Y_holdout - (yhat > decisionVariable); 
    end
    pctCorr_targetXVal_vis = mean(ysub_targetHoldouts == 0);
    pctCorr_targetXVal_vis_locZero = mean(ysub_targetHoldouts == 0);
        
    pctCorrTarget_visXori_XVal(iexp,end) = pctCorr_targetXVal_vis;
    pctCorrTarget_visXori_XVal_locZero(iexp,end) = pctCorr_targetXVal_vis_locZero;
    
    pctCorrXValOri_target = nan(1,nOri);  
    pctCorrOri_target = nan(1,nOri);
    pctCorrXValOri_target_locZero = nan(1,nOri);  
    pctCorrOri_target_locZero = nan(1,nOri);
   for iori = 1:nOri
        trInd = find(tOri == orientations(iori));
        X = X_targets(trInd,:);
                
        yhat = glmval([B_intc; B_best],X,'logit') > decisionVariable;
        pctCorrOri_target(iori) = mean((Y_targets_vis(trInd) - yhat) == 0);  
       
        yhat = glmval([B_intc_locZero; B_locZero],X,'logit') > decisionVariable;
        pctCorrOri_target_locZero(iori) = mean((Y_targets_vis(trInd) - yhat) == 0);  
        
        nTrials = length(trInd);
        ysub_targetHoldouts = nan(1,nTrials);
        ysub_targetHoldouts_locZero = nan(1,nTrials);
        for i = 1:nTrials
            X_holdout = X(i,:);
            oriTrInd = trInd(i);
            othersInd = [1:(oriTrInd-1),(oriTrInd+1):nTrials_all];
            X_others = X_targets(othersInd,:);
            
            [B_vis_xval,stats_vis_xval]=lassoglm(...
                X_others,Y_targets_vis(othersInd),'binomial');
            [~,loc]=min(stats_vis.Deviance);
            locZero = find(any(B_vis_xval == 0,1),1);

            B_best_xval = B_vis_xval(:,loc);
            B_intc_xval = stats_vis_xval.Intercept(loc);
            
            B_locZero_xval = B_vis_xval(:,locZero);
            B_intcLocZero_xval = stats_vis_xval.Intercept(locZero);
            
            yhat_targetHoldouts = glmval(...
                [B_intc_xval; B_best_xval],X_holdout,'logit') > decisionVariable;
            Y_holdout = Y_targets_vis(oriTrInd);
            ysub_targetHoldouts(i) = Y_holdout - yhat_targetHoldouts;
            
            yhat_targetHoldouts = glmval(...
                [B_intcLocZero_xval; B_locZero_xval],X_holdout,'logit') > decisionVariable;
            ysub_targetHoldouts_locZero(i) = Y_holdout - yhat_targetHoldouts;
        end
        pctCorrXValOri_target(iori) = mean(ysub_targetHoldouts == 0);
        pctCorrXValOri_target_locZero(iori) = mean(ysub_targetHoldouts_locZero == 0);
   end
    
   ind = ismember(allOris,orientations);
    pctCorrTarget_visXori(iexp,ind) = pctCorrOri_target;
   pctCorrTarget_visXori_XVal(iexp,ind) = pctCorrXValOri_target;
    pctCorrTarget_visXori_locZero(iexp,ind) = pctCorrOri_target_locZero;
   pctCorrTarget_visXori_XVal_locZero(iexp,ind) = pctCorrXValOri_target_locZero;
end

%%

%%
figure
suptitle('Visual Trials: Target LR Weights')
subplot 221
scatter(targetWeightsAll,targetWeightsAll_lasso,'ko')
hold on
plot(weightLim,weightLim,'k--')
figXAxis([],'No Reg',weightLim)
figYAxis([],'Regularized',weightLim)
figAxForm
title('Minimum Deviance')

subplot 222
scatter(targetWeightsAll,targetWeightsAll_lassoLocZero,'ko')
hold on
plot(weightLim,weightLim,'k--')
figXAxis([],'No Reg',weightLim)
figYAxis([],'Regularized',weightLim)
figAxForm
title('Lowest Weight Minimized')

subplot 223
x = [allOris 100];
y = pctCorrTarget_visXori_XVal.*100;
h = errorbar(x,nanmean(y,1),ste(y,1),'ko');
h.MarkerFaceColor = 'k';
hold on
for i = 1:nexp
    ind = ~isnan(y(i,:));
    if ~all(ind == 0)
        h = plot(x(ind), y(i,ind),'-');
        h.Color = [0.5 0.5 0.5];
    end
end
figXAxis([],'Orientation',[-10 110])
figYAxis([],'% Correct',[0 110])
figAxForm
title('XVal: Minimum Deviance')

subplot 224
x = [allOris 100];
y = pctCorrTarget_visXori_XVal_locZero.*100;
h = errorbar(x,nanmean(y,1),ste(y,1),'ko');
h.MarkerFaceColor = 'k';
hold on
for i = 1:nexp
    ind = ~isnan(y(i,:));
    if ~all(ind == 0)
        h = plot(x(ind), y(i,ind),'-');
        h.Color = [0.5 0.5 0.5];
    end
end
figXAxis([],'Orientation',[-10 110])
figYAxis([],'% Correct',[0 110])
figAxForm
title('XVal: Lowest Weight Minimized')


print(fullfile(fnout,'LR_targetweights_testReg'),'-dpdf','-fillpage')