close all
clear all

FSAV_V1_decode

load('X:\home\ashley\Analysis\FSAV Choice\FSAV_decodeData.mat')
load('X:\home\ashley\Analysis\FSAV Summaries\FSAV_attentionV1\attentionV1_startAlign_FSAV_attnData')
fnout = 'X:\home\ashley\Analysis\FSAV Choice';
nexp = length(dcExpt);

doExptPlots = 0;
%%
outcomecolmat = {'k';'r'};
avcolmat = {'k','c'};

ampBinsFine = [0 exp(linspace(log(0.0001),log(1),7))];
ampBins = [0 0.0001 0.1 1];
weightsDiffZeroAlpha = 0.05/2;
nPC = 15;
weightLim = [-2.2 5];
weightSubLim = [-2.5 2.5];
corrLim = [-1 1];
siLim = [-12 12];
adaptLim = [-1.5 1.5];

minTrN = 17;
theta90Threshold = 11.25;
decisionVariable = 0.5;
detectThresholdHR = 0.8;
detectThresholdOri = 45;

oriBinsFine = [0 1 16 22.5 32 45 64 90];
oriBins = [0 1 32 90];
%%

amplitudes = [];
orientations = [];
for iexp = 1:nexp
    tOri = dcExpt(iexp).trialOrientation;
    tOri = binAndRelabelTrialOrientation(tOri,oriBins);
    orientations = cat(2,orientations,unique(tOri));
    tAmp = dcExpt(iexp).audTrialAmp;
    tAmp = binAndRelabelTrialOrientation(tAmp,ampBinsFine);
    amplitudes = cat(2,amplitudes,unique(tAmp));
end
allOris = unique(orientations);
allAmps = unique(amplitudes);
%% weights, predictions, cross validation; Visual & Auditory Trials
si = [];
avMod = [];
firstRespCells = [];
targetRespCells = [];
pctVisAdapt = [];
pctAudAdapt = [];

respCells_targetOnly = [];
respCells_baseOnly = [];
respCells_targetAndBase = [];
respCells_lateBase = [];
suppCells = [];

targetWeightsAll_vis = [];
detectWeightsAll_vis = [];

targetCorrsAll_vis = [];
detectCorrsAll_vis = [];

targetCorrsAll_catch = [];
detectCorrsAll_catch = [];

targetWeightsAll_aud = [];
detectWeightsAll_aud = [];

targetCorrsAll_aud = [];
detectCorrsAll_aud = [];

nCells = zeros(1,nexp);
pcNumber = [];
exptName = cell(1,nexp);

pctCorrectTarget_Xval_vis = nan(1,nexp);
pctCorrectDetect_Xval_vis = nan(1,nexp);
pctCorrectTarget_catch = nan(1,nexp);
pctCorrectDetect_catch = nan(1,nexp);
pctCorrectTarget_aud = nan(1,nexp);
pctCorrectDetect_aud = nan(1,nexp);

for iexp=1:nexp
    
    if ~strcmp(dcExpt(iexp).mouse,attnInfoExpt(iexp).ms)
        error('Mouse name from attention data does not match decode data')
    end
    
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
    tOri = tOri(matchTrialsInd);

   X_targets=dcExpt(iexp).stimResp';
   
   cellIdx = dcExpt(iexp).signifResponsiveCells' & ...
       dcExpt(iexp).oriFitTheta90 < theta90Threshold;
   nCells(iexp) = sum(cellIdx);

   X_targets=X_targets(matchTrialsInd,cellIdx);
   Y_yeses_vis=Y_yeses_vis(matchTrialsInd);   
   Y_targets_vis=Y_targets_vis(matchTrialsInd); 
   
   X_targets = zscore(X_targets);
   
   si = cat(2,si,attnInfoExpt(iexp).si(cellIdx));
   avMod = cat(2,avMod,attnInfoExpt(iexp).avModTest(cellIdx));
   
   tc = dcExpt(iexp).targetResponsiveCells(cellIdx);
   fbrc = dcExpt(iexp).firstBaseResponsiveCells(cellIdx);
   
   respCells_targetOnly = cat(1,respCells_targetOnly,tc & ~fbrc);
   respCells_baseOnly = cat(1,respCells_baseOnly,~tc & fbrc);
   respCells_targetAndBase = cat(1,respCells_targetAndBase, tc & fbrc);
   respCells_lateBase = cat(1,respCells_lateBase,...
       attnInfoExpt(iexp).lateBaseRespCells(cellIdx)');   
   suppCells = cat(1,suppCells,attnInfoExpt(iexp).lateBaseSuppCells(cellIdx)');
   
   baseResponsiveCells = dcExpt(iexp).firstBaseResponsiveCells(cellIdx);
   firstRespCells = cat(1,firstRespCells,baseResponsiveCells);
   targetRespCells = cat(1,targetRespCells, ...
       dcExpt(iexp).targetResponsiveCells(cellIdx));
   ad_vis = attnInfoExpt(iexp).visAdapt(cellIdx);
   ad_aud = attnInfoExpt(iexp).audAdapt(cellIdx);
   ad_vis(~baseResponsiveCells) = nan;
   ad_aud(~baseResponsiveCells) = nan;
   
   pctVisAdapt = cat(2,pctVisAdapt,ad_vis);
   pctAudAdapt = cat(2,pctAudAdapt,ad_aud);
   
   idx=find(~isnan(sum(X_targets,2)));
   Y_yeses_vis=Y_yeses_vis(idx,1);   
   Y_targets_vis=Y_targets_vis(idx,1);   
   X_targets=X_targets(idx,:);
   
%    [c,score,latent] = pca(X_targets);
%    if sum(cellIdx) < nPC
%        npc = sum(cellIdx);
%    else
%        npc = nPC;
%    end
%    X_prime_vis = score(:,1:npc);
%    X_prime_vis = zscore(X_prime_vis);
%    pcNumber = cat(2,pcNumber,1:npc);
%     
   detectCorr = corr(Y_yeses_vis,X_targets)';
   targetCorr = corr(Y_targets_vis,X_targets)';
   
   detectCorrsAll_vis = cat(1,detectCorrsAll_vis,detectCorr);
   targetCorrsAll_vis = cat(1,targetCorrsAll_vis,targetCorr);
   
   C=eye(size(X_targets,2));
   p=1;
   [temp1,dev1,targetGLM_vis]=glmfit(X_targets,Y_targets_vis,'binomial');
   [temp2,dev2,detectGLM_vis]=glmfit(X_targets,Y_yeses_vis,'binomial');
   idx=find(targetGLM_vis.p>p|detectGLM_vis.p>p);
   temp1(idx)=0;
   temp2(idx)=0;
   B_targets = C*temp1(2:end,1);
   B_yeses = C*temp2(2:end,1);
   
   targetWeightsAll_vis = cat(1,targetWeightsAll_vis,B_targets);
   detectWeightsAll_vis = cat(1,detectWeightsAll_vis,B_yeses);
   
   dv_target = mean(Y_targets_vis);
   dv_detect = mean(Y_yeses_vis);
   
   
   if ~isempty(dcExpt(iexp).catchResp)
   %***test holdout trials on matched trials
       tCatchOri = dcExpt(iexp).catchOrientation;
       tCatchOri = binAndRelabelTrialOrientation(tCatchOri,oriBins);
       catchOri = round(unique(tCatchOri),2,'significant');
       catchOriID = findgroups(tCatchOri);
       nCatchOri = histcounts(catchOriID,1:max(catchOriID)+1);
       nCatch = length(catchOri);
       matchedValidVisTrialInd = [];
       for i = 1:nCatch
           ind1 = find(tOri == catchOri(i));
           n = nCatchOri(i);
           ind2 = randsample(ind1,n);
           matchedValidVisTrialInd = sort(cat(2,matchedValidVisTrialInd,ind2));
       end

       X_catchMatched = X_targets(matchedValidVisTrialInd,:);

       nt_match = length(matchedValidVisTrialInd);
       nt_all = size(X_targets,1);
       ysub_target = nan(1,nt_match);
       ysub_detect = nan(1,nt_match);
       for i = 1:nt_match
           ind = matchedValidVisTrialInd(i);
            X_holdout = X_targets(ind,:);
            othersInd = [1:(ind-1),(ind+1):nt_all];
            X_others = X_targets(othersInd,:);

           [~,~,targetGLM_holdout]=glmfit(X_others,Y_targets_vis(othersInd),'binomial');
           [~,~,detectGLM_holdout]=glmfit(X_others,Y_yeses_vis(othersInd),'binomial');

           y_holdout = Y_targets_vis(ind);
           yhat = glmval(targetGLM_holdout.beta,X_holdout,'logit');
           ysub_target(i) = y_holdout - (yhat > dv_target);

           y_holdout = Y_yeses_vis(ind);
           yhat = glmval(detectGLM_holdout.beta,X_holdout,'logit');
           ysub_detect(i) = y_holdout - (yhat > dv_detect);
       end
       pctCorrectTarget_Xval_vis(iexp) = mean(ysub_target == 0);
       pctCorrectDetect_Xval_vis(iexp) = mean(ysub_detect == 0);

   %***CATCH TRIALS
       X_targets = dcExpt(iexp).catchResp';
       X_targets = X_targets(:,cellIdx);
       trOut=dcExpt(iexp).catchOutcome;
        [Y_yeses_catch, Y_targets_catch] = getStimAndBehaviorYs(trOut);

       detectCorr_catch = corr(Y_yeses_catch,X_targets);
       targetCorr_catch = corr(Y_targets_catch,X_targets);
   
        yhat = glmval(targetGLM_vis.beta,X_targets,'logit') > dv_target;
        pctCorrectTarget_catch(iexp) = mean((Y_yeses_catch - yhat) == 0);
        
        yhat = glmval(detectGLM_vis.beta,X_targets,'logit') > dv_detect;
        pctCorrectDetect_catch(iexp) = mean((Y_targets_catch - yhat) == 0);        
   else
       detectCorr_catch = nan(1,sum(cellIdx));
       targetCorr_catch = nan(1,sum(cellIdx));
   end
    targetCorrsAll_catch = cat(2,targetCorrsAll_catch,targetCorr_catch);
    detectCorrsAll_catch =  cat(2,detectCorrsAll_catch,detectCorr_catch);
   
   %***AUDITORY TRIALS
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

   X_targets=dcExpt(iexp).audStimResp';
   X_targets=X_targets(matchTrialsInd,cellIdx);
   Y_yeses_aud=Y_yeses_aud(matchTrialsInd);   
   Y_targets_aud=Y_targets_aud(matchTrialsInd);  
   
   X_targets = zscore(X_targets);
      
   detectCorr = corr(Y_yeses_aud,X_targets)';
   targetCorr = corr(Y_targets_aud,X_targets)';
   
   detectCorrsAll_aud = cat(1,detectCorrsAll_aud,detectCorr);
   targetCorrsAll_aud = cat(1,targetCorrsAll_aud,targetCorr);
   
   C=eye(size(X_targets,2));
   p=1;
   [temp1,dev1,targetGLM_aud]=glmfit(X_targets,Y_targets_aud,'binomial');
   [temp2,dev2,detectGLM_aud]=glmfit(X_targets,Y_yeses_aud,'binomial');
   idx=find(targetGLM_aud.p>p|detectGLM_aud.p>p);
   temp1(idx)=0;
   temp2(idx)=0;
   B_targets = C*temp1(2:end,1);
   B_yeses = C*temp2(2:end,1);
   
   targetWeightsAll_aud = cat(1,targetWeightsAll_aud,B_targets);
   detectWeightsAll_aud = cat(1,detectWeightsAll_aud,B_yeses);
   
    yhat = glmval(targetGLM_vis.beta,X_targets,'logit') > dv_target;
    pctCorrectTarget_aud(iexp) = mean((Y_yeses_aud - yhat) == 0);

    yhat = glmval(detectGLM_vis.beta,X_targets,'logit') > dv_detect;
    pctCorrectDetect_aud(iexp) = mean((Y_targets_aud - yhat) == 0); 
  
   
end
respCells_targetOnly = logical(respCells_targetOnly);
respCells_firstBase = logical(respCells_baseOnly);
respCells_targetAndBase = logical(respCells_targetAndBase);
respCells_lateBaseOnly = logical(respCells_lateBase &...
    ~respCells_baseOnly & ~respCells_targetAndBase);
suppCells = logical(suppCells);
%%
targetWeights_sub = targetWeightsAll_vis - targetWeightsAll_aud;
detectWeights_sub = detectWeightsAll_vis - detectWeightsAll_aud;
%%
setFigParams4Print('portrait')
set(0,'defaultAxesFontSize',10)
figure
suptitle('Visual Trials Model, Significantly Modulated Cells in Red')
ind = logical(avMod);
subplot 321
scatter(si,targetCorrsAll_vis,'ko')
hold on
scatter(si(ind),targetCorrsAll_vis(ind),'ro')
figXAxis([],'V-A Selectivity',siLim)
figYAxis([],'Correlation',corrLim)
figAxForm
title('Target')
subplot 322
scatter(si,detectCorrsAll_vis,'ko')
hold on
scatter(si(ind),detectCorrsAll_vis(ind),'ro')
figXAxis([],'V-A Selectivity',siLim)
figYAxis([],'Correlation',corrLim)
figAxForm
title('Detect')
subplot 323
scatter(si,targetWeightsAll_vis,'ko')
hold on
scatter(si(ind),targetWeightsAll_vis(ind),'ro')
figXAxis([],'V-A Selectivity',siLim)
figYAxis([],'Weight',weightLim)
figAxForm
title('Target')
subplot 324
scatter(si,detectWeightsAll_vis,'ko')
hold on
scatter(si(ind),detectWeightsAll_vis(ind),'ro')
figXAxis([],'V-A Selectivity',siLim)
figYAxis([],'Weight',weightLim)
figAxForm
title('Detect')

subplot 325
ind = si < 0;
scatter(targetWeightsAll_vis(ind),detectWeightsAll_vis(ind),'bo');
hold on
ind = si > 0;
scatter(targetWeightsAll_vis(ind),detectWeightsAll_vis(ind),'mo');
figXAxis([],'Target Weight',weightLim)
figYAxis([],'Detect Weight',weightLim)
figAxForm
title('V-A Selectivity')
legend({'si < 0';'si > 0'},'location','northwest')
subplot 326
scatter(targetWeightsAll_vis,detectWeightsAll_vis,'ko');
hold on
ind = logical(avMod);
scatter(targetWeightsAll_vis(ind),detectWeightsAll_vis(ind),'ro');
hold on
figXAxis([],'Target Weight',weightLim)
figYAxis([],'Detect Weight',weightLim)
figAxForm
title('Signifcantly Modulated Cells')

print(fullfile(fnout,'LRmodelXattn'),'-dpdf','-fillpage')

%%
setFigParams4Print('portrait')
set(0,'defaultAxesFontSize',12)
figure
suptitle('Visual Trials Model, Significantly Modulated Cells in Blue(target resp only) and Red(first stim resp only) and Magenta (both)')
ind1 = firstRespCells' & ~targetRespCells';
ind2 = ~firstRespCells' & targetRespCells';
ind3 = firstRespCells' & targetRespCells';
subplot 221
scatter(si,targetCorrsAll_vis,'ko')
hold on
scatter(si(ind1),targetCorrsAll_vis(ind1),'ro')
scatter(si(ind2),targetCorrsAll_vis(ind2),'bo')
scatter(si(ind3),targetCorrsAll_vis(ind3),'mo')
figXAxis([],'V-A Selectivity',siLim)
figYAxis([],'Correlation',corrLim)
figAxForm
title('Target')
subplot 222
scatter(si,detectCorrsAll_vis,'ko')
hold on
scatter(si(ind1),detectCorrsAll_vis(ind1),'ro')
scatter(si(ind2),detectCorrsAll_vis(ind2),'bo')
scatter(si(ind3),detectCorrsAll_vis(ind3),'mo')
figXAxis([],'V-A Selectivity',siLim)
figYAxis([],'Correlation',corrLim)
figAxForm
title('Detect')
subplot 223
scatter(si,targetWeightsAll_vis,'ko')
hold on
scatter(si(ind1),targetWeightsAll_vis(ind1),'ro')
scatter(si(ind2),targetWeightsAll_vis(ind2),'bo')
scatter(si(ind3),detectCorrsAll_vis(ind3),'mo')
figXAxis([],'V-A Selectivity',siLim)
figYAxis([],'Weight',weightLim)
figAxForm
title('Target')
subplot 224
scatter(si,detectWeightsAll_vis,'ko')
hold on
scatter(si(ind1),detectWeightsAll_vis(ind1),'ro')
scatter(si(ind2),detectWeightsAll_vis(ind2),'bo')
scatter(si(ind3),detectCorrsAll_vis(ind3),'mo')
figXAxis([],'V-A Selectivity',siLim)
figYAxis([],'Weight',weightLim)
figAxForm
title('Detect')

figure
scatter(targetWeightsAll_vis(ind1),detectWeightsAll_vis(ind1),'ro')
hold on
scatter(targetWeightsAll_vis(ind2),detectWeightsAll_vis(ind2),'bo')
scatter(targetWeightsAll_vis(ind3),detectCorrsAll_vis(ind3),'mo')
figXAxis([],'Target',weightLim)
figYAxis([],'Detect',weightLim)
figAxForm
title('LR Weights')
leg = legend({'first stim only';'target only';'both'});
title(leg,'Responds to:')

print(fullfile(fnout,'LRmodelXcellType'),'-dpdf')
%%

setFigParams4Print('portrait')
set(0,'defaultAxesFontSize',12)
figure
suptitle('Visual Trials Model')
% ind1 = avMod & firstRespCells' & ~targetRespCells';
% ind2 = avMod & ~firstRespCells' & targetRespCells';
% ind3 = avMod & firstRespCells' & targetRespCells';
ind1 = logical(avMod');
% ind1 = avMod' & respCells_targetOnly & ~respCells_lateBaseOnly;
% ind2 = avMod' & (respCells_firstBase | respCells_lateBaseOnly);
% ind3 = avMod' & (respCells_targetAndBase |...
%     (respCells_lateBaseOnly & respCells_targetOnly));
% ind4 = avMod' & suppCells & ~(ind1 | ind2 | ind3);
subplot 221
h=plot(si,targetCorrsAll_vis,'o');
h.Color = [0.5 0.5 0.5];
h.MarkerFaceColor = [0.5 0.5 0.5];
hold on
h=plot(si(ind1),targetCorrsAll_vis(ind1),'ko');
h.MarkerFaceColor = 'k';
% h=plot(si(ind2),targetCorrsAll_vis(ind2),'bo')
% h=plot(si(ind3),targetCorrsAll_vis(ind3),'mo')
% scatter(si(ind4),targetCorrsAll_vis(ind4),'go')
figXAxis([],'V-A Selectivity',siLim)
figYAxis([],'Correlation',corrLim)
figAxForm
title('Target')
subplot 222
h=plot(si,detectCorrsAll_vis,'o');
h.Color = [0.5 0.5 0.5];
h.MarkerFaceColor = [0.5 0.5 0.5];
hold on
h=plot(si(ind1),detectCorrsAll_vis(ind1),'ko');
h.MarkerFaceColor = 'k';
% scatter(si(ind2),detectCorrsAll_vis(ind2),'bo')
% scatter(si(ind3),detectCorrsAll_vis(ind3),'mo')
% scatter(si(ind4),detectCorrsAll_vis(ind4),'go')
figXAxis([],'V-A Selectivity',siLim)
figYAxis([],'Correlation',corrLim)
figAxForm
title('Detect')
subplot 223
h=plot(si,targetWeightsAll_vis,'o');
h.Color = [0.5 0.5 0.5];
h.MarkerFaceColor = [0.5 0.5 0.5];
hold on
h=plot(si(ind1),targetWeightsAll_vis(ind1),'ko');
h.MarkerFaceColor = 'k';
% scatter(si(ind2),targetWeightsAll_vis(ind2),'bo')
% scatter(si(ind3),targetWeightsAll_vis(ind3),'mo')
% scatter(si(ind4),targetWeightsAll_vis(ind4),'go')
figXAxis([],'V-A Selectivity',siLim)
figYAxis([],'Weight',weightLim)
figAxForm
title('Target')
subplot 224
h=plot(si,detectWeightsAll_vis,'ko');
h.Color = [0.5 0.5 0.5];
h.MarkerFaceColor = [0.5 0.5 0.5];
hold on
h=plot(si(ind1),detectWeightsAll_vis(ind1),'ko');
h.MarkerFaceColor = 'k';
% scatter(si(ind2),detectWeightsAll_vis(ind2),'bo')
% scatter(si(ind3),detectWeightsAll_vis(ind3),'mo')
% scatter(si(ind4),detectWeightsAll_vis(ind4),'go')
figXAxis([],'V-A Selectivity',siLim)
figYAxis([],'Weight',weightLim)
figAxForm
title('Detect')
% legend({'All';'Target Only';'Base Only';'Base & Target';'Base Supp Only'},...
%     'location','northeast')
legend({'All';'Selective Only'},...
    'location','northeast')

print(fullfile(fnout,'LRmodelXattn2'),'-dpdf','-fillpage')


figure
suptitle('Visual Trials Model, Significantly Modulated Cells in Blue(target resp only) and Red(first stim resp only) and Magenta (both)')
subplot 221
x = si(ind1);
y = targetCorrsAll_vis(ind1);
errorbar(mean(x), mean(y), ste(y,1), ste(y,1), ste(x,2),ste(x,2),'ro');
hold on
x = si(ind2);
y =targetCorrsAll_vis(ind2); 
errorbar(mean(x), mean(y), ste(y,1), ste(y,1), ste(x,2),ste(x,2),'bo');
x = si(ind3);
y = targetCorrsAll_vis(ind3);
errorbar(mean(x), mean(y), ste(y,1), ste(y,1), ste(x,2),ste(x,2),'mo');
figXAxis([],'V-A Selectivity',siLim)
figYAxis([],'Correlation',corrLim)
figAxForm
title('Target')
subplot 222
x = si(ind1);
y = detectCorrsAll_vis(ind1);
errorbar(mean(x), mean(y), ste(y,1), ste(y,1), ste(x,2),ste(x,2),'ro');
hold on
x = si(ind2);
y =detectCorrsAll_vis(ind2); 
errorbar(mean(x), mean(y), ste(y,1), ste(y,1), ste(x,2),ste(x,2),'bo');
x = si(ind3);
y = detectCorrsAll_vis(ind3);
errorbar(mean(x), mean(y), ste(y,1), ste(y,1), ste(x,2),ste(x,2),'mo');
figXAxis([],'V-A Selectivity',siLim)
figYAxis([],'Correlation',corrLim)
figAxForm
title('Detect')
subplot 223
x = si(ind1);
y = targetWeightsAll_vis(ind1);
errorbar(mean(x), mean(y), ste(y,1), ste(y,1), ste(x,2),ste(x,2),'ro');
hold on
x = si(ind2);
y =targetWeightsAll_vis(ind2); 
errorbar(mean(x), mean(y), ste(y,1), ste(y,1), ste(x,2),ste(x,2),'bo');
x = si(ind3);
y = targetWeightsAll_vis(ind3);
errorbar(mean(x), mean(y), ste(y,1), ste(y,1), ste(x,2),ste(x,2),'mo');
figXAxis([],'V-A Selectivity',siLim)
figYAxis([],'Weight',weightLim)
figAxForm
title('Target')
subplot 224
x = si(ind1);
y = detectWeightsAll_vis(ind1);
errorbar(mean(x), mean(y), ste(y,1), ste(y,1), ste(x,2),ste(x,2),'ro');
hold on
x = si(ind2);
y =detectWeightsAll_vis(ind2); 
errorbar(mean(x), mean(y), ste(y,1), ste(y,1), ste(x,2),ste(x,2),'bo');
x = si(ind3);
y = detectWeightsAll_vis(ind3);
errorbar(mean(x), mean(y), ste(y,1), ste(y,1), ste(x,2),ste(x,2),'mo');
figXAxis([],'V-A Selectivity',siLim)
figYAxis([],'Weight',weightLim)
figAxForm
title('Detect')

print(fullfile(fnout,'LRmodelXattn2_summary'),'-dpdf','-fillpage')

%%
setFigParams4Print('landscape')
set(0,'defaultAxesFontSize',10)
figure
suptitle('Visual Trials Model, % Adapt = mean(last 4 cycles)/first')
subplot 241
scatter(pctVisAdapt,targetCorrsAll_vis,'ko')
figXAxis([],'% Adapt - Visual',adaptLim)
figYAxis([],'Correlation',corrLim)
figAxForm
title('Target')
subplot 242
scatter(pctVisAdapt,detectCorrsAll_vis,'ko')
figXAxis([],'% Adapt - Visual',adaptLim)
figYAxis([],'Correlation',corrLim)
figAxForm
title('Detect')
subplot 243
scatter(pctVisAdapt,targetWeightsAll_vis,'ko')
figXAxis([],'% Adapt - Visual',adaptLim)
figYAxis([],'Weight',weightLim)
figAxForm
title('Target')
subplot 244
scatter(pctVisAdapt,detectWeightsAll_vis,'ko')
figXAxis([],'% Adapt - Visual',adaptLim)
figYAxis([],'Weight',weightLim)
figAxForm
title('Detect')

subplot 245
scatter(pctAudAdapt,targetCorrsAll_vis,'ko')
figXAxis([],'% Adapt - Auditory',adaptLim)
figYAxis([],'Correlation',corrLim)
figAxForm
title('Target')
subplot 246
scatter(pctAudAdapt,detectCorrsAll_vis,'ko')
figXAxis([],'% Adapt - Auditory',adaptLim)
figYAxis([],'Correlation',corrLim)
figAxForm
title('Detect')
subplot 247
scatter(pctAudAdapt,targetWeightsAll_vis,'ko')
figXAxis([],'% Adapt - Auditory',adaptLim)
figYAxis([],'Weight',weightLim)
figAxForm
title('Target')
subplot 248
scatter(pctAudAdapt,detectWeightsAll_vis,'ko')
figXAxis([],'% Adapt - Auditory',adaptLim)
figYAxis([],'Weight',weightLim)
figAxForm
title('Detect')

print(fullfile(fnout,'LRmodelXadaptation'),'-dpdf','-fillpage')

figure;
subplot 131
ind = ~isnan(pctVisAdapt) & ~isnan(pctAudAdapt);
a = pctVisAdapt(ind);
b = pctAudAdapt(ind);
abar = mean(a);
aerr = ste(a,2);
bbar = mean(b);
berr = ste(b,2);
h = bar(1,abar);
h.BarWidth = 1;
h.EdgeColor = 'none';
h.FaceColor = avcolmat{1};
hold on
errorbar(1,abar,aerr,avcolmat{1});
h = bar(2,bbar);
h.BarWidth = 1;
h.EdgeColor = 'none';
h.FaceColor = avcolmat{2};
errorbar(2,bbar,berr,avcolmat{2});
hold on
figXAxis([],'',[0 3],1:2,{'Visual','Auditory'})
figYAxis([],'% Adapt',[0 0.5])
figAxForm([],0)

subplot 132
scatter(a,b,'ko')
hold on
plot(adaptLim,adaptLim,'k--')
errorbar(mean(a),mean(b),ste(b,2),ste(b,2),ste(a,2),ste(a,2),'ro')
figXAxis([],'Visual',adaptLim)
figYAxis([],'Auditory',adaptLim)
figAxForm
title('% Adapt')

subplot 133
ind = pctVisAdapt < pctAudAdapt;
scatter(targetWeightsAll_vis(ind),detectWeightsAll_vis(ind),'bo');
hold on
ind = pctVisAdapt > pctAudAdapt;
scatter(targetWeightsAll_vis(ind),detectWeightsAll_vis(ind),'ro');
figXAxis([],'Target Weight',weightLim)
figYAxis([],'Detect Weight',weightLim)
figAxForm
legend({'% Vis < % Aud';'% Vis > % Aud'},'location','northwest')

print(fullfile(fnout,'LRmodelXadaptation2'),'-dpdf','-fillpage')
%%
setFigParams4Print('landscape')
set(0,'defaultAxesFontSize',16)

figure
subplot 121
x = si;
y = targetWeights_sub;
ybar = mean(y);
yerr = ste(y,2);
scatter(x,y,'ko')
hold on
ind = avMod' & firstRespCells & ~targetRespCells;
scatter(x(ind),y(ind),'ro')
ind = avMod' & ~firstRespCells & targetRespCells;
scatter(x(ind),y(ind),'bo')
ind = avMod' & firstRespCells & targetRespCells;
scatter(x(ind),y(ind),'mo')

figXAxis([],'Selectivity',siLim)
figYAxis([],'B_V - B_A',weightSubLim)
figAxForm
hline(0,'k--')
vline(0,'k--')
title('Target')
subplot 122
y = detectWeights_sub;
ind = logical(avMod);
ybar = mean(y);
yerr = ste(y,2);
scatter(x,y,'ko')
hold on
ind = avMod' & firstRespCells & ~targetRespCells;
scatter(x(ind),y(ind),'ro')
ind = avMod' & ~firstRespCells & targetRespCells;
scatter(x(ind),y(ind),'bo')
ind = avMod' & firstRespCells & targetRespCells;
scatter(x(ind),y(ind),'mo')
figXAxis([],'Selectivity',siLim)
figYAxis([],'B_V - B_A',weightSubLim)
figAxForm
hline(0,'k--')
vline(0,'k--')
title('Detect')

%% catch trial correlations
setFigParams4Print('portrait')
set(0,'defaultAxesFontSize',12)

figure
ind1 = avMod == 1 & si > 0; 
ind2 = avMod == 1 & si < 0;
subplot 221
x = targetCorrsAll_vis;
y = targetCorrsAll_catch;
scatter(x,y,'ko');
hold on
scatter(x(ind1),y(ind1),'ro');
scatter(x(ind2),y(ind2),'bo');
plot(corrLim,corrLim,'k--')
figXAxis([],'Valid Visual',corrLim)
figYAxis([],'Invalid Visual',corrLim)
figAxForm
title('Target Correlation')
subplot 222
x = targetCorrsAll_aud;
y = targetCorrsAll_catch;
scatter(x,y,'ko');
hold on
scatter(x(ind1),y(ind1),'ro');
scatter(x(ind2),y(ind2),'bo');
plot(corrLim,corrLim,'k--')
figXAxis([],'Valid Auditory',corrLim)
figYAxis([],'Invalid Visual',corrLim)
figAxForm
title('Target Correlation')
subplot 223
x = detectCorrsAll_vis;
y = detectCorrsAll_catch;
scatter(x,y,'ko');
hold on
scatter(x(ind1),y(ind1),'ro');
scatter(x(ind2),y(ind2),'bo');
plot(corrLim,corrLim,'k--')
figXAxis([],'Valid Visual',corrLim)
figYAxis([],'Invalid Visual',corrLim)
figAxForm
title('Detect Correlation')
subplot 224
x = detectCorrsAll_aud;
y = detectCorrsAll_catch;
scatter(x,y,'ko');
hold on
scatter(x(ind1),y(ind1),'ro');
scatter(x(ind2),y(ind2),'bo');
plot(corrLim,corrLim,'k--')
figXAxis([],'Valid Auditory',corrLim)
figYAxis([],'Invalid Visual',corrLim)
figAxForm
title('Detect Correlation')

print(fullfile(fnout,'corrXattnAll'),'-dpdf','-fillpage')

figure
ind1 = avMod == 1 & si > 0;
ind2 = avMod == 1 & si < 0;
subplot 221
x = targetCorrsAll_vis;
y = targetCorrsAll_catch;
errorbar(nanmean(x),nanmean(y),ste(y,2),ste(y,2),ste(x,1),ste(x,1),'ko');
hold on
errorbar(nanmean(x(ind1)),nanmean(y(ind1)),ste(y(ind1),2),ste(y(ind1),2),ste(x(ind1),1),ste(x(ind1),1),'ro');
errorbar(nanmean(x(ind2)),nanmean(y(ind2)),ste(y(ind2),2),ste(y(ind2),2),ste(x(ind2),1),ste(x(ind2),1),'bo');
plot(corrLim,corrLim,'k--')
figXAxis([],'Valid Visual',corrLim)
figYAxis([],'Invalid Visual',corrLim)
figAxForm
title('Target Correlation')
subplot 222
x = targetCorrsAll_aud;
y = targetCorrsAll_catch;
errorbar(nanmean(x),nanmean(y),ste(y,2),ste(y,2),ste(x,1),ste(x,1),'ko');
hold on
errorbar(nanmean(x(ind1)),nanmean(y(ind1)),ste(y(ind1),2),ste(y(ind1),2),ste(x(ind1),1),ste(x(ind1),1),'ro');
errorbar(nanmean(x(ind2)),nanmean(y(ind2)),ste(y(ind2),2),ste(y(ind2),2),ste(x(ind2),1),ste(x(ind2),1),'bo');
plot(corrLim,corrLim,'k--')
figXAxis([],'Valid Auditory',corrLim)
figYAxis([],'Invalid Visual',corrLim)
figAxForm
title('Target Correlation')
subplot 223
x = detectCorrsAll_vis;
y = detectCorrsAll_catch;
errorbar(nanmean(x),nanmean(y),ste(y,2),ste(y,2),ste(x,1),ste(x,1),'ko');
hold on
errorbar(nanmean(x(ind1)),nanmean(y(ind1)),ste(y(ind1),2),ste(y(ind1),2),ste(x(ind1),1),ste(x(ind1),1),'ro');
errorbar(nanmean(x(ind2)),nanmean(y(ind2)),ste(y(ind2),2),ste(y(ind2),2),ste(x(ind2),1),ste(x(ind2),1),'bo');
plot(corrLim,corrLim,'k--')
figXAxis([],'Valid Visual',corrLim)
figYAxis([],'Invalid Visual',corrLim)
figAxForm
title('Detect Correlation')
subplot 224
x = detectCorrsAll_aud;
y = detectCorrsAll_catch;
errorbar(nanmean(x),nanmean(y),ste(y,2),ste(y,2),ste(x,1),ste(x,1),'ko');
hold on
errorbar(nanmean(x(ind1)),nanmean(y(ind1)),ste(y(ind1),2),ste(y(ind1),2),ste(x(ind1),1),ste(x(ind1),1),'ro');
errorbar(nanmean(x(ind2)),nanmean(y(ind2)),ste(y(ind2),2),ste(y(ind2),2),ste(x(ind2),1),ste(x(ind2),1),'bo');
plot(corrLim,corrLim,'k--')
figXAxis([],'Valid Auditory',corrLim)
figYAxis([],'Invalid Visual',corrLim)
figAxForm
title('Detect Correlation')
legend({'All';'SI > 0';'SI < 0'},'location','southwest')

print(fullfile(fnout,'corrXattnMean'),'-dpdf','-fillpage')

figure
suptitle('Visual Trials Model, Significantly Modulated Cells in Blue(target resp only) and Red(first stim resp only) and Magenta (both)')
ind1 = avMod & firstRespCells' & ~targetRespCells';
ind2 = avMod & ~firstRespCells' & targetRespCells';
ind3 = avMod & firstRespCells' & targetRespCells';
subplot 221
scatter(si,targetCorrsAll_catch,'ko')
hold on
scatter(si(ind1),targetCorrsAll_catch(ind1),'ro')
scatter(si(ind2),targetCorrsAll_catch(ind2),'bo')
scatter(si(ind3),targetCorrsAll_catch(ind3),'mo')
figXAxis([],'V-A Selectivity',siLim)
figYAxis([],'Correlation',corrLim)
figAxForm
title('Target')
subplot 222
scatter(si,detectCorrsAll_catch,'ko')
hold on
scatter(si(ind1),detectCorrsAll_catch(ind1),'ro')
scatter(si(ind2),detectCorrsAll_catch(ind2),'bo')
scatter(si(ind3),detectCorrsAll_catch(ind3),'mo')
figXAxis([],'V-A Selectivity',siLim)
figYAxis([],'Correlation',corrLim)
figAxForm
title('Detect')

print(fullfile(fnout,'InvCorrXselectiivty'),'-dpdf','-fillpage')
%% test model on catch trials
figure
subplot 121
ind = ~isnan(pctCorrectTarget_catch) & pctCorrectTarget_Xval_vis > 0;

for i = 1:nexp
    if ind(i)
        hold on
        x = 1:3;
        y = [pctCorrectTarget_Xval_vis(i) ...
            pctCorrectTarget_catch(i) pctCorrectTarget_aud(i)];
        h = plot(x,y,'-');
        h.Color = [0.5 0.5 0.5];
    end
end
n = sum(ind);
x = ones(n,1);
y = pctCorrectTarget_Xval_vis(ind);
errorbar(mean(x),mean(y),ste(y,2),'ko')
x = ones(n,1).*2;
y = pctCorrectTarget_catch(ind);
errorbar(mean(x),mean(y),ste(y,2),'ko')
x = ones(n,1).*3;
y = pctCorrectTarget_aud(ind);
errorbar(mean(x),mean(y),ste(y,2),'ko')
figXAxis([],'Test',[0 4],1:3,{'Val Vis HO';'Inv Vis';'Aud'})
figYAxis([],'% Correct',[0 1])
figAxForm
title('Target')

subplot 122
ind = ~isnan(pctCorrectDetect_catch) & pctCorrectDetect_Xval_vis > 0;
for i = 1:nexp
    if ind(i)
        hold on
        x = 1:3;
        y = [pctCorrectDetect_Xval_vis(i) ...
            pctCorrectDetect_catch(i) pctCorrectDetect_aud(i)];
        h = plot(x,y,'-');
        h.Color = [0.5 0.5 0.5];
    end
end
n = sum(ind);
x = ones(n,1);
y = pctCorrectDetect_Xval_vis(ind);
errorbar(mean(x),mean(y),ste(y,2),'ko')
x = ones(n,1).*2;
y = pctCorrectDetect_catch(ind);
errorbar(mean(x),mean(y),ste(y,2),'ko')
x = ones(n,1).*3;
y = pctCorrectDetect_aud(ind);
errorbar(mean(x),mean(y),ste(y,2),'ko')
figXAxis([],'',[0 4],1:3,{'Val Vis HO';'Inv Vis';'Aud'})
figYAxis([],'% Correct',[0 1])
figAxForm
title('Detect')

print(fullfile(fnout,'modelTestInvalid'),'-dpdf','-fillpage')

%%
setFigParams4Print('portrait')
set(0,'defaultAxesFontSize',12)
figure
suptitle('Visual Trials Model')
ind1 = respCells_targetOnly;
ind2 = (respCells_firstBase | respCells_targetAndBase);
ind3 = respCells_lateBaseOnly;
ind4 = suppCells;
% ind1 = avMod' & respCells_targetOnly;
% ind2 = avMod' & (respCells_baseOnly | respCells_targetAndBase);
% ind3 = avMod' & (respCells_lateBase & ~ind2);
% ind4 = avMod' & suppCells;
subplot 221
scatter(si,targetCorrsAll_vis,'ko')
hold on
scatter(si(ind1),targetCorrsAll_vis(ind1),'ro')
scatter(si(ind2),targetCorrsAll_vis(ind2),'bo')
scatter(si(ind3),targetCorrsAll_vis(ind3),'mo')
scatter(si(ind4),targetCorrsAll_vis(ind4),'go')
figXAxis([],'V-A Selectivity',siLim)
figYAxis([],'Correlation',corrLim)
figAxForm
title('Target')
legend({'All';'Target';'First Base';'Late Base';'Base Supp'})

subplot 222
scatter(si,detectCorrsAll_vis,'ko')
hold on
scatter(si(ind1),detectCorrsAll_vis(ind1),'ro')
scatter(si(ind2),detectCorrsAll_vis(ind2),'bo')
scatter(si(ind3),detectCorrsAll_vis(ind3),'mo')
scatter(si(ind4),detectCorrsAll_vis(ind4),'go')
figXAxis([],'V-A Selectivity',siLim)
figYAxis([],'Correlation',corrLim)
figAxForm
title('Detect')
subplot 223
scatter(si,targetWeightsAll_vis,'ko')
hold on
scatter(si(ind1),targetWeightsAll_vis(ind1),'ro')
scatter(si(ind2),targetWeightsAll_vis(ind2),'bo')
scatter(si(ind3),detectCorrsAll_vis(ind3),'mo')
scatter(si(ind4),detectCorrsAll_vis(ind4),'go')
figXAxis([],'V-A Selectivity',siLim)
figYAxis([],'Weight',weightLim)
figAxForm
title('Target')
subplot 224
scatter(si,detectWeightsAll_vis,'ko')
hold on
scatter(si(ind1),detectWeightsAll_vis(ind1),'ro')
scatter(si(ind2),detectWeightsAll_vis(ind2),'bo')
scatter(si(ind3),detectCorrsAll_vis(ind3),'mo')
scatter(si(ind4),detectCorrsAll_vis(ind4),'go')
figXAxis([],'V-A Selectivity',siLim)
figYAxis([],'Weight',weightLim)
figAxForm
title('Detect')

figure
subplot 121
scatter(targetWeightsAll_vis(ind1),detectWeightsAll_vis(ind1),'ro')
hold on
scatter(targetWeightsAll_vis(ind2),detectWeightsAll_vis(ind2),'bo')
scatter(targetWeightsAll_vis(ind3),detectWeightsAll_vis(ind3),'mo')
scatter(targetWeightsAll_vis(ind4),detectWeightsAll_vis(ind4),'go')
figXAxis([],'Target',weightLim)
figYAxis([],'Detect',weightLim)
figAxForm
title('LR Weights')
leg = legend({'Target';'First Base';'Late Base';'Base Supp'},...
    'location','northwest');
title(leg,'Responds to:')
subplot 122
x = targetWeightsAll_vis;
y = detectWeightsAll_vis;
errorbar(mean(x(ind1)),mean(y(ind1)),ste(y(ind1),1),ste(y(ind1),1),...
    ste(x(ind1),1),ste(x(ind1),1),'ro')
hold on
errorbar(mean(x(ind2)),mean(y(ind2)),ste(y(ind2),1),ste(y(ind2),1),...
    ste(x(ind2),1),ste(x(ind2),1),'bo')
errorbar(mean(x(ind3)),mean(y(ind3)),ste(y(ind3),1),ste(y(ind3),1),...
    ste(x(ind3),1),ste(x(ind3),1),'mo')
errorbar(mean(x(ind4)),mean(y(ind4)),ste(y(ind4),1),ste(y(ind4),1),...
    ste(x(ind4),1),ste(x(ind4),1),'go')
figXAxis([],'Target',weightLim)
figYAxis([],'Detect',weightLim)
figAxForm

% print(fullfile(fnout,'LRmodelXcellType2'),'-dpdf')