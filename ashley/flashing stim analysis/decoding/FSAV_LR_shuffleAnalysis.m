close all
clear all

load('Z:\Analysis\FSAV Choice\FSAV_decodeData_exampleExpts.mat')
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
pcv = cell(1,nexp);
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
   
   pcv{iexp} = cumsum(laten)./sum(latent);
   exptName{iexp} = [dcExpt(iexp).mouse '-' dcExpt(iexp).date];
   
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
figure
for iexp = 1:nexp
    h = plot(pcv{iexp});
    hold on
    
end
hline(0.95,'k--')
figXAxis([],'N Components',[])
figYAxis([],'% Var Expl.',[])
figAxForm
print(fullfile(fnout,'pca_pctVarExplained'),'-dpdf','-fillpage')
%%
setFigParams4Print('portrait')
set(0,'defaultAxesFontSize',12)
weightLim = [-2 3];
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
hline(0,'k:')
vline(0,'k:')
figXAxis([],'Target LR Weight',weightLim)
figYAxis([],'Target Correlation',corrLim)
figAxForm
subplot 325   
s = scatter(targetWeightsAll_vis,detectWeightsAll_vis,'bo');
s.MarkerFaceColor = [1 1 1];
hold on
hline(0,'k:')
vline(0,'k:')
figXAxis([],'Target LR Weight',weightLim)
figYAxis([],'Detect LR Weight',weightLim)
figAxForm

subplot 322
scatter(detectWeightsAll_shuf_vis,detectCorrsAll_shuf_vis,'bo')
hold on
hline(0,'k:')
vline(0,'k:')
figXAxis([],'Detect LR Weight',weightLim)
figYAxis([],'Detect Correlation',corrLim)
figAxForm
title('Shuffled Trials')
subplot 324
scatter(targetWeightsAll_shuf_vis,targetCorrsAll_shuf_vis,'bo')
hold on
hline(0,'k:')
vline(0,'k:')
figXAxis([],'Target LR Weight',weightLim)
figYAxis([],'Target Correlation',corrLim)
figAxForm
subplot 326   
s = scatter(targetWeightsAll_shuf_vis,detectWeightsAll_shuf_vis,'bo');
s.MarkerFaceColor = [1 1 1];
hold on
hline(0,'k:')
vline(0,'k:')
figXAxis([],'Target LR Weight',weightLim)
figYAxis([],'Detect LR Weight',weightLim)
figAxForm

print(fullfile(fnout,'corrXWeight_shuffle_vis'),'-dpdf','-fillpage')

figure
suptitle('Auditory Trials')
subplot 321
scatter(detectWeightsAll_aud,detectCorrsAll_aud,'bo')
hold on
hline(0,'k:')
vline(0,'k:')
figXAxis([],'Detect LR Weight',weightLim)
figYAxis([],'Detect Correlation',corrLim)
figAxForm
title('All Trials')
subplot 323
scatter(targetWeightsAll_aud,targetCorrsAll_aud,'bo')
hold on
hline(0,'k:')
vline(0,'k:')
figXAxis([],'Target LR Weight',weightLim)
figYAxis([],'Target Correlation',corrLim)
figAxForm
subplot 325   
s = scatter(targetWeightsAll_aud,detectWeightsAll_aud,'bo');
s.MarkerFaceColor = [1 1 1];
hold on
hline(0,'k:')
vline(0,'k:')
figXAxis([],'Target LR Weight',weightLim)
figYAxis([],'Detect LR Weight',weightLim)
figAxForm

subplot 322
scatter(detectWeightsAll_shuf_aud,detectCorrsAll_shuf_aud,'bo')
hold on
hline(0,'k:')
vline(0,'k:')
figXAxis([],'Detect LR Weight',weightLim)
figYAxis([],'Detect Correlation',corrLim)
figAxForm
title('Shuffled Trials')
subplot 324
scatter(targetWeightsAll_shuf_aud,targetCorrsAll_shuf_aud,'bo')
hold on
hline(0,'k:')
vline(0,'k:')
figXAxis([],'Target LR Weight',weightLim)
figYAxis([],'Target Correlation',corrLim)
figAxForm
subplot 326   
s = scatter(targetWeightsAll_shuf_aud,detectWeightsAll_shuf_aud,'bo');
s.MarkerFaceColor = [1 1 1];
hold on
hline(0,'k:')
vline(0,'k:')
figXAxis([],'Target LR Weight',weightLim)
figYAxis([],'Detect LR Weight',weightLim)
figAxForm

print(fullfile(fnout,'corrXWeight_shuffle_aud'),'-dpdf','-fillpage')