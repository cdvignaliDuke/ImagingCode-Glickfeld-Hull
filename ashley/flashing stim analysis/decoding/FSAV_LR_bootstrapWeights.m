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
weightLim = [-3 4];
weightSubLim = [-2.5 2.5];
corrLim = [-1 1];
siLim = [-12 12];
adaptLim = [-1.5 1.5];
maxCellN = 10;
nBoot = 1000;

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
tuned = [];
pctVisAdapt = [];
pctAudAdapt = [];

respCells_targetOnly = [];
respCells_baseOnly = [];
respCells_targetAndBase = [];
respCells_lateBase = [];
suppCells = [];

targetWeightsAll_vis_boot = [];
detectWeightsAll_vis_boot = [];

targetWeightsAll_aud_boot = [];
detectWeightsAll_aud_boot = [];

testsPerCell = [];

targetCorrsAll_vis = [];
detectCorrsAll_vis = [];
targetCorrsAll_aud = [];
detectCorrsAll_aud = [];
    

% nCells = zeros(1,nexp);
pcNumber = [];
exptName = cell(1,nexp);

for iexp=1:nexp
    
    if ~strcmp(dcExpt(iexp).mouse,attnInfoExpt(iexp).ms)
        error('Mouse name from attention data does not match decode data')
    end
    
   disp([dcExpt(iexp).mouse '-' dcExpt(iexp).date])
    %****VISUAL TRIALS
   trOut=dcExpt(iexp).trialOutcome;
    [Y_yeses_vis, Y_targets_vis] = getStimAndBehaviorYs(trOut);
   trOut=dcExpt(iexp).audTrialOutcome;
    [Y_yeses_aud, Y_targets_aud] = getStimAndBehaviorYs(trOut);
    
    
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
    
    X_targets = dcExpt(iexp).stimResp';
   X_targets=X_targets(matchTrialsInd,:);
   nCells = size(X_targets,2);
   
   Y_yeses_vis=Y_yeses_vis(matchTrialsInd);   
   Y_targets_vis=Y_targets_vis(matchTrialsInd); 

   detectCorr = corr(Y_yeses_vis,X_targets)';
   targetCorr = corr(Y_targets_vis,X_targets)';
    
    targetCorrsAll_vis = cat(1,targetCorrsAll_vis,targetCorr);
    detectCorrsAll_vis = cat(1,detectCorrsAll_vis,detectCorr);
    
    X_targets = dcExpt(iexp).audStimResp';
    detectCorr = corr(Y_yeses_aud,X_targets)';
   targetCorr = corr(Y_targets_aud,X_targets)';
    
    targetCorrsAll_aud = cat(1,targetCorrsAll_aud,targetCorr);
    detectCorrsAll_aud = cat(1,detectCorrsAll_aud,detectCorr);
    
   targetWeights_boot = cell(nBoot,nCells);
   detectWeights_boot = cell(nBoot,nCells);
   targetWeights_aud_boot = cell(nBoot,nCells);
   detectWeights_aud_boot = cell(nBoot,nCells);
    for iboot = 1:nBoot
       X_targets=dcExpt(iexp).stimResp';

       cellIdx = false(nCells,1);
       cellIdx(randperm(numel(cellIdx),maxCellN)) = true;

       X_targets=X_targets(matchTrialsInd,cellIdx);
       X_targets = zscore(X_targets);
       
       C=eye(size(X_targets,2));
       p=1;
       [temp1,dev1,targetGLM_vis]=glmfit(X_targets,Y_targets_vis,'binomial');
       [temp2,dev2,detectGLM_vis]=glmfit(X_targets,Y_yeses_vis,'binomial');
       idx=find(targetGLM_vis.p>p|detectGLM_vis.p>p);
       temp1(idx)=0;
       temp2(idx)=0;
       B_targets = C*temp1(2:end,1);
       B_yeses = C*temp2(2:end,1);
       
       targetWeights_boot(iboot,cellIdx) = num2cell(B_targets);
       detectWeights_boot(iboot,cellIdx) = num2cell(B_yeses);
       
       X_targets = dcExpt(iexp).audStimResp';
       X_targets=X_targets(:,cellIdx);
       X_targets = zscore(X_targets);
       
       C=eye(size(X_targets,2));
       p=1;
       [temp1,dev1,targetGLM_aud]=glmfit(X_targets,Y_targets_aud,'binomial');
       [temp2,dev2,detectGLM_aud]=glmfit(X_targets,Y_yeses_aud,'binomial');
       idx=find(targetGLM_aud.p>p|detectGLM_aud.p>p);
       temp1(idx)=0;
       temp2(idx)=0;
       B_targets = C*temp1(2:end,1);
       B_yeses = C*temp2(2:end,1);
       
       targetWeights_aud_boot(iboot,cellIdx) = num2cell(B_targets);
       detectWeights_aud_boot(iboot,cellIdx) = num2cell(B_yeses);
       
    end   
    testsPerCell = cat(2,testsPerCell,...
        sum(cellfun(@(x) ~isempty(x),targetWeights_boot),1));

    targetWeightsAll_vis_boot = cat(2,targetWeightsAll_vis_boot,...
        nanmean(celleqel2mat_padded(targetWeights_boot),1));
    detectWeightsAll_vis_boot = cat(2,detectWeightsAll_vis_boot,...
        nanmean(celleqel2mat_padded(detectWeights_boot),1));
    targetWeightsAll_aud_boot = cat(2,targetWeightsAll_aud_boot,...
        nanmean(celleqel2mat_padded(targetWeights_aud_boot),1));
    detectWeightsAll_aud_boot = cat(2,detectWeightsAll_aud_boot,...
        nanmean(celleqel2mat_padded(detectWeights_aud_boot),1));
    
   si = cat(2,si,attnInfoExpt(iexp).si);
   avMod = cat(2,avMod,attnInfoExpt(iexp).avModTest);
   tuned = cat(2,tuned,dcExpt(iexp).oriFitTheta90 < theta90Threshold);
   
   
   tc = dcExpt(iexp).targetResponsiveCells;
   fbrc = dcExpt(iexp).firstBaseResponsiveCells;
   
   respCells_targetOnly = cat(1,respCells_targetOnly,tc & ~fbrc);
   respCells_baseOnly = cat(1,respCells_baseOnly,~tc & fbrc);
   respCells_targetAndBase = cat(1,respCells_targetAndBase, tc & fbrc);
   respCells_lateBase = cat(1,respCells_lateBase,...
       attnInfoExpt(iexp).lateBaseRespCells');   
   suppCells = cat(1,suppCells,attnInfoExpt(iexp).lateBaseSuppCells');
   
   baseResponsiveCells = dcExpt(iexp).firstBaseResponsiveCells;
   ad_vis = attnInfoExpt(iexp).visAdapt;
   ad_aud = attnInfoExpt(iexp).audAdapt;
   ad_vis(~baseResponsiveCells) = nan;
   ad_aud(~baseResponsiveCells) = nan;
   
   pctVisAdapt = cat(2,pctVisAdapt,ad_vis);
   pctAudAdapt = cat(2,pctAudAdapt,ad_aud);
   
   
   
end
respCells_targetOnly = logical(respCells_targetOnly);
respCells_firstBase = logical(respCells_baseOnly);
respCells_targetAndBase = logical(respCells_targetAndBase);
respCells_lateBaseOnly = logical(respCells_lateBase &...
    ~respCells_baseOnly & ~respCells_targetAndBase);
suppCells = logical(suppCells);

%% boot strapped weights for visual or auditory trials
nbins = 100;
figure;
histogram(targetWeightsAll_vis_boot,nbins)
hold on
histogram(detectWeightsAll_vis_boot,nbins)

figure
suptitle('Visual Trials')
respCells_all = respCells_targetOnly | respCells_firstBase | respCells_targetAndBase;
subplot 221
plot(targetCorrsAll_vis(respCells_all),targetWeightsAll_vis_boot(respCells_all),'ko')
figXAxis([],'Correlation',corrLim)
figYAxis([],'Weight',weightLim)
figAxForm
title('Target')

subplot 222
plot(detectCorrsAll_vis(respCells_all),detectWeightsAll_vis_boot(respCells_all),'ko')
figXAxis([],'Correlation',corrLim)
figYAxis([],'Weight',weightLim)
figAxForm
title('Detect')

subplot 223
plot(targetWeightsAll_vis_boot(respCells_all),detectWeightsAll_vis_boot(respCells_all),'ko')
figXAxis([],'Target',weightLim)
figYAxis([],'Detect',weightLim)
figAxForm

subplot 224
ind = respCells_all' & targetWeightsAll_vis_boot < 0;
x = targetWeightsAll_vis_boot(ind);
y = detectWeightsAll_vis_boot(ind);
[~,p] = ttest(y);
plot(ones(sum(ind),1).*mean(x),y,'k.')
hold on
h = errorbar(mean(x),mean(y),ste(y,2),'ko');
if p < 0.025
    h.Color = 'r';
end
ind = respCells_all' & targetWeightsAll_vis_boot > 0;
x = targetWeightsAll_vis_boot(ind);
y = detectWeightsAll_vis_boot(ind);
[~,p] = ttest(y);
plot(ones(sum(ind),1).*mean(x),y,'k.')
hold on
h = errorbar(mean(x),mean(y),ste(y,2),'ko');
if p < 0.025
    h.Color = 'r';
end
figXAxis([],'Target',weightLim)
figYAxis([],'Detect',weightLim)
figAxForm
hline(0,'k--')
vline(0,'k--')
plot(weightLim,weightLim,'k--')

print(fullfile(fnout,'LRweight_vis_bootWeight'),'-dpdf','-fillpage')

figure
suptitle('Auditory Trials')
subplot 221
plot(targetCorrsAll_aud(respCells_all),targetWeightsAll_aud_boot(respCells_all),'ko')
figXAxis([],'Correlation',corrLim)
figYAxis([],'Weight',weightLim)
figAxForm
title('Target')

subplot 222
plot(detectCorrsAll_aud(respCells_all),detectWeightsAll_aud_boot(respCells_all),'ko')
figXAxis([],'Correlation',corrLim)
figYAxis([],'Weight',weightLim)
figAxForm
title('Detect')

subplot 223
plot(targetWeightsAll_aud_boot(respCells_all),detectWeightsAll_aud_boot(respCells_all),'ko')
figXAxis([],'Target',weightLim)
figYAxis([],'Detect',weightLim)
figAxForm

subplot 224
ind = respCells_all' & targetWeightsAll_aud_boot < 0;
x = targetWeightsAll_aud_boot(ind);
y = detectWeightsAll_aud_boot(ind);
[~,p] = ttest(y);
plot(ones(sum(ind),1).*mean(x),y,'k.')
hold on
h = errorbar(mean(x),mean(y),ste(y,2),'ko');
if p < 0.025
    h.Color = 'r';
end
ind = respCells_all' & targetWeightsAll_aud_boot > 0;
x = targetWeightsAll_aud_boot(ind);
y = detectWeightsAll_aud_boot(ind);
[~,p] = ttest(y);
plot(ones(sum(ind),1).*mean(x),y,'k.')
hold on
h = errorbar(mean(x),mean(y),ste(y,2),'ko');
if p < 0.025
    h.Color = 'r';
end
figXAxis([],'Target',weightLim)
figYAxis([],'Detect',weightLim)
figAxForm
hline(0,'k--')
vline(0,'k--')
plot(weightLim,weightLim,'k--')

print(fullfile(fnout,'LRweight_aud_bootWeight'),'-dpdf','-fillpage')
%% compare visual and auditory bootstrapped weights for each cell
figure
subplot 221
x = targetCorrsAll_vis;
y = targetCorrsAll_aud;
linFitCoeffs = polyfit(x,y,1);
yfit = polyval(linFitCoeffs,x);
Rsquared = 1 - (sum((y-yfit).^2)/sum((y - mean(y)).^2));
plot(x,y,'ko')
hold on
plot(x,yfit,'r-');
plot(corrLim,corrLim,'k--')
figXAxis([],'Visual',corrLim);
figYAxis([],'Auditory',corrLim);
figAxForm
title(sprintf('Target Correlation: R^2 = %s',num2str(Rsquared)))
subplot 222
x = detectCorrsAll_vis;
y = detectCorrsAll_aud;
linFitCoeffs = polyfit(x,y,1);
yfit = polyval(linFitCoeffs,x);
Rsquared = 1 - (sum((y-yfit).^2)/sum((y - mean(y)).^2));
plot(x,y,'ko')
hold on
plot(x,yfit,'r-');
plot(corrLim,corrLim,'k--')
figXAxis([],'Visual',corrLim);
figYAxis([],'Auditory',corrLim);
figAxForm
title(sprintf('Detect Correlation: R^2 = %s',num2str(Rsquared)))
subplot 223
x = targetWeightsAll_vis_boot;
y = targetWeightsAll_aud_boot;
linFitCoeffs = polyfit(x,y,1);
yfit = polyval(linFitCoeffs,x);
Rsquared = 1 - (sum((y-yfit).^2)/sum((y - mean(y)).^2));
plot(x,y,'ko')
hold on
plot(x,yfit,'r-');
plot(weightLim,weightLim,'k--')
figXAxis([],'Visual',weightLim);
figYAxis([],'Auditory',weightLim);
figAxForm
title(sprintf('Target Weight: R^2 = %s',num2str(Rsquared)));
subplot 224
x = detectWeightsAll_vis_boot;
y = detectWeightsAll_aud_boot;
linFitCoeffs = polyfit(x,y,1);
yfit = polyval(linFitCoeffs,x);
Rsquared = 1 - (sum((y-yfit).^2)/sum((y - mean(y)).^2));
plot(x,y,'ko')
hold on
plot(x,yfit,'r-');
plot(weightLim,weightLim,'k--')
figXAxis([],'Visual',weightLim);
figYAxis([],'Auditory',weightLim);
figAxForm
title(sprintf('Detect Weight: R^2 = %s',num2str(Rsquared)));

print(fullfile(fnout,'VisVsAud_bootWeight'),'-dpdf','-fillpage')
%% attention vs bootstrapped weights
figure
suptitle('Visual Trials, task-responsive neurons')
ind0 = respCells_targetOnly | respCells_targetAndBase | ...
    respCells_firstBase | respCells_lateBaseOnly | suppCells;
ind1 = respCells_targetOnly & avMod';
ind2 = respCells_targetAndBase & avMod';
ind3 = (respCells_firstBase | respCells_lateBaseOnly) & avMod';
% ind3 = (respCells_firstBase) & avMod';
ind4 = suppCells & ~(ind1 | ind2 | ind3) & avMod';
x = si;
subplot 221
y = targetCorrsAll_vis;
plot(x,y,'k.')
hold on
plot(x(ind1),y(ind1),'ro');
plot(x(ind2),y(ind2),'mo');
plot(x(ind3),y(ind3),'bo');
plot(x(ind4),y(ind4),'go');
figXAxis([],'V-A Selectivity',siLim);
figYAxis([],'Correlation',corrLim);
figAxForm
title('Target')
subplot 222
y = detectCorrsAll_vis;
plot(x,y,'k.')
hold on
plot(x(ind1),y(ind1),'ro');
plot(x(ind2),y(ind2),'mo');
plot(x(ind3),y(ind3),'bo');
plot(x(ind4),y(ind4),'go');
figXAxis([],'V-A Selectivity',siLim);
figYAxis([],'Correlation',corrLim);
figAxForm
title('Detect')
subplot 223
y = targetWeightsAll_vis_boot;
plot(x,y,'k.')
hold on
plot(x(ind1),y(ind1),'ro');
plot(x(ind2),y(ind2),'mo');
plot(x(ind3),y(ind3),'bo');
plot(x(ind4),y(ind4),'go');
figXAxis([],'V-A Selectivity',siLim);
figYAxis([],'Weight',weightLim);
figAxForm
title('Target')
subplot 224
y = detectWeightsAll_vis_boot;
plot(x,y,'k.')
hold on
plot(x(ind1),y(ind1),'ro');
plot(x(ind2),y(ind2),'mo');
plot(x(ind3),y(ind3),'bo');
plot(x(ind4),y(ind4),'go');
figXAxis([],'V-A Selectivity',siLim);
figYAxis([],'Weight',weightLim);
figAxForm
title('Detect')
legend({'All';'Target Only';'Base & Target';'Base Only';'Base Supp Only'},...
    'location','northeast')


print(fullfile(fnout,'LRmodelXattn_bootWeight'),'-dpdf','-fillpage')

%%
