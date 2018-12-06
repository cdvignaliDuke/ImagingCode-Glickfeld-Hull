close all
clear all

load('X:\home\ashley\Analysis\FSAV Choice\FSAV_decodeData.mat')
load('X:\home\ashley\Analysis\FSAV Summaries\FSAV_attentionV1\attentionV1_startAlign_FSAV_attnData')
fnout = 'X:\home\ashley\Analysis\FSAV Choice';
nexp = length(dcExpt);

doExptPlots = 0;
%%
minTrN = 17;
theta90Threshold = 11.25;
decisionVariable = 0.5;
detectThresholdHR = 0.8;
detectThresholdOri = 45;
oriBinsFine = [0 1 16 22.5 32 45 64 90];
oriBins = [0 1 32 90];
% oriPrefBins = [0 22.5:22.5:157.5 181];
oriPrefBins = [0 22.5:45:157.5 181];
ampBins = [0 0.0001 0.1 1];
weightsDiffZeroAlpha = 0.05/2;
%%
targetWeightsAll_vis = [];
detectWeightsAll_vis = [];

targetCorrsAll_vis = [];
detectCorrsAll_vis = [];

targetWeightsAll_aud = [];
detectWeightsAll_aud = [];

targetCorrsAll_aud = [];
detectCorrsAll_aud = [];

respCells_targetOnly = [];
respCells_baseOnly = [];
respCells_targetAndBase = [];
respCells_lateBase = [];
suppCells = [];

nCells = zeros(1,nexp);
oriPrefAll = [];
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
   
   tc = dcExpt(iexp).targetResponsiveCells(cellIdx);
   fbrc = dcExpt(iexp).firstBaseResponsiveCells(cellIdx);
   
   respCells_targetOnly = cat(1,respCells_targetOnly,tc & ~fbrc);
   respCells_baseOnly = cat(1,respCells_baseOnly,~tc & fbrc);
   respCells_targetAndBase = cat(1,respCells_targetAndBase, tc & fbrc);
   respCells_lateBase = cat(1,respCells_lateBase,...
       attnInfoExpt(iexp).lateBaseRespCells(cellIdx)');   
   suppCells = cat(1,suppCells,attnInfoExpt(iexp).lateBaseSuppCells(cellIdx)');
   
   X_targets=X_targets(matchTrialsInd,cellIdx);
   
   Y_yeses=Y_yeses(matchTrialsInd);   
   Y_targets=Y_targets(matchTrialsInd); 
   
   detectCorr = corr(Y_yeses,X_targets)';
   targetCorr = corr(Y_targets,X_targets)';
   
   detectCorrsAll_vis = cat(1,detectCorrsAll_vis,detectCorr);
   targetCorrsAll_vis = cat(1,targetCorrsAll_vis,targetCorr);
  
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
 
   targetWeightsAll_vis = cat(1,targetWeightsAll_vis,B_targets);
   detectWeightsAll_vis = cat(1,detectWeightsAll_vis,B_yeses);

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
end
respCells_targetOnly = logical(respCells_targetOnly);
respCells_firstBase = logical(respCells_baseOnly);
respCells_targetAndBase = logical(respCells_targetAndBase);
respCells_lateBaseOnly = logical(respCells_lateBase &...
    ~respCells_baseOnly & ~respCells_targetAndBase);
suppCells = logical(suppCells);
%%
[n,~,oriPrefBinID] = histcounts(oriPrefAll,oriPrefBins);
oriPrefBinID(oriPrefBinID == max(oriPrefBinID)) = 1;
n(1) = n(1)+n(end);
n = n(1:end-1);
maxBin = max(oriPrefBinID);

detectWeightOriBin_vis = nan(1,maxBin);
detectWeightOriBinErr_vis = nan(1,maxBin);
targetWeightOriBin_vis = nan(1,maxBin);
targetWeightOriBinErr_vis = nan(1,maxBin);
detectWeightOriBin_aud = nan(1,maxBin);
detectWeightOriBinErr_aud = nan(1,maxBin);
targetWeightOriBin_aud = nan(1,maxBin);
targetWeightOriBinErr_aud = nan(1,maxBin);
for i = 1:maxBin
    ind = oriPrefBinID == i;
    
    detectWeightOriBin_vis(i) = mean(detectWeightsAll_vis(ind));
    detectWeightOriBinErr_vis(i) = ste(detectWeightsAll_vis(ind),1);
    targetWeightOriBin_vis(i) = mean(targetWeightsAll_vis(ind));
    targetWeightOriBinErr_vis(i) = ste(targetWeightsAll_vis(ind),1);
    
    detectWeightOriBin_aud(i) = mean(detectWeightsAll_aud(ind));
    detectWeightOriBinErr_aud(i) = ste(detectWeightsAll_aud(ind),1);
    targetWeightOriBin_aud(i) = mean(targetWeightsAll_aud(ind));
    targetWeightOriBinErr_aud(i) = ste(targetWeightsAll_aud(ind),1);
end

[p_anova_detect_vis,~,stats_detect_vis] = anova1(detectWeightsAll_vis,oriPrefBinID);
[p_anova_detect_aud,~,stats_detect_aud] = anova1(detectWeightsAll_aud,oriPrefBinID);
[p_anova_target_vis,~,stats_target_vis] = anova1(targetWeightsAll_vis,oriPrefBinID);
[p_anova_target_aud,~,stats_target_aud] = anova1(targetWeightsAll_aud,oriPrefBinID);

%%
weightLim = [-3 4];
audWeightLim = [-2 2];
polarWeightLim = [0 4];
polarAudWeightLim = [0 2];

%%

setFigParams4Print('portrait')
set(0,'defaultAxesFontSize',12)

figure
suptitle('Visual Trials')
subplot 321
bar(1:maxBin,n)
figXAxis([],'Ori Pref Bin',[0 maxBin+1],[1:maxBin],oriPrefBins(2:end-1))
figYAxis([],'N Cells',[]);
figAxForm
subplot 322
errorbar(1:maxBin,targetWeightOriBin_vis,targetWeightOriBinErr_vis,'ko-');
hold on
errorbar(1:maxBin,detectWeightOriBin_vis,detectWeightOriBinErr_vis,'bo-');
figXAxis([],'Ori Pref Bin',[0 maxBin+1],[1:maxBin],oriPrefBins(2:end-1))
figYAxis([],'Weight',[]);
figAxForm
legendStr = strcat({'Target ';'Detect '},...
    {num2str(round(p_anova_target_vis,2,'significant'));...
    num2str(round(p_anova_detect_vis,2,'significant'))});
legend(legendStr,'location','northeast')

subplot 323
h = plot(oriPrefAll,targetWeightsAll_vis,'o');
h.Color = 'k';
figXAxis([],'Orientation Pref (deg)',[-10 190],0:45:180,0:45:180)
figYAxis([],'Weight',weightLim)
figAxForm
title('Target')
subplot 324
h = plot(oriPrefAll,detectWeightsAll_vis,'o');
h.Color = 'k';
figXAxis([],'Orientation Pref (deg)',[-10 190],0:45:180,0:45:180)
figYAxis([],'Weight',weightLim)
figAxForm
title('Detect')

subplot 325
h = polarplot(deg2rad(oriPrefAll),targetWeightsAll_vis','bo');
h.Parent.ThetaTick = [0:45:360];
h.Parent.RLim = polarWeightLim;
h.Parent.RTick = [0:2:4];
title('Target Weight')
subplot 326
h = polarplot(deg2rad(oriPrefAll),detectWeightsAll_vis','bo');
h.Parent.ThetaTick = [0:45:360];
h.Parent.RLim = polarWeightLim;
h.Parent.RTick = [0:2:4];
title('Detect Weight')

print(fullfile(fnout,'tuningXweight_matched_vis'),'-dpdf','-fillpage')

%%
respCellType = {respCells_targetOnly;respCells_firstBase;respCells_targetAndBase};
respCellColor = {'r';'b';'m'};
respCellName = {'Target Only'; 'Base Only'; 'Both'};

detectWeightOriBin_vis_respCells = cell(1,3);
targetWeightOriBin_vis_respCells = cell(1,3);
detectWeightOriBin_aud_respCells = cell(1,3);
targetWeightOriBin_aud_respCells = cell(1,3);

detectWeightOriBin_vis_respCells(:) = {nan(1,maxBin)};
targetWeightOriBin_vis_respCells(:) = {nan(1,maxBin)};
detectWeightOriBin_aud_respCells(:) = {nan(1,maxBin)};
targetWeightOriBin_aud_respCells(:) = {nan(1,maxBin)};
for icell = 1:3
    respCellInd = respCellType{icell}';
    for i = 1:maxBin
        ind = oriPrefBinID == i & respCellInd;

        detectWeightOriBin_vis_respCells{icell}(i) = mean(detectWeightsAll_vis(ind));
        targetWeightOriBin_vis_respCells{icell}(i) = mean(targetWeightsAll_vis(ind));

        detectWeightOriBin_aud_respCells{icell}(i) = mean(detectWeightsAll_aud(ind));
        targetWeightOriBin_aud_respCells{icell}(i) = mean(targetWeightsAll_aud(ind));
    end
end


%%
figure

subplot 221
h = plot(oriPrefAll,targetWeightsAll_vis,'o');
h.Color = 'k';
hold on
h = plot(oriPrefAll(respCells_firstBase),...
    targetWeightsAll_vis(respCells_firstBase),'o');
h.Color = 'b';
h = plot(oriPrefAll(respCells_targetAndBase),...
    targetWeightsAll_vis(respCells_targetAndBase),'o');
h.Color = 'm';
h = plot(oriPrefAll(respCells_targetOnly),...
    targetWeightsAll_vis(respCells_targetOnly),'o');
h.Color = 'r';
figXAxis([],'Orientation Pref (deg)',[-10 190],0:45:180,0:45:180)
figYAxis([],'Weight',weightLim)
figAxForm
title('Target')

subplot 222
h = plot(oriPrefAll,detectWeightsAll_vis,'o');
h.Color = 'k';
hold on
h = plot(oriPrefAll(respCells_firstBase),...
    detectWeightsAll_vis(respCells_firstBase),'o');
h.Color = 'b';
h = plot(oriPrefAll(respCells_targetAndBase),...
    detectWeightsAll_vis(respCells_targetAndBase),'o');
h.Color = 'm';
h = plot(oriPrefAll(respCells_targetOnly),...
    detectWeightsAll_vis(respCells_targetOnly),'o');
h.Color = 'r';
figXAxis([],'Orientation Pref (deg)',[-10 190],0:45:180,0:45:180)
figYAxis([],'Weight',weightLim)
figAxForm
title('Detect')

subplot 223
for i = 1:3
    x = 1:maxBin;
    y = targetWeightOriBin_vis_respCells{i};
    hold on
    h = plot(x,y,'o-');
    h.Color = respCellColor{i};
end
hline(0,'k--')
figXAxis([],'Orientation Pref (deg)',[0 maxBin+1],1:maxBin,0:45:180)
figYAxis([],'Weight',weightLim)
figAxForm
title('Target')
subplot 224
for i = 1:3
    x = 1:maxBin;
    y = detectWeightOriBin_vis_respCells{i};
    hold on
    h = plot(x,y,'o-');
    h.Color = respCellColor{i};
end
hline(0,'k--')
figXAxis([],'Orientation Pref (deg)',[0 maxBin+1],1:maxBin,0:45:180)
figYAxis([],'Weight',weightLim)
figAxForm
title('Target')
legend(respCellName)

print(fullfile(fnout,'tuningXweight_respCellTypes_vis'),'-dpdf','-fillpage')


%%

figure
suptitle('Auditory Trials')
subplot 321
bar(1:maxBin,n)
figXAxis([],'Ori Pref Bin',[0 maxBin+1],[1:maxBin],oriPrefBins(1:end-1))
figYAxis([],'N Cells',[]);
figAxForm
subplot 322
errorbar(1:maxBin,targetWeightOriBin_aud,targetWeightOriBinErr_aud,'ko-');
hold on
errorbar(1:maxBin,detectWeightOriBin_aud,detectWeightOriBinErr_aud,'bo-');
figXAxis([],'Ori Pref Bin',[0 maxBin+1],[1:maxBin],oriPrefBins(1:end-1))
figYAxis([],'Weight',[]);
figAxForm
legendStr = strcat({'Target ';'Detect '},...
    {num2str(round(p_anova_target_aud,2,'significant'));...
    num2str(round(p_anova_detect_aud,2,'significant'))});
legend(legendStr,'location','northeast')

subplot 323
h = plot(oriPrefAll,targetWeightsAll_aud,'o');
h.Color = 'k';
figXAxis([],'Orientation Pref (deg)',[-10 190],0:45:180,0:45:180)
figYAxis([],'Weight',audWeightLim)
figAxForm
title('Target')
subplot 324
h = plot(oriPrefAll,detectWeightsAll_aud,'o');
h.Color = 'k';
figXAxis([],'Orientation Pref (deg)',[-10 190],0:45:180,0:45:180)
figYAxis([],'Weight',audWeightLim)
figAxForm
title('Detect')

subplot 325
h = polarplot(deg2rad(oriPrefAll),targetWeightsAll_aud','bo');
h.Parent.ThetaTick = [0:45:360];
h.Parent.RLim = polarAudWeightLim;
h.Parent.RTick = [0:2:4];
title('Target Weight')
subplot 326
h = polarplot(deg2rad(oriPrefAll),detectWeightsAll_aud','bo');
h.Parent.ThetaTick = [0:45:360];
h.Parent.RLim = polarAudWeightLim;
h.Parent.RTick = [0:2:4];
title('Detect Weight')

print(fullfile(fnout,'tuningXweight_matched_aud'),'-dpdf','-fillpage')

%%
figure
subplot 121
errorbar(1:maxBin,targetWeightOriBin_vis,targetWeightOriBinErr_vis,'ko-');
hold on
errorbar(1:maxBin,targetWeightOriBin_aud,targetWeightOriBinErr_aud,'bo-');
figXAxis([],'Ori Pref Bin',[0 maxBin+1],[1:maxBin],oriPrefBins(2:end-1))
figYAxis([],'Weight',[]);
figAxForm
title('Target')

subplot 122
errorbar(1:maxBin,detectWeightOriBin_vis,detectWeightOriBinErr_vis,'ko-');
hold on
errorbar(1:maxBin,detectWeightOriBin_aud,detectWeightOriBinErr_aud,'bo-');
figXAxis([],'Ori Pref Bin',[0 maxBin+1],[1:maxBin],oriPrefBins(2:end-1))
figYAxis([],'Weight',[]);
figAxForm
legend({'Visual';'Auditory'},'location','northeast')
title('Detect')

print(fullfile(fnout,'tuningXweight_matched_compareVisAud'),'-dpdf','-fillpage')

%%
detectWeightOriBin_vis_respCells = cell(1,3);
targetWeightOriBin_vis_respCells = cell(1,3);
detectWeightOriBin_aud_respCells = cell(1,3);
targetWeightOriBin_aud_respCells = cell(1,3);

detectWeightOriBin_vis_respCells(:) = {nan(1,maxBin)};
targetWeightOriBin_vis_respCells(:) = {nan(1,maxBin)};
detectWeightOriBin_aud_respCells(:) = {nan(1,maxBin)};
detectWeightOriBin_aud_respCells(:) = {nan(1,maxBin)};
for icell = 1:3
    respCellInd = respCellType{icell}';
    for i = 1:maxBin
        ind = oriPrefBinID == i & respCellInd;

        detectWeightOriBin_vis_respCells{icell}(i) = mean(detectWeightsAll_vis(ind));
        targetWeightOriBin_vis_respCells{icell}(i) = mean(targetWeightsAll_vis(ind));

        detectWeightOriBin_aud_respCells{icell}(i) = mean(detectWeightsAll_aud(ind));
        targetWeightOriBin_aud_respCells{icell}(i) = mean(targetWeightsAll_aud(ind));
    end
end


%%
figure

subplot 221
h = plot(oriPrefAll,targetWeightsAll_aud,'o');
h.Color = 'k';
hold on
h = plot(oriPrefAll(respCells_firstBase),...
    targetWeightsAll_aud(respCells_firstBase),'o');
h.Color = 'b';
h = plot(oriPrefAll(respCells_targetAndBase),...
    targetWeightsAll_aud(respCells_targetAndBase),'o');
h.Color = 'm';
h = plot(oriPrefAll(respCells_targetOnly),...
    targetWeightsAll_aud(respCells_targetOnly),'o');
h.Color = 'r';
figXAxis([],'Orientation Pref (deg)',[-10 190],0:45:180,0:45:180)
figYAxis([],'Weight',weightLim)
figAxForm
title('Target')

subplot 222
h = plot(oriPrefAll,detectWeightsAll_aud,'o');
h.Color = 'k';
hold on
h = plot(oriPrefAll(respCells_firstBase),...
    detectWeightsAll_aud(respCells_firstBase),'o');
h.Color = 'b';
h = plot(oriPrefAll(respCells_targetAndBase),...
    detectWeightsAll_aud(respCells_targetAndBase),'o');
h.Color = 'm';
h = plot(oriPrefAll(respCells_targetOnly),...
    detectWeightsAll_aud(respCells_targetOnly),'o');
h.Color = 'r';
figXAxis([],'Orientation Pref (deg)',[-10 190],0:45:180,0:45:180)
figYAxis([],'Weight',weightLim)
figAxForm
title('Detect')

subplot 223
for i = 1:3
    x = 1:maxBin;
    y = targetWeightOriBin_aud_respCells{i};
    hold on
    h = plot(x,y,'o-');
    h.Color = respCellColor{i};
end
hline(0,'k--')
figXAxis([],'Orientation Pref (deg)',[0 maxBin+1],1:maxBin,0:45:180)
figYAxis([],'Weight',weightLim)
figAxForm
title('Target')
subplot 224
for i = 1:3
    x = 1:maxBin;
    y = detectWeightOriBin_aud_respCells{i};
    hold on
    h = plot(x,y,'o-');
    h.Color = respCellColor{i};
end
hline(0,'k--')
figXAxis([],'Orientation Pref (deg)',[0 maxBin+1],1:maxBin,0:45:180)
figYAxis([],'Weight',weightLim)
figAxForm
title('Target')
legend(respCellName)

print(fullfile(fnout,'tuningXweight_respCellTypes_aud'),'-dpdf','-fillpage')

%%
figure
subplot 221
h = plot(oriPrefAll,targetWeightsAll_vis,'o');
h.Color = 'k';
hold on
h = plot(oriPrefAll(suppCells),...
    targetWeightsAll_vis(suppCells),'o');
h.Color = 'b';
h = plot(oriPrefAll(respCells_lateBaseOnly),...
    targetWeightsAll_vis(respCells_lateBaseOnly),'o');
h.Color = 'm';
h = plot(oriPrefAll(respCells_firstBase),...
    targetWeightsAll_vis(respCells_firstBase),'o');
h.Color = 'r';
figXAxis([],'Orientation Pref (deg)',[-10 190],0:45:180,0:45:180)
figYAxis([],'Weight',weightLim)
figAxForm
title('Target')

subplot 222
h = plot(oriPrefAll,detectWeightsAll_vis,'o');
h.Color = 'k';
hold on
h = plot(oriPrefAll(suppCells),...
    detectWeightsAll_vis(suppCells),'o');
h.Color = 'b';
h = plot(oriPrefAll(respCells_lateBaseOnly),...
    detectWeightsAll_vis(respCells_lateBaseOnly),'o');
h.Color = 'm';
h = plot(oriPrefAll(respCells_firstBase),...
    detectWeightsAll_vis(respCells_firstBase),'o');
h.Color = 'r';
figXAxis([],'Orientation Pref (deg)',[-10 190],0:45:180,0:45:180)
figYAxis([],'Weight',weightLim)
figAxForm
title('Detect')

%%
