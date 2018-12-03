% need to run FSAV_attnAnticipationAnalysis and FSAV_trialOutcomeAnalysis
% and FSAV_decodeDataStruct before running this analysis
%%
close all
clear all

ds = 'FSAV_attentionV1';

load('Z:\home\ashley\Analysis\FSAV Choice\FSAV_decodeData.mat')
load(['Z:\home\ashley\Analysis\FSAV Summaries\' ds '\attentionV1_startAlign_FSAV_attnData'])
fnout = ['Z:\home\ashley\Analysis\FSAV Summaries\FSAV Choice'];
nexp = length(dcExpt);

doExptPlots = 0;
%%
nBoot = 100;
minTrN = 15;
theta90Threshold = 11.25;
taskTuningFitCutoff = 0.666;
decisionVariable = 0.5;
detectThresholdHR = 0.8;
detectThresholdOri = 45;
% oriBinsFine = [0 1 16 22.5 32 45 64 90];
oriBins = [0 1 32 90];
orientationEdges = [0 32 90];
nBinnedOris = length(orientationEdges);
oriPrefBins = [0 22.5:45:157.5 181];
oriPrefBinCenter = [0:45:(180-45)];
ampBins = [0 0.0001 0.1 1];
ampEdges = [0 0.1 1];
nBinnedAmps = length(ampEdges);
weightsDiffZeroAlpha = 0.05/2;
%% calculated weights for detect and target LR models
% trials are matched for orientation of target - first binned into har
% targets, easy targets, and distractors, by what the variable oriBins is
% set to. Then trial numbers are balanced by selecting a random number of
% target trials from the easy and hard bins so that the sum total of
% targets is equal to the number of distractor trials.
targetWeightsAll_vis = [];
detectWeightsAll_vis = [];
targetWeightsAll_shuf_vis = [];
detectWeightsAll_shuf_vis = [];

targetCorrsAll_vis = [];
detectCorrsAll_vis = [];
targetCorrsAll_shuf_vis = [];
detectCorrsAll_shuf_vis = [];

targetWeightsAll_aud = [];
detectWeightsAll_aud = [];
targetWeightsAll_shuf_aud = [];
detectWeightsAll_shuf_aud = [];

targetCorrsAll_aud = [];
detectCorrsAll_aud = [];
targetCorrsAll_shuf_aud = [];
detectCorrsAll_shuf_aud = [];

pctCorrTarget_visXori = nan(nexp,nBinnedOris+1);
pctCorrDetect_visXori = nan(nexp,nBinnedOris+1);
pctCorrTarget_visXori_XVal = nan(nexp,nBinnedOris+1);
pctCorrDetect_visXori_XVal = nan(nexp,nBinnedOris+1);

pctCorrTarget_audXamp = nan(nexp,nBinnedAmps+1);
pctCorrDetect_audXamp = nan(nexp,nBinnedAmps+1);
pctCorrTarget_audXamp_XVal = nan(nexp,nBinnedAmps+1);
pctCorrDetect_audXamp_XVal = nan(nexp,nBinnedAmps+1);

pctCorr_target_testAudXamp = nan(nexp,nBinnedAmps+1);
pctCorr_detect_testAudXamp = nan(nexp,nBinnedAmps+1);
pctCorr_target_testVisXori = nan(nexp,nBinnedOris+1);
pctCorr_detect_testVisXori = nan(nexp,nBinnedOris+1);

respCells_targetOnly = [];
respCells_baseOnly = [];
respCells_targetAndBase = [];
respCells_lateBase = [];
suppCells = [];
respCells_dist = [];
respCells_targetAndDist = [];
respCells_target = [];

osi = [];
si_resp = [];
si_win = [];
avMod = [];
lateVisResp = [];
lateAudResp = [];
lateSubResp = [];

taskTuningFit = [];

targetAuroc = [];
detectAuroc = [];

nCells = zeros(1,nexp);
taskOriRespAll = [];
taskOriRespAll_firstStim0 = [];
oriPrefAll = [];
exptName = cell(1,nexp);
for iexp=1:nexp
    disp([dcExpt(iexp).mouse '-' dcExpt(iexp).date])
    %***VISUAL TRIALS
   trOut=dcExpt(iexp).trialOutcome;
    [Y_yeses_vis, Y_targets_vis] = getStimAndBehaviorYs(trOut);
    
   X_targets_vis=dcExpt(iexp).stimResp';
   
   cellIdx = dcExpt(iexp).signifResponsiveCells' & ...
       dcExpt(iexp).oriFitTheta90 < theta90Threshold;
   nCells(iexp) = sum(cellIdx);
    
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
    taskOriResp = nan(nBinnedOris,nCells(iexp));
    taskOriResp_firstStim0 = nan(nBinnedOris,nCells(iexp));
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
                taskOriResp_firstStim0(iori,:) = dcExpt(iexp).firstStimResp(cellIdx);
            else
                indSample = randsample(ind,minBinN);
                matchTrialsInd = cat(2,matchTrialsInd,indSample);
                oriInd = orientationEdges == orientations(iori);
                taskOriResp_firstStim0(oriInd,:) = mean(X_targets_vis(indSample,cellIdx),1);
            end
            oriInd = orientationEdges == orientations(iori);
            taskOriResp(oriInd,:) = mean(X_targets_vis(indSample,cellIdx),1);
        end
    end
    tOri = tOri(matchTrialsInd);
    taskOriRespAll = cat(2,taskOriRespAll,taskOriResp);
    taskOriRespAll_firstStim0 = cat(2,taskOriRespAll_firstStim0,taskOriResp_firstStim0);
    
   
   si_resp = cat(2,si_resp,attnInfoExpt(iexp).si_resp(cellIdx));
   si_win = cat(2,si_win,attnInfoExpt(iexp).si_win(cellIdx));
   avMod = cat(2,avMod,attnInfoExpt(iexp).avModTest(cellIdx));
    lateVisResp = cat(2,lateVisResp,attnInfoExpt(iexp).vLateResp(cellIdx));
    lateAudResp = cat(2,lateAudResp,attnInfoExpt(iexp).aLateResp(cellIdx));
    lateSubResp = cat(2,lateSubResp,...
        attnInfoExpt(iexp).vLateResp(cellIdx) - ...
        attnInfoExpt(iexp).aLateResp(cellIdx));
   
   tc = dcExpt(iexp).targetResponsiveCells(cellIdx);
   fbrc = dcExpt(iexp).firstBaseResponsiveCells(cellIdx);
   
   respCells_targetOnly = cat(1,respCells_targetOnly,tc & ~fbrc);
   respCells_baseOnly = cat(1,respCells_baseOnly,~tc & fbrc);
   respCells_targetAndBase = cat(1,respCells_targetAndBase, tc & fbrc);
   respCells_lateBase = cat(1,respCells_lateBase,...
       attnInfoExpt(iexp).lateBaseRespCells(cellIdx)');   
   suppCells = cat(1,suppCells,attnInfoExpt(iexp).lateBaseSuppCells(cellIdx)');
   
   respCells_dist = cat(1,respCells_dist,fbrc);
   respCells_targetAndDist = cat(1,respCells_targetAndDist,fbrc & tc);
   respCells_target = cat(1,respCells_target,tc);
   
   X_targets_vis=X_targets_vis(matchTrialsInd,cellIdx);
   
   taskfit = getTuningReliability(tOri,X_targets_vis',taskOriResp,nBoot);
   taskTuningFit = cat(1,taskTuningFit,taskfit); 
   
   Y_yeses_vis=Y_yeses_vis(matchTrialsInd);   
   Y_targets_vis=Y_targets_vis(matchTrialsInd); 
   
    target_auroc = nan(1,nCells(iexp));
    detect_auroc = nan(1,nCells(iexp));
    for i = 1:nCells(iexp)
        detect_auroc(i) = roc_gh(X_targets_vis(Y_yeses_vis == 0,i),...
            X_targets_vis(Y_yeses_vis == 1,i));
        target_auroc(i) = roc_gh(X_targets_vis(Y_targets_vis == 0,i),...
            X_targets_vis(Y_targets_vis == 1,i));
    end
    targetAuroc = cat(2,targetAuroc,target_auroc);
    detectAuroc = cat(2,detectAuroc,detect_auroc);
    
   detectCorr = corr(Y_yeses_vis,X_targets_vis)';
   targetCorr = corr(Y_targets_vis,X_targets_vis)';
   
   detectCorrsAll_vis = cat(1,detectCorrsAll_vis,detectCorr);
   targetCorrsAll_vis = cat(1,targetCorrsAll_vis,targetCorr);
    
   X_targets_vis=bsxfun(@plus,X_targets_vis,-mean(X_targets_vis));
   X_targets_vis=bsxfun(@times,X_targets_vis,1./std(X_targets_vis));
   C=eye(size(X_targets_vis,2));%C=inv(sqrtm(cov(X)));
   X_targets_vis=X_targets_vis*C;
   p=1;
   [temp1,dev1,targetGLM_vis]=glmfit(X_targets_vis,Y_targets_vis,'binomial');
   [temp2,dev2,detectGLM_vis]=glmfit(X_targets_vis,Y_yeses_vis,'binomial');
   idx=find(targetGLM_vis.p>p|detectGLM_vis.p>p);
   temp1(idx)=0;
   temp2(idx)=0;
   B_targets = C*temp1(2:end,1);
   B_yeses = C*temp2(2:end,1);
 
   targetWeightsAll_vis = cat(1,targetWeightsAll_vis,B_targets);
   detectWeightsAll_vis = cat(1,detectWeightsAll_vis,B_yeses);
   
   % test vis model with vis trials
   dv_detect_vis = mean(Y_yeses_vis);
   dv_target_vis = mean(Y_targets_vis);
   
   pctCorrTarget_visXori(iexp,end) = getPctCorr_trainData(...
       targetGLM_vis,X_targets_vis,Y_targets_vis,dv_target_vis);
   pctCorrDetect_visXori(iexp,end) = getPctCorr_trainData(...
       detectGLM_vis,X_targets_vis,Y_yeses_vis,dv_detect_vis);
   pctCorrTarget_visXori_XVal(iexp,end) = getPctCorr_hoData(...
       X_targets_vis,Y_targets_vis,dv_target_vis);
   pctCorrDetect_visXori_XVal(iexp,end) = getPctCorr_hoData(...
       X_targets_vis,Y_yeses_vis,dv_detect_vis);
   for iori = 1:nBinnedOris
       exptOriInd = find(orientations == orientationEdges(iori));
       if ~isempty(exptOriInd)
           tOriInd = tOri == orientationEdges(iori);
           X_ori = X_targets_vis(tOriInd,:);
           Y_targets_ori = Y_targets_vis(tOriInd);
           Y_yeses_ori = Y_yeses_vis(tOriInd);
       
           pctCorrTarget_visXori(iexp,exptOriInd) = getPctCorr_trainData(...
               targetGLM_vis,X_ori,Y_targets_ori,dv_target_vis);
           pctCorrDetect_visXori(iexp,exptOriInd) = getPctCorr_trainData(...
               detectGLM_vis,X_ori,Y_yeses_ori,dv_detect_vis);
           pctCorrTarget_visXori_XVal(iexp,exptOriInd) = getPctCorr_hoData_subGroup(...
               X_targets_vis,Y_targets_vis,find(tOriInd),dv_target_vis);
           pctCorrDetect_visXori_XVal(iexp,exptOriInd) = getPctCorr_hoData_subGroup(...
               X_targets_vis,Y_yeses_vis,find(tOriInd),dv_detect_vis);
       end
   end
   
   % shuffled data
    Y_yeses_shuf = Y_yeses_vis(randperm(length(Y_yeses_vis)));
    Y_targets_shuf = Y_targets_vis(randperm(length(Y_targets_vis)));
    
   detectCorr_shuf = corr(Y_yeses_shuf,X_targets_vis)';
   targetCorr_shuf = corr(Y_targets_shuf,X_targets_vis)';
   detectCorrsAll_shuf_vis = cat(1,detectCorrsAll_shuf_vis,detectCorr_shuf);
   targetCorrsAll_shuf_vis = cat(1,targetCorrsAll_shuf_vis,targetCorr_shuf);
   
   [temp1,dev1,targetGLM_vis]=glmfit(X_targets_vis,Y_targets_shuf,'binomial');
   [temp2,dev2,detectGLM_vis]=glmfit(X_targets_vis,Y_yeses_shuf,'binomial');
   idx=find(targetGLM_vis.p>p|detectGLM_vis.p>p);
   temp1(idx)=0;
   temp2(idx)=0;
   B_targets_shuf = C*temp1(2:end,1);
   B_yeses_shuf = C*temp2(2:end,1);
   targetWeightsAll_shuf_vis = cat(1,targetWeightsAll_shuf_vis,B_targets_shuf);
   detectWeightsAll_shuf_vis = cat(1,detectWeightsAll_shuf_vis,B_yeses_shuf);
   
    % orientaiton preference
   oriPref = dcExpt(iexp).oriPref(cellIdx);
   oriPrefAll = cat(2,oriPrefAll,oriPref);
   osi = cat(2,osi,dcExpt(iexp).osi(cellIdx));
   
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
    tAmp = tAmp(matchTrialsInd);

   X_targets_aud=dcExpt(iexp).audStimResp';

   X_targets_aud=X_targets_aud(matchTrialsInd,cellIdx);
   
   Y_yeses_aud=Y_yeses_aud(matchTrialsInd);   
   Y_targets_aud=Y_targets_aud(matchTrialsInd);   

   detectCorr = corr(Y_yeses_aud,X_targets_aud)';
   targetCorr = corr(Y_targets_aud,X_targets_aud)';

   detectCorrsAll_aud = cat(1,detectCorrsAll_aud,detectCorr);
   targetCorrsAll_aud = cat(1,targetCorrsAll_aud,targetCorr);

   X_targets_aud=bsxfun(@plus,X_targets_aud,-mean(X_targets_aud));
   X_targets_aud=bsxfun(@times,X_targets_aud,1./std(X_targets_aud));
   C=eye(size(X_targets_aud,2));%C=inv(sqrtm(cov(X)));
   X_targets_aud=X_targets_aud*C;
   p=1;
   [temp1,dev1,targetGLM_aud]=glmfit(X_targets_aud,Y_targets_aud,'binomial');
   [temp2,dev2,detectGLM_aud]=glmfit(X_targets_aud,Y_yeses_aud,'binomial');
   idx=find(targetGLM_aud.p>p|detectGLM_aud.p>p);
   temp1(idx)=0;
   temp2(idx)=0;
   B_targets = C*temp1(2:end,1);
   B_yeses = C*temp2(2:end,1);

   targetWeightsAll_aud = cat(1,targetWeightsAll_aud,B_targets);
   detectWeightsAll_aud = cat(1,detectWeightsAll_aud,B_yeses);
   
   % test vis model with vis trials
   dv_detect_aud = mean(Y_yeses_aud);
   dv_target_aud = mean(Y_targets_aud);
   
   pctCorrTarget_audXamp(iexp,end) = getPctCorr_trainData(...
       targetGLM_aud,X_targets_aud,Y_targets_aud,dv_target_aud);
   pctCorrDetect_audXamp(iexp,end) = getPctCorr_trainData(...
       detectGLM_aud,X_targets_aud,Y_yeses_aud,dv_detect_aud);
   pctCorrTarget_audXamp_XVal(iexp,end) = getPctCorr_hoData(...
       X_targets_aud,Y_targets_aud,dv_target_aud);
   pctCorrDetect_audXamp_XVal(iexp,end) = getPctCorr_hoData(...
       X_targets_aud,Y_yeses_aud,dv_detect_aud);
   for iamp = 1:nBinnedAmps
       exptAmpInd = find(amplitudes == ampEdges(iamp));
       if ~isempty(exptOriInd)
           tAmpInd = tAmp == ampEdges(iamp);
           X_amp = X_targets_aud(tAmpInd,:);
           Y_targets_amp = Y_targets_aud(tAmpInd);
           Y_yeses_amp = Y_yeses_aud(tAmpInd);
       
           pctCorrTarget_audXamp(iexp,exptAmpInd) = getPctCorr_trainData(...
               targetGLM_aud,X_amp,Y_targets_amp,dv_target_aud);
           pctCorrDetect_audXamp(iexp,exptAmpInd) = getPctCorr_trainData(...
               detectGLM_aud,X_amp,Y_yeses_amp,dv_detect_aud);
           pctCorrTarget_audXamp_XVal(iexp,exptAmpInd) = getPctCorr_hoData_subGroup(...
               X_targets_aud,Y_targets_aud,find(tAmpInd),dv_target_aud);
           pctCorrDetect_audXamp_XVal(iexp,exptAmpInd) = getPctCorr_hoData_subGroup(...
               X_targets_aud,Y_yeses_aud,find(tAmpInd),dv_detect_aud);
       end
   end
   
   % shuffled data
    Y_yeses_shuf = Y_yeses_aud(randperm(length(Y_yeses_aud)));
    Y_targets_shuf = Y_targets_aud(randperm(length(Y_targets_aud)));
    
   detectCorr_shuf = corr(Y_yeses_shuf,X_targets_aud)';
   targetCorr_shuf = corr(Y_targets_shuf,X_targets_aud)';
   detectCorrsAll_shuf_aud = cat(1,detectCorrsAll_shuf_aud,detectCorr_shuf);
   targetCorrsAll_shuf_aud = cat(1,targetCorrsAll_shuf_aud,targetCorr_shuf);
   
   [temp1,dev1,targetGLM_aud]=glmfit(X_targets_aud,Y_targets_shuf,'binomial');
   [temp2,dev2,detectGLM_aud]=glmfit(X_targets_aud,Y_yeses_shuf,'binomial');
   idx=find(targetGLM_aud.p>p|detectGLM_aud.p>p);
   temp1(idx)=0;
   temp2(idx)=0;
   B_targets_shuf = C*temp1(2:end,1);
   B_yeses_shuf = C*temp2(2:end,1);
   targetWeightsAll_shuf_aud = cat(1,targetWeightsAll_shuf_aud,B_targets_shuf);
   detectWeightsAll_shuf_aud = cat(1,detectWeightsAll_shuf_aud,B_yeses_shuf);
   
   % ***Test Auditory with Visual Model
   
   pctCorr_target_testAudXamp(iexp,end) = getPctCorr_trainData(...
       targetGLM_vis,X_targets_aud,Y_targets_aud,dv_target_aud);
   pctCorr_detect_testAudXamp(iexp,end) = getPctCorr_trainData(...
       detectGLM_vis,X_targets_aud,Y_yeses_aud,dv_detect_aud);
   for iamp = 1:nBinnedAmps
       exptAmpInd = find(amplitudes == ampEdges(iamp));
       if ~isempty(exptOriInd)
           tAmpInd = tAmp == ampEdges(iamp);
           X_amp = X_targets_aud(tAmpInd,:);
           Y_targets_amp = Y_targets_aud(tAmpInd);
           Y_yeses_amp = Y_yeses_aud(tAmpInd);
       
           pctCorr_target_testAudXamp(iexp,exptAmpInd) = getPctCorr_trainData(...
               targetGLM_vis,X_amp,Y_targets_amp,dv_target_aud);
           pctCorr_detect_testAudXamp(iexp,exptAmpInd) = getPctCorr_trainData(...
               detectGLM_vis,X_amp,Y_yeses_amp,dv_detect_aud);
       end
   end
   
   % ***Test Visual with Auditory Model
   
   pctCorr_target_testVisXori(iexp,end) = getPctCorr_trainData(...
       targetGLM_aud,X_targets_vis,Y_targets_vis,dv_target_vis);
   pctCorr_detect_testVisXori(iexp,end) = getPctCorr_trainData(...
       detectGLM_aud,X_targets_vis,Y_yeses_vis,dv_detect_vis);
   for iori = 1:nBinnedOris
       exptOriInd = find(orientations == orientationEdges(iori));
       if ~isempty(exptOriInd)
           tOriInd = tOri == orientationEdges(iori);
           X_ori = X_targets_vis(tOriInd,:);
           Y_targets_ori = Y_targets_vis(tOriInd);
           Y_yeses_ori = Y_yeses_vis(tOriInd);
       
           pctCorr_target_testVisXori(iexp,exptOriInd) = getPctCorr_trainData(...
               targetGLM_aud,X_ori,Y_targets_ori,dv_target_vis);
           pctCorr_detect_testVisXori(iexp,exptOriInd) = getPctCorr_trainData(...
               detectGLM_aud,X_ori,Y_yeses_ori,dv_detect_vis);
       end
   end
   
   
end
respCells_targetOnly = logical(respCells_targetOnly);
respCells_firstBase = logical(respCells_baseOnly);
respCells_targetAndBase = logical(respCells_targetAndBase);
respCells_lateBaseOnly = logical(respCells_lateBase &...
    ~respCells_baseOnly & ~respCells_targetAndBase);
suppCells = logical(suppCells);

respCells_dist = logical(respCells_dist);
respCells_targetAndDist = logical(respCells_targetAndDist);
respCells_target = logical(respCells_target);

[~,taskOriPrefID] = max(taskOriRespAll);
taskOriPref = orientationEdges(taskOriPrefID);
isTaskTuned = taskTuningFit > taskTuningFitCutoff;
%% bin weights by orienation preference, test effect of orientation by ANOVA
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

%% bin weights by ori pref and cell responsivity type
respCellType = {respCells_target & ~respCells_dist;respCells_dist & ~respCells_target;respCells_targetAndDist};
respCellColor = {'r';'b';'m'};
respCellName = {'Target Only'; 'Base Only'; 'Both'};

ind1 = respCells_target & ~respCells_dist;
ind2 = respCells_targetAndDist;
ind3 = respCells_dist & ~respCells_target;

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

%weights binned by task tuning
detectWeightTaskOriBin_vis_respCells = cell(1,3);
targetWeightTaskOriBin_vis_respCells = cell(1,3);
detectWeightTaskOriBin_aud_respCells = cell(1,3);
targetWeightTaskOriBin_aud_respCells = cell(1,3);

detectWeightTaskOriBin_err_vis_respCells = cell(1,3);
targetWeightTaskOriBin_err_vis_respCells = cell(1,3);
detectWeightTaskOriBin_err_aud_respCells = cell(1,3);
targetWeightTaskOriBin_err_aud_respCells = cell(1,3);

detectWeightTaskOriBin_vis_respCells(:) = {nan(1,nBinnedOris)};
targetWeightTaskOriBin_vis_respCells(:) = {nan(1,nBinnedOris)};
detectWeightTaskOriBin_aud_respCells(:) = {nan(1,nBinnedOris)};
targetWeightTaskOriBin_aud_respCells(:) = {nan(1,nBinnedOris)};

detectWeightTaskOriBin_err_vis_respCells(:) = {nan(1,nBinnedOris)};
targetWeightTaskOriBin_err_vis_respCells(:) = {nan(1,nBinnedOris)};
detectWeightTaskOriBin_err_aud_respCells(:) = {nan(1,nBinnedOris)};
targetWeightTaskOriBin_err_aud_respCells(:) = {nan(1,nBinnedOris)};
for icell = 1:3
    respCellInd = respCellType{icell}';
    for i = 1:nBinnedOris
        ind = taskOriPrefID == i & respCellInd & isTaskTuned';

        detectWeightTaskOriBin_vis_respCells{icell}(i) = mean(detectWeightsAll_vis(ind));
        targetWeightTaskOriBin_vis_respCells{icell}(i) = mean(targetWeightsAll_vis(ind));
        detectWeightTaskOriBin_err_vis_respCells{icell}(i) = ste(detectWeightsAll_vis(ind),1);
        targetWeightTaskOriBin_err_vis_respCells{icell}(i) = ste(targetWeightsAll_vis(ind),1);

        detectWeightTaskOriBin_aud_respCells{icell}(i) = mean(detectWeightsAll_aud(ind));
        targetWeightTaskOriBin_aud_respCells{icell}(i) = mean(targetWeightsAll_aud(ind));
        detectWeightTaskOriBin_err_aud_respCells{icell}(i) = ste(detectWeightsAll_aud(ind),1);
        targetWeightTaskOriBin_err_aud_respCells{icell}(i) = ste(targetWeightsAll_aud(ind),1);
    end
end

%% bin detect weights by target weights > or < 0

ind = targetWeightsAll_vis < 0;
targetWeightsBinned_vis{1} = targetWeightsAll_vis(ind);
detectWeightsBinned_vis{1} = detectWeightsAll_vis(ind);
targetWeightsBinned_vis{2} = targetWeightsAll_vis(~ind);
detectWeightsBinned_vis{2} = detectWeightsAll_vis(~ind);

ind = targetWeightsAll_aud < 0;
targetWeightsBinned_aud{1} = targetWeightsAll_aud(ind);
detectWeightsBinned_aud{1} = detectWeightsAll_aud(ind);
targetWeightsBinned_aud{2} = targetWeightsAll_aud(~ind);
detectWeightsBinned_aud{2} = detectWeightsAll_aud(~ind);

targetWeightsBinned_shuff_vis = cell(1,2);
detectWeightsBinned_shuff_vis = cell(1,2);
ind = targetWeightsAll_shuf_vis < 0;
targetWeightsBinned_shuff_vis{1} = targetWeightsAll_shuf_vis(ind);
detectWeightsBinned_shuff_vis{1} = detectWeightsAll_shuf_vis(ind);
targetWeightsBinned_shuff_vis{2} = targetWeightsAll_shuf_vis(~ind);
detectWeightsBinned_shuff_vis{2} = detectWeightsAll_shuf_vis(~ind);

targetWeightsBinned_shuff_aud = cell(1,2);
detectWeightsBinned_shuff_aud = cell(1,2);
ind = targetWeightsAll_shuf_aud < 0;
targetWeightsBinned_shuff_aud{1} = targetWeightsAll_shuf_aud(ind);
detectWeightsBinned_shuff_aud{1} = detectWeightsAll_shuf_aud(ind);
targetWeightsBinned_shuff_aud{2} = targetWeightsAll_shuf_aud(~ind);
detectWeightsBinned_shuff_aud{2} = detectWeightsAll_shuf_aud(~ind);
%%
weightLim = [-2 4];
% audWeightLim = [-2 2];
polarWeightLim = [0 3];
% polarAudWeightLim = [0 2];
siLim = [-12 12];
respLim = [-0.05 0.1];
respLimCropped = [-0.03 0.03];
corrLim = [-1 1];
%% plot weight and selectivity by late dF/F response

figure
suptitle('Visual Trials Model')

subplot 221
h=plot(lateVisResp,targetWeightsAll_vis,'k.');
h.MarkerFaceColor = 'k';
figXAxis([],'Late Visual Response (dF/F)',respLim)
figYAxis([],'Target Weight',weightLim)
figAxForm

subplot 222
h=plot(lateVisResp,detectWeightsAll_vis,'k.');
h.MarkerFaceColor = 'k';
figXAxis([],'Late Visual Response (dF/F)',respLim)
figYAxis([],'Detect Weight',weightLim)
figAxForm

subplot 223
h=plot(lateVisResp,targetWeightsAll_vis,'k.');
h.MarkerFaceColor = 'k';
figXAxis([],'Late Visual Response (dF/F)',respLimCropped)
figYAxis([],'Target Weight',weightLim)
figAxForm

subplot 224
h=plot(lateVisResp,detectWeightsAll_vis,'k.');
h.MarkerFaceColor = 'k';
figXAxis([],'Late Visual Response (dF/F)',respLimCropped)
figYAxis([],'Detect Weight',weightLim)
figAxForm

print(fullfile(fnout,'LRmodelXlateVisResp'),'-dpdf','-fillpage')
%% plot weight by selectivity

setFigParams4Print('portrait')
set(0,'defaultAxesFontSize',12)
figure
suptitle('Visual Trials Model')
ind1 = logical(avMod');
subplot 221
h=plot(si_resp,targetCorrsAll_vis,'o');
h.Color = [0.5 0.5 0.5];
h.MarkerFaceColor = [0.5 0.5 0.5];
hold on
h=plot(si_resp(ind1),targetCorrsAll_vis(ind1),'ko');
h.MarkerFaceColor = 'k';
figXAxis([],'V-A Selectivity',siLim)
figYAxis([],'Correlation',corrLim)
figAxForm
title('Target')
subplot 222
h=plot(si_resp,detectCorrsAll_vis,'o');
h.Color = [0.5 0.5 0.5];
h.MarkerFaceColor = [0.5 0.5 0.5];
hold on
h=plot(si_resp(ind1),detectCorrsAll_vis(ind1),'ko');
h.MarkerFaceColor = 'k';
figXAxis([],'V-A Selectivity',siLim)
figYAxis([],'Correlation',corrLim)
figAxForm
title('Detect')
subplot 223
h=plot(si_resp,targetWeightsAll_vis,'o');
h.Color = [0.5 0.5 0.5];
h.MarkerFaceColor = [0.5 0.5 0.5];
hold on
h=plot(si_resp(ind1),targetWeightsAll_vis(ind1),'ko');
h.MarkerFaceColor = 'k';
figXAxis([],'V-A Selectivity',siLim)
figYAxis([],'Weight',weightLim)
figAxForm
title('Target')
subplot 224
h=plot(si_resp,detectWeightsAll_vis,'ko');
h.Color = [0.5 0.5 0.5];
h.MarkerFaceColor = [0.5 0.5 0.5];
hold on
h=plot(si_resp(ind1),detectWeightsAll_vis(ind1),'ko');
h.MarkerFaceColor = 'k';
figXAxis([],'V-A Selectivity',siLim)
figYAxis([],'Weight',weightLim)
figAxForm
title('Detect')
legend({'All';'Selective Only'},...
    'location','northeast')

print(fullfile(fnout,'LRmodelXattn'),'-dpdf','-fillpage')

figure
ind1 = respCells_target & ~respCells_dist;
ind2 = respCells_targetAndDist;
ind3 = respCells_dist & ~respCells_target;
subplot 221
h=plot(si_resp(ind1),targetWeightsAll_vis(ind1),'ro');
h.MarkerFaceColor = h.Color;
hold on
h=plot(si_resp(ind2),targetWeightsAll_vis(ind2),'mo');
h.MarkerFaceColor = h.Color;
h=plot(si_resp(ind3),targetWeightsAll_vis(ind3),'bo');
h.MarkerFaceColor = h.Color;
% h=plot(si_resp(ind1),targetWeightsAll_vis(ind1),'ko');
% h.MarkerFaceColor = 'k';
figXAxis([],'V-A Resp SI',siLim)
figYAxis([],'Weight',weightLim)
figAxForm
hline(0,'k:')
vline(0,'k:')
title('Target')
subplot 222
h=plot(si_resp(ind1),detectWeightsAll_vis(ind1),'ro');
h.MarkerFaceColor = h.Color;
hold on
h=plot(si_resp(ind2),detectWeightsAll_vis(ind2),'mo');
h.MarkerFaceColor = h.Color;
h=plot(si_resp(ind3),detectWeightsAll_vis(ind3),'bo');
h.MarkerFaceColor = h.Color;
% h=plot(si_resp(ind1),detectWeightsAll_vis(ind1),'ko');
% h.MarkerFaceColor = 'k';
figXAxis([],'V-A Resp SI',siLim)
figYAxis([],'Weight',weightLim)
figAxForm
hline(0,'k:')
vline(0,'k:')
title('Detect')
subplot 223
h=plot(si_win(ind1),targetWeightsAll_vis(ind1),'ro');
h.MarkerFaceColor = h.Color;
hold on
h=plot(si_win(ind2),targetWeightsAll_vis(ind2),'mo');
h.MarkerFaceColor = h.Color;
h=plot(si_win(ind3),targetWeightsAll_vis(ind3),'bo');
h.MarkerFaceColor = h.Color;
% h=plot(si_resp(ind1),targetWeightsAll_vis(ind1),'ko');
% h.MarkerFaceColor = 'k';
figXAxis([],'V-A Win SI',siLim)
figYAxis([],'Weight',weightLim)
figAxForm
hline(0,'k:')
vline(0,'k:')
title('Target')
subplot 224
h=plot(si_win(ind1),detectWeightsAll_vis(ind1),'ro');
h.MarkerFaceColor = h.Color;
hold on
h=plot(si_win(ind2),detectWeightsAll_vis(ind2),'mo');
h.MarkerFaceColor = h.Color;
h=plot(si_win(ind3),detectWeightsAll_vis(ind3),'bo');
h.MarkerFaceColor = h.Color;
% h=plot(si_resp(ind1),detectWeightsAll_vis(ind1),'ko');
% h.MarkerFaceColor = 'k';
figXAxis([],'V-A Win SI',siLim)
figYAxis([],'Weight',weightLim)
figAxForm
hline(0,'k:')
vline(0,'k:')
title('Detect')

print(fullfile(fnout,'LRmodelXattn_resp&win'),'-dpdf','-fillpage')

figure
ind1 = respCells_target & ~respCells_dist;
ind2 = respCells_targetAndDist;
ind3 = respCells_dist & ~respCells_target;
subplot 221
h=errorbar(mean(si_resp(ind1)),mean(targetWeightsAll_vis(ind1)),...
    ste(targetWeightsAll_vis(ind1),1),ste(targetWeightsAll_vis(ind1),1),...
    ste(si_resp(ind1),2),ste(si_resp(ind1),2),'ro');
hold on
h=errorbar(mean(si_resp(ind2)),mean(targetWeightsAll_vis(ind2)),...
    ste(targetWeightsAll_vis(ind2),1),ste(targetWeightsAll_vis(ind2),1),...
    ste(si_resp(ind2),2),ste(si_resp(ind2),2),'mo');
h=errorbar(mean(si_resp(ind3)),mean(targetWeightsAll_vis(ind3)),...
    ste(targetWeightsAll_vis(ind3),1),ste(targetWeightsAll_vis(ind3),1),...
    ste(si_resp(ind3),2),ste(si_resp(ind3),2),'bo');
% h=plot(si_resp(ind2),targetWeightsAll_vis(ind1),'ko');
% h.MarkerFaceColor = 'k';
figXAxis([],'V-A Resp SI',[-4 4])
figYAxis([],'Weight',[-1 1])
figAxForm
hline(0,'k:')
vline(0,'k:')
title('Target')
subplot 222
h=errorbar(mean(si_resp(ind1)),mean(detectWeightsAll_vis(ind1)),...
    ste(detectWeightsAll_vis(ind1),1),ste(detectWeightsAll_vis(ind1),1),...
    ste(si_resp(ind1),2),ste(si_resp(ind1),2),'ro');
hold on
h=errorbar(mean(si_resp(ind2)),mean(detectWeightsAll_vis(ind2)),...
    ste(detectWeightsAll_vis(ind2),1),ste(detectWeightsAll_vis(ind2),1),...
    ste(si_resp(ind2),2),ste(si_resp(ind2),2),'mo');
h=errorbar(mean(si_resp(ind3)),mean(detectWeightsAll_vis(ind3)),...
    ste(detectWeightsAll_vis(ind3),1),ste(detectWeightsAll_vis(ind3),1),...
    ste(si_resp(ind3),2),ste(si_resp(ind3),2),'bo');
% h=plot(si_resp(ind1),detectWeightsAll_vis(ind1),'ko');
% h.MarkerFaceColor = 'k';
figXAxis([],'V-A Resp SI',[-4 4])
figYAxis([],'Weight',[-1 1])
figAxForm
hline(0,'k:')
vline(0,'k:')
title('Detect')
subplot 223
h=errorbar(mean(si_win(ind1)),mean(targetWeightsAll_vis(ind1)),...
    ste(targetWeightsAll_vis(ind1),1),ste(targetWeightsAll_vis(ind1),1),...
    ste(si_win(ind1),2),ste(si_win(ind1),2),'ro');
hold on
h=errorbar(mean(si_win(ind2)),mean(targetWeightsAll_vis(ind2)),...
    ste(targetWeightsAll_vis(ind2),1),ste(targetWeightsAll_vis(ind2),1),...
    ste(si_win(ind2),2),ste(si_win(ind2),2),'mo');
h=errorbar(mean(si_win(ind3)),mean(targetWeightsAll_vis(ind3)),...
    ste(targetWeightsAll_vis(ind3),1),ste(targetWeightsAll_vis(ind3),1),...
    ste(si_win(ind3),2),ste(si_win(ind3),2),'bo');
% h=plot(si_resp(ind1),targetWeightsAll_vis(ind1),'ko');
% h.MarkerFaceColor = 'k';
figXAxis([],'V-A Win SI',[-2 2])
figYAxis([],'Weight',[-.75 .75])
figAxForm
hline(0,'k:')
vline(0,'k:')
title('Target')
subplot 224
h=errorbar(mean(si_win(ind1)),mean(detectWeightsAll_vis(ind1)),...
    ste(detectWeightsAll_vis(ind1),1),ste(detectWeightsAll_vis(ind1),1),...
    ste(si_win(ind1),2),ste(si_win(ind1),2),'ro');
hold on
h=errorbar(mean(si_win(ind2)),mean(detectWeightsAll_vis(ind2)),...
    ste(detectWeightsAll_vis(ind2),1),ste(detectWeightsAll_vis(ind2),1),...
    ste(si_win(ind2),2),ste(si_win(ind2),2),'mo');
h=errorbar(mean(si_win(ind3)),mean(detectWeightsAll_vis(ind3)),...
    ste(detectWeightsAll_vis(ind3),1),ste(detectWeightsAll_vis(ind3),1),...
    ste(si_win(ind3),2),ste(si_win(ind3),2),'bo');
% h=plot(si_resp(ind1),detectWeightsAll_vis(ind1),'ko');
% h.MarkerFaceColor = 'k';
figXAxis([],'V-A Win SI',[-2 2])
figYAxis([],'Weight',[-.75 .75])
figAxForm
hline(0,'k:')
vline(0,'k:')
title('Detect')

print(fullfile(fnout,'LRmodelXattn_resp&win_sum'),'-dpdf','-fillpage')

figure
plot(si_resp,si_win,'ko')
hold on
plot([-10 10],[-10 10],'k--')
figXAxis([],'V-A Resp SI',[-10 10])
figYAxis([],'V-A Win SI',[-10 10])
figAxForm
hline(0,'k:')
vline(0,'k:')
title('Selectivity Correlation of Model Neurons')

print(fullfile(fnout,'siRespXsiWin_modelCells'),'-dpdf','-fillpage')

setFigParams4Print('portrait')
set(0,'defaultAxesFontSize',12)
figure
suptitle('Visual Trials Model')
ind1 = respCells_targetOnly;
ind2 = respCells_targetAndBase;
ind3 = respCells_firstBase;
subplot 221
h=plot(si_resp,targetCorrsAll_vis,'o');
h.Color = [0.5 0.5 0.5];
h.MarkerFaceColor = [0.5 0.5 0.5];
hold on
h=plot(si_resp(ind1),targetCorrsAll_vis(ind1),'ro');
h.MarkerFaceColor = 'r';
h=plot(si_resp(ind2),targetCorrsAll_vis(ind2),'mo');
h.MarkerFaceColor = 'm';
h=plot(si_resp(ind3),targetCorrsAll_vis(ind3),'bo');
h.MarkerFaceColor = 'b';
figXAxis([],'V-A Selectivity',siLim)
figYAxis([],'Correlation',corrLim)
figAxForm
title('Target')
subplot 222
h=plot(si_resp,detectCorrsAll_vis,'o');
h.Color = [0.5 0.5 0.5];
h.MarkerFaceColor = [0.5 0.5 0.5];
hold on
h=plot(si_resp(ind1),detectCorrsAll_vis(ind1),'ro');
h.MarkerFaceColor = 'r';
h=plot(si_resp(ind2),detectCorrsAll_vis(ind2),'mo');
h.MarkerFaceColor = 'm';
h=plot(si_resp(ind3),detectCorrsAll_vis(ind3),'bo');
h.MarkerFaceColor = 'b';
figXAxis([],'V-A Selectivity',siLim)
figYAxis([],'Correlation',corrLim)
figAxForm
title('Detect')
subplot 223
h=plot(si_resp,targetWeightsAll_vis,'o');
h.Color = [0.5 0.5 0.5];
h.MarkerFaceColor = [0.5 0.5 0.5];
hold on
h=plot(si_resp(ind1),targetWeightsAll_vis(ind1),'ro');
h.MarkerFaceColor = 'r';
h=plot(si_resp(ind2),targetWeightsAll_vis(ind2),'mo');
h.MarkerFaceColor = 'm';
h=plot(si_resp(ind3),targetWeightsAll_vis(ind3),'bo');
h.MarkerFaceColor = 'b';
figXAxis([],'V-A Selectivity',siLim)
figYAxis([],'Weight',weightLim)
figAxForm
title('Target')
subplot 224
h=plot(si_resp,detectWeightsAll_vis,'ko');
h.Color = [0.5 0.5 0.5];
h.MarkerFaceColor = [0.5 0.5 0.5];
hold on
h=plot(si_resp(ind1),detectWeightsAll_vis(ind1),'ro');
h.MarkerFaceColor = 'r';
h=plot(si_resp(ind2),detectWeightsAll_vis(ind2),'mo');
h.MarkerFaceColor = 'm';
h=plot(si_resp(ind3),detectWeightsAll_vis(ind3),'bo');
h.MarkerFaceColor = 'b';
figXAxis([],'V-A Selectivity',siLim)
figYAxis([],'Weight',weightLim)
figAxForm
title('Detect')
legend({'All';'Target Only';'Target & Distractor';'Distractor Only'},...
    'location','northeast')

print(fullfile(fnout,'LRmodelXattn2'),'-dpdf','-fillpage')

%% bin selectivity by SI and task responsiveness
siBins = [-12 -6.25:2.5:6.25 12];
[nSiPerBin,~,siBinID] = histcounts(discretize(si_resp,siBins,'IncludedEdge','right'));
nSiBins = length(nSiPerBin);
siBinLabel = [8 siBins(2:end-2)+1.25 8];

setFigParams4Print('portrait')
set(0,'defaultAxesFontSize',12)
figure
suptitle('Visual Trials Model')
ind1 = respCells_targetOnly;
ind2 = respCells_targetAndBase;
ind3 = respCells_firstBase;
subplot 221
x = 1:nSiBins;
y = nan(1,nSiBins);
yerr = nan(1,nSiBins);
for i = 1:nSiBins
    ind_si = siBinID' == i;
    y(i) = mean(targetWeightsAll_vis(ind1 & ind_si));
    yerr(i) = ste(targetWeightsAll_vis(ind1 & ind_si),1);
end
h=plot(x,y,'ro-');
h.MarkerFaceColor = [1 1 1];
hold on
y = nan(1,nSiBins);
yerr = nan(1,nSiBins);
for i = 1:nSiBins
    ind_si = siBinID' == i;
    y(i) = mean(targetWeightsAll_vis(ind2 & ind_si));
    yerr(i) = ste(targetWeightsAll_vis(ind2 & ind_si),1);
end
h=plot(x,y,'mo-');
h.MarkerFaceColor = [1 1 1];
y = nan(1,nSiBins);
yerr = nan(1,nSiBins);
for i = 1:nSiBins
    ind_si = siBinID' == i;
    y(i) = mean(targetWeightsAll_vis(ind3 & ind_si));
    yerr(i) = ste(targetWeightsAll_vis(ind3 & ind_si),1);
end
h=plot(x,y,'bo-');
h.MarkerFaceColor = [1 1 1];
figXAxis([],'V-A Selectivity',[0 nSiBins+1],x,siBinLabel)
figYAxis([],'Weight',weightLim)
figAxForm
title('Target')
subplot 222
x = 1:nSiBins;
y = nan(1,nSiBins);
yerr = nan(1,nSiBins);
for i = 1:nSiBins
    ind_si = siBinID' == i;
    y(i) = mean(detectWeightsAll_vis(ind1 & ind_si));
    yerr(i) = ste(detectWeightsAll_vis(ind1 & ind_si),1);
end
h=plot(x,y,'ro-');
h.MarkerFaceColor = [1 1 1];
hold on
y = nan(1,nSiBins);
yerr = nan(1,nSiBins);
for i = 1:nSiBins
    ind_si = siBinID' == i;
    y(i) = mean(detectWeightsAll_vis(ind2 & ind_si));
    yerr(i) = ste(detectWeightsAll_vis(ind2 & ind_si),1);
end
h=plot(x,y,'mo-');
h.MarkerFaceColor = [1 1 1];
y = nan(1,nSiBins);
yerr = nan(1,nSiBins);
for i = 1:nSiBins
    ind_si = siBinID' == i;
    y(i) = mean(detectWeightsAll_vis(ind3 & ind_si));
    yerr(i) = ste(detectWeightsAll_vis(ind3 & ind_si),1);
end
h=plot(x,y,'bo-');
h.MarkerFaceColor = [1 1 1];
figXAxis([],'V-A Selectivity',[0 nSiBins+1],x,siBinLabel)
figYAxis([],'Weight',weightLim)
figAxForm
title('Detect')

subplot 223
h=plot(si_resp,targetWeightsAll_vis,'o');
h.Color = [0.5 0.5 0.5];
h.MarkerFaceColor = [0.5 0.5 0.5];
hold on
h=plot(si_resp(ind1),targetWeightsAll_vis(ind1),'ro');
h.MarkerFaceColor = 'r';
h=plot(si_resp(ind2),targetWeightsAll_vis(ind2),'mo');
h.MarkerFaceColor = 'm';
h=plot(si_resp(ind3),targetWeightsAll_vis(ind3),'bo');
h.MarkerFaceColor = 'b';
figXAxis([],'V-A Selectivity',siLim)
figYAxis([],'Weight',weightLim)
figAxForm
title('Target')
subplot 224
h=plot(si_resp,detectWeightsAll_vis,'ko');
h.Color = [0.5 0.5 0.5];
h.MarkerFaceColor = [0.5 0.5 0.5];
hold on
h=plot(si_resp(ind1),detectWeightsAll_vis(ind1),'ro');
h.MarkerFaceColor = 'r';
h=plot(si_resp(ind2),detectWeightsAll_vis(ind2),'mo');
h.MarkerFaceColor = 'm';
h=plot(si_resp(ind3),detectWeightsAll_vis(ind3),'bo');
h.MarkerFaceColor = 'b';
figXAxis([],'V-A Selectivity',siLim)
figYAxis([],'Weight',weightLim)
figAxForm
title('Detect')
legend({'All';'Target Only';'Target & Distractor';'Distractor Only'},...
    'location','northeast')

print(fullfile(fnout,'LRmodelXattnSummary'),'-dpdf','-fillpage')

%% plot weights by orientation
weightLim_ori = [-1.3 1.3];

setFigParams4Print('portrait')
set(0,'defaultAxesFontSize',12)

figure
suptitle('Visual Trials')
subplot 321
bar(1:maxBin,n)
figXAxis([],'Ori Pref Bin',[0 maxBin+1],[1:maxBin],oriPrefBinCenter)
figYAxis([],'N Cells',[]);
figAxForm
subplot 322
errorbar(1:maxBin,targetWeightOriBin_vis,targetWeightOriBinErr_vis,'ko-');
hold on
errorbar(1:maxBin,detectWeightOriBin_vis,detectWeightOriBinErr_vis,'bo-');
figXAxis([],'Ori Pref Bin',[0 maxBin+1],[1:maxBin],oriPrefBinCenter)
figYAxis([],'Weight',weightLim_ori);
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

%% plot weights, correlations, shuffles, and quantification of detect vs. target weight

setFigParams4Print('portrait')
set(0,'defaultAxesFontSize',12)
%***VISUAL TRIALS
figure
suptitle('Visual Trials')
subplot 321
h=plot(detectWeightsAll_vis,detectCorrsAll_vis,'ko');
h.MarkerFaceColor = 'k';
hold on
hline(0,'k:')
vline(0,'k:')
figXAxis([],'Detect LR Weight',weightLim)
figYAxis([],'Detect Correlation',corrLim)
figAxForm
title('All Trials')
subplot 323
h=plot(targetWeightsAll_vis,targetCorrsAll_vis,'ko');
h.MarkerFaceColor = 'k';
hold on
figXAxis([],'Target LR Weight',weightLim)
figYAxis([],'Target Correlation',corrLim)
hline(0,'k:')
vline(0,'k:')
figAxForm
subplot 325   
h=plot(targetWeightsAll_vis,detectWeightsAll_vis,'ko');
h.MarkerFaceColor = 'k';
hold on
figXAxis([],'Target LR Weight',weightLim)
figYAxis([],'Detect LR Weight',weightLim)
hline(0,'k:')
vline(0,'k:')
figAxForm
subplot 326   
h=plot(targetWeightsAll_shuf_vis,detectWeightsAll_shuf_vis,'ko');
h.MarkerFaceColor = 'k';
hold on
figXAxis([],'Target LR Weight',weightLim)
figYAxis([],'Detect LR Weight',weightLim)
hline(0,'k:')
vline(0,'k:')
figAxForm

subplot 322
p = nan(1,2);
for i = 1:2
    x = mean(targetWeightsBinned_vis{i});
    xerr = ste(targetWeightsBinned_vis{i},1);
    y = detectWeightsBinned_vis{i};
    yerr = ste(y,1);
    [ttest_pass,p(i)] = ttest(y,[],'alpha',0.025);
    plot(x,y,'k.');
    hold on
    h = errorbar(x,mean(y),yerr,yerr,xerr,xerr,'ko');
    if ttest_pass
        h.Color = 'r';
    end
end
figXAxis([],'Target LR Weight',weightLim)
figYAxis([],'Detect LR Weight',weightLim)
hline(0,'k:')
vline(0,'k:')
figAxForm
title(sprintf('hi/lo target, p = %s/%s',round(p,2,'significant')))

subplot 324
p = nan(1,2);
for i = 1:2
    x = mean(targetWeightsBinned_shuff_vis{i});
    xerr = ste(targetWeightsBinned_shuff_vis{i},1);
    y = detectWeightsBinned_shuff_vis{i};
    yerr = ste(y,1);
    [ttest_pass,p(i)] = ttest(y,[],'alpha',0.025);
    plot(x,y,'k.');
    hold on
    h = errorbar(x,mean(y),yerr,yerr,xerr,xerr,'ko');
    if ttest_pass
        h.Color = 'r';
    end
end
figXAxis([],'Target LR Weight',weightLim)
figYAxis([],'Detect LR Weight',weightLim)
hline(0,'k:')
vline(0,'k:')
figAxForm
title(sprintf('SHUFFLE: hi/lo target, p = %s/%s',round(p,2,'significant')))

print(fullfile(fnout,'corrXWeight_shuffle_vis'),'-dpdf','-fillpage')

%***AUDITORY TRIALS
figure
suptitle('Auditory Trials')
subplot 321
h = plot(detectWeightsAll_aud,detectCorrsAll_aud,'ko');
h.MarkerFaceColor = 'k';
hold on
figXAxis([],'Detect LR Weight',weightLim)
figYAxis([],'Detect Correlation',corrLim)
figAxForm
hline(0,'k:')
vline(0,'k:')
title('All Trials')
subplot 323
h = plot(targetWeightsAll_aud,targetCorrsAll_aud,'ko')
h.MarkerFaceColor = 'k';
hold on
figXAxis([],'Target LR Weight',weightLim)
figYAxis([],'Target Correlation',corrLim)
hline(0,'k:')
vline(0,'k:')
figAxForm
subplot 325   
h = plot(targetWeightsAll_aud,detectWeightsAll_aud,'ko');
h.MarkerFaceColor = 'k';
hold on
figXAxis([],'Target LR Weight',weightLim)
figYAxis([],'Detect LR Weight',weightLim)
hline(0,'k:')
vline(0,'k:')
figAxForm
subplot 326   
h = plot(targetWeightsAll_shuf_aud,detectWeightsAll_shuf_aud,'ko');
h.MarkerFaceColor = 'k';
hold on
figXAxis([],'Target LR Weight',weightLim)
figYAxis([],'Detect LR Weight',weightLim)
hline(0,'k:')
vline(0,'k:')
figAxForm


subplot 322
p = nan(1,2);
for i = 1:2
    x = mean(targetWeightsBinned_aud{i});
    xerr = ste(targetWeightsBinned_aud{i},1);
    y = detectWeightsBinned_aud{i};
    yerr = ste(y,1);
    [ttest_pass,p(i)] = ttest(y,[],'alpha',0.025);
    plot(x,y,'k.');
    hold on
    h = errorbar(x,mean(y),yerr,yerr,xerr,xerr,'ko');
    if ttest_pass
        h.Color = 'r';
    end
end
figXAxis([],'Target LR Weight',weightLim)
figYAxis([],'Detect LR Weight',weightLim)
hline(0,'k:')
vline(0,'k:')
figAxForm
title(sprintf('hi/lo target, p = %s/%s',round(p,2,'significant')))

subplot 324
p = nan(1,2);
for i = 1:2
    x = mean(targetWeightsBinned_shuff_aud{i});
    xerr = ste(targetWeightsBinned_shuff_aud{i},1);
    y = detectWeightsBinned_shuff_aud{i};
    yerr = ste(y,1);
    [ttest_pass,p(i)] = ttest(y,[],'alpha',0.025);
    plot(x,y,'k.');
    hold on
    h = errorbar(x,mean(y),yerr,yerr,xerr,xerr,'ko');
    if ttest_pass
        h.Color = 'r';
    end
end
figXAxis([],'Target LR Weight',weightLim)
figYAxis([],'Detect LR Weight',weightLim)
hline(0,'k:')
vline(0,'k:')
figAxForm
title(sprintf('SHUFFLE: hi/lo target, p = %s/%s',round(p,2,'significant')))

print(fullfile(fnout,'corrXWeight_shuffle_aud'),'-dpdf','-fillpage')

%% plot orientation ...
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
figYAxis([],'Weight',weightLim_ori)
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
figYAxis([],'Weight',weightLim_ori)
figAxForm
title('Detect')
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
figYAxis([],'Weight',weightLim_ori);
figAxForm
legendStr = strcat({'Target ';'Detect '},...
    {num2str(round(p_anova_target_aud,2,'significant'));...
    num2str(round(p_anova_detect_aud,2,'significant'))});
legend(legendStr,'location','northeast')

subplot 323
h = plot(oriPrefAll,targetWeightsAll_aud,'o');
h.Color = 'k';
figXAxis([],'Orientation Pref (deg)',[-10 190],0:45:180,0:45:180)
figYAxis([],'Weight',weightLim)
figAxForm
title('Target')
subplot 324
h = plot(oriPrefAll,detectWeightsAll_aud,'o');
h.Color = 'k';
figXAxis([],'Orientation Pref (deg)',[-10 190],0:45:180,0:45:180)
figYAxis([],'Weight',weightLim)
figAxForm
title('Detect')

subplot 325
h = polarplot(deg2rad(oriPrefAll),targetWeightsAll_aud','bo');
h.Parent.ThetaTick = [0:45:360];
h.Parent.RLim = polarWeightLim;
h.Parent.RTick = [0:2:4];
title('Target Weight')
subplot 326
h = polarplot(deg2rad(oriPrefAll),detectWeightsAll_aud','bo');
h.Parent.ThetaTick = [0:45:360];
h.Parent.RLim = polarWeightLim;
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
figYAxis([],'Weight',weightLim_ori)
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
figYAxis([],'Weight',weightLim_ori)
figAxForm
title('Target')
legend(respCellName)

print(fullfile(fnout,'tuningXweight_respCellTypes_aud'),'-dpdf','-fillpage')

%% task selectivity, orientation tuning, signal analysis

targetSN = abs(targetAuroc - 0.5);
detectSN = abs(detectAuroc - 0.5);

setFigParams4Print('landscape')
set(0,'defaultAxesFontSize',12)
figure
suptitle(sprintf(...
    'Task Preferences are sorted into zero, hard (<%s), and easy targets',...
    num2str(orientationEdges(2))))
for iori = 1:nBinnedOris
    subplot(3,nBinnedOris,iori)
    ind = taskOriPref == orientationEdges(iori);
    hold on
    histogram(oriPrefAll(ind),[0:20:180]);
    title(sprintf('Task Ori Pref: = %s',num2str(orientationEdges(iori))))
    figXAxis([],'Passive Ori Pref',[-1 181])
    figYAxis([], 'N Cells',[])
    figAxForm
end

subplot(3,nBinnedOris,4)
h = plot(oriPrefAll(respCells_firstBase),...
    targetSN(respCells_firstBase),'bo');
hold on
h = plot(oriPrefAll(respCells_targetAndBase),...
    targetSN(respCells_targetAndBase),'mo');
h = plot(oriPrefAll(respCells_targetOnly),...
    targetSN(respCells_targetOnly),'ro');
title('Target vs Distractor')
figXAxis([],'Passive Ori Pref',[-1 181])
figYAxis([], 'Target Sig:Noise',[0 0.5])
figAxForm
subplot(3,nBinnedOris,5)
h = plot(oriPrefAll(respCells_firstBase),...
    detectSN(respCells_firstBase),'bo');
hold on
h = plot(oriPrefAll(respCells_targetAndBase),...
    detectSN(respCells_targetAndBase),'mo');
h = plot(oriPrefAll(respCells_targetOnly),...
    detectSN(respCells_targetOnly),'ro');
title('Yes vs. No')
figXAxis([],'Passive Ori Pref',[-1 181])
figYAxis([], 'Choice Sig:Noise',[0 0.5])
figAxForm

subplot(3,nBinnedOris,7)
h = plot(taskOriPref(respCells_firstBase),...
    targetSN(respCells_firstBase),'bo');
hold on
h = plot(taskOriPref(respCells_targetAndBase)+3,...
    targetSN(respCells_targetAndBase),'mo');
h = plot(taskOriPref(respCells_targetOnly)+6,...
    targetSN(respCells_targetOnly),'ro');
title('Target vs Distractor')
figXAxis([],'Task Ori Pref',[-1 181])
figYAxis([], 'Target Sig:Noise',[0 0.5])
figAxForm
subplot(3,nBinnedOris,8)
h = plot(taskOriPref(respCells_firstBase),...
    detectSN(respCells_firstBase),'bo');
hold on
h = plot(taskOriPref(respCells_targetAndBase)+3,...
    detectSN(respCells_targetAndBase),'mo');
h = plot(taskOriPref(respCells_targetOnly)+6,...
    detectSN(respCells_targetOnly),'ro');
title('Yes vs. No')
figXAxis([],'Task Ori Pref',[-1 181])
figYAxis([], 'Choice Sig:Noise',[0 0.5])
figAxForm

subplot(3,nBinnedOris,6)
h = plot(oriPrefAll(respCells_firstBase),...
    si_resp(respCells_firstBase),'bo');
hold on
h = plot(oriPrefAll(respCells_targetAndBase),...
    si_resp(respCells_targetAndBase),'mo');
h = plot(oriPrefAll(respCells_targetOnly),...
    si_resp(respCells_targetOnly),'ro');
title('Anticipation Selectivity')
figXAxis([],'Passive Ori Pref',[0 181])
figYAxis([], 'V-A Selectivity',[-12 12])
figAxForm
subplot(3,nBinnedOris,9)
h = plot(taskOriPref(respCells_firstBase),...
    si_resp(respCells_firstBase),'bo');
hold on
h = plot(taskOriPref(respCells_targetAndBase)+3,...
    si_resp(respCells_targetAndBase),'mo');
h = plot(taskOriPref(respCells_targetOnly)+6,...
    si_resp(respCells_targetOnly),'ro');
title('Anticipation Selectivity')
figXAxis([],'Task Ori Pref',[0 181])
figYAxis([], 'V-A Selectivity',[-12 12])
figAxForm

print(fullfile(fnout,'tuningXsnXsiXtaskresp_vis'),'-dpdf','-fillpage')

%% weight by task tuning
setFigParams4Print('portrait')
set(0,'defaultAxesFontSize',12)
figure

subplot 221
% h = plot(taskOriPref,targetWeightsAll_vis,'o');
% h.Color = 'k';
hold on
h = plot(taskOriPref(respCells_firstBase),...
    targetWeightsAll_vis(respCells_firstBase),'.');
h.Color = 'b';
h = plot(taskOriPref(respCells_targetAndBase),...
    targetWeightsAll_vis(respCells_targetAndBase),'.');
h.Color = 'm';
h = plot(taskOriPref(respCells_targetOnly),...
    targetWeightsAll_vis(respCells_targetOnly),'.');
h.Color = 'r';
figXAxis([],'Orientation Pref (deg)',[-10 190],orientationEdges,orientationEdges)
figYAxis([],'Weight',weightLim)
figAxForm
title('Target')

subplot 222
% h = plot(taskOriPref,detectWeightsAll_vis,'o');
% h.Color = 'k';
hold on
h = plot(taskOriPref(respCells_firstBase),...
    detectWeightsAll_vis(respCells_firstBase),'.');
h.Color = 'b';
h = plot(taskOriPref(respCells_targetAndBase),...
    detectWeightsAll_vis(respCells_targetAndBase),'.');
h.Color = 'm';
h = plot(taskOriPref(respCells_targetOnly),...
    detectWeightsAll_vis(respCells_targetOnly),'.');
h.Color = 'r';
figXAxis([],'Orientation Pref (deg)',[-10 190],orientationEdges,orientationEdges)
figYAxis([],'Weight',weightLim)
figAxForm
title('Detect')

subplot 223
% h = plot(taskOriPref,targetWeightsAll_vis,'o');
% h.Color = 'k';
hold on
h = plot(ones(sum(respCells_firstBase),1),...
    targetWeightsAll_vis(respCells_firstBase),'.');
h.Color = 'b';
h = plot(ones(sum(respCells_targetAndBase),1)*2,...
    targetWeightsAll_vis(respCells_targetAndBase),'.');
h.Color = 'm';
h = plot(ones(sum(respCells_targetOnly),1)*3,...
    targetWeightsAll_vis(respCells_targetOnly),'.');
h.Color = 'r';
% figXAxis([],'Orientation Pref (deg)',[-10 190],orientationEdges,orientationEdges)
% figYAxis([],'Weight',weightLim)
% figAxForm
% title('Target')
for i = 1:3
    x = 1:nBinnedOris;
    y = targetWeightTaskOriBin_vis_respCells{i};
    hold on
    h = plot(x,y,'o-');
    h.Color = respCellColor{i};
end
hline(0,'k--')
figXAxis([],'Orientation Pref (deg)',[0 nBinnedOris+1],1:nBinnedOris,orientationEdges)
figYAxis([],'Weight',weightLim)
figAxForm
title('Target')
subplot 224
hold on
h = plot(ones(sum(respCells_firstBase),1),...
    detectWeightsAll_vis(respCells_firstBase),'.');
h.Color = 'b';
h = plot(ones(sum(respCells_targetAndBase),1)*2,...
    detectWeightsAll_vis(respCells_targetAndBase),'.');
h.Color = 'm';
h = plot(ones(sum(respCells_targetOnly),1)*3,...
    detectWeightsAll_vis(respCells_targetOnly),'.');
h.Color = 'r';
for i = 1:3
    x = 1:nBinnedOris;
    y = detectWeightTaskOriBin_vis_respCells{i};
    hold on
    h = plot(x,y,'o-');
    h.Color = respCellColor{i};
end
hline(0,'k--')
figXAxis([],'Orientation Pref (deg)',[0 nBinnedOris+1],1:nBinnedOris,orientationEdges)
figYAxis([],'Weight',weightLim)
figAxForm
title('Detect')
legend(respCellName)

print(fullfile(fnout,'tasktuningXweight_vis'),'-dpdf','-fillpage')

%%
setFigParams4Print('portrait')
set(0,'defaultAxesFontSize',12)
figure
subplot 221
title('Train Visual, Target')
hold on
for i = 1:nexp
    x = [1,2];
    y = [pctCorrTarget_visXori(i,end),pctCorr_target_testAudXamp(i,end)].*100;
    h = plot(x,y,'-');
    h.Color = [0.5 0.5 0.5];
end
y = mean([pctCorrTarget_visXori(:,end),pctCorr_target_testAudXamp(:,end)].*100,1);
yerr = ste([pctCorrTarget_visXori(:,end),pctCorr_target_testAudXamp(:,end)].*100,1);
h = errorbar(x,y,yerr,'ko-');
h.MarkerFaceColor = [1 1 1];
figXAxis([],'Test',[0 3],x,{'Visual';'Auditory'})
figYAxis([],'% Correct',[0 110])
figAxForm
hline(50,'k--')
subplot 222
title('Train Visual, Detect')
hold on
for i = 1:nexp
    x = [1,2];
    y = [pctCorrDetect_visXori(i,end),pctCorr_detect_testAudXamp(i,end)]*100;
    h = plot(x,y,'-');
    h.Color = [0.5 0.5 0.5];
end
y = mean([pctCorrDetect_visXori(:,end),pctCorr_detect_testAudXamp(:,end)].*100,1);
yerr = ste([pctCorrDetect_visXori(:,end),pctCorr_detect_testAudXamp(:,end)].*100,1);
h = errorbar(x,y,yerr,'ko-');
h.MarkerFaceColor = [1 1 1];
figXAxis([],'Test',[0 3],x,{'Visual';'Auditory'})
figYAxis([],'% Correct',[0 110])
figAxForm
hline(50,'k--')

subplot 223
title('Train Auditory, Target')
hold on
for i = 1:nexp
    x = [1,2];
    y = [pctCorr_target_testVisXori(i,end),pctCorrTarget_audXamp(i,end)].*100;
    h = plot(x,y,'-');
    h.Color = [0.5 0.5 0.5];
end
y = mean([pctCorr_target_testVisXori(:,end),pctCorrTarget_audXamp(:,end)].*100,1);
yerr = ste([pctCorr_target_testVisXori(:,end),pctCorrTarget_audXamp(:,end)].*100,1);
h = errorbar(x,y,yerr,'ko-');
h.MarkerFaceColor = [1 1 1];
figXAxis([],'Test',[0 3],x,{'Visual';'Auditory'})
figYAxis([],'% Correct',[0 110])
figAxForm
hline(50,'k--')
subplot 224
title('Train Auditory, Detect')
hold on
for i = 1:nexp
    x = [1,2];
    y = [pctCorr_detect_testVisXori(i,end),pctCorrDetect_audXamp(i,end)]*100;
    h = plot(x,y,'-');
    h.Color = [0.5 0.5 0.5];
end
y = mean([pctCorr_detect_testVisXori(:,end),pctCorrDetect_audXamp(:,end)].*100,1);
yerr = ste([pctCorr_detect_testVisXori(:,end),pctCorrDetect_audXamp(:,end)].*100,1);
h = errorbar(x,y,yerr,'ko-');
h.MarkerFaceColor = [1 1 1];
figXAxis([],'Test',[0 3],x,{'Visual';'Auditory'})
figYAxis([],'% Correct',[0 110])
figAxForm
hline(50,'k--')

print(fullfile(fnout,'crossModalTest'),'-dpdf','-fillpage')

%%
figure
suptitle('Weights and corrlations done with PCs')

subplot 221
hold on
x = targetCorrsAll_vis;
y = targetCorrsAll_aud;
h = scatter(x,y,'ko');
h.MarkerFaceColor = 'k';
linFitCoeffs = polyfit(x,y,1);
yfit = polyval(linFitCoeffs,x);
R_squared = 1 - (sum((y-yfit).^2)/sum((y - mean(y)).^2));
h = plot(x,yfit,'r-');
plot(corrLim,corrLim,'k--')
figXAxis([],'Visual',corrLim);
figYAxis([],'Auditory',corrLim);
figAxForm
title(sprintf('Target Correlation,R^2 = %s',num2str(R_squared)))

subplot 222
hold on
x = detectCorrsAll_vis;
y = detectCorrsAll_aud;
h = scatter(x,y,'ko');
h.MarkerFaceColor = 'k';
linFitCoeffs = polyfit(x,y,1);
yfit = polyval(linFitCoeffs,x);
R_squared = 1 - (sum((y-yfit).^2)/sum((y - mean(y)).^2));
h = plot(x,yfit,'r-');
plot(corrLim,corrLim,'k--')
figXAxis([],'Visual',corrLim);
figYAxis([],'Auditory',corrLim);
figAxForm
title(sprintf('Detect Correlation, R^2 = %s',num2str(R_squared)))

subplot 223
hold on
x = targetWeightsAll_vis;
y = targetWeightsAll_aud;
h = scatter(x,y,'ko');
h.MarkerFaceColor = 'k';
linFitCoeffs = polyfit(x,y,1);
yfit = polyval(linFitCoeffs,x);
R_squared = 1 - (sum((y-yfit).^2)/sum((y - mean(y)).^2));
[R,p]= corrcoef(x,y);
R = round(R(1,2),2,'significant');
p = round(p(1,2),2,'significant');
h = plot(x,yfit,'r-');
plot(weightLim,weightLim,'k--')
figXAxis([],'Visual',weightLim);
figYAxis([],'Auditory',weightLim);
figAxForm
title(sprintf('Target Weight, R^2 = %s, R=%s,p=%s',num2str(R_squared),num2str(R),num2str(p)))

subplot 224
hold on
x = detectWeightsAll_vis;
y = detectWeightsAll_aud;
h = scatter(x,y,'ko');
h.MarkerFaceColor = 'k';
linFitCoeffs = polyfit(x,y,1);
yfit = polyval(linFitCoeffs,x);
R_squared = 1 - (sum((y-yfit).^2)/sum((y - mean(y)).^2));
[R,p]= corrcoef(x,y);
R = round(R(1,2),2,'significant');
p = round(p(1,2),2,'significant');
h = plot(x,yfit,'r-');
plot(weightLim,weightLim,'k--')
figXAxis([],'Visual',weightLim);
figYAxis([],'Auditory',weightLim);
figAxForm
title(sprintf('Detect Weight, R^2 = %s, R=%s,p=%s',num2str(R_squared),num2str(R),num2str(p)))

print(fullfile(fnout,'crossModalCorrelations'),'-dpdf','-fillpage')

figure
weightLim_shuf = [-1 1];
subplot 223
hold on
x = targetWeightsAll_shuf_vis;
y = targetWeightsAll_shuf_aud;
h = scatter(x,y,'ko');
h.MarkerFaceColor = 'k';
linFitCoeffs = polyfit(x,y,1);
yfit = polyval(linFitCoeffs,x);
R_squared = 1 - (sum((y-yfit).^2)/sum((y - mean(y)).^2));
[R,p]= corrcoef(x,y);
R = round(R(1,2),2,'significant');
p = round(p(1,2),2,'significant');
h = plot(x,yfit,'r-');
plot(weightLim_shuf,weightLim_shuf,'k--')
figXAxis([],'Visual',weightLim_shuf);
figYAxis([],'Auditory',weightLim_shuf);
figAxForm
title(sprintf('Target Weight, R^2 = %s, R=%s,p=%s',num2str(R_squared),num2str(R),num2str(p)))

subplot 224
hold on
x = detectWeightsAll_shuf_vis;
y = detectWeightsAll_shuf_aud;
h = scatter(x,y,'ko');
h.MarkerFaceColor = 'k';
linFitCoeffs = polyfit(x,y,1);
yfit = polyval(linFitCoeffs,x);
R_squared = 1 - (sum((y-yfit).^2)/sum((y - mean(y)).^2));
[R,p]= corrcoef(x,y);
R = round(R(1,2),2,'significant');
p = round(p(1,2),2,'significant');
h = plot(x,yfit,'r-');
plot(weightLim_shuf,weightLim_shuf,'k--')
figXAxis([],'Visual',weightLim_shuf);
figYAxis([],'Auditory',weightLim_shuf);
figAxForm
title(sprintf('Dectect Weight, R^2 = %s, R=%s,p=%s',num2str(R_squared),num2str(R),num2str(p)))

print(fullfile(fnout,'crossModalCorrelations_shuffle'),'-dpdf','-fillpage')
%% percent correct and cross valiation by orienation or amplitude

pctCorrOri_detect_sub = pctCorrDetect_visXori - ...
    pctCorrDetect_visXori_XVal;
pctCorrOri_target_sub = pctCorrTarget_visXori - ...
    pctCorrTarget_visXori_XVal;

set(0,'defaultAxesFontSize',16)
figure
suptitle('Summary Across Experiments')
subplot 231
for i = 1:nexp
    x = 1:4;
    y = pctCorrDetect_visXori(i,:).*100;
    h = plot(x,y,'-');
    h.Color = [0.5 0.5 0.5];
    hold on
end
h = errorbar(x,nanmean(pctCorrDetect_visXori.*100,1),...
    ste(pctCorrDetect_visXori.*100,1),'ko-');
h.MarkerFaceColor = [1 1 1];
xLabel = [cellfun(@(a)...
    num2str(round(a,2,'significant')),num2cell(orientationEdges),'unif',0),...
    {'All Trials'}];
figXAxis([],'Orientation (deg)',[0 5],1:4,xLabel)
figYAxis([],'% Correct',[0 110])
figAxForm
hline(50,'k--')
title('Visual Detect - Training')

subplot 232
for i = 1:nexp
    x = 1:4;
    y = pctCorrDetect_visXori_XVal(i,:).*100;
    h = plot(x,y,'-');
    h.Color = [0.5 0.5 0.5];
    hold on
end
h = errorbar(x,nanmean(pctCorrDetect_visXori_XVal.*100,1),...
    ste(pctCorrDetect_visXori_XVal.*100,1),'ko-');
h.MarkerFaceColor = [1 1 1];
xLabel = [cellfun(@(a)...
    num2str(round(a,2,'significant')),num2cell(orientationEdges),'unif',0),...
    {'All Trials'}];
figXAxis([],'Orientation (deg)',[0 5],1:4,xLabel)
figYAxis([],'% Correct',[0 110])
figAxForm
hline(50,'k--')
title('Visual Detect - Hold-out')

subplot 233
for i = 1:nexp
    x = 1:4;
    y = pctCorrOri_detect_sub(i,:).*100;
    h = plot(x,y,'-');
    h.Color = [0.5 0.5 0.5];
    hold on
end
h = errorbar(x,nanmean(pctCorrOri_detect_sub.*100,1),...
    ste(pctCorrOri_detect_sub.*100,1),'ko-');
h.MarkerFaceColor = [1 1 1];
xLabel = [cellfun(@(a)...
    num2str(round(a,2,'significant')),num2cell(orientationEdges),'unif',0),...
    {'All Trials'}];
figXAxis([],'Orientation (deg)',[0 5],1:4,xLabel)
figYAxis([],'Train - Test (% Correct)',[-20 20])
figAxForm
hline(0,'k--')
title('Visual Detect - Diff')

subplot 234
% allExpt = nan(length(allOris),nexp);
for i = 1:nexp
    x = 1:4;
    y = pctCorrTarget_visXori(i,:).*100;
    h = plot(x,y,'-');
    h.Color = [0.5 0.5 0.5];
    hold on
end
h = errorbar(x,nanmean(pctCorrTarget_visXori.*100,1),...
    ste(pctCorrTarget_visXori.*100,1),'ko-');
h.MarkerFaceColor = [1 1 1];
xLabel = [cellfun(@(a)...
    num2str(round(a,2,'significant')),num2cell(orientationEdges),'unif',0),...
    {'All Trials'}];
figXAxis([],'Orientation (deg)',[0 5],1:4,xLabel)
figYAxis([],'% Correct',[0 110])
figAxForm
hline(50,'k--')
title('Visual Target Model')

subplot 235
for i = 1:nexp
    x = 1:4;
    y = pctCorrTarget_visXori_XVal(i,:).*100;
    h = plot(x,y,'-');
    h.Color = [0.5 0.5 0.5];
    hold on
end
h = errorbar(x,nanmean(pctCorrTarget_visXori_XVal.*100,1),...
    ste(pctCorrTarget_visXori_XVal.*100,1),'ko-');
h.MarkerFaceColor = [1 1 1];
xLabel = [cellfun(@(a)...
    num2str(round(a,2,'significant')),num2cell(orientationEdges),'unif',0),...
    {'All Trials'}];
figXAxis([],'Orientation (deg)',[0 5],1:4,xLabel)
figYAxis([],'% Correct',[0 110])
figAxForm
hline(50,'k--')
title('Visual Target - Hold-out')

subplot 236
for i = 1:nexp
    x = 1:4;
    y = pctCorrOri_target_sub(i,:).*100;
    h = plot(x,y,'-');
    h.Color = [0.5 0.5 0.5];
    hold on
end
h = errorbar(x,nanmean(pctCorrOri_target_sub.*100,1),...
    ste(pctCorrOri_target_sub.*100,1),'ko-');
h.MarkerFaceColor = [1 1 1];
xLabel = [cellfun(@(a)...
    num2str(round(a,2,'significant')),num2cell(orientationEdges),'unif',0),...
    {'All Trials'}];
figXAxis([],'Orientation (deg)',[0 5],1:4,xLabel)
figYAxis([],'Train - Test (% Correct)',[-20 20])
figAxForm
hline(0,'k--')
title('Visual Target - Diff')

print(fullfile(fnout,'pctCorr_visModel'),'-dpdf','-fillpage')


pctCorrOri_detect_sub = pctCorrDetect_audXamp- ...
    pctCorrDetect_audXamp_XVal;
pctCorrOri_target_sub = pctCorrTarget_audXamp - ...
    pctCorrTarget_audXamp_XVal;

set(0,'defaultAxesFontSize',16)
figure
suptitle('Summary Across Experiments')
subplot 231
for i = 1:nexp
    x = 1:4;
    y = pctCorrDetect_audXamp(i,:).*100;
    ind = ~isnan(y);
    h = plot(x(ind),y(ind),'-');
    h.Color = [0.5 0.5 0.5];
    hold on
end
h = errorbar(x,nanmean(pctCorrDetect_audXamp.*100,1),...
    ste(pctCorrDetect_audXamp.*100,1),'ko-');
h.MarkerFaceColor = [1 1 1];
xLabel = [cellfun(@(a)...
    num2str(round(a,2,'significant')),num2cell(ampEdges),'unif',0),...
    {'All Trials'}];
figXAxis([],'Volume (% Max)',[0 5],1:4,xLabel)
figYAxis([],'% Correct',[0 110])
figAxForm
hline(50,'k--')
title('Auditory Detect - Training')

subplot 232
for i = 1:nexp
    x = 1:4;
    y = pctCorrDetect_audXamp_XVal(i,:).*100;
    ind = ~isnan(y);
    h = plot(x(ind),y(ind),'-');
    h.Color = [0.5 0.5 0.5];
    hold on
end
h = errorbar(x,nanmean(pctCorrDetect_audXamp_XVal.*100,1),...
    ste(pctCorrDetect_audXamp_XVal.*100,1),'ko-');
h.MarkerFaceColor = [1 1 1];
xLabel = [cellfun(@(a)...
    num2str(round(a,2,'significant')),num2cell(ampEdges),'unif',0),...
    {'All Trials'}];
figXAxis([],'Volume (% Max)',[0 5],1:4,xLabel)
figYAxis([],'% Correct',[0 110])
figAxForm
hline(50,'k--')
title('Auditory Detect - Hold-out')

subplot 233
for i = 1:nexp
    x = 1:4;
    y = pctCorrOri_detect_sub(i,:).*100;
    ind = ~isnan(y);
    h = plot(x(ind),y(ind),'-');
    h.Color = [0.5 0.5 0.5];
    hold on
end
h = errorbar(x,nanmean(pctCorrOri_detect_sub.*100,1),...
    ste(pctCorrOri_detect_sub.*100,1),'ko-');
h.MarkerFaceColor = [1 1 1];
xLabel = [cellfun(@(a)...
    num2str(round(a,2,'significant')),num2cell(ampEdges),'unif',0),...
    {'All Trials'}];
figXAxis([],'Volume (% Max)',[0 5],1:4,xLabel)
figYAxis([],'Train - Test (% Correct)',[-20 20])
figAxForm
hline(0,'k--')
title('Auditory Detect - Diff')

subplot 234
% allExpt = nan(length(allOris),nexp);
for i = 1:nexp
    x = 1:4;
    y = pctCorrTarget_audXamp(i,:).*100;
    ind = ~isnan(y);
    h = plot(x(ind),y(ind),'-');
    h.Color = [0.5 0.5 0.5];
    hold on
end
h = errorbar(x,nanmean(pctCorrTarget_audXamp.*100,1),...
    ste(pctCorrTarget_audXamp.*100,1),'ko-');
h.MarkerFaceColor = [1 1 1];
xLabel = [cellfun(@(a)...
    num2str(round(a,2,'significant')),num2cell(ampEdges),'unif',0),...
    {'All Trials'}];
figXAxis([],'Volume (% Max)',[0 5],1:4,xLabel)
figYAxis([],'% Correct',[0 110])
figAxForm
hline(50,'k--')
title('Auditory Target Model')

subplot 235
for i = 1:nexp
    x = 1:4;
    y = pctCorrTarget_audXamp_XVal(i,:).*100;
    ind = ~isnan(y);
    h = plot(x(ind),y(ind),'-');
    h.Color = [0.5 0.5 0.5];
    hold on
end
h = errorbar(x,nanmean(pctCorrTarget_audXamp_XVal.*100,1),...
    ste(pctCorrTarget_audXamp_XVal.*100,1),'ko-');
h.MarkerFaceColor = [1 1 1];
xLabel = [cellfun(@(a)...
    num2str(round(a,2,'significant')),num2cell(ampEdges),'unif',0),...
    {'All Trials'}];
figXAxis([],'Volume (% Max)',[0 5],1:4,xLabel)
figYAxis([],'% Correct',[0 110])
figAxForm
hline(50,'k--')
title('Auditory Target - Hold-out')

subplot 236
for i = 1:nexp
    x = 1:4;
    y = pctCorrOri_detect_sub(i,:).*100;
    ind = ~isnan(y);
    h = plot(x(ind),y(ind),'-');
    h.Color = [0.5 0.5 0.5];
    hold on
end
h = errorbar(x,nanmean(pctCorrOri_detect_sub.*100,1),...
    ste(pctCorrOri_detect_sub.*100,1),'ko-');
h.MarkerFaceColor = [1 1 1];
xLabel = [cellfun(@(a)...
    num2str(round(a,2,'significant')),num2cell(ampEdges),'unif',0),...
    {'All Trials'}];
figXAxis([],'Volume (% Max)',[0 5],1:4,xLabel)
figYAxis([],'Train - Test (% Correct)',[-20 20])
figAxForm
hline(0,'k--')
title('Auditory Target - Diff')

print(fullfile(fnout,'pctCorr_audModel'),'-dpdf','-fillpage')