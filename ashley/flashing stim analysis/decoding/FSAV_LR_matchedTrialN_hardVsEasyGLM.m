close all
clear all

bxParams_FSAV

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
oriThreshold = 45;

oriBinsFine = [0 1 16 22.5 32 45 64 90];
oriBins = [0 1 32 90];
ampBinsFine = [0 exp(linspace(log(0.0001),log(1),7))];
% ampBins = [0 exp(linspace(log(0.0001),log(1),3))];
ampBins = [0 0.0001 0.1 1];
weightsDiffZeroAlpha = 0.05/2;
nPC = 10;
maxCells = 15;
weightLim = [-2.2 5];
corrLim = [-1 1];
HRLim = [-10 110];
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

% targetWeightsAll_aud = [];
% detectWeightsAll_aud = [];
% 
% targetCorrsAll_aud = [];
% detectCorrsAll_aud = [];

nCells = zeros(1,nexp);
nTrialsExpt = zeros(2,nexp);
% nTrialsExpt_aud = zeros(1,nexp);
% pcNumber = [];
exptName = cell(1,nexp);

% pctCorrTarget_visXori = nan(nexp,length(allOris)+1);
% pctCorrDetect_visXori = nan(nexp,length(allOris)+1);
% pctCorrTarget_visXori_XVal = nan(nexp,length(allOris)+1);
% pctCorrDetect_visXori_XVal = nan(nexp,length(allOris)+1);

pctCorrTarget = nan(nexp,2);
pctCorrDetect = nan(nexp,2);
pctCorrTarget_XValGroup = nan(nexp,2);
pctCorrDetect_XValGroup = nan(nexp,2);
pctCorrTarget_XValHO = nan(nexp,2);
pctCorrDetect_XValHO = nan(nexp,2);

probEaTrial_target = cell(1,2);
probEaTrial_detect = cell(1,2);
oriEaTrial = cell(1,2);

% pctCorrTarget_audXamp = nan(nexp,length(allAmps)+1);
% pctCorrDetect_audXamp = nan(nexp,length(allAmps)+1);
% pctCorrTarget_audXamp_XVal = nan(nexp,length(allAmps)+1);
% pctCorrDetect_audXamp_XVal = nan(nexp,length(allAmps)+1);

for iexp=1:nexp
    
   disp([dcExpt(iexp).mouse '-' dcExpt(iexp).date])
    %****VISUAL TRIALS
    
   cellIdx = dcExpt(iexp).signifResponsiveCells' & ...
       dcExpt(iexp).oriFitTheta90 < theta90Threshold;
    nCells(iexp) = sum(cellIdx);
    
    if nCells(iexp) > maxCells
        ind = find(cellIdx);
        cellSample = randsample(ind,maxCells);
        ind = false(1,length(cellIdx));
        ind(cellSample) = true;
        cellIdx = ind;
        nCells(iexp) = sum(cellIdx);
    end
    
%    trOut=dcExpt(iexp).trialOutcome;
    tOri = round(dcExpt(iexp).trialOrientation,2,'significant');
%     tOri = binAndRelabelTrialOrientation(tOri,oriBinsFine);
    orientations = unique(tOri);
    nOri = length(orientations);
%     nt = length(tOri);
    
%     ind1 = strcmp(trOut,'fa');
%     ind2 = strcmp(trOut,'cr');    
%     ind3 = find(ind1 |ind2);
%     randZeroTrialsInd1 = [randsample(find(ind1), round(sum(ind1)./2)) ...
%         randsample(find(ind2), round(sum(ind2)./2))];
%     randZeroTrialsInd2 = ind3(~ismember(find(ind3), randZeroTrialsInd1));
%     easyInd = false(1,nt);
%     easyInd(randZeroTrialsInd1) = true;
    easyInd = tOri >= oriThreshold | tOri == 0;
%     hardInd = false(1,nt);
%     hardInd(randZeroTrialsInd2) = true;
    hardInd = tOri < oriThreshold | tOri == 0;
   nTrialsExpt(1,iexp) = sum(easyInd);
   nTrialsExpt(2,iexp) = sum(hardInd);
    
    %***Easy Trials***
   trOut=dcExpt(iexp).trialOutcome;
    trOut = trOut(easyInd);
    [Y_yeses, Y_targets] = getStimAndBehaviorYs(trOut);
   
   X_targets=dcExpt(iexp).stimResp(cellIdx,easyInd)';
    
   X_targets = zscore(X_targets);
    
   detectCorr = corr(Y_yeses,X_targets)';
   targetCorr = corr(Y_targets,X_targets)';
   
   detectCorrsAll_vis = cat(1,detectCorrsAll_vis,detectCorr);
   targetCorrsAll_vis = cat(1,targetCorrsAll_vis,targetCorr);
   
   C=eye(size(X_targets,2));
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
   
   dv_target = mean(Y_targets);
   dv_detect = mean(Y_yeses);

    yhat_target = glmval(targetGLM.beta,X_targets,'logit') > dv_target;
    pctCorr_target_vis = mean((Y_targets-yhat_target) == 0);
   
    yhat_detect = glmval(detectGLM.beta,X_targets,'logit') > dv_detect;
    pctCorr_detect_vis = mean((Y_yeses-yhat_detect) == 0);
    
    pctCorrTarget(iexp,1) = pctCorr_target_vis.*100;
    pctCorrDetect(iexp,1) = pctCorr_detect_vis.*100;
    
    probEaTrial_target{1} = cat(1,probEaTrial_target{1},...
        glmval(targetGLM.beta,X_targets,'logit'));
    probEaTrial_detect{1} = cat(1,probEaTrial_detect{1},...
        glmval(detectGLM.beta,X_targets,'logit'));
    oriEaTrial{1} = cat(2,oriEaTrial{1},tOri(easyInd));
    
   nt = length(trOut);
   ysub_target_ho = nan(1,nt);
   ysub_detect_ho = nan(1,nt);
   for i = 1:nt
       X_holdout = X_targets(i,:);
       othersInd = [1:(i-1),(i+1):nt];
       X_others = X_targets(othersInd,:);
       
       [~,~,targetGLM_ho]=glmfit(X_others,Y_targets(othersInd),'binomial');
       [~,~,detectGLM_ho]=glmfit(X_others,Y_yeses(othersInd),'binomial');
       
       y_holdout = Y_targets(i);
       yhat = glmval(targetGLM_ho.beta,X_holdout,'logit');
       ysub_target_ho(i) = y_holdout - (yhat > dv_target);
       
       y_holdout = Y_yeses(i);
       yhat = glmval(detectGLM_ho.beta,X_holdout,'logit');
       ysub_detect_ho(i) = y_holdout - (yhat > dv_detect);
   end   
    pctCorrTarget_XValHO(iexp,1) = mean(ysub_target_ho == 0).*100;
    pctCorrDetect_XValHO(iexp,1) = mean(ysub_detect_ho == 0).*100;
   
    %test other
    ind = hardInd & tOri ~= 0;
    trOut=dcExpt(iexp).trialOutcome;
    trOut = trOut(ind);
    [Y_yeses, Y_targets] = getStimAndBehaviorYs(trOut);
   
    X_targets=dcExpt(iexp).stimResp(cellIdx,ind)';
   
   yhat_target = glmval(targetGLM.beta,X_targets,'logit') > dv_target;
    pctCorr_target_vis = mean((Y_targets-yhat_target) == 0);
   
    yhat_detect = glmval(detectGLM.beta,X_targets,'logit') > dv_detect;
    pctCorr_detect_vis = mean((Y_yeses-yhat_detect) == 0); 
    
   pctCorrTarget_XValGroup(iexp,1) = pctCorr_target_vis.*100;
   pctCorrDetect_XValGroup(iexp,1) = pctCorr_detect_vis.*100;
   
   
   %***Hard Trials***
   trOut=dcExpt(iexp).trialOutcome;
    trOut = trOut(hardInd);
    [Y_yeses, Y_targets] = getStimAndBehaviorYs(trOut);
   
   X_targets=dcExpt(iexp).stimResp(cellIdx,hardInd)';
    
   X_targets = zscore(X_targets);
    
   detectCorr = corr(Y_yeses,X_targets)';
   targetCorr = corr(Y_targets,X_targets)';
   
   detectCorrsAll_vis = cat(1,detectCorrsAll_vis,detectCorr);
   targetCorrsAll_vis = cat(1,targetCorrsAll_vis,targetCorr);
   
   C=eye(size(X_targets,2));
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
   
   dv_target = mean(Y_targets);
   dv_detect = mean(Y_yeses);

    yhat_target = glmval(targetGLM.beta,X_targets,'logit') > dv_target;
    pctCorr_target_vis = mean((Y_targets-yhat_target) == 0);
   
    yhat_detect = glmval(detectGLM.beta,X_targets,'logit') > dv_detect;
    pctCorr_detect_vis = mean((Y_yeses-yhat_detect) == 0);
    
    pctCorrTarget(iexp,2) = pctCorr_target_vis.*100;
    pctCorrDetect(iexp,2) = pctCorr_detect_vis.*100;
    
    probEaTrial_target{2} = cat(1,probEaTrial_target{2},...
        glmval(targetGLM.beta,X_targets,'logit'));
    probEaTrial_detect{2} = cat(1,probEaTrial_detect{2},...
        glmval(detectGLM.beta,X_targets,'logit'));
    oriEaTrial{2} = cat(2,oriEaTrial{2},tOri(hardInd));
    
    nt = length(trOut);
   ysub_target_ho = nan(1,nt);
   ysub_detect_ho = nan(1,nt);
   for i = 1:nt
       X_holdout = X_targets(i,:);
       othersInd = [1:(i-1),(i+1):nt];
       X_others = X_targets(othersInd,:);
       
       [~,~,targetGLM_ho]=glmfit(X_others,Y_targets(othersInd),'binomial');
       [~,~,detectGLM_ho]=glmfit(X_others,Y_yeses(othersInd),'binomial');
       
       y_holdout = Y_targets(i);
       yhat = glmval(targetGLM_ho.beta,X_holdout,'logit');
       ysub_target_ho(i) = y_holdout - (yhat > dv_target);
       
       y_holdout = Y_yeses(i);
       yhat = glmval(detectGLM_ho.beta,X_holdout,'logit');
       ysub_detect_ho(i) = y_holdout - (yhat > dv_detect);
   end   
    pctCorrTarget_XValHO(iexp,2) = mean(ysub_target_ho == 0).*100;
    pctCorrDetect_XValHO(iexp,2) = mean(ysub_detect_ho == 0).*100;
    %test other
    ind = easyInd & tOri ~= 0;
    trOut=dcExpt(iexp).trialOutcome;
    trOut = trOut(ind);
    [Y_yeses, Y_targets] = getStimAndBehaviorYs(trOut);
   
    X_targets=dcExpt(iexp).stimResp(cellIdx,ind)';
   
   yhat_target = glmval(targetGLM.beta,X_targets,'logit') > dv_target;
    pctCorr_target_vis = mean((Y_targets-yhat_target) == 0);
   
    yhat_detect = glmval(detectGLM.beta,X_targets,'logit') > dv_detect;
    pctCorr_detect_vis = mean((Y_yeses-yhat_detect) == 0); 
    
   pctCorrTarget_XValGroup(iexp,2) = pctCorr_target_vis.*100;
   pctCorrDetect_XValGroup(iexp,2) = pctCorr_detect_vis.*100;
   
end

%%
setFigParams4Print('portrait')
set(0,'defaultAxesFontSize',12)
figure
suptitle(sprintf('Visual Trials; Easy trials: target >= %s deg,',num2str(oriThreshold)))
subplot 321
plot(pctCorrTarget(:,1),pctCorrTarget(:,2),'ko')
hold on
plot(HRLim,HRLim,'--')
figXAxis([],'Easy Trials (% corr)',HRLim)
figYAxis([],'Hard Trials (% corr)',HRLim)
figAxForm
title('Target: Train & test same data')
subplot 322
plot(pctCorrDetect(:,1),pctCorrDetect(:,2),'ko')
hold on
plot(HRLim,HRLim,'--')
figXAxis([],'Easy Trials (% corr)',HRLim)
figYAxis([],'Hard Trials (% corr)',HRLim)
figAxForm
title('Detect: Train & test same data')

subplot 323
plot(pctCorrTarget_XValHO(:,1),pctCorrTarget(:,1),'ko')
hold on
plot(pctCorrTarget_XValHO(:,2),pctCorrTarget(:,2),'bo')
plot(HRLim,HRLim,'--')
figXAxis([],'Test Holdout Data (% corr)',HRLim)
figYAxis([],'Test Training Data (% corr)',HRLim)
figAxForm
title('Target')
legend({'Easy';'Hard';'Unity'},'location','southeast')
subplot 324
plot(pctCorrDetect_XValHO(:,1),pctCorrDetect(:,1),'ko')
hold on
plot(pctCorrDetect_XValHO(:,2),pctCorrDetect(:,2),'bo')
plot(HRLim,HRLim,'--')
figXAxis([],'Test Holdout Data (% corr)',HRLim)
figYAxis([],'Test Training Data (% corr)',HRLim)
figAxForm
title('Detect')

subplot 325
plot(pctCorrTarget_XValGroup(:,1),pctCorrTarget_XValGroup(:,2),'ko')
hold on
plot(HRLim,HRLim,'--')
figXAxis([],'Train Easy Trials (% corr)',HRLim)
figYAxis([],'Train Hard Trials (% corr)',HRLim)
figAxForm
title('Target: Test opposite data')
subplot 326
plot(pctCorrDetect_XValGroup(:,1),pctCorrDetect_XValGroup(:,2),'ko')
hold on
plot(HRLim,HRLim,'--')
figXAxis([],'Train Easy Trials (% corr)',HRLim)
figYAxis([],'Train Hard Trials (% corr)',HRLim)
figAxForm
title('Detect: Test opposite data')
print(fullfile(fnout,'pctCorr_compareEasyHardTrialTrain'),'-dpdf','-fillpage')