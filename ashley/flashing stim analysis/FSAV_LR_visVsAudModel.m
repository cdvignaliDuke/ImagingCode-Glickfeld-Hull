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
weightsDiffZeroAlpha = 0.05/2;
%%
visTargetGLM = cell(1,nexp);
visDetectGLM = cell(1,nexp);
pctCorr_target_audTest = nan(1,nexp);
pctCorr_detect_audTest = nan(1,nexp);
pctCorr_target_visTest = nan(1,nexp);
pctCorr_detect_visTest = nan(1,nexp);
for iexp=1:nexp
   trOut=dcExpt(iexp).trialOutcome;
    [Y_yeses, Y_targets] = getStimAndBehaviorYs(trOut);

   X_vis=dcExpt(iexp).stimResp';
   
   cellIdx = dcExpt(iexp).signifResponsiveCells' & ...
       dcExpt(iexp).oriFitTheta90 < theta90Threshold;

   X_vis=X_vis(:,cellIdx);
   
   idx=find(~isnan(sum(X_vis,2)));
   Y_yeses=Y_yeses(idx,1);   
   Y_targets=Y_targets(idx,1);   
   X_vis=X_vis(idx,:);

   X_vis=bsxfun(@plus,X_vis,-mean(X_vis));
   X_vis=bsxfun(@times,X_vis,1./std(X_vis));
   
   X_aud = dcExpt(iexp).audStimResp';
   X_aud = X_aud(:,cellIdx);
   idx=find(~isnan(sum(X_aud,2)));
   X_aud = X_aud(idx,:);
   
   nTrialsVis = size(X_vis,1);
   nTrialsAud = size(X_aud,1);
 
   C=eye(size(X_vis,2));%C=inv(sqrtm(cov(X)));
   X_visTransform=X_vis*C;
   [~,~,visTargetGLM{iexp}]=glmfit(X_visTransform,Y_targets,'binomial');
   [~,~,visDetectGLM{iexp}]=glmfit(X_visTransform,Y_yeses,'binomial');
   
   % test auditory trials
   yhat_target = glmval(visTargetGLM{iexp}.beta,X_aud,'logit');
   pctCorr_target_audTest(iexp) = sum(yhat_target > decisionVariable)./nTrialsAud;
   
   yhat_detect = glmval(visDetectGLM{iexp}.beta,X_aud,'logit');
   pctCorr_detect_audTest(iexp) = sum(yhat_detect > decisionVariable)./nTrialsAud;
   
   
   trOut=dcExpt(iexp).audTrialOutcome;
    [Y_yeses, Y_targets] = getStimAndBehaviorYs(trOut);
   
   idx=find(~isnan(sum(X_aud,2)));
   Y_yeses=Y_yeses(idx,1);   
   Y_targets=Y_targets(idx,1);
   
   C=eye(size(X_aud,2));%C=inv(sqrtm(cov(X)));
   X_audTransform=X_aud*C;
   [~,~,audTargetGLM{iexp}]=glmfit(X_audTransform,Y_targets,'binomial');
   [~,~,audDetectGLM{iexp}]=glmfit(X_audTransform,Y_yeses,'binomial');

   
   % test visual trials
   yhat_target = glmval(audTargetGLM{iexp}.beta,X_vis,'logit');
   pctCorr_target_visTest(iexp) = sum(yhat_target > decisionVariable)./nTrialsVis;
   
   yhat_detect = glmval(audDetectGLM{iexp}.beta,X_vis,'logit');
   pctCorr_detect_visTest(iexp) = sum(yhat_detect > decisionVariable)./nTrialsVis;
   
end

%%
setFigParams4Print('landscape')
figure
suptitle('Cross-modality test of GLM, each expt')
subplot 121
x = ones(1,nexp);
y = pctCorr_target_audTest.*100;
h = scatter(x,y,'ko');
hold on
h = errorbar(1,mean(y),ste(y,2),'ko');
h.MarkerFaceColor = 'k';
x = ones(1,nexp).*2;
y = pctCorr_target_visTest.*100;
h = scatter(x,y,'ko');
h = errorbar(2,mean(y),ste(y,2),'ko');
h.MarkerFaceColor = 'k';
figXAxis([],'',[0 3],1:2,{'visGLM, audTest';'audGLM, visTest'})
figYAxis([],'% Target Present',[0 110])
figAxForm
title('Target GLM')

subplot 122
x = ones(1,nexp);
y = pctCorr_detect_audTest.*100;
h = scatter(x,y,'ko');
hold on
h = errorbar(1,mean(y),ste(y,2),'ko');
h.MarkerFaceColor = 'k';
x = ones(1,nexp).*2;
y = pctCorr_detect_visTest.*100;
h = scatter(x,y,'ko');
h = errorbar(2,mean(y),ste(y,2),'ko');
h.MarkerFaceColor = 'k';
figXAxis([],'',[0 3],1:2,{'visGLM, audTest';'audGLM, visTest'})
figYAxis([],'% Target Present',[0 110])
figAxForm
title('Detect GLM')

print(fullfile(fnout,'crossModalTest'),'-dpdf','-fillpage')