close all
clear all

load('Z:\Analysis\FSAV Choice\FSAV_respwinData.mat')
fnout = 'Z:\Analysis\FSAV Choice';
nexp = length(dcExpt);

doExptPlots = 0;
%%
theta90Threshold = 11.25;
% decisionVariable = 0.5;
detectThresholdHR = 0.8;
detectThresholdOri = 45;
oriBinsFine = [0 1 16 22.5 32 45 64 90];
oriBins = [0 1 45 90];
weightsDiffZeroAlpha = 0.05/2;
nRespwin = size(dcExpt(1).respwin,2);
optimalVisDelayInd = 6;
visRTRange = [250 570];
audRTRange = [150 450];
frameRateHz = 30;
nFrMonitorDelay_target = 2; %using cTargetOn
visStimOffsetMs = (nFrMonitorDelay_target/frameRateHz).*1000;


respwinMs = cell(1,nRespwin);
for rw = 1:nRespwin
    respwinMs{rw} = dcExpt(1).respwin(rw).winMs;
end
respwinMs = cell2mat(respwinMs');

dv_detect_vis = nan(1,nexp);
dv_target_vis = nan(1,nexp);
dv_detect_aud = nan(1,nexp);
dv_target_aud = nan(1,nexp);
for iexp = 1:nexp
    trOut=dcExpt(iexp).respwin(rw).trialOutcome;
    [Y_yeses, Y_targets] = getStimAndBehaviorYs(trOut);
    dv_detect_vis(iexp) = mean(Y_yeses);
    dv_target_vis(iexp) = mean(Y_targets);
    
    trOut=dcExpt(iexp).respwin(rw).audTrialOutcome;
    [Y_yeses, Y_targets] = getStimAndBehaviorYs(trOut);
    dv_detect_aud(iexp) = mean(Y_yeses);
    dv_target_aud(iexp) = mean(Y_targets);
end

figure
x = 1:2;
y = cat(2,mean(dv_detect_vis),mean(dv_detect_aud));
yerr = cat(2,ste(dv_detect_vis,2),ste(dv_detect_aud,2));
subplot 121
h = errorbar(x,y,yerr,'ko');
h.MarkerFaceColor = 'k';
figXAxis([],'',[0 3],x,{'Vis';'Aud'})
figYAxis([],'Fraction of Trials',[0 1])
figAxForm([],0)
title('Yes')

x = 1:2;
y = cat(2,mean(dv_target_vis),mean(dv_target_aud));
yerr = cat(2,ste(dv_target_vis,2),ste(dv_target_aud,2));
subplot 122
h = errorbar(x,y,yerr,'ko');
h.MarkerFaceColor = 'k';
figXAxis([],'',[0 3],x,{'Vis';'Aud'})
figYAxis([],'Fraction of Trials',[0 1])
figAxForm([],0)
title('Target Present')

print(fullfile(fnout,'decisionVariables'),'-dpdf')
%% visual trials
pctCorr_target_vis = nan(nexp,nRespwin);
pctCorr_detect_vis = nan(nexp,nRespwin);

pctCorr_target_audTest = nan(nexp,nRespwin);
pctCorr_detect_audTest = nan(nexp,nRespwin);

for iexp = 1:nexp
   for rw = 1:nRespwin
       trOut=dcExpt(iexp).respwin(rw).trialOutcome;
        [Y_yeses, Y_targets] = getStimAndBehaviorYs(trOut);

       X_targets=dcExpt(iexp).respwin(rw).stimResp';

       cellIdx = dcExpt(iexp).signifResponsiveCells' & ...
           dcExpt(iexp).oriFitTheta90 < theta90Threshold;

       X_targets=X_targets(:,cellIdx);

       idx=find(~isnan(sum(X_targets,2)));
       Y_yeses=Y_yeses(idx,1);   
       Y_targets=Y_targets(idx,1);   
       X_targets=X_targets(idx,:);
       
       trOut=dcExpt(iexp).respwin(rw).audTrialOutcome;
        [Y_yeses_aud, Y_targets_aud] = getStimAndBehaviorYs(trOut);
       X_targets_aud = dcExpt(iexp).respwin(rw).audStimResp';
       X_targets_aud = X_targets_aud(:,cellIdx);

        decisionVariable = dv_target_vis(iexp);
        [~,~,targetGLM] = glmfit(X_targets,Y_targets,'binomial');
        yhat_target = glmval(targetGLM.beta,X_targets,'logit') > decisionVariable;
        pctCorr_target_vis(iexp,rw) = mean((Y_targets-yhat_target) == 0);
        
        yhat_target = glmval(targetGLM.beta,X_targets_aud,'logit') > decisionVariable;
        pctCorr_target_audTest(iexp,rw) = mean((Y_targets_aud-yhat_target) == 0);
    
        decisionVariable = dv_detect_vis(iexp);
        [~,~,detectGLM] = glmfit(X_targets,Y_yeses,'binomial');
        yhat_detect = glmval(detectGLM.beta,X_targets,'logit') > decisionVariable;
        pctCorr_detect_vis(iexp,rw) = mean((Y_yeses-yhat_detect) == 0);
        
        yhat_detect = glmval(detectGLM.beta,X_targets_aud,'logit') > decisionVariable;
        pctCorr_detect_audTest(iexp,rw) = mean((Y_yeses_aud-yhat_detect) == 0);
   end   
end

%% auditory trials
pctCorr_target_aud = nan(nexp,nRespwin);
pctCorr_detect_aud = nan(nexp,nRespwin);

pctCorr_target_visTest = nan(nexp,nRespwin);
pctCorr_detect_visTest = nan(nexp,nRespwin);

for iexp = 1:nexp
   for rw = 1:nRespwin
       trOut=dcExpt(iexp).respwin(rw).audTrialOutcome;
        [Y_yeses, Y_targets] = getStimAndBehaviorYs(trOut);

       X_targets=dcExpt(iexp).respwin(rw).audStimResp';

       cellIdx = dcExpt(iexp).signifResponsiveCells' & ...
           dcExpt(iexp).oriFitTheta90 < theta90Threshold;

       X_targets=X_targets(:,cellIdx);

       idx=find(~isnan(sum(X_targets,2)));
       Y_yeses=Y_yeses(idx,1);   
       Y_targets=Y_targets(idx,1);   
       X_targets=X_targets(idx,:);

       trOut=dcExpt(iexp).respwin(rw).trialOutcome;
        [Y_yeses_vis, Y_targets_vis] = getStimAndBehaviorYs(trOut);
       X_targets_vis = dcExpt(iexp).respwin(rw).stimResp';
       X_targets_vis = X_targets_vis(:,cellIdx);
       
        decisionVariable = dv_target_aud(iexp);
        [~,~,targetGLM] = glmfit(X_targets,Y_targets,'binomial');
        yhat_target = glmval(targetGLM.beta,X_targets,'logit') > decisionVariable;
        pctCorr_target_aud(iexp,rw) = mean((Y_targets-yhat_target) == 0);
        
        yhat_target = glmval(targetGLM.beta,X_targets_vis,'logit') > decisionVariable;
        pctCorr_target_visTest(iexp,rw) = mean((Y_targets_vis-yhat_target) == 0);
        
        decisionVariable = dv_detect_aud(iexp);
        [~,~,detectGLM] = glmfit(X_targets,Y_yeses,'binomial');
        yhat_detect = glmval(detectGLM.beta,X_targets,'logit') > decisionVariable;
        pctCorr_detect_aud(iexp,rw) = mean((Y_yeses-yhat_detect) == 0);
        
        yhat_detect = glmval(detectGLM.beta,X_targets_vis,'logit') > decisionVariable;
        pctCorr_detect_visTest(iexp,rw) = mean((Y_yeses_vis-yhat_detect) == 0);
   end   
end


%% plot % correct by time
winMs = diff(respwinMs(1,:));
tt = respwinMs(:,1)-visStimOffsetMs;
optimalWin = respwinMs(optimalVisDelayInd,:) - visStimOffsetMs;
tooFastTimeMs = 100 - visStimOffsetMs;

setFigParams4Print('landscape')
set(0,'defaultAxesFontSize',16)
figure
suptitle('Visual Trials')
subplot 121
for i = 1:nexp
    h = plot(tt, pctCorr_target_vis(i,:).*100,'-');
    h.Color = [0.5 0.5 0.5];
    hold on
end
y = mean(pctCorr_target_vis,1).*100;
yerr = ste(pctCorr_target_vis,1).*100;
h = errorbar(tt,y,yerr,'k-');
h.LineWidth = 1;
figXAxis([],'Resp Win Start',[tt(1) - winMs tt(end)+winMs],...
    tt(1:2:end),round(tt(1:2:end)))
figYAxis([],'% Correct',[0 110])
figAxForm
title({'Target GLM';sprintf('win size: %s',num2str(winMs))})
vline(tooFastTimeMs,'k--')
vline(optimalWin,'r--')


subplot 122
for i = 1:nexp
    h = plot(tt, pctCorr_detect_vis(i,:).*100,'-');
    h.Color = [0.5 0.5 0.5];
    hold on
end
y = mean(pctCorr_detect_vis,1).*100;
yerr = ste(pctCorr_detect_vis,1).*100;
h = errorbar(tt,y,yerr,'k-');
h.LineWidth = 1;
figXAxis([],'Resp Win Start',[tt(1) - winMs tt(end)+winMs],...
    tt(1:2:end),round(tt(1:2:end)))
figYAxis([],'% Correct',[0 110])
figAxForm
title({'Detect GLM';sprintf('win size: %s',num2str(winMs))})
vline(tooFastTimeMs,'k--')
vline(optimalWin,'r--')

print(fullfile(fnout,'pctCorr_byRespwin_vis'),'-dpdf','-fillpage')

figure
suptitle('Auditory Trials')
subplot 121
for i = 1:nexp
    h = plot(tt, pctCorr_target_aud(i,:).*100,'-');
    h.Color = [0.5 0.5 0.5];
    hold on
end
y = mean(pctCorr_target_aud,1).*100;
yerr = ste(pctCorr_target_aud,1).*100;
h = errorbar(tt,y,yerr,'k-');
h.LineWidth = 1;
figXAxis([],'Resp Win Start',[tt(1) - winMs tt(end)+winMs],...
    tt(1:2:end),round(tt(1:2:end)))
figYAxis([],'% Correct',[0 110])
figAxForm
title({'Target GLM';sprintf('win size: %s',num2str(winMs))})
vline(tooFastTimeMs,'k--')
vline(optimalWin,'r--')


subplot 122
for i = 1:nexp
    h = plot(tt, pctCorr_detect_aud(i,:).*100,'-');
    h.Color = [0.5 0.5 0.5];
    hold on
end
y = mean(pctCorr_detect_aud,1).*100;
yerr = ste(pctCorr_detect_aud,1).*100;
h = errorbar(tt,y,yerr,'k-');
h.LineWidth = 1;
figXAxis([],'Resp Win Start',[tt(1) - winMs tt(end)+winMs],...
    tt(1:2:end),round(tt(1:2:end)))
figYAxis([],'% Correct',[0 110])
figAxForm
title({'Detect GLM';sprintf('win size: %s',num2str(winMs))})
vline(tooFastTimeMs,'k--')
vline(optimalWin,'r--')

print(fullfile(fnout,'pctCorr_byRespwin_aud'),'-dpdf','-fillpage')

%% plot cross-modality test
setFigParams4Print('landscape')
set(0,'defaultAxesFontSize',16)
figure
suptitle('Visual GLM,Test Auditory')
subplot 121
for i = 1:nexp
    h = plot(tt, pctCorr_target_audTest(i,:).*100,'-');
    h.Color = [0.5 0.5 0.5];
    hold on
end
y = mean(pctCorr_target_audTest,1).*100;
yerr = ste(pctCorr_target_audTest,1).*100;
h = errorbar(tt,y,yerr,'k-');
h.LineWidth = 1;
figXAxis([],'Resp Win Start',[tt(1) - winMs tt(end)+winMs],...
    tt(1:2:end),round(tt(1:2:end)))
figYAxis([],'% Correct',[0 110])
figAxForm
title('Target')
vline(tooFastTimeMs,'k--')
vline(optimalWin,'r--')


subplot 122
for i = 1:nexp
    h = plot(tt, pctCorr_detect_audTest(i,:).*100,'-');
    h.Color = [0.5 0.5 0.5];
    hold on
end
y = mean(pctCorr_detect_audTest,1).*100;
yerr = ste(pctCorr_detect_audTest,1).*100;
h = errorbar(tt,y,yerr,'k-');
h.LineWidth = 1;
figXAxis([],'Resp Win Start',[tt(1) - winMs tt(end)+winMs],...
    tt(1:2:end),round(tt(1:2:end)))
figYAxis([],'% Correct',[0 110])
figAxForm
title('Detect')
vline(tooFastTimeMs,'k--')
vline(optimalWin,'r--')

print(fullfile(fnout,'pctCorr_byRespwin_visGLMaudTest'),'-dpdf','-fillpage')

figure
suptitle('Auditory GLM, Test Visual')
subplot 121
for i = 1:nexp
    h = plot(tt, pctCorr_target_visTest(i,:).*100,'-');
    h.Color = [0.5 0.5 0.5];
    hold on
end
y = mean(pctCorr_target_visTest,1).*100;
yerr = ste(pctCorr_target_visTest,1).*100;
h = errorbar(tt,y,yerr,'k-');
h.LineWidth = 1;
figXAxis([],'Resp Win Start',[tt(1) - winMs tt(end)+winMs],...
    tt(1:2:end),round(tt(1:2:end)))
figYAxis([],'% Correct',[0 110])
figAxForm
title('Target')
vline(tooFastTimeMs,'k--')
vline(optimalWin,'r--')


subplot 122
for i = 1:nexp
    h = plot(tt, pctCorr_detect_visTest(i,:).*100,'-');
    h.Color = [0.5 0.5 0.5];
    hold on
end
y = mean(pctCorr_detect_visTest,1).*100;
yerr = ste(pctCorr_detect_visTest,1).*100;
h = errorbar(tt,y,yerr,'k-');
h.LineWidth = 1;
figXAxis([],'Resp Win Start',[tt(1) - winMs tt(end)+winMs],...
    tt(1:2:end),round(tt(1:2:end)))
figYAxis([],'% Correct',[0 110])
figAxForm
title('Detect')
vline(tooFastTimeMs,'k--')
vline(optimalWin,'r--')

print(fullfile(fnout,'pctCorr_byRespwin_audGLMvisTest'),'-dpdf','-fillpage')


