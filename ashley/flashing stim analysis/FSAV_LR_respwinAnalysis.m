close all
clear all

load('Z:\Analysis\FSAV Choice\FSAV_respwinData.mat')
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
nRespwin = size(dcExpt(1).respwin,2);
optimalVisDelayInd = 2;
visRTRange = [250 570];
audRTRange = [150 450];

respwinMs = cell(1,3);
for rw = 1:nRespwin
    respwinMs{rw} = dcExpt(1).respwin(rw).winMs;
end
%%
targetWeightsAll = cell(1,nRespwin);
detectWeightsAll = cell(1,nRespwin);
targetWeightsVar = cell(1,nRespwin);
detectWeightsVar = cell(1,nRespwin);
nCells = cell(1,nRespwin);
targetCorrsAll = cell(1,nRespwin);
detectCorrsAll = cell(1,nRespwin);
oriPrefAll = cell(1,nRespwin);
for iexp=1:nexp
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

       detectCorr = corr(Y_yeses,X_targets)';
       targetCorr = corr(Y_targets,X_targets)';

       detectCorrsAll{rw} = cat(1,detectCorrsAll{rw},detectCorr);
       targetCorrsAll{rw} = cat(1,targetCorrsAll{rw},targetCorr);

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

       targetWeightsAll{rw} = cat(1,targetWeightsAll{rw},B_targets);
       detectWeightsAll{rw} = cat(1,detectWeightsAll{rw},B_yeses);

       targetWeightsVar{rw} = cat(1,targetWeightsVar{rw},var(B_targets));
       detectWeightsVar{rw} = cat(1,detectWeightsVar{rw},var(B_yeses));
       nCells{rw} = cat(1,nCells{rw},length(B_targets));

       oriPref = dcExpt(iexp).oriPref(cellIdx);
       oriPrefAll{rw} = cat(2,oriPrefAll{rw},oriPref);

       if doExptPlots
           if rw == optimalVisDelayInd
               rtLim = [0 600];
               yLim = [0 55];
               rtBins = 0:25:600;
            visRT = dcExpt(iexp).visRT;
            audRT = dcExpt(iexp).audRT;

            figure
            suptitle([dcExpt(iexp).mouse '-' dcExpt(iexp).date])
            subplot 121
            h = histogram(visRT,rtBins);
            hold on 
            for irw = 1:nRespwin
                h = plot(dcExpt(iexp).respwin(irw).winMs,[irw irw]+50,'b-');
                h.LineWidth = 2;
            end
            title('Visual Trials')
            figXAxis([],'Reaction Time',rtLim)
            figYAxis([],'N Trials',yLim)
            figAxForm
            vline(visRTRange,'k:')
            subplot 122
            h = histogram(audRT,rtBins);
            hold on 
            for irw = 1:nRespwin
                h = plot(dcExpt(iexp).respwin(irw).winMs,[irw irw]+50,'b-');
                h.LineWidth = 2;
            end
            title('Auditory Trials')
            figXAxis([],'Reaction Time',rtLim)
            figYAxis([],'N Trials',yLim)
            figAxForm
            vline(audRTRange,'k:')

            print(fullfile(fnout,...
                [dcExpt(iexp).mouse '_' dcExpt(iexp).date '_rt']),...
                '-dpdf','-fillpage')
           end
%         figure
%         suptitle({[dcExpt(iexp).mouse '-' dcExpt(iexp).date],...
%             sprintf('Response Window %s',num2str(rw))})
%         subplot 221
%         scatter(B_yeses,detectCorr,'bo')
%         hold on
%         hline(0,'k:')
%         vline(0,'k:')
%         figXAxis([],'Detect LR Weight',[])
%         figYAxis([],'Detect Correlation',[])
%         figAxForm
%         subplot 222
%         scatter(B_targets,targetCorr,'bo')
%         hold on
%         hline(0,'k:')
%         vline(0,'k:')
%         figXAxis([],'Target LR Weight',[])
%         figYAxis([],'Target Correlation',[])
%         figAxForm
%         subplot 223   
%         s = scatter(B_targets,B_yeses,'bo');
%         s.MarkerFaceColor = [1 1 1];
%         hold on
%         hline(0,'k:')
%         vline(0,'k:')
%         figXAxis([],'Target LR Weight',[])
%         figYAxis([],'Detect LR Weight',[])
%         figAxForm
%         subplot 224
%         axis off
%         text(0.5,1,sprintf('N Cells = %s',num2str(sum(cellIdx))))
%         hold on
%         text(0.5, 0.5,sprintf('N Trials = %s',num2str(length(Y_yeses))))
%         figXAxis([],'',[0 2])
%         figYAxis([],'',[0 2])

%         print(fullfile(fnout,...
%             [dcExpt(iexp).mouse '_' dcExpt(iexp).date '_corrXweight']),...
%             '-dpdf','-fillpage')

%         figure
%         suptitle({[dcExpt(iexp).mouse '-' dcExpt(iexp).date],...
%             sprintf('Response Window %s',num2str(rw))})
%         subplot 221
%         scatter(oriPref,B_targets,'bo')
%         hold on
%         hline(0,'k:')
%         vline(0,'k:')
%         figXAxis([],'Orientation Preference',[-10 190],[0:45:180],[0:45:180])
%         figYAxis([],'Target LR Weight',[])
%         figAxForm
%         subplot 222
%         scatter(oriPref,targetCorr,'bo')
%         hold on
%         hline(0,'k:')
%         vline(0,'k:')
%         figXAxis([],'Orientation Preference',[-10 190],[0:45:180],[0:45:180])
%         figYAxis([],'Target Correlation',[])
%         figAxForm
%         subplot 223
%         scatter(oriPref,B_yeses,'bo')
%         hold on
%         hline(0,'k:')
%         vline(0,'k:')
%         figXAxis([],'Orientation Preference',[-10 190],[0:45:180],[0:45:180])
%         figYAxis([],'Detect LR Weight',[])
%         figAxForm
%         subplot 224
%         scatter(oriPref,detectCorr,'bo')
%         hold on
%         hline(0,'k:')
%         vline(0,'k:')
%         figXAxis([],'Orientation Preference',[-10 190],[0:45:180],[0:45:180])
%         figYAxis([],'Detect Correlation',[])
%         figAxForm

%         print(fullfile(fnout,...
%             [dcExpt(iexp).mouse '_' dcExpt(iexp).date '_oriPrefXcorrXweight']),...
%             '-dpdf','-fillpage')
       end
    end
end

%%
weightLim = [-1 3];
corrLim = [-1 1];

figure
subplot 331
scatter(detectWeightsAll{1},detectCorrsAll{1},'bo')
hold on
hline(0,'k:')
vline(0,'k:')
figXAxis([],'Detect LR Weight',weightLim)
figYAxis([],'Detect Correlation',corrLim)
figAxForm
title('Resp Win 1')
subplot 332
scatter(targetWeightsAll{1},targetCorrsAll{1},'bo')
hold on
hline(0,'k:')
vline(0,'k:')
figXAxis([],'Target LR Weight',weightLim)
figYAxis([],'Target Correlation',corrLim)
figAxForm
subplot 333   
s = scatter(targetWeightsAll{1},detectWeightsAll{1},'bo');
s.MarkerFaceColor = [1 1 1];
hold on
hline(0,'k:')
vline(0,'k:')
figXAxis([],'Target LR Weight',weightLim)
figYAxis([],'Detect LR Weight',weightLim)
figAxForm

subplot 334
scatter(detectWeightsAll{2},detectCorrsAll{2},'bo')
hold on
hline(0,'k:')
vline(0,'k:')
figXAxis([],'Detect LR Weight',weightLim)
figYAxis([],'Detect Correlation',corrLim)
figAxForm
title('Resp Win 2')
subplot 335
scatter(targetWeightsAll{2},targetCorrsAll{2},'bo')
hold on
hline(0,'k:')
vline(0,'k:')
figXAxis([],'Target LR Weight',weightLim)
figYAxis([],'Target Correlation',corrLim)
figAxForm
subplot 336   
s = scatter(targetWeightsAll{2},detectWeightsAll{2},'bo');
s.MarkerFaceColor = [1 1 1];
hold on
hline(0,'k:')
vline(0,'k:')
figXAxis([],'Target LR Weight',weightLim)
figYAxis([],'Detect LR Weight',weightLim)
figAxForm

subplot 337
scatter(detectWeightsAll{3},detectCorrsAll{3},'bo')
hold on
hline(0,'k:')
vline(0,'k:')
figXAxis([],'Detect LR Weight',weightLim)
figYAxis([],'Detect Correlation',corrLim)
figAxForm
title('Resp Win 3')
subplot 338
scatter(targetWeightsAll{3},targetCorrsAll{3},'bo')
hold on
hline(0,'k:')
vline(0,'k:')
figXAxis([],'Target LR Weight',weightLim)
figYAxis([],'Target Correlation',corrLim)
figAxForm
subplot 339   
s = scatter(targetWeightsAll{3},detectWeightsAll{3},'bo');
s.MarkerFaceColor = [1 1 1];
hold on
hline(0,'k:')
vline(0,'k:')
figXAxis([],'Target LR Weight',weightLim)
figYAxis([],'Detect LR Weight',weightLim)
figAxForm

print(fullfile(fnout,'corrXweight_respwin'),'-dpdf','-fillpage')

figure
suptitle('Resp Win 1 vs Resp Win 2')
subplot 221
s = scatter(targetCorrsAll{1},targetCorrsAll{2},'bo');
s.MarkerFaceColor = [1 1 1];
hold on
plot(corrLim,corrLim,'k--')
hline(0,'k:')
vline(0,'k:')
figXAxis([],'Resp Win 1',corrLim)
figYAxis([],'Resp Win 2',corrLim)
figAxForm
title('Target Correlation')
subplot 222
s = scatter(detectCorrsAll{1},detectCorrsAll{2},'bo');
s.MarkerFaceColor = [1 1 1];
hold on
plot(corrLim,corrLim,'k--')
hline(0,'k:')
vline(0,'k:')
figXAxis([],'Resp Win 1',corrLim)
figYAxis([],'Resp Win 2',corrLim)
figAxForm
title('Detect Correlation')
subplot 223
s = scatter(targetWeightsAll{1},targetWeightsAll{2},'bo');
s.MarkerFaceColor = [1 1 1];
hold on
plot(weightLim,weightLim,'k--')
hline(0,'k:')
vline(0,'k:')
figXAxis([],'Resp Win 1',weightLim)
figYAxis([],'Resp Win 2',weightLim)
figAxForm
title('Target LR Weights')
subplot 224
s = scatter(detectWeightsAll{1},detectWeightsAll{2},'bo');
s.MarkerFaceColor = [1 1 1];
hold on
plot(weightLim,weightLim,'k--')
hline(0,'k:')
vline(0,'k:')
figXAxis([],'Resp Win 1',weightLim)
figYAxis([],'Resp Win 2',weightLim)
figAxForm
title('Detect LR Weights')

print(fullfile(fnout,'compareRespwin1and2'),'-dpdf','-fillpage')

figure
suptitle('Resp Win 3 vs Resp Win 2')
subplot 221
s = scatter(targetCorrsAll{3},targetCorrsAll{2},'bo');
s.MarkerFaceColor = [1 1 1];
hold on
plot(corrLim,corrLim,'k--')
hline(0,'k:')
vline(0,'k:')
figXAxis([],'Resp Win 3',corrLim)
figYAxis([],'Resp Win 2',corrLim)
figAxForm
title('Target Correlation')
subplot 222
s = scatter(detectCorrsAll{3},detectCorrsAll{2},'bo');
s.MarkerFaceColor = [1 1 1];
hold on
plot(corrLim,corrLim,'k--')
hline(0,'k:')
vline(0,'k:')
figXAxis([],'Resp Win 3',corrLim)
figYAxis([],'Resp Win 2',corrLim)
figAxForm
title('Detect Correlation')
subplot 223
s = scatter(targetWeightsAll{3},targetWeightsAll{2},'bo');
s.MarkerFaceColor = [1 1 1];
hold on
plot(weightLim,weightLim,'k--')
hline(0,'k:')
vline(0,'k:')
figXAxis([],'Resp Win 3',weightLim)
figYAxis([],'Resp Win 2',weightLim)
figAxForm
title('Target LR Weights')
subplot 224
s = scatter(detectWeightsAll{3},detectWeightsAll{2},'bo');
s.MarkerFaceColor = [1 1 1];
hold on
plot(weightLim,weightLim,'k--')
hline(0,'k:')
vline(0,'k:')
figXAxis([],'Resp Win 3',weightLim)
figYAxis([],'Resp Win 2',weightLim)
figAxForm
title('Detect LR Weights')

print(fullfile(fnout,'compareRespwin3and2'),'-dpdf','-fillpage')
%%
targetWeightsAll_aud = cell(1,nRespwin);
detectWeightsAll_aud = cell(1,nRespwin);
targetCorrsAll_aud = cell(1,nRespwin);
detectCorrsAll_aud = cell(1,nRespwin);
for iexp=1:nexp
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

       detectCorr = corr(Y_yeses,X_targets)';
       targetCorr = corr(Y_targets,X_targets)';

       detectCorrsAll_aud{rw} = cat(1,detectCorrsAll_aud{rw},detectCorr);
       targetCorrsAll_aud{rw} = cat(1,targetCorrsAll_aud{rw},targetCorr);

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

       targetWeightsAll_aud{rw} = cat(1,targetWeightsAll_aud{rw},B_targets);
       detectWeightsAll_aud{rw} = cat(1,detectWeightsAll_aud{rw},B_yeses);
    end
end

%%
weightLim = [-1.1 2];
corrLim = [-1 1];

figure
suptitle('Auditory Trials')
subplot 331
scatter(detectWeightsAll_aud{1},detectCorrsAll_aud{1},'bo')
hold on
hline(0,'k:')
vline(0,'k:')
figXAxis([],'Detect LR Weight',weightLim)
figYAxis([],'Detect Correlation',corrLim)
figAxForm
title('Resp Win 1')
subplot 332
scatter(targetWeightsAll_aud{1},targetCorrsAll_aud{1},'bo')
hold on
hline(0,'k:')
vline(0,'k:')
figXAxis([],'Target LR Weight',weightLim)
figYAxis([],'Target Correlation',corrLim)
figAxForm
subplot 333   
s = scatter(targetWeightsAll_aud{1},detectWeightsAll_aud{1},'bo');
s.MarkerFaceColor = [1 1 1];
hold on
hline(0,'k:')
vline(0,'k:')
figXAxis([],'Target LR Weight',weightLim)
figYAxis([],'Detect LR Weight',weightLim)
figAxForm

subplot 334
scatter(detectWeightsAll_aud{2},detectCorrsAll_aud{2},'bo')
hold on
hline(0,'k:')
vline(0,'k:')
figXAxis([],'Detect LR Weight',weightLim)
figYAxis([],'Detect Correlation',corrLim)
figAxForm
title('Resp Win 2')
subplot 335
scatter(targetWeightsAll_aud{2},targetCorrsAll_aud{2},'bo')
hold on
hline(0,'k:')
vline(0,'k:')
figXAxis([],'Target LR Weight',weightLim)
figYAxis([],'Target Correlation',corrLim)
figAxForm
subplot 336   
s = scatter(targetWeightsAll_aud{2},detectWeightsAll_aud{2},'bo');
s.MarkerFaceColor = [1 1 1];
hold on
hline(0,'k:')
vline(0,'k:')
figXAxis([],'Target LR Weight',weightLim)
figYAxis([],'Detect LR Weight',weightLim)
figAxForm

subplot 337
scatter(detectWeightsAll_aud{3},detectCorrsAll_aud{3},'bo')
hold on
hline(0,'k:')
vline(0,'k:')
figXAxis([],'Detect LR Weight',weightLim)
figYAxis([],'Detect Correlation',corrLim)
figAxForm
title('Resp Win 3')
subplot 338
scatter(targetWeightsAll_aud{3},targetCorrsAll_aud{3},'bo')
hold on
hline(0,'k:')
vline(0,'k:')
figXAxis([],'Target LR Weight',weightLim)
figYAxis([],'Target Correlation',corrLim)
figAxForm
subplot 339   
s = scatter(targetWeightsAll_aud{3},detectWeightsAll_aud{3},'bo');
s.MarkerFaceColor = [1 1 1];
hold on
hline(0,'k:')
vline(0,'k:')
figXAxis([],'Target LR Weight',weightLim)
figYAxis([],'Detect LR Weight',weightLim)
figAxForm

print(fullfile(fnout,'corrXweight_respwin_audTrials'),'-dpdf','-fillpage')

figure
suptitle('Resp Win 1 vs Resp Win 2')
subplot 221
s = scatter(targetCorrsAll_aud{1},targetCorrsAll_aud{2},'bo');
s.MarkerFaceColor = [1 1 1];
hold on
plot(corrLim,corrLim,'k--')
hline(0,'k:')
vline(0,'k:')
figXAxis([],'Resp Win 1',corrLim)
figYAxis([],'Resp Win 2',corrLim)
figAxForm
title('Target Correlation')
subplot 222
s = scatter(detectCorrsAll_aud{1},detectCorrsAll_aud{2},'bo');
s.MarkerFaceColor = [1 1 1];
hold on
plot(corrLim,corrLim,'k--')
hline(0,'k:')
vline(0,'k:')
figXAxis([],'Resp Win 1',corrLim)
figYAxis([],'Resp Win 2',corrLim)
figAxForm
title('Detect Correlation')
subplot 223
s = scatter(targetWeightsAll_aud{1},targetWeightsAll_aud{2},'bo');
s.MarkerFaceColor = [1 1 1];
hold on
plot(weightLim,weightLim,'k--')
hline(0,'k:')
vline(0,'k:')
figXAxis([],'Resp Win 1',weightLim)
figYAxis([],'Resp Win 2',weightLim)
figAxForm
title('Target LR Weights')
subplot 224
s = scatter(detectWeightsAll_aud{1},detectWeightsAll_aud{2},'bo');
s.MarkerFaceColor = [1 1 1];
hold on
plot(weightLim,weightLim,'k--')
hline(0,'k:')
vline(0,'k:')
figXAxis([],'Resp Win 1',weightLim)
figYAxis([],'Resp Win 2',weightLim)
figAxForm
title('Detect LR Weights')

print(fullfile(fnout,'compareRespwin1and2_audTrials'),'-dpdf','-fillpage')

figure
suptitle('Resp Win 3 vs Resp Win 2')
subplot 221
s = scatter(targetCorrsAll_aud{3},targetCorrsAll_aud{2},'bo');
s.MarkerFaceColor = [1 1 1];
hold on
plot(corrLim,corrLim,'k--')
hline(0,'k:')
vline(0,'k:')
figXAxis([],'Resp Win 3',corrLim)
figYAxis([],'Resp Win 2',corrLim)
figAxForm
title('Target Correlation')
subplot 222
s = scatter(detectCorrsAll_aud{3},detectCorrsAll_aud{2},'bo');
s.MarkerFaceColor = [1 1 1];
hold on
plot(corrLim,corrLim,'k--')
hline(0,'k:')
vline(0,'k:')
figXAxis([],'Resp Win 3',corrLim)
figYAxis([],'Resp Win 2',corrLim)
figAxForm
title('Detect Correlation')
subplot 223
s = scatter(targetWeightsAll_aud{3},targetWeightsAll_aud{2},'bo');
s.MarkerFaceColor = [1 1 1];
hold on
plot(weightLim,weightLim,'k--')
hline(0,'k:')
vline(0,'k:')
figXAxis([],'Resp Win 3',weightLim)
figYAxis([],'Resp Win 2',weightLim)
figAxForm
title('Target LR Weights')
subplot 224
s = scatter(detectWeightsAll_aud{3},detectWeightsAll_aud{2},'bo');
s.MarkerFaceColor = [1 1 1];
hold on
plot(weightLim,weightLim,'k--')
hline(0,'k:')
vline(0,'k:')
figXAxis([],'Resp Win 3',weightLim)
figYAxis([],'Resp Win 2',weightLim)
figAxForm
title('Detect LR Weights')

print(fullfile(fnout,'compareRespwin3and2_audTrials'),'-dpdf','-fillpage')

%%

weightLim = [-1.1 3];
corrLim = [-1 1];
corrAlpha = 0.05/4;

for rw = 1:nRespwin
    
    figure
    suptitle(sprintf('Response Window %s',num2str(rw)))
    subplot 221
    x = targetCorrsAll{rw};
    y = targetCorrsAll_aud{rw};
    [r,p] = corr(x,y,'Type','Pearson');
    P = polyfit(x,y,1);
    yFit = polyval(P,x);
    scatter(x,y,'ko')
    hold on
    if p < corrAlpha
        h = plot(x,yFit,'r-');
        h.LineWidth = 1;
    else
        plot(x,yFit,'k-')
        h.LineWidth = 1;
    end
    plot(corrLim,corrLim,'k--')
    figXAxis([],'Visual Trials',corrLim)
    figYAxis([],'Auditory Trials',corrLim)
    figAxForm
    vline(0,'k:')
    hline(0,'k:')
    title({'Target Correlation';...
        sprintf('r = %s, p = %s',num2str(r),num2str(p))})
    subplot 222
    x = detectCorrsAll{rw};
    y = detectCorrsAll_aud{rw};
    [r,p] = corr(x,y,'Type','Pearson');
    P = polyfit(x,y,1);
    yFit = polyval(P,x);
    scatter(x,y,'ko')
    hold on
    if p < corrAlpha
        plot(x,yFit,'r-')
        h.LineWidth = 1;
    else
        plot(x,yFit,'k-')
        h.LineWidth = 1;
    end
    plot(corrLim,corrLim,'k--')
    figXAxis([],'Visual Trials',corrLim)
    figYAxis([],'Auditory Trials',corrLim)
    figAxForm
    vline(0,'k:')
    hline(0,'k:')
    title({'Detect Correlation';...
        sprintf('r = %s, p = %s',num2str(r),num2str(p))})
    subplot 223
    x = targetWeightsAll{rw};
    y = targetWeightsAll_aud{rw};
    [r,p] = corr(x,y,'Type','Pearson');
    P = polyfit(x,y,1);
    yFit = polyval(P,x);
    scatter(x,y,'ko')
    hold on
    if p < corrAlpha
        plot(x,yFit,'r-')
        h.LineWidth = 1;
    else
        plot(x,yFit,'k-')
        h.LineWidth = 1;
    end
    plot(weightLim,weightLim,'k--')
    figXAxis([],'Visual Trials',weightLim)
    figYAxis([],'Auditory Trials',weightLim)
    figAxForm
    vline(0,'k:')
    hline(0,'k:')
    title({'Target Weight';...
        sprintf('r = %s, p = %s',num2str(r),num2str(p))})
    subplot 224
    x = detectWeightsAll{rw};
    y = detectWeightsAll_aud{rw};
    [r,p] = corr(x,y,'Type','Pearson');
    P = polyfit(x,y,1);
    yFit = polyval(P,x);
    scatter(x,y,'ko')
    hold on
    if p < corrAlpha
        plot(x,yFit,'r-')
        h.LineWidth = 1;
    else
        plot(x,yFit,'k-')
        h.LineWidth = 1;
    end
    plot(weightLim,weightLim,'k--')
    figXAxis([],'Visual Trials',weightLim)
    figYAxis([],'Auditory Trials',weightLim)
    figAxForm
    vline(0,'k:')
    hline(0,'k:')
    title({'Detect Weight';...
        sprintf('r = %s, p = %s',num2str(r),num2str(p))})    
    
    saveName = sprintf('compareVisAudCorrsWeightsRespwin%s',num2str(rw));
    print(fullfile(fnout,saveName),'-dpdf','-fillpage')
end