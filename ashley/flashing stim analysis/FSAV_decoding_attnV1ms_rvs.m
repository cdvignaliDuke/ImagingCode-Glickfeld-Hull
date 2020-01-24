clear all
close all
ds = 'FSAV_attentionV1_noAttn';
analysisDate = '200122';
titleStr = ds(6:end);
rc = behavConstsAV;
imgParams_FSAV
bxParams_FSAV_attnV1ms
fn = fullfile(rc.ashley, 'Manuscripts','Attention V1','Matlab Figs',...
    [titleStr '_' datestr(now,'yymmdd') '_']); 
load([fn(1:end-7) analysisDate '_imgAnalysisData'])
dc_noAttn = decodeAnalysis;
mouseName_noAttn = {respCellsExpt.mouse};

ds = 'FSAV_V1_naive_GCaMP6m';
analysisDate = '200117';
titleStr = ds(6:end);
rc = behavConstsAV;
imgParams_FSAV
bxParams_FSAV_attnV1ms
fn = fullfile(rc.ashley, 'Manuscripts','Attention V1','Matlab Figs',...
    [titleStr '_' datestr(now,'yymmdd') '_']); 
load([fn(1:end-7) analysisDate '_imgAnalysisData'])
dc_naive = decodeAnalysis;
mouseName_naive = {respCellsExpt.mouse};

ds = 'FSAV_attentionV1';
analysisDate = '200122';
titleStr = ds(6:end);
rc = behavConstsAV;
imgParams_FSAV
bxParams_FSAV_attnV1ms
fn = fullfile(rc.ashley, 'Manuscripts','Attention V1','Matlab Figs',...
    [titleStr '_' datestr(now,'yymmdd') '_']); 
load([fn(1:end-7) analysisDate '_imgAnalysisData'])
dc_attn = decodeAnalysis;
mouseName_attn = {respCellsExpt.mouse};

fnout = fullfile(rc.ashley, 'Manuscripts','Attention V1','Matlab Figs',...
    [titleStr '_' datestr(now,'yymmdd') '_decoding_']);

%%
load(fullfile(rc.ashley, 'Manuscripts','Attention V1','Matlab Figs','bxStats.mat'))
    
%% indexing and naming variables
avName = {'Vis','Aud'};
targetName = {'D','HT','ET'};
modelName = {'Stimulus','Choice'};

nexp_attn = size(dc_attn,2);
nexp_noAttn = size(dc_noAttn,2);
nexp_naive = size(dc_naive,2);

stimModInd = 1;
choiceModInd = 2;
%% performance all trials
%model performance - attention mice
mdl_attn = struct;
for iav = 1:2
    mdl_attn(iav).name = avName{iav};
    mdl_attn(iav).mdl(stimModInd).name = 'Stimulus';
    mdl_attn(iav).mdl(choiceModInd).name = 'Choice';
    mdl_attn(iav).mdl(stimModInd).pctCorr_all = nan(1,nexp_attn);
    mdl_attn(iav).mdl(choiceModInd).pctCorr_all = nan(1,nexp_attn);
    for iexp = 1:nexp_attn
        mdl_attn(iav).mdl(stimModInd).pctCorr_all(iexp) = ...
            dc_attn(iexp).av(iav).pctCorrectAllTarget_holdout;
        mdl_attn(iav).mdl(choiceModInd).pctCorr_all(iexp) = ...
            dc_attn(iexp).av(iav).pctCorrectAllDetect_holdout;
    end
end
goodPerfInd_attn = mdl_attn(visualTrials).mdl(stimModInd).pctCorr_all...
    > pctCorrThresh;
%model performance - no attention mice
mdl_noAttn = struct;
for iav = 1:2
    mdl_noAttn(iav).name = avName{iav};
    mdl_noAttn(iav).mdl(stimModInd).name = 'Stimulus';
    mdl_noAttn(iav).mdl(choiceModInd).name = 'Choice';
    mdl_noAttn(iav).mdl(stimModInd).pctCorr_all = nan(1,nexp_noAttn);
    mdl_noAttn(iav).mdl(choiceModInd).pctCorr_all = nan(1,nexp_noAttn);
    for iexp = 1:nexp_noAttn
        mdl_noAttn(iav).mdl(stimModInd).pctCorr_all(iexp) = ...
            dc_noAttn(iexp).av(iav).pctCorrectAllTarget_holdout;
        mdl_noAttn(iav).mdl(choiceModInd).pctCorr_all(iexp) = ...
            dc_noAttn(iexp).av(iav).pctCorrectAllDetect_holdout;
    end
end
goodPerfInd_noAttn = mdl_noAttn(visualTrials).mdl(stimModInd).pctCorr_all...
    > pctCorrThresh;

%model performance - naive mice
mdl_naive = struct;
for iav = 1:2
    mdl_naive(iav).name = avName{iav};
    mdl_naive(iav).mdl(stimModInd).name = 'Stimulus';
    mdl_naive(iav).mdl(stimModInd).pctCorr_all = nan(1,nexp_naive);
    for iexp = 1:nexp_naive
        mdl_naive(iav).mdl(stimModInd).pctCorr_all(iexp) = ...
            dc_naive(iexp).av(iav).pctCorrectAllTarget_holdout;
    end
end
goodPerfInd_naive = mdl_naive(visualTrials).mdl(stimModInd).pctCorr_all...
    > pctCorrThresh;

mdlPerf_fig = figure;
suptitle(sprintf(...
    'All Trials Performance (held-out data), pct corr must be > %s in vis stim model',...
    num2str(pctCorrThresh)))
for iav = 1:2
    for imod = 1:2
        if iav == 1
            iplot = imod;
        else
            iplot = imod+2;
        end
        subplot(2,2,iplot)
        y1 = mdl_attn(iav).mdl(imod).pctCorr_all(goodPerfInd_attn);
        plot(ones(1,length(y1)),y1,'k.','MarkerSize',10)
        hold on
        y2 = mdl_noAttn(iav).mdl(imod).pctCorr_all(goodPerfInd_noAttn);
        plot(ones(1,length(y2)).*2,y2,'k.','MarkerSize',10)
        if imod == 1
            x = 1:3;
            y3 = mdl_naive(iav).mdl(imod).pctCorr_all(goodPerfInd_naive);
            plot(ones(1,length(y3))*3,y3,'k.','MarkerSize',10);
            y_all = [mean(y1),mean(y2),mean(y3)];
            y_all_ste = [ste(y1,2),ste(y2,2),ste(y3,2)];
            figXAxis([],'',[0 4],1:3,{'Attn','No Attn','Naive'})
        else
            x = 1:2;
            y_all = [mean(y1),mean(y2)];
            y_all_ste = [ste(y1,2),ste(y2,2)];
            figXAxis([],'',[0 3],1:2,{'Attn','No Attn'})
        end
        errorbar(x,y_all,y_all_ste,'.','MarkerSize',20)
        figYAxis([],'% Correct',[0 1],0:0.2:1)
        figAxForm([],0)
        hline(pctCorrThresh,'k:')
        hline(0.5,'k-')
        [~,p] = ttest2(y1,y2);
        title(sprintf('%s %s Model, p=%s',mdl_attn(iav).name,...
            mdl_attn(iav).mdl(imod).name,num2str(round(p,2,'significant'))))
    end
end

print([fnout 'pctCorrect'],'-dpdf','-fillpage')
%% performance across stim types
for iav = 1:2
    mdl_attn(iav).mdl(stimModInd).pctCorr_xstim = nan(3,nexp_attn);
    mdl_attn(iav).mdl(choiceModInd).pctCorr_xstim = nan(3,nexp_attn);
    for iexp = 1:nexp_attn
        mdl_attn(iav).mdl(stimModInd).pctCorr_xstim(:,iexp) = ...
            dc_attn(iexp).av(iav).pctCorrectXStimTarget_holdout;
        mdl_attn(iav).mdl(choiceModInd).pctCorr_xstim(:,iexp) = ...
            dc_attn(iexp).av(iav).pctCorrectXStimDetect_holdout;
    end
    mdl_noAttn(iav).mdl(stimModInd).pctCorr_xstim = nan(3,nexp_noAttn);
    mdl_noAttn(iav).mdl(choiceModInd).pctCorr_xstim = nan(3,nexp_noAttn);
    for iexp = 1:nexp_noAttn
        mdl_noAttn(iav).mdl(stimModInd).pctCorr_xstim(:,iexp) = ...
            dc_noAttn(iexp).av(iav).pctCorrectXStimTarget_holdout;
        mdl_noAttn(iav).mdl(choiceModInd).pctCorr_xstim(:,iexp) = ...
            dc_noAttn(iexp).av(iav).pctCorrectXStimDetect_holdout;
    end
    
    mdl_naive(iav).mdl(stimModInd).pctCorr_xstim = nan(3,nexp_naive);
    for iexp = 1:nexp_naive
        mdl_naive(iav).mdl(stimModInd).pctCorr_xstim(:,iexp) = ...
            dc_naive(iexp).av(iav).pctCorrectXStimTarget_holdout;
    end    
end

mdlPerfxStim_fig = figure;
suptitle(sprintf(...
    'Ea. Stim Performance (held-out data), pct corr must be > %s in vis stim model',...
    num2str(pctCorrThresh)))
x=1:3;
for iav = 1:2
    for imod = 1:2
        if iav == 1
            iplot = imod;
        else
            iplot = imod+2;
        end
        subplot(3,4,iplot)
        hold on
        y = mdl_attn(iav).mdl(imod).pctCorr_xstim(:,goodPerfInd_attn);
        ind = sum(isnan(y),2)<minTrN;
        for i = 1:size(y,2)
            nanInd = ~isnan(y(:,i)) & ind;
            plot(x(nanInd),y(nanInd,i),'k-')
        end
        errorbar(x(ind),nanmean(y(ind,:),2),ste(y(ind,:),2),'.-','MarkerSize',20)
        figXAxis([],'Test',[0 4],x,targetName)
        figYAxis([],'% Correct',[0 1],0:0.2:1)
        figAxForm
        hline(pctCorrThresh,'k:')
        hline(0.5,'k-')
        [~,pChance] = ttest(y',0.5);
        anovaTest = anova1(y',[],'off');
        title({'Attn Mice';sprintf('%s %s Model, anovaP=%s',mdl_attn(iav).name,...
            mdl_attn(iav).mdl(imod).name,num2str(round(anovaTest,2,'significant')))})
        for i = 1:3
            text(i,0.2,num2str(round(pChance(i),2,'significant')))
        end
        text(0,0.2,'pChnc=')
        
        if iav == 1
            iplot = imod+4;
        else
            iplot = imod+6;
        end
        subplot(3,4,iplot)
        hold on
        y = mdl_noAttn(iav).mdl(imod).pctCorr_xstim(:,goodPerfInd_noAttn);
        ind = sum(isnan(y),2)<minTrN;
        for i = 1:size(y,2)
            nanInd = ~isnan(y(:,i)) & ind;
            plot(x(nanInd),y(nanInd,i),'k-')
        end
        errorbar(x(ind),nanmean(y(ind,:),2),ste(y(ind,:),2),'.-','MarkerSize',20)
        figXAxis([],'Test',[0 4],x,targetName)
        figYAxis([],'% Correct',[0 1],0:0.2:1)
        figAxForm
        hline(pctCorrThresh,'k:')
        hline(0.5,'k-')
        [~,pChance] = ttest(y',0.5);
        anovaTest = anova1(y',[],'off');
        title({'No Attn Mice';sprintf('%s %s Model, anovaP=%s',mdl_noAttn(iav).name,...
            mdl_noAttn(iav).mdl(imod).name,num2str(round(anovaTest,2,'significant')))})
        for i = 1:3
            text(i,0.2,num2str(round(pChance(i),2,'significant')))
        end
        text(0,0.2,'pChnc=')
        if imod == 1
            if iav == 1
                iplot = imod+8;
            else
                iplot = imod+10;
            end
            subplot(3,4,iplot)
            hold on
            y = mdl_naive(iav).mdl(imod).pctCorr_xstim(:,goodPerfInd_naive);
            ind = sum(isnan(y),2)<minTrN;
            for i = 1:size(y,2)
                nanInd = ~isnan(y(:,i)) & ind;
                plot(x(nanInd),y(nanInd,i),'k-')
            end
            errorbar(x(ind),nanmean(y(ind,:),2),ste(y(ind,:),2),'.-','MarkerSize',20)
            figXAxis([],'Test',[0 4],x,targetName)
            figYAxis([],'% Correct',[0 1],0:0.2:1)
            figAxForm
            hline(pctCorrThresh,'k:')
            hline(0.5,'k-')
            [~,pChance] = ttest(y',0.5);
            anovaTest = anova1(y',[],'off');
            title({'Naive Mice';sprintf('%s %s Model, anovaP=%s',mdl_naive(iav).name,...
                mdl_naive(iav).mdl(imod).name,num2str(round(anovaTest,2,'significant')))})
            for i = 1:3
                text(i,0.2,num2str(round(pChance(i),2,'significant')))
            end
            text(0,0.2,'pChnc=')
        end
    end 
end

print([fnout 'pctCorrectXStim'],'-dpdf','-fillpage')
%% cross-modality tested model performance
for iav = 1:2
    mdl_attn(iav).mdl(stimModInd).pctCorr_otherAV = nan(1,nexp_attn);
    mdl_attn(iav).mdl(choiceModInd).pctCorr_otherAV = nan(1,nexp_attn);
    for iexp = 1:nexp_attn
        mdl_attn(iav).mdl(stimModInd).pctCorr_otherAV(iexp) = ...
            dc_attn(iexp).av(iav).pctCorrectTarget_otherAV;
        mdl_attn(iav).mdl(choiceModInd).pctCorr_otherAV(iexp) = ...
            dc_attn(iexp).av(iav).pctCorrectDetect_otherAV;
    end
    mdl_noAttn(iav).mdl(stimModInd).pctCorr_otherAV = nan(1,nexp_noAttn);
    mdl_noAttn(iav).mdl(choiceModInd).pctCorr_otherAV = nan(1,nexp_noAttn);
    for iexp = 1:nexp_noAttn
        mdl_noAttn(iav).mdl(stimModInd).pctCorr_otherAV(iexp) = ...
            dc_noAttn(iexp).av(iav).pctCorrectTarget_otherAV;
        mdl_noAttn(iav).mdl(choiceModInd).pctCorr_otherAV(iexp) = ...
            dc_noAttn(iexp).av(iav).pctCorrectDetect_otherAV;
    end
    
    mdl_naive(iav).mdl(stimModInd).pctCorr_otherAV = nan(1,nexp_naive);
    for iexp = 1:nexp_naive
        mdl_naive(iav).mdl(stimModInd).pctCorr_otherAV(iexp) = ...
            dc_naive(iexp).av(iav).pctCorrectTarget_otherAV;
    end
end

xModalMdlPerf_fig = figure;
suptitle(sprintf(...
    'Cross-modality model test, pct corr must be > %s in vis stim model',...
    num2str(pctCorrThresh)))
for iav = 1:2
    for imod = 1:2
        if iav == 1
            iplot = imod;
            otherAVInd = 2;
        else
            iplot = imod+2;
            otherAVInd = 1;
        end
        subplot(2,4,iplot)
        ind = mdl_attn(iav).mdl(imod).pctCorr_all...
            > pctCorrThresh;
        y1 = mdl_attn(iav).mdl(imod).pctCorr_all(ind);
        y2 = mdl_attn(iav).mdl(imod).pctCorr_otherAV(ind);
        plot(1:2,[y1',y2'],'k-')
        hold on
        errorbar(1:2,[mean(y1),mean(y2)],[ste(y1,2),ste(y2,2)],...
            '.','MarkerSize',20)
        figXAxis([],'Test',[0 3],1:2,...
            {mdl_attn(iav).name,mdl_attn(otherAVInd).name})
        figYAxis([],'% Correct',[0 1],0:0.2:1)
        figAxForm([],0)
        hline(pctCorrThresh,'k:')
        hline(0.5,'k-')
        [~,p] = ttest2(y1,y2);
        title({'Attn Mice';sprintf('Train %s, %s Model, p=%s',mdl_attn(iav).name,...
            mdl_attn(iav).mdl(imod).name,num2str(round(p,2,'significant')))})
        [~,p_chance] = ttest([y1',y2'],0.5);
        text([0:2],ones(1,3)*.2,{'pCh',num2str(round(p_chance(1),2,'significant')),...
            num2str(round(p_chance(2),2,'significant'))})
        
        if iav == 1
            iplot = imod+4;
        else
            iplot = imod+6;
        end
        subplot(2,4,iplot)
        ind = mdl_noAttn(iav).mdl(imod).pctCorr_all...
            > pctCorrThresh;
        y1 = mdl_noAttn(iav).mdl(imod).pctCorr_all(ind);
        y2 = mdl_noAttn(iav).mdl(imod).pctCorr_otherAV(ind);
        plot(1:2,[y1',y2'],'k-')
        hold on
        errorbar(1:2,[mean(y1),mean(y2)],[ste(y1,2),ste(y2,2)],...
            '.','MarkerSize',20)
        figXAxis([],'Test',[0 3],1:2,...
            {mdl_noAttn(iav).name,mdl_noAttn(otherAVInd).name})
        figYAxis([],'% Correct',[0 1],0:0.2:1)
        figAxForm([],0)
        hline(pctCorrThresh,'k:')
        hline(0.5,'k-')
        [~,p] = ttest2(y1,y2);
        title({'No Attn Mice';sprintf('Train %s, %s Model, p=%s',mdl_noAttn(iav).name,...
            mdl_noAttn(iav).mdl(imod).name,num2str(round(p,2,'significant')))})
        [~,p_chance] = ttest([y1',y2'],0.5);
        text([0:2],ones(1,3)*.2,{'pCh',num2str(round(p_chance(1),2,'significant')),...
            num2str(round(p_chance(2),2,'significant'))})
    end 
end

print([fnout 'pctCorrectOtherAVModel'],'-dpdf','-fillpage')

figure
for iav = 1:2
    if iav == 1
        otherAVInd = 2;
    else
        otherAVInd = 1;
    end
    subplot(1,2,iav)
    ind = mdl_naive(iav).mdl(stimModInd).pctCorr_all...
            > pctCorrThresh;
    y1 = mdl_naive(iav).mdl(stimModInd).pctCorr_all(ind);
    y2 = mdl_naive(iav).mdl(stimModInd).pctCorr_otherAV(ind);
    plot(1:2,[y1',y2'],'k-')
    hold on
    errorbar(1:2,[mean(y1),mean(y2)],[ste(y1,2),ste(y2,2)],...
        '.','MarkerSize',20)
    figXAxis([],'Test',[0 3],1:2,...
        {mdl_naive(iav).name,mdl_naive(otherAVInd).name})
    figYAxis([],'% Correct',[0 1],0:0.2:1)
    figAxForm([],0)
    hline(pctCorrThresh,'k:')
    hline(0.5,'k-')
    [~,p] = ttest2(y1,y2);
    title({'Naive Mice';sprintf('Train %s, %s Model, p=%s',mdl_naive(iav).name,...
        mdl_naive(iav).mdl(stimModInd).name,num2str(round(p,2,'significant')))})
    [~,p_chance] = ttest([y1',y2'],0.5);
    text([0:2],ones(1,3)*.2,{'pCh',num2str(round(p_chance(1),2,'significant')),...
        num2str(round(p_chance(2),2,'significant'))})
end

print([fnout 'pctCorrectOtherAVModel_naive'],'-dpdf','-fillpage')
%% weights
for iav = 1:2
    mdl_attn(iav).mdl(stimModInd).weights = cell(1,nexp_attn);
    mdl_attn(iav).mdl(choiceModInd).weights = cell(1,nexp_attn);
    for iexp = 1:nexp_attn
        mdl_attn(iav).mdl(stimModInd).weights{iexp} = ...
            dc_attn(iexp).av(iav).weightTarget;
        mdl_attn(iav).mdl(choiceModInd).weights{iexp} = ...
            dc_attn(iexp).av(iav).weightDetect;
    end
    mdl_attn(iav).mdl(stimModInd).normWeights = ...
        cellfun(@(x) x./max(abs(x)), mdl_attn(iav).mdl(stimModInd).weights,'unif',0);
    mdl_attn(iav).mdl(choiceModInd).normWeights = ...
        cellfun(@(x) x./max(abs(x)), mdl_attn(iav).mdl(choiceModInd).weights,'unif',0);
    
    mdl_noAttn(iav).mdl(stimModInd).weights = cell(1,nexp_noAttn);
    mdl_noAttn(iav).mdl(choiceModInd).weights = cell(1,nexp_noAttn);
    for iexp = 1:nexp_noAttn
        mdl_noAttn(iav).mdl(stimModInd).weights{iexp} = ...
            dc_noAttn(iexp).av(iav).weightTarget;
        mdl_noAttn(iav).mdl(choiceModInd).weights{iexp} = ...
            dc_noAttn(iexp).av(iav).weightDetect;
    end
    mdl_noAttn(iav).mdl(stimModInd).normWeights = ...
        cellfun(@(x) x./max(abs(x)), mdl_noAttn(iav).mdl(stimModInd).weights,'unif',0);
    mdl_noAttn(iav).mdl(choiceModInd).normWeights = ...
        cellfun(@(x) x./max(abs(x)), mdl_noAttn(iav).mdl(choiceModInd).weights,'unif',0);
end

goodPerfInd_av_attn = mdl_attn(visualTrials).mdl(stimModInd).pctCorr_all...
    > pctCorrThresh & ...
    mdl_attn(auditoryTrials).mdl(stimModInd).pctCorr_all > pctCorrThresh;
goodPerfInd_av_noAttn = mdl_noAttn(visualTrials).mdl(stimModInd).pctCorr_all...
    > pctCorrThresh & ...
    mdl_noAttn(auditoryTrials).mdl(stimModInd).pctCorr_all > pctCorrThresh;

weights_AVcorr_attn = cell(1,2);
weights_AVcorr_noAttn = cell(1,2);
for imod = 1:2
    weights_AVcorr_attn{imod} = cellfun(@(x,y) corr(x,y),...
        mdl_attn(visualTrials).mdl(imod).normWeights(goodPerfInd_av_attn),...
        mdl_attn(auditoryTrials).mdl(imod).normWeights(goodPerfInd_av_attn));
    weights_AVcorr_noAttn{imod} = cellfun(@(x,y) corr(x,y),...
        mdl_noAttn(visualTrials).mdl(imod).normWeights(goodPerfInd_av_noAttn),...
        mdl_noAttn(auditoryTrials).mdl(imod).normWeights(goodPerfInd_av_noAttn));
end
[~,p_wAVCorr]=cellfun(@(x,y) ttest2(x,y),weights_AVcorr_attn,weights_AVcorr_noAttn);

mdlWeights_fig = figure;
suptitle({'Weights';...
    sprintf('pct corr must be > %s in vis & aud stim models',...
    num2str(pctCorrThresh))})
for imod = 1:2
    subplot(2,4,imod)
    y1 = cell2mat(mdl_attn(visualTrials).mdl(imod).normWeights(goodPerfInd_av_attn)');
    y2 = cell2mat(mdl_attn(auditoryTrials).mdl(imod).normWeights(goodPerfInd_av_attn)');
    plot(y1,y2,'k.','MarkerSize',10)
    hold on
    mdl_fitline = fitlm(y1,y2);
    rsq = mdl_fitline.Rsquared.Ordinary;
    yfit = predict(mdl_fitline,y1);
    plot(y1,yfit,'-')
    figXAxis([],'Visual Weight',[-1 1])
    figYAxis([],'Auditory Weight',[-1 1])
    figAxForm
    hline(0,'k:')
    vline(0,'k:')
    r = corr(y1,y2);
    title({'Attn Mice';sprintf('%s Model, r=%s',...
        mdl_attn(iav).mdl(imod).name,num2str(round(r,2,'significant')))})
    
    
    subplot(2,4,imod+4)
    oriEdges = [0:15:90];
    weightAVangle = rad2deg(abs(atan(y2./y1)));
    histogram(weightAVangle,oriEdges,'Normalization','probability')
    figXAxis([],'Angle of A/V Weight',[-10 100],0:15:90)
    figYAxis([],'Fraction of Cells',[0 .3])
    figAxForm
    s = skewness(weightAVangle);
    title({'Attn Mice';sprintf('%s Model, skew=%s',...
        mdl_attn(iav).mdl(imod).name,num2str(round(s,2,'significant')))})
    
    
    subplot(2,4,imod+2)
    y1 = cell2mat(mdl_noAttn(visualTrials).mdl(imod).normWeights(goodPerfInd_av_noAttn)');
    y2 = cell2mat(mdl_noAttn(auditoryTrials).mdl(imod).normWeights(goodPerfInd_av_noAttn)');
    plot(y1,y2,'k.','MarkerSize',10)
    hold on
    mdl_fitline = fitlm(y1,y2);
    rsq = mdl_fitline.Rsquared.Ordinary;
    yfit = predict(mdl_fitline,y1);
    plot(y1,yfit,'-')
    figXAxis([],'Visual Weight',[-1 1])
    figYAxis([],'Auditory Weight',[-1 1])
    figXAxis([],'Visual Weight',[])
    figYAxis([],'Auditory Weight',[])
    figAxForm
    hline(0,'k:')
    vline(0,'k:')
    r = corr(y1,y2);
    title({'No Attn Mice';sprintf('%s Model, r=%s',...
        mdl_noAttn(iav).mdl(imod).name,num2str(round(r,2,'significant')))})
    
    subplot(2,4,imod+6)
    oriEdges = [0:15:90];
    weightAVangle = rad2deg(abs(atan(y2./y1)));
    histogram(weightAVangle,oriEdges,'Normalization','probability')
    figXAxis([],'Angle of A/V Weight',[-10 100],0:15:90)
    figYAxis([],'Fraction of Cells',[0 .3])
    figAxForm
    s = skewness(weightAVangle);
    title({'No Attn Mice';sprintf('%s Model, skew=%s',...
        mdl_attn(iav).mdl(imod).name,num2str(round(s,2,'significant')))})
    
    if imod == 2
        figure
        ind = find(round(weightAVangle)==60,1);
        plot(y1(ind),y2(ind),'.','MarkerSize',20)
        hold on
        plot([0 y1(ind)],[0 y2(ind)],'-')
        figXAxis([],'Visual Weight',[-.1 .6])
        figYAxis([],'Auditory Weight',[-.1 .6])
        figAxForm
        hline(0,'k:')
        vline(0,'k:')
        title({'No Attn Mice';sprintf('%s Model, AVangle=%s',...
            mdl_attn(iav).mdl(imod).name,...
            num2str(round(weightAVangle(ind),2,'significant')))})
        print([fnout 'modelWeights_exCellAngle'],'-dpdf')
        figure(mdlWeights_fig)
    end
    
end 
print([fnout 'modelWeights'],'-dpdf','-fillpage')

%% attention effect on model and behavior
weights_modCorr_attn = cell(1,2);
weights_modCorr_noAttn = cell(1,2);
pctCorrectDiff_bwModels_attn = cell(1,2);
pctCorrectDiff_bwModels_noAttn = cell(1,2);
pctCorrectDiff_xMod_attn = cell(2,2);
pctCorrectDiff_xMod_noAttn = cell(2,2);
av_mod_gpi_attn = cell(2,2);
av_mod_gpi_noAttn = cell(2,2);
for iav = 1:2
    weights_modCorr_attn{iav} = cellfun(@(x,y) corr(x,y),...
        mdl_attn(visualTrials).mdl(stimModInd).normWeights(goodPerfInd_av_attn),...
        mdl_attn(auditoryTrials).mdl(choiceModInd).normWeights(goodPerfInd_av_attn));
    weights_modCorr_noAttn{iav} = cellfun(@(x,y) corr(x,y),...
        mdl_noAttn(iav).mdl(stimModInd).normWeights(goodPerfInd_av_noAttn),...
        mdl_noAttn(iav).mdl(choiceModInd).normWeights(goodPerfInd_av_noAttn));
    pctCorrectDiff_bwModels_attn{iav} = abs(mdl_attn(iav).mdl(stimModInd).pctCorr_all - ...
        mdl_attn(iav).mdl(choiceModInd).pctCorr_all);
    pctCorrectDiff_bwModels_noAttn{iav} = abs(mdl_noAttn(iav).mdl(stimModInd).pctCorr_all - ...
        mdl_noAttn(iav).mdl(choiceModInd).pctCorr_all);
    for imod = 1:2
        pctCorrectDiff_xMod_attn{iav,imod} = mdl_attn(iav).mdl(imod).pctCorr_all - ...
            mdl_attn(iav).mdl(imod).pctCorr_otherAV;
        av_mod_gpi_attn{iav,imod} = ...
            mdl_attn(iav).mdl(imod).pctCorr_all>pctCorrThresh;
        pctCorrectDiff_xMod_noAttn{iav,imod} = mdl_noAttn(iav).mdl(imod).pctCorr_all - ...
            mdl_noAttn(iav).mdl(imod).pctCorr_otherAV;
        av_mod_gpi_noAttn{iav,imod} = ...
            mdl_noAttn(iav).mdl(imod).pctCorr_all>pctCorrThresh;
    end
end
[~,p_wModCorr]=cellfun(@(x,y) ttest2(x,y),weights_modCorr_attn,weights_modCorr_noAttn);

ind = ismember(bxStats.mouseNames,unique(mouseName_attn));
[~, bxDataSortInd] = sort(bxStats.mouseNames(ind));
hrDiffAll = bxStats.av(visualTrials).matchedHRDiff(ind);
hrDiff_attn = hrDiffAll(bxDataSortInd);
ind = ismember(bxStats.mouseNames,unique(mouseName_noAttn));
[~, bxDataSortInd] = sort(bxStats.mouseNames(ind));
hrDiffAll = bxStats.av(visualTrials).matchedHRDiff(ind);
hrDiff_noAttn = hrDiffAll(bxDataSortInd);

mice_attn = unique(mouseName_attn);
pctCorrDiff_bwModels_eaMouse_attn = nan(length(mice_attn),2);
pctCorrDiff_xMod_eaMouse_attn = cell(2,2);
for im = 1:length(mice_attn)
    ind = strcmp(mouseName_attn,mice_attn(im)) & goodPerfInd_attn;
    pctCorrDiff_bwModels_eaMouse_attn(im,:) = cellfun(@(x) mean(x(ind)),...
        pctCorrectDiff_bwModels_attn);
    pctCorrDiff_xMod_eaMouse_attn = cellfun(@(x,y,z)...
        cat(1,x,mean(y(z & strcmp(mouseName_attn,mice_attn(im))))),...
        pctCorrDiff_xMod_eaMouse_attn,pctCorrectDiff_xMod_attn,...
        av_mod_gpi_attn,'unif',0);
end
mice_noAttn = unique(mouseName_noAttn);
pctCorrDiff_bwModels_eaMouse_noAttn = nan(length(mice_noAttn),2);
pctCorrDiff_xMod_eaMouse_noAttn = cell(2,2);
for im = 1:length(mice_noAttn)
    ind = strcmp(mouseName_noAttn,mice_noAttn(im)) & goodPerfInd_noAttn;
    pctCorrDiff_bwModels_eaMouse_noAttn(im,:) = cellfun(@(x) mean(x(ind)),...
        pctCorrectDiff_bwModels_noAttn);
    pctCorrDiff_xMod_eaMouse_noAttn = cellfun(@(x,y,z)...
        cat(1,x,mean(y(z & strcmp(mouseName_noAttn,mice_noAttn(im))))),...
        pctCorrDiff_xMod_eaMouse_noAttn,pctCorrectDiff_xMod_noAttn,...
        av_mod_gpi_noAttn,'unif',0);
end

figure
for iav = 1:2
    subplot(2,2,iav)
    plot(weights_modCorr_attn{iav},pctCorrectDiff_bwModels_attn{iav}(goodPerfInd_av_attn),'.','MarkerSize',20)
    hold on
    plot(weights_modCorr_noAttn{iav},pctCorrectDiff_bwModels_noAttn{iav}(goodPerfInd_av_noAttn),'.','MarkerSize',20)
    figXAxis([],'Stimuls to Choice Weight Correlation',[-0.5 1])
    figYAxis([],'Stimulus-Choice Pct Corect (abs val)',[0 0.2])
    figAxForm
    title(avName{iav})
    
    subplot(2,2,iav+2)
    plot(hrDiff_attn,pctCorrDiff_bwModels_eaMouse_attn(:,iav),'.','MarkerSize',20)
    hold on
    plot(hrDiff_noAttn,pctCorrDiff_bwModels_eaMouse_noAttn(:,iav),'.','MarkerSize',20)
    figXAxis([],'Val-Inv. Hit Rate',[0 0.5])
    figYAxis([],'Stimulus-Choice Pct Corect (abs val)',[0 0.2])
    figAxForm
end

figure
for iav = 1:2
    for imod = 1:2
        subplot(2,2,iav+((imod-1)*2))
        plot(hrDiff_attn,pctCorrDiff_xMod_eaMouse_attn{iav,imod},'.','MarkerSize',20)
        hold on
        plot(hrDiff_noAttn,pctCorrDiff_xMod_eaMouse_noAttn{iav,imod},'.','MarkerSize',20)
        x = cat(1,hrDiff_attn,hrDiff_noAttn);
        y = cat(1,pctCorrDiff_xMod_eaMouse_attn{iav,imod},...
            pctCorrDiff_xMod_eaMouse_noAttn{iav,imod});
        ind = ~isnan(x)&~isnan(y);
        [r,p] = corr(x(ind),y(ind));
        figXAxis([],'Val-Inv. Hit Rate',[-0.2 0.5])
        figYAxis([],'Within - Cross Modality Pct Corr',[-0.05 0.35])
        figAxForm
        title(sprintf('%s %s Model, r=%s, p=%s',avName{iav},modelName{imod},...
            sigfigString(r),sigfigString(p)))
    end
end
print([fnout 'xModalModelPerfbyAttnBehavior'],'-dpdf','-fillpage')

%% example experiment PCs
exExpt = 11;
npcs = dc_attn(exExpt).nPCs;
stimExPCs_fig = figure;
colormap(brewermap([],'Blues'))
choiceExPCs_fig = figure;
colormap(brewermap([],'Reds'))
weightsExPCs_fig = figure;
colormap gray
pcNum = cellfun(@num2str,num2cell(1:npcs),'unif',0)';
for iav = 1:2
    pcs = dc_attn(exExpt).av(iav).respAllCells(:,1:npcs);
    pc_wChoice = dc_attn(exExpt).av(iav).weightDetect_pcs(1:npcs);
    pc_wStim = dc_attn(exExpt).av(iav).weightTarget_pcs(1:npcs);
    trOut = dc_attn(exExpt).av(iav).trOut;
    [choiceTrOut, stimTrOut] = getStimAndBehaviorYs(trOut);
    choiceResp = nan(npcs,2);
    stimResp = nan(npcs,2);
    for i = 1:2
        choiceResp(:,i) = mean(pcs(choiceTrOut==(i-1),:),1);
        stimResp(:,i) = mean(pcs(stimTrOut==(i-1),:),1);
    end
    if iav == 1
%         [~,sortInd_choice] = sort(pc_wChoice);
        [~,sortInd_stim] = sort(pc_wStim);
    end
    
    pcLabel_stim = cellfun(@(x,y) [x '-' y],...
        cellfun(@num2str,num2cell(round(pc_wStim(sortInd_stim),2,'significant')),'unif',0),...
        pcNum(sortInd_stim),'unif',0);
    pcLabel_choice = cellfun(@(x,y) [x '-' y],...
        cellfun(@num2str,num2cell(round(pc_wChoice(sortInd_stim),2,'significant')),'unif',0),...
        pcNum(sortInd_stim),'unif',0);
    
    figure(stimExPCs_fig)
    subplot(1,2,iav)
    imagesc(stimResp(sortInd_stim,:))
	clim([-0.6 1.2])
    colorbar
    figXAxis([],'Stimulus',[],1:2,{'Dist','Tar'})
    figYAxis([],'PC #',[],1:npcs,pcLabel_stim)
    figAxForm([],0)
    title(avName{iav})
    
    figure(choiceExPCs_fig)
    subplot(1,2,iav)
    imagesc(choiceResp(sortInd_stim,:))
    clim([-0.6 1.2])
    colorbar
    figXAxis([],'Choice',[],1:2,{'No','Yes'})
    figYAxis([],'PC #',[],1:npcs,pcLabel_choice)
    figAxForm([],0)
    title(avName{iav})
    
    figure(weightsExPCs_fig)
    subplot(2,2,iav)
    imagesc(pc_wStim(sortInd_stim))
    clim([-1 3])
    colorbar
    figXAxis([],'',[])
    figYAxis([],'Stimulus Weight',[])
    figAxForm([],0)
    title(avName{iav})
    subplot(2,2,iav+2)
    imagesc(pc_wChoice(sortInd_stim))
    clim([-1 3])
    colorbar
    figXAxis([],'',[])
    figYAxis([],'Choice Weight',[])
    figAxForm([],0)
    title(avName{iav})
end
figure(stimExPCs_fig)
print([fnout 'exExpt_pcHM_stim'],'-dpdf')
figure(choiceExPCs_fig)
print([fnout 'exExpt_pcHM_choice'],'-dpdf')
figure(weightsExPCs_fig)
print([fnout 'exExpt_pcHM_weights'],'-dpdf')

%% catch trial performance
catchExptInd_attn = true(1,nexp_attn);
for imod = 1:2
    mdl_attn(visualTrials).mdl(imod).pctCorrect_inv = nan(1,nexp_attn);
    mdl_attn(visualTrials).mdl(imod).pctCorrect_valMatched = nan(1,nexp_attn);
    mdl_attn(visualTrials).mdl(imod).pctCorrect_invAudTest = nan(1,nexp_attn);
    mdl_attn(visualTrials).mdl(imod).pctCorrect_valMatchedAudTest = nan(1,nexp_attn);
    mdl_attn(visualTrials).mdl(imod).catchTrialResult_inv = cell(1,nexp_attn);
    mdl_attn(visualTrials).mdl(imod).catchTrialResult_valMatched = cell(1,nexp_attn);
    mdl_attn(visualTrials).mdl(imod).catchTrialResult_invAudTest = cell(1,nexp_attn);
    mdl_attn(visualTrials).mdl(imod).catchTrialResult_valMatchedAudTest = cell(1,nexp_attn);
    
    for iexp = 1:nexp_attn
        if isempty(dc_attn(iexp).av(visualTrials).invalidPctCorrectDetect)
            catchExptInd_attn(iexp) = false;
            continue
        end
        if imod == 1
            modType = 'Target';
        else
            modType = 'Detect';
        end
        mdl_attn(visualTrials).mdl(imod).pctCorrect_inv(iexp) = ...
            eval(['dc_attn(iexp).av(visualTrials).invalidPctCorrect' modType]);
        mdl_attn(visualTrials).mdl(imod).pctCorrect_valMatched(iexp) = ...
            eval(['dc_attn(iexp).av(visualTrials).validMatchedPctCorrect'...
            modType '_holdout']);            
        mdl_attn(visualTrials).mdl(imod).pctCorrect_invAudTest(iexp) = ...
            eval(['dc_attn(iexp).av(visualTrials).invalidPctCorrect'...
            modType '_testAudModel']);  
        mdl_attn(visualTrials).mdl(imod).pctCorrect_valMatchedAudTest(iexp) = ...
            eval(['dc_attn(iexp).av(visualTrials).validMatchedPctCorrect'...
            modType '_testAudModel']);  
        
        mdl_attn(visualTrials).mdl(imod).catchTrialResult_inv{iexp} = ...
            eval(['dc_attn(iexp).av(visualTrials).invalidCorrectTrials' modType])';
        mdl_attn(visualTrials).mdl(imod).catchTrialResult_valMatched{iexp} = ...
            eval(['dc_attn(iexp).av(visualTrials).validMatchedCorrectTrials' modType])';            
        mdl_attn(visualTrials).mdl(imod).catchTrialResult_invAudTest{iexp} = ...
            eval(['dc_attn(iexp).av(visualTrials).invalidCorrectTrials'...
            modType '_testAud'])';  
        mdl_attn(visualTrials).mdl(imod).catchTrialResult_valMatchedAudTest{iexp} = ...
            eval(['dc_attn(iexp).av(visualTrials).validMatchedCorrectTrials'...
            modType '_testAud'])';  
    end
end

catchExptInd_noAttn = true(1,nexp_noAttn);
for imod = 1:2
    mdl_noAttn(visualTrials).mdl(imod).pctCorrect_inv = nan(1,nexp_noAttn);
    mdl_noAttn(visualTrials).mdl(imod).pctCorrect_valMatched = nan(1,nexp_noAttn);
    mdl_noAttn(visualTrials).mdl(imod).pctCorrect_invAudTest = nan(1,nexp_noAttn);
    mdl_noAttn(visualTrials).mdl(imod).pctCorrect_valMatchedAudTest = nan(1,nexp_noAttn);
    mdl_noAttn(visualTrials).mdl(imod).catchTrialResult_inv = cell(1,nexp_noAttn);
    mdl_noAttn(visualTrials).mdl(imod).catchTrialResult_valMatched = cell(1,nexp_noAttn);
    mdl_noAttn(visualTrials).mdl(imod).catchTrialResult_invAudTest = cell(1,nexp_noAttn);
    mdl_noAttn(visualTrials).mdl(imod).catchTrialResult_valMatchedAudTest = cell(1,nexp_noAttn);
    
    for iexp = 1:nexp_noAttn
        if isempty(dc_noAttn(iexp).av(visualTrials).invalidPctCorrectDetect)
            catchExptInd_noAttn(iexp) = false;
            continue
        end
        if imod == 1
            modType = 'Target';
        else
            modType = 'Detect';
        end
        mdl_noAttn(visualTrials).mdl(imod).pctCorrect_inv(iexp) = ...
            eval(['dc_noAttn(iexp).av(visualTrials).invalidPctCorrect' modType]);
        mdl_noAttn(visualTrials).mdl(imod).pctCorrect_valMatched(iexp) = ...
            eval(['dc_noAttn(iexp).av(visualTrials).validMatchedPctCorrect'...
            modType '_holdout']);            
        mdl_noAttn(visualTrials).mdl(imod).pctCorrect_invAudTest(iexp) = ...
            eval(['dc_noAttn(iexp).av(visualTrials).invalidPctCorrect'...
            modType '_testAudModel']);  
        mdl_noAttn(visualTrials).mdl(imod).pctCorrect_valMatchedAudTest(iexp) = ...
            eval(['dc_noAttn(iexp).av(visualTrials).validMatchedPctCorrect'...
            modType '_testAudModel']);  
        
        mdl_noAttn(visualTrials).mdl(imod).catchTrialResult_inv{iexp} = ...
            eval(['dc_noAttn(iexp).av(visualTrials).invalidCorrectTrials' modType])';
        mdl_noAttn(visualTrials).mdl(imod).catchTrialResult_valMatched{iexp} = ...
            eval(['dc_noAttn(iexp).av(visualTrials).validMatchedCorrectTrials' modType])';            
        mdl_noAttn(visualTrials).mdl(imod).catchTrialResult_invAudTest{iexp} = ...
            eval(['dc_noAttn(iexp).av(visualTrials).invalidCorrectTrials'...
            modType '_testAud'])';  
        mdl_noAttn(visualTrials).mdl(imod).catchTrialResult_valMatchedAudTest{iexp} = ...
            eval(['dc_noAttn(iexp).av(visualTrials).validMatchedCorrectTrials'...
            modType '_testAud'])';  
    end
end

figure
for imod = 1:2
    ind = mdl_attn(visualTrials).mdl(imod).pctCorrect_valMatched > pctCorrThresh ...
        & catchExptInd_attn;
    y = cat(1,mdl_attn(visualTrials).mdl(imod).pctCorrect_valMatched(ind),...
        mdl_attn(visualTrials).mdl(imod).pctCorrect_inv(ind),...
        mdl_attn(visualTrials).mdl(imod).pctCorrect_invAudTest(ind));
    
    subplot(2,2,imod)
    plot(1:3,y','k.-','MarkerSize',10)
    hold on
    errorbar(1:3,mean(y',1),ste(y',1),'.','MarkerSize',20)
    figXAxis([],'Train-Test',[0 4],1:3,{'Vis-ValVis','Vis-InvVis','Aud-InvVis'})
    figYAxis([],'Pct Correct',[0 1])
    figAxForm
    title(sprintf('Attn, %s Model',modelName{imod}))
    hline(0.5,'k:')
    
    ind = mdl_noAttn(visualTrials).mdl(imod).pctCorrect_valMatched > pctCorrThresh ...
        & catchExptInd_noAttn;
    y = cat(1,mdl_noAttn(visualTrials).mdl(imod).pctCorrect_valMatched(ind),...
        mdl_noAttn(visualTrials).mdl(imod).pctCorrect_inv(ind),...
        mdl_noAttn(visualTrials).mdl(imod).pctCorrect_invAudTest(ind));
    
    subplot(2,2,imod+2)
    plot(1:3,y','k.-','MarkerSize',10)
    hold on
    errorbar(1:3,mean(y',1),ste(y',1),'.','MarkerSize',20)
    figXAxis([],'Train-Test',[0 4],1:3,{'Vis-ValVis','Vis-InvVis','Aud-InvVis'})
    figYAxis([],'Pct Correct',[0 1])
    figAxForm
    title(sprintf('No Attn, %s Model',modelName{imod}))
    hline(0.5,'k:')
end
print([fnout 'invalidTrialPerf'],'-dpdf','-fillpage')

figure
for imod = 1:2
    subplot(2,2,imod)
    ind = cellfun(@(x) ~isempty(x),mdl_attn(visualTrials).mdl(imod).catchTrialResult_inv) & ...
        mdl_attn(visualTrials).mdl(imod).pctCorr_all > pctCorrThresh;
    invTrialResult = cell2mat(mdl_attn(visualTrials).mdl(imod).catchTrialResult_inv(ind));
    [invHR,invCI] = binofit(sum(invTrialResult),length(invTrialResult));
    valTrialResult = cell2mat(mdl_attn(visualTrials).mdl(imod).catchTrialResult_valMatched(ind));
    [valHR,valCI] = binofit(sum(valTrialResult),length(valTrialResult));
    invTrialResult_audModel = cell2mat(mdl_attn(visualTrials).mdl(imod).catchTrialResult_invAudTest(ind));
    [invHR_audModel,invCI_audModel] = binofit(sum(invTrialResult_audModel),length(invTrialResult_audModel));
    valTrialResult_audModel = cell2mat(mdl_attn(visualTrials).mdl(imod).catchTrialResult_valMatchedAudTest(ind));
    [valHR_audModel,valCI_audModel] = binofit(sum(valTrialResult_audModel),length(valTrialResult_audModel));
    
    y = [valHR,invHR,valHR_audModel,invHR_audModel];
    y_lerr = y - [valCI(1),invCI(1),valCI_audModel(1),invCI_audModel(1)];
    y_uerr = [valCI(2),invCI(2),valCI_audModel(2),invCI_audModel(2)]-y;
    errorbar(1:4,y, y_lerr,y_uerr,'.','MarkerSize',20)
    figXAxis([],'Train-Test',[0 5],1:4,{'Vis-Val','Vis-Inv','Aud-Val','Aud-Inv'})
    figYAxis([],'Pct Correct',[0 1])
    figAxForm
    hline(0.5,'k:')
    title(sprintf('Attn; %s Model, Vis Trials',modelName{imod}))
    
    subplot(2,2,imod+2)
    ind = cellfun(@(x) ~isempty(x),mdl_noAttn(visualTrials).mdl(imod).catchTrialResult_inv) & ...
        mdl_noAttn(visualTrials).mdl(imod).pctCorr_all > pctCorrThresh;
    invTrialResult = cell2mat(mdl_noAttn(visualTrials).mdl(imod).catchTrialResult_inv(ind));
    [invHR,invCI] = binofit(sum(invTrialResult),length(invTrialResult));
    valTrialResult = cell2mat(mdl_noAttn(visualTrials).mdl(imod).catchTrialResult_valMatched(ind));
    [valHR,valCI] = binofit(sum(valTrialResult),length(valTrialResult));
    invTrialResult_audModel = cell2mat(mdl_noAttn(visualTrials).mdl(imod).catchTrialResult_invAudTest(ind));
    [invHR_audModel,invCI_audModel] = binofit(sum(invTrialResult_audModel),length(invTrialResult_audModel));
    valTrialResult_audModel = cell2mat(mdl_noAttn(visualTrials).mdl(imod).catchTrialResult_valMatchedAudTest(ind));
    [valHR_audModel,valCI_audModel] = binofit(sum(valTrialResult_audModel),length(valTrialResult_audModel));
    
    y = [valHR,invHR,valHR_audModel,invHR_audModel];
    y_lerr = y - [valCI(1),invCI(1),valCI_audModel(1),invCI_audModel(1)];
    y_uerr = [valCI(2),invCI(2),valCI_audModel(2),invCI_audModel(2)]-y;
    errorbar(1:4,y, y_lerr,y_uerr,'.','MarkerSize',20)
    figXAxis([],'Train-Test',[0 5],1:4,{'Vis-Val','Vis-Inv','Aud-Val','Aud-Inv'})
    figYAxis([],'Pct Correct',[0 1])
    figAxForm
    hline(0.5,'k:')
    title(sprintf('No Attn; %s Model, Vis Trials',modelName{imod}))
end
print([fnout 'invalidTrialPerf_combExpt'],'-dpdf','-fillpage')

figure
ind = goodPerfInd_attn & catchExptInd_attn;
x = mdl_attn(visualTrials).mdl(stimModInd).pctCorrect_valMatched(ind) -...
    mdl_attn(visualTrials).mdl(stimModInd).pctCorrect_inv(ind);
y = mdl_attn(visualTrials).mdl(choiceModInd).pctCorrect_valMatched(ind) -...
    mdl_attn(visualTrials).mdl(choiceModInd).pctCorrect_inv(ind);
subplot(2,2,1)
plot(x,y,'.','MarkerSize',20)
hold on
plot([-0.5 0.7],[-0.5 0.7],'k:')
figXAxis([],'Stim Model, Val-Inv',[-0.2 0.7])
figYAxis([],'Choice Model, Val-Inv',[-0.2 0.7])
figAxForm
vline(0,'k:')
hline(0,'k:')
title('No Attn')

valInvDiff_stim_eaMouse = nan(1,length(mice_attn));
valInvDiff_choice_eaMouse = nan(1,length(mice_attn));
for im = 1:length(mice_attn)
    ind = goodPerfInd_attn & catchExptInd_attn & ...
        strcmp(mouseName_attn,mice_attn(im));
    x = mdl_attn(visualTrials).mdl(stimModInd).pctCorrect_valMatched(ind) -...
        mdl_attn(visualTrials).mdl(stimModInd).pctCorrect_inv(ind);
    y = mdl_attn(visualTrials).mdl(choiceModInd).pctCorrect_valMatched(ind) -...
        mdl_attn(visualTrials).mdl(choiceModInd).pctCorrect_inv(ind);
    valInvDiff_stim_eaMouse(im) = mean(x);
    valInvDiff_choice_eaMouse(im) = mean(y);
end

subplot(2,2,3)
plot(hrDiff_attn,valInvDiff_stim_eaMouse,'.','MarkerSize',20)
hold on
figXAxis([],'HR Val-Inv',[-0 0.5])
figYAxis([],'Stim Model, Val-Inv',[-0.2 0.7])
figAxForm
vline(0,'k:')
hline(0,'k:')
subplot(2,2,4)
plot(hrDiff_attn,valInvDiff_choice_eaMouse,'.','MarkerSize',20)
hold on
figXAxis([],'HR Val-Inv',[-0 0.5])
figYAxis([],'Stim Model, Val-Inv',[-0.2 0.7])
figAxForm
vline(0,'k:')
hline(0,'k:')

ind = goodPerfInd_noAttn & catchExptInd_noAttn;
x = mdl_noAttn(visualTrials).mdl(stimModInd).pctCorrect_valMatched(ind) -...
    mdl_noAttn(visualTrials).mdl(stimModInd).pctCorrect_inv(ind);
y = mdl_noAttn(visualTrials).mdl(choiceModInd).pctCorrect_valMatched(ind) -...
    mdl_noAttn(visualTrials).mdl(choiceModInd).pctCorrect_inv(ind);
subplot(2,2,2)
plot(x,y,'.','MarkerSize',20)
hold on
plot([-0.5 0.7],[-0.5 0.7],'k:')
figXAxis([],'Stim Model, Val-Inv',[-0.2 0.7])
figYAxis([],'Choice Model, Val-Inv',[-0.2 0.7])
figAxForm
vline(0,'k:')
hline(0,'k:')
title('No Attn')

valInvDiff_stim_eaMouse = nan(1,length(mice_noAttn));
valInvDiff_choice_eaMouse = nan(1,length(mice_noAttn));
for im = 1:length(mice_noAttn)
    ind = goodPerfInd_noAttn & catchExptInd_noAttn & ...
        strcmp(mouseName_noAttn,mice_noAttn(im));
    x = mdl_noAttn(visualTrials).mdl(stimModInd).pctCorrect_valMatched(ind) -...
        mdl_noAttn(visualTrials).mdl(stimModInd).pctCorrect_inv(ind);
    y = mdl_noAttn(visualTrials).mdl(choiceModInd).pctCorrect_valMatched(ind) -...
        mdl_noAttn(visualTrials).mdl(choiceModInd).pctCorrect_inv(ind);
    valInvDiff_stim_eaMouse(im) = mean(x);
    valInvDiff_choice_eaMouse(im) = mean(y);
end

subplot(2,2,3)
plot(hrDiff_noAttn,valInvDiff_stim_eaMouse,'.','MarkerSize',20)
hold on
figXAxis([],'HR Val-Inv',[-0 0.5])
figYAxis([],'Stim Model, Val-Inv',[-0.2 0.7])
figAxForm
vline(0,'k:')
hline(0,'k:')
legend({'Attn','No Attn'})
subplot(2,2,4)
plot(hrDiff_noAttn,valInvDiff_choice_eaMouse,'.','MarkerSize',20)
hold on
figXAxis([],'HR Val-Inv',[-0 0.5])
figYAxis([],'Stim Model, Val-Inv',[-0.2 0.7])
figAxForm
vline(0,'k:')
hline(0,'k:')
legend({'Attn','No Attn'})
print([fnout 'invalidTrialPerfxBehavior'],'-dpdf','-fillpage')

%% distractor only choice model performance
for iav = 1:2
    mdl_attn(iav).mdl(choiceModInd).pctCorr_distOnly = nan(1,nexp_attn);
    mdl_attn(iav).mdl(choiceModInd).pctCorr_distOnly_otherAV = nan(1,nexp_attn);
    for iexp = 1:nexp_attn
        mdl_attn(iav).mdl(choiceModInd).pctCorr_distOnly(iexp) = ...
            dc_attn(iexp).av(iav).pctCorrectDetect_distOnly_holdout;
        mdl_attn(iav).mdl(choiceModInd).pctCorr_distOnly_otherAV(iexp) = ...
            dc_attn(iexp).av(iav).pctCorrectDetect_distOnly_otherAV;
    end
    
    mdl_noAttn(iav).mdl(choiceModInd).pctCorr_distOnly = nan(1,nexp_noAttn);
    mdl_noAttn(iav).mdl(choiceModInd).pctCorr_distOnly_otherAV = nan(1,nexp_noAttn);
    for iexp = 1:nexp_noAttn
        mdl_noAttn(iav).mdl(choiceModInd).pctCorr_distOnly(iexp) = ...
            dc_noAttn(iexp).av(iav).pctCorrectDetect_distOnly_holdout;
        mdl_noAttn(iav).mdl(choiceModInd).pctCorr_distOnly_otherAV(iexp) = ...
            dc_noAttn(iexp).av(iav).pctCorrectDetect_distOnly_otherAV;
    end
end

figure
suptitle('Distractor Only Trained')
for iav = 1:2
    if iav == 1
        otherAV = 2;
    else
        otherAV = 1;
    end
    subplot(2,2,iav)
    y = cat(1,mdl_attn(iav).mdl(choiceModInd).pctCorr_distOnly,...
        mdl_attn(iav).mdl(choiceModInd).pctCorr_distOnly_otherAV);
    ind = y(1,:)>pctCorrThresh;
    plot(1:2,y(:,ind),'k-')
    hold on
    errorbar(1:2,mean(y(:,ind),2),ste(y(:,ind),2),'.','MarkerSize',20);
    [~,p] = ttest(y(1,:),y(2,:));
    [~,p_chance] = ttest(y');
    figXAxis([],'Test',[0 3],1:2,{avName{iav},avName{otherAV}})
    figYAxis([],'Pct Correct',[0 1])
    figAxForm
    hline(0.5,'k-')
    hline(pctCorrThresh,'k:')
    title(sprintf('Attn, %s Choice Model, p=%s',avName{iav},num2str(round(p,2,'significant'))))    
    text([0:2],ones(1,3)*.2,{'pCh',num2str(round(p_chance(1),2,'significant')),...
        num2str(round(p_chance(2),2,'significant'))})
    
    subplot(2,2,iav+2)
    y = cat(1,mdl_noAttn(iav).mdl(choiceModInd).pctCorr_distOnly,...
        mdl_noAttn(iav).mdl(choiceModInd).pctCorr_distOnly_otherAV);
    ind = y(1,:)>pctCorrThresh;
    plot(1:2,y(:,ind),'k-')
    hold on
    errorbar(1:2,mean(y(:,ind),2),ste(y(:,ind),2),'.','MarkerSize',20);
    [~,p] = ttest(y(1,:),y(2,:));
    [~,p_chance] = ttest(y');
    figXAxis([],'Test',[0 3],1:2,{avName{iav},avName{otherAV}})
    figYAxis([],'Pct Correct',[0 1])
    figAxForm
    hline(0.5,'k-')
    hline(pctCorrThresh,'k:')
    title(sprintf('No Attn, %s Choice Model, p=%s',avName{iav},num2str(round(p,2,'significant'))))    
    text([0:2],ones(1,3)*.2,{'pCh',num2str(round(p_chance(1),2,'significant')),...
        num2str(round(p_chance(2),2,'significant'))})
end
print([fnout 'pctCorrectOtherAVModel_distOnly'],'-dpdf','-fillpage')

figure
suptitle('Distractor Only Trained')
for iav = 1:2
    subplot(1,2,iav)
    y1 = mdl_attn(iav).mdl(choiceModInd).pctCorr_distOnly(goodPerfInd_attn);
    plot(ones(1,sum(goodPerfInd_attn)),y1,'k.','MarkerSize',10)
    hold on
    errorbar(1,mean(y1),ste(y1,2),'.','MarkerSize',20)
    y2 = mdl_noAttn(iav).mdl(choiceModInd).pctCorr_distOnly(goodPerfInd_noAttn);
    plot(ones(1,sum(goodPerfInd_noAttn)).*2,y2,'k.','MarkerSize',10)
    hold on
    errorbar(2,mean(y2),ste(y2,2),'.','MarkerSize',20)
    [~,p] = ttest2(y1,y2);
    figXAxis([],'Test',[0 3],1:2,{'Attn','No Attn'})
    figYAxis([],'Pct Correct',[0 1])
    figAxForm
    hline(0.5,'k-')
    hline(pctCorrThresh,'k:')
    title(sprintf('%s Choice Model',avName{iav}))     
end
print([fnout 'pctCorrect_distOnly'],'-dpdf','-fillpage')    

figure
suptitle('Distractor Only Training')
for iav = 1:2
    y = mdl_attn(iav).mdl(choiceModInd).pctCorr_distOnly - ...
        mdl_attn(iav).mdl(choiceModInd).pctCorr_distOnly_otherAV;
    y_ms_attn = nan(1,length(mice_attn));
    for im = 1:length(mice_attn)
        ind = mdl_attn(iav).mdl(choiceModInd).pctCorr_distOnly >pctCorrThresh & ...
            strcmp(mouseName_attn,mice_attn(im));
        y_ms_attn(im) = mean(y(ind));
    end    
    subplot(1,2,iav)
    plot(hrDiff_attn,y_ms_attn,'.','MarkerSize',20)
    hold on
    
    y = mdl_noAttn(iav).mdl(choiceModInd).pctCorr_distOnly - ...
        mdl_noAttn(iav).mdl(choiceModInd).pctCorr_distOnly_otherAV;
    y_ms_noAttn = nan(1,length(mice_noAttn));
    for im = 1:length(mice_noAttn)
        ind = mdl_noAttn(iav).mdl(choiceModInd).pctCorr_distOnly >pctCorrThresh & ...
            strcmp(mouseName_noAttn,mice_noAttn(im));
        y_ms_noAttn(im) = mean(y(ind));
    end    
    plot(hrDiff_noAttn,y_ms_noAttn,'.','MarkerSize',20)
    
    r = corr([hrDiff_attn;hrDiff_noAttn],[y_ms_attn';y_ms_noAttn']);
    
    figXAxis([],'HR Val-Inv',[-0.2 0.5])
    figYAxis([],'Within-Other Pct Correct',[])
    figAxForm
    title(sprintf('%s Choice Model, r=%s',avName{iav},...
        num2str(round(r,2,'significant'))))
    legend({'Attn','No Attn'})
end
print([fnout 'pctCorrectXbehavior_distOnly'],'-dpdf','-fillpage')   
%% performance across response windows
nBaselineFr = 30;
nwins=15;
winInd = [1:2,6:nwins];
movWinLabelFr = 30:(30+nwins-1);
movWinLabelMs = (movWinLabelFr - (nBaselineFr+nVisDelayFr)).*(1000/frameRateHz);
respWinTT = ([respwin(1) respwin(end)] - (nBaselineFr+nVisDelayFr))...
    .*(1000/frameRateHz);
baseWinTT = (...
    [basewin_0(1) basewin_0(end)] - (nBaselineFr+nVisDelayFr))...
    .*(1000/frameRateHz);

for iav = 1:2
    mdl_attn(iav).mdl(stimModInd).pctCorr_respWin = nan(nwins,nexp_attn);
    mdl_attn(iav).mdl(choiceModInd).pctCorr_respWin = nan(nwins,nexp_attn);
    for iexp = 1:nexp_attn
        mdl_attn(iav).mdl(stimModInd).pctCorr_respWin(:,iexp) = ...
            dc_attn(iexp).av(iav).pctCorrectTargetMovRespWin_holdout;        
        mdl_attn(iav).mdl(choiceModInd).pctCorr_respWin(:,iexp) = ...
            dc_attn(iexp).av(iav).pctCorrectDetectMovRespWin_holdout;
    end
    mdl_noAttn(iav).mdl(stimModInd).pctCorr_respWin = nan(nwins,nexp_noAttn);
    mdl_noAttn(iav).mdl(choiceModInd).pctCorr_respWin = nan(nwins,nexp_noAttn);
    for iexp = 1:nexp_noAttn
        mdl_noAttn(iav).mdl(stimModInd).pctCorr_respWin(:,iexp) = ...
            dc_noAttn(iexp).av(iav).pctCorrectTargetMovRespWin_holdout;        
        mdl_noAttn(iav).mdl(choiceModInd).pctCorr_respWin(:,iexp) = ...
            dc_noAttn(iexp).av(iav).pctCorrectDetectMovRespWin_holdout;
    end
end

figure
x = movWinLabelMs(winInd);
targetTimes = x(x>0);
for iav = 1:2
    for imod = 1:2
        if imod == 1
            iplot = 0;
        else
            iplot = 2;
        end
        subplot(4,2,iplot+iav)
        hold on
        y = mdl_attn(iav).mdl(imod).pctCorr_respWin(winInd,:);
        ind = mdl_attn(iav).mdl(imod).pctCorr_all > pctCorrThresh;
        h = plot(x,y(:,ind),'k-');
        errorbar(x,mean(y(:,ind),2),ste(y(:,ind),2),'.-')
        firstTimes = nan(1,nexp_attn);
        for iexp = 1:nexp_attn
            firstTimes(iexp) = targetTimes(find(y(x>0,iexp)>pctCorrThresh,1));
        end
        errorbar(mean(firstTimes(ind)),0.9,[],[],...
            ste(firstTimes(ind),2),ste(firstTimes(ind),2),'.')
        figXAxis([],'Resp. Win. Center (Relative to Stim)',[x(1) x(end)],x(1:2:length(x)),round(x(1:2:length(x))))
        hline(0.5,'k-')
        figYAxis([],'Frac. Correct',[0 1])
        vline(0,'k:')
        figAxForm
        title(sprintf('Attn %s %s Model',avName{iav},modelName{imod}))
        vline(respWinTT,'r--')
        vline(baseWinTT,'k--')
        hline(pctCorrThresh,'k:')
        if iav == 1
            vline(visRTwindow(1),'b--')
        else
            vline(audRTwindow(1),'b--')
        end
        
        if imod == 1
            iplot = 4;
        else
            iplot = 6;
        end
        subplot(4,2,iplot+iav)
        hold on
        y = mdl_noAttn(iav).mdl(imod).pctCorr_respWin(winInd,:);
        ind = mdl_noAttn(iav).mdl(imod).pctCorr_all > pctCorrThresh;
        h = plot(x,y(:,ind),'k-');
        errorbar(x,mean(y(:,ind),2),ste(y(:,ind),2),'.-')
        firstTimes = nan(1,nexp_noAttn);
        for iexp = 1:nexp_noAttn
            firstTimes(iexp) = targetTimes(find(y(x>0,iexp)>pctCorrThresh,1));
        end
        errorbar(mean(firstTimes(ind)),0.9,[],[],...
            ste(firstTimes(ind),2),ste(firstTimes(ind),2),'.')
        figXAxis([],'Resp. Win. Center (Relative to Stim)',[x(1) x(end)],x(1:2:length(x)),round(x(1:2:length(x))))
        hline(0.5,'k-')
        figYAxis([],'Frac. Correct',[0 1])
        vline(0,'k:')
        figAxForm
        title(sprintf('No Attn %s %s Model',avName{iav},modelName{imod}))
        vline(respWinTT,'r--')
        vline(baseWinTT,'k--')
        hline(pctCorrThresh,'k:')
        if iav == 1
            vline(visRTwindow(1),'b--')
        else
            vline(audRTwindow(1),'b--')
        end
    end
end

print([fnout 'pctCorrectXrespWin'],'-dpdf','-fillpage')   

%% stim and choice weight correlation

figure
for iav = 1:2
    ind = mdl_attn(iav).mdl(stimModInd).pctCorr_all > pctCorrThresh;
    subplot(2,4,iav)
    x = cell2mat(mdl_attn(iav).mdl(stimModInd).normWeights(ind)');
    y = cell2mat(mdl_attn(iav).mdl(choiceModInd).normWeights(ind)');
    hold on
    plot(x,y,'.','MarkerSize',10)
    
    mdl_fitline = fitlm(x,y);
    rsq = mdl_fitline.Rsquared.Ordinary;
    yfit = predict(mdl_fitline,x);
    plot(x,yfit,'-')
    [r,p] = corr(x,y);
    title(sprintf('Attn, %s Models,r=%s,p=%s',avName{iav},...
        num2str(round(r,2,'significant')),num2str(round(p,2,'significant'))))    
    figXAxis([],'Stimulus Weight',[-1 1])
    figYAxis([],'Choice Weight',[-1 1])
    figAxForm
    
    ind = mdl_noAttn(iav).mdl(stimModInd).pctCorr_all > pctCorrThresh;
    subplot(2,4,iav+2)
    x = cell2mat(mdl_noAttn(iav).mdl(stimModInd).normWeights(ind)');
    y = cell2mat(mdl_noAttn(iav).mdl(choiceModInd).normWeights(ind)');
    hold on
    plot(x,y,'.','MarkerSize',10)
    
    mdl_fitline = fitlm(x,y);
    rsq = mdl_fitline.Rsquared.Ordinary;
    yfit = predict(mdl_fitline,x);
    plot(x,yfit,'-')
    [r,p] = corr(x,y);
    title(sprintf('No Attn, %s Models,r=%s,p=%s',avName{iav},...
        num2str(round(r,2,'significant')),num2str(round(p,2,'significant'))))    
    figXAxis([],'Stimulus Weight',[-1 1])
    figYAxis([],'Choice Weight',[-1 1])
    figAxForm
end
subplot 245
ind = mdl_attn(visualTrials).mdl(stimModInd).pctCorr_all > pctCorrThresh & ...
    mdl_attn(auditoryTrials).mdl(stimModInd).pctCorr_all > pctCorrThresh;
x = cellfun(@corr,mdl_attn(visualTrials).mdl(stimModInd).normWeights(ind),...
    mdl_attn(visualTrials).mdl(choiceModInd).normWeights(ind));
y = cellfun(@corr,mdl_attn(auditoryTrials).mdl(stimModInd).normWeights(ind),...
    mdl_attn(auditoryTrials).mdl(choiceModInd).normWeights(ind));
plot(x,y,'.','MarkerSize',20)
hold on
plot(0:1,0:1,'k--')
figXAxis([],'Vis Stim W:Choice W Corr',[0 1])
figYAxis([],'Aud Stim W:Choice W Corr',[0 1])
figAxForm
title('Attn Mice')

subplot 246
plot(1:2,[x(ind)',y(ind)'],'k-','MarkerSize',20)
hold on
errorbar(1,mean(x(ind)),ste(x(ind),2),'.','MarkerSize',20)
% plot(ones(1,sum(ind)).*2,y(ind),'.','MarkerSize',20)
errorbar(2,mean(y(ind)),ste(y(ind),2),'.','MarkerSize',20)
figXAxis([],'',[0 3],1:2,avName)
figYAxis([],'Stim:Choice Weight Corr',[0 1])
figAxForm([],0)

subplot 247
ind = mdl_noAttn(visualTrials).mdl(stimModInd).pctCorr_all > pctCorrThresh & ...
    mdl_noAttn(auditoryTrials).mdl(stimModInd).pctCorr_all > pctCorrThresh;
x = cellfun(@corr,mdl_noAttn(visualTrials).mdl(stimModInd).normWeights,...
    mdl_noAttn(visualTrials).mdl(choiceModInd).weights);
y = cellfun(@corr,mdl_noAttn(auditoryTrials).mdl(stimModInd).normWeights,...
    mdl_noAttn(auditoryTrials).mdl(choiceModInd).weights);
plot(x,y,'.','MarkerSize',20)
hold on
plot(0:1,0:1,'k--')
figXAxis([],'Vis Stim W:Choice W Corr',[0 1])
figYAxis([],'Aud Stim W:Choice W Corr',[0 1])
figAxForm
title('No Attn Mice')

subplot 248
plot(1:2,[x(ind)',y(ind)'],'k-','MarkerSize',20)
hold on
errorbar(1,mean(x(ind)),ste(x(ind),2),'.','MarkerSize',20)
plot(ones(1,sum(ind)).*2,y(ind),'.','MarkerSize',20)
errorbar(2,mean(y(ind)),ste(y(ind),2),'.','MarkerSize',20)
figXAxis([],'',[0 3],1:2,avName)
figYAxis([],'Stim:Choice Weight Corr',[0 1])
figAxForm([],0)

print([fnout 'stimVsChoiceWeights'],'-dpdf','-fillpage')

%% combo model performance
mdl_attn(visualTrials).mdl(stimModInd).comboWeight = cell(1,nexp_attn);
mdl_attn(visualTrials).mdl(choiceModInd).comboWeight = cell(1,nexp_attn);
mdl_noAttn(visualTrials).mdl(stimModInd).comboWeight = cell(1,nexp_noAttn);
mdl_noAttn(visualTrials).mdl(choiceModInd).comboWeight = cell(1,nexp_noAttn);
for iav = 1:2
    mdl_attn(iav).mdl(stimModInd).pctCorr_combo = nan(1,nexp_attn);
    mdl_attn(iav).mdl(choiceModInd).pctCorr_combo = nan(1,nexp_attn);
    for iexp = 1:nexp_attn
        mdl_attn(iav).mdl(stimModInd).pctCorr_combo(iexp) = ...
            dc_attn(iexp).av(iav).pctCorrectTarget_comboTrain;
        mdl_attn(iav).mdl(choiceModInd).pctCorr_combo(iexp) = ...
            dc_attn(iexp).av(iav).pctCorrectDetect_comboTrain;
        if iav==1
            mdl_attn(visualTrials).mdl(stimModInd).comboWeight{iexp} = ...
                dc_attn(iexp).comboTrainWeightTarget;
            mdl_attn(visualTrials).mdl(choiceModInd).comboWeight{iexp} = ...
                dc_attn(iexp).comboTrainWeightDetect;
        end
    end
    mdl_noAttn(iav).mdl(stimModInd).pctCorr_combo = nan(1,nexp_noAttn);
    mdl_noAttn(iav).mdl(choiceModInd).pctCorr_combo = nan(1,nexp_noAttn);
    for iexp = 1:nexp_noAttn
        mdl_noAttn(iav).mdl(stimModInd).pctCorr_combo(iexp) = ...
            dc_noAttn(iexp).av(iav).pctCorrectTarget_comboTrain;
        mdl_noAttn(iav).mdl(choiceModInd).pctCorr_combo(iexp) = ...
            dc_noAttn(iexp).av(iav).pctCorrectDetect_comboTrain;
        if iav==1
            mdl_noAttn(visualTrials).mdl(stimModInd).comboWeight{iexp} = ...
                dc_noAttn(iexp).comboTrainWeightTarget;
            mdl_noAttn(visualTrials).mdl(choiceModInd).comboWeight{iexp} = ...
                dc_noAttn(iexp).comboTrainWeightDetect;
        end
    end
end

for imod = 1:2
    mdl_attn(visualTrials).mdl(imod).normComboWeight = cellfun(@(x) ...
        x./max(abs(x)),mdl_attn(visualTrials).mdl(imod).comboWeight,'unif',0);
    mdl_noAttn(visualTrials).mdl(imod).normComboWeight = cellfun(@(x) ...
        x./max(abs(x)),mdl_noAttn(visualTrials).mdl(imod).comboWeight,'unif',0);
end

figure
for iav = 1:2
    for imod = 1:2
        if imod == 1
            iplot = 0;
        else
            iplot = 2;
        end
        subplot(4,2,iav+iplot)
        ind = mdl_attn(iav).mdl(imod).pctCorr_all > pctCorrThresh;
        y = cat(2,mdl_attn(iav).mdl(imod).pctCorr_all(ind)',...
            mdl_attn(iav).mdl(imod).pctCorr_combo(ind)');
        [~,p] = ttest(y(:,1),y(:,2));
        plot(1:2,y,'k-')
        hold on
        errorbar(1:2,mean(y,1),ste(y,1),'.')
        figXAxis([],'Train',[0 3],1:2,{avName{iav},'Vis+Aud'})
        figYAxis([],'Frac. Correct',[0 1])
        figAxForm
        hline(0.5,'k-')
        hline(pctCorrThresh,'k:')
        title(sprintf('Attn; Test %s, %s Model,p=%s',avName{iav},modelName{imod},...
            num2str(round(p,2,'significant'))))
        
        if imod == 1
            iplot = 4;
        else
            iplot = 6;
        end
        subplot(4,2,iav+iplot)
        ind = mdl_noAttn(iav).mdl(imod).pctCorr_all > pctCorrThresh;
        y = cat(2,mdl_noAttn(iav).mdl(imod).pctCorr_all(ind)',...
            mdl_noAttn(iav).mdl(imod).pctCorr_combo(ind)');
        [~,p] = ttest(y(:,1),y(:,2));
        plot(1:2,y,'k-')
        hold on
        errorbar(1:2,mean(y,1),ste(y,1),'.')
        figXAxis([],'Train',[0 3],1:2,{avName{iav},'Vis+Aud'})
        figYAxis([],'Frac. Correct',[0 1])
        figAxForm
        hline(0.5,'k-')
        hline(pctCorrThresh,'k:')
        title(sprintf('No Attn; Test %s, %s Model,p=%s',avName{iav},modelName{imod},...
            num2str(round(p,2,'significant'))))
    end
end
print([fnout 'comboModelPerf'],'-dpdf','-fillpage')


figure
for imod = 1:2
    subplot(2,2,imod)
    ind = mdl_attn(visualTrials).mdl(imod).pctCorr_all > pctCorrThresh & ...
        mdl_attn(auditoryTrials).mdl(imod).pctCorr_all > pctCorrThresh;
    maxWeight = cellfun(@(x,y,z) max(abs(cat(1,x,y,z))),...
        mdl_attn(visualTrials).mdl(imod).weights,...
        mdl_attn(auditoryTrials).mdl(imod).weights,...
        mdl_attn(visualTrials).mdl(imod).comboWeight,'unif',0);
    y = [cell2mat(cellfun(@(x,m) abs(x./m),...
        mdl_attn(visualTrials).mdl(imod).weights(ind),maxWeight(ind),'unif',0)'),...
        cell2mat(cellfun(@(x,m) abs(x./m),...
        mdl_attn(auditoryTrials).mdl(imod).weights(ind),maxWeight(ind),'unif',0)'),...
        cell2mat(cellfun(@(x,m) abs(x./m),...
        mdl_attn(visualTrials).mdl(imod).comboWeight(ind),maxWeight(ind),'unif',0)')];
    [p,~,stats] = anova1(y,[],'off');
%     y = abs([cell2mat(mdl_attn(visualTrials).mdl(imod).normWeights(ind)'),...
%         cell2mat(mdl_attn(auditoryTrials).mdl(imod).normWeights(ind)'),...
%         cell2mat(mdl_attn(visualTrials).mdl(imod).normComboWeight(ind)')]);
    tbl = multcompare(stats,'display','off');
    fprintf('Attn %s Model:\n',modelName{imod})
    disp(tbl(:,[1:2,6]))
    hold on
    errorbar(1:3,mean(y,1),ste(y,1),'.')
    figXAxis([],'',[0 4],1:3,cat(2,avName,{'Vis+Aud'}))
    figYAxis([],'Abs Weight (norm across conditions)',[0 0.5])
    figAxForm
    title(sprintf('Attn %s Model,p=%s',modelName{imod},sigfigString(p)))

    subplot(2,2,imod+2)
    ind = mdl_noAttn(visualTrials).mdl(imod).pctCorr_all > pctCorrThresh & ...
        mdl_noAttn(auditoryTrials).mdl(imod).pctCorr_all > pctCorrThresh;
    maxWeight = cellfun(@(x,y,z) max(abs(cat(1,x,y,z))),...
        mdl_noAttn(visualTrials).mdl(imod).weights,...
        mdl_noAttn(auditoryTrials).mdl(imod).weights,...
        mdl_noAttn(visualTrials).mdl(imod).comboWeight,'unif',0);
    y = [cell2mat(cellfun(@(x,m) abs(x./m),...
        mdl_noAttn(visualTrials).mdl(imod).weights(ind),maxWeight(ind),'unif',0)'),...
        cell2mat(cellfun(@(x,m) abs(x./m),...
        mdl_noAttn(auditoryTrials).mdl(imod).weights(ind),maxWeight(ind),'unif',0)'),...
        cell2mat(cellfun(@(x,m) abs(x./m),...
        mdl_noAttn(visualTrials).mdl(imod).comboWeight(ind),maxWeight(ind),'unif',0)')];
    [p,~,stats] = anova1(y,[],'off');
%     y = abs([cell2mat(mdl_noAttn(visualTrials).mdl(imod).normWeights(ind)'),...
%         cell2mat(mdl_noAttn(auditoryTrials).mdl(imod).normWeights(ind)'),...
%         cell2mat(mdl_noAttn(visualTrials).mdl(imod).normComboWeight(ind)')]);
    tbl = multcompare(stats,'display','off');
    fprintf('No Attn %s Model:\n',modelName{imod})
    disp(tbl(:,[1:2,6]))
    hold on
    errorbar(1:3,mean(y,1),ste(y,1),'.')
    figXAxis([],'',[0 4],1:3,cat(2,avName,{'Vis+Aud'}))
    figYAxis([],'Abs Weight (norm across conditions)',[0 0.5])
    figAxForm
    title(sprintf('No Attn %s Model,p=%s',modelName{imod},sigfigString(p)))
end

print([fnout 'comboModelWeights'],'-dpdf','-fillpage')