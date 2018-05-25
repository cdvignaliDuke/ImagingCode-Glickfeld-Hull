clear all
close all
ds = 'FSAV_V1_SOM';
cellsOrDendrites = 1;
doRawData = 0;
%%
rc = behavConstsAV;
if strcmp(rc.name,'ashle') & isempty(ds)
    dataGroup = ['awFSAVdatasets' ds];
elseif strcmp(rc.name,'ashle') & strcmp(ds(1:3),'FSA')
    dataGroup = ds;
elseif strcmp(rc.name,'ashle')
    dataGroup = ['awFSAVdatasets' ds];
else
    dataGroup = [];
end
eval(dataGroup)
titleStr = ds;
if strcmp(titleStr, '')
    titleStr = 'V1_100ms';
elseif strcmp(titleStr(1),'_')
    titleStr = titleStr(2:end);    
elseif strcmp(titleStr(1:4), 'FSAV')
    titleStr = titleStr(6:end);
end
str = unique({expt.SubNum});
mouse_str = ['i' strjoin(str,'_i')];
if doRawData == 1
    mouse_str = [mouse_str '_rawDataTC'];
    titleStr = [titleStr '_rawDataTC'];
end

if isempty(ds)
    load(fullfile(rc.caOutputDir,dataGroup,[mouse_str '_CaSummary' ds '.mat']));
elseif strcmp(ds(1:3),'FSA')
    if cellsOrDendrites == 1
        load(fullfile(rc.caOutputDir,dataGroup,[mouse_str '_CaSummary_cells' ds(5:end) '.mat']));    
    elseif cellsOrDendrites == 2
        load(fullfile(rc.caOutputDir,dataGroup,[mouse_str '_CaSummary_dendrites' ds(5:end) '.mat']));    
    end        
else
    load(fullfile(rc.caOutputDir,dataGroup,[mouse_str '_CaSummary' ds '.mat']));
end
% load(fullfile(rc.caOutputDir,dataGroup,[titleStr '_' mouse_str '_modCells.mat']));
if cellsOrDendrites == 1
    fnout = fullfile(rc.caOutputDir, dataGroup,[titleStr '_startAlign']); 
elseif cellsOrDendrites == 2
    fnout = fullfile(rc.caOutputDir, dataGroup,[titleStr '_den_startAlign']); 
end
    

%%
outcomecolmat = {'k';'r'};
avcolmat = {'k','c'};
avmat = {'visual';'auditory'};
outcomemat = {'hit';'miss'};
pressAlign = 1;
targetAlign = 2;
visualTrials = 1;
auditoryTrials = 2;
prevTrialVisual = 0;
prevTrialAuditory = 1;
nav = 2;
hitTrials = 1;
missTrials = 2;
nout = 2;
cycLengthFr = mouse(1).expt(1).info.cyc_time;
frRateHz = expt(1).frame_rate;
onTimeFr = 0.1*frRateHz;
nMonitorDelayFr = 2;
nBaselineFr = mouse(1).expt(1).pre_event_frames;
timestamp_1cyc = (-nBaselineFr:cycLengthFr-1);
plotTimeRange_1cyc = [-5 chop(timestamp_1cyc(end),2)];
nexp = 0;
for imouse = 1:size(mouse,2)
    nexp = nexp+size(mouse(imouse).expt,2);
end
exptDates = {expt.date};
exptMice = {expt.SubNum};

basewin = 1:34;
basewin_0 = 32:34;
respwin = 36:38;
basewinTarget = 31:33;
respwinTarget = 35:37;
nCycles = 8;
cellGroupsAlpha = 0.05;
lateCycles = 5:nCycles;
lateWinFrames = cycLengthFr*length(lateCycles);
tuningReliabilityThresh = 45;
minRespThreshold = 0.002;

%%
exptName = cell(1,nexp);
exptSN = nan(1,nexp);
responsiveCellsAV = cell(nexp,1);
lateRespCellsAV = cell(nexp,1);
lateSuppCellsAV = cell(nexp,1);
longTrialTCAV = cell(1,nexp);
longTrialTCAV_sem = cell(1,nexp);
longTrialTCEaExpt = cell(nexp,2);
trialTCEaCycEaCell = cell(nCycles,nexp,nav);
tcEaCellEaCycAllOutcome = cell(nCycles,nexp,nav);
lateWinRespAVEaExpt = cell(nexp,nav);
lateWinSIEaExpt = cell(nexp,1);
lateSIEaExpt = cell(nexp,1);
lateCycWinSIEaExpt = cell(nexp,1);
lateCycRespEaExpt = cell(nexp,nav);
lateCycMeanEaExpt = cell(nexp,nav);
ttestLateAVEaExpt = cell(nexp,1);
ttestLateWinAVEaExpt = cell(nexp,1);
ttestLateCycWinAVEaExpt = cell(nexp,1);
oriTuningFitEaExpt = cell(nexp,1);
oriTuningFitReliability = cell(nexp,1);
for imouse = 1:size(mouse,2)
    for iexp = 1:size(mouse(imouse).expt,2)
        if imouse == 1 && iexp == 1
            exptN = 1;
        else
            exptN = exptN+1;
        end
        exptName{exptN} = [mouse(imouse).expt(iexp).mouse_name '-' ...
            mouse(imouse).expt(iexp).date];
        exptSN(exptN) = str2num(mouse(imouse).expt(iexp).mouse_name);
        dAV_cycle1 = cat(3,...
            mouse(imouse).expt(iexp).align(pressAlign).av(visualTrials).outcome(hitTrials).cmlvCycResp{1},...
            mouse(imouse).expt(iexp).align(pressAlign).av(auditoryTrials).outcome(hitTrials).cmlvCycResp{1});
        responsiveCellsAV{exptN} = ttest(...
            squeeze(mean(dAV_cycle1(respwin,:,:),1)),...
            squeeze(mean(dAV_cycle1(basewin,:,:),1)),...
            'dim',2,'tail','right','alpha',cellGroupsAlpha)' ... 
            & mean(mean(dAV_cycle1(respwin,:,:),3),1) > minRespThreshold;
        
        dAV_nCycles = cat(3,mouse(imouse).expt(iexp).align(pressAlign).av(visualTrials).outcome(hitTrials).cmlvCycResp{nCycles},...
            mouse(imouse).expt(iexp).align(pressAlign).av(auditoryTrials).outcome(hitTrials).cmlvCycResp{nCycles});
        longTrialTCAV{exptN} = mean(dAV_nCycles,3);
        longTrialTCAV_sem{exptN} = ste(dAV_nCycles,3);
        
        dV_nCycles = mouse(imouse).expt(iexp).align(pressAlign).av(...
            visualTrials).outcome(hitTrials).cmlvCycResp{nCycles};
        dA_nCycles = mouse(imouse).expt(iexp).align(pressAlign).av(...
            auditoryTrials).outcome(hitTrials).cmlvCycResp{nCycles};
        lateRespCellsAV{exptN} = logical(ttest(...
            squeeze(mean(dAV_nCycles((basewin(end)+1):end,:,:),1)),...
            squeeze(mean(dAV_nCycles(basewin,:,:),1)),...
            'dim',2,'tail','right','alpha',cellGroupsAlpha)');
        lateSuppCellsAV{exptN} = logical(ttest(...
            squeeze(mean(dAV_nCycles((basewin(end)+1):end,:,:),1)),...
            squeeze(mean(dAV_nCycles(basewin,:,:),1)),...
            'dim',2,'tail','left','alpha',cellGroupsAlpha)');
        
        responseAllTrials = cell(nCycles,nav);
        cycMeanAllTrials = cell(nCycles,nav);
        for iav = 1:nav
            d = mouse(imouse).expt(iexp).align(pressAlign).av(iav);        
            for icyc = 1:nCycles
                tcThisCycAllOutcome = cat(3,d.outcome(hitTrials).cmlvCycResp{icyc});
                extraFrames = cycLengthFr*(icyc-1);
%                 cycWin = (extraFrames:(cycLengthFr*icyc))+30;
                trialTCEaCycEaCell{icyc,exptN,iav} = mean(tcThisCycAllOutcome,3);
                baselineAllTrials = mean(tcThisCycAllOutcome(basewin_0+extraFrames,:,:),1);
                responseAllTrials{icyc,iav} = mean(tcThisCycAllOutcome(...
                    respwin+extraFrames,:,:),1) - baselineAllTrials;
                cycMeanAllTrials{icyc,iav} = mean(tcThisCycAllOutcome(...
                    (respwin(1)+extraFrames):end,:,:),1);
                tcEaCellEaCycAllOutcome{icyc,exptN,iav} = mean(bsxfun(@minus,tcThisCycAllOutcome(...
                    end-(cycLengthFr+nBaselineFr)+1:end,:,:),...
                    baselineAllTrials),3);
            end
        end
        longTrialTCnCyc = cell(1,nav);
        longTrialTCnCyc{visualTrials} = mouse(imouse).expt(iexp).align(...
            pressAlign).av(visualTrials).outcome(hitTrials).cmlvCycResp{nCycles};
        longTrialTCnCyc{auditoryTrials} = mouse(imouse).expt(iexp).align(...
            pressAlign).av(auditoryTrials).outcome(hitTrials).cmlvCycResp{nCycles};
        lateWinRespAV = cellfun(@(x) squeeze(mean(x(end-lateWinFrames:end,:,:),1))',...
            longTrialTCnCyc,'unif',0);
        lateWinRespAVEaExpt(exptN,:) = cellfun(@(x) mean(x,1),lateWinRespAV,'unif',0);
        lateWinSIEaExpt{exptN} = getSelectivityIndex(...
            lateWinRespAV{visualTrials},lateWinRespAV{auditoryTrials});
        ttestLateWinAVEaExpt{exptN} = ttest2(...
            lateWinRespAV{visualTrials},lateWinRespAV{auditoryTrials},...
            'dim',1,'alpha',0.05);
        longTrialTCEaExpt(exptN,:) = longTrialTCnCyc;
        
        responseAllTrials = cellfun(@(x)squeeze(x)',responseAllTrials,'unif',0);
        lateCycRespEaExpt{exptN,visualTrials} = ...
            mean(cell2mat(responseAllTrials(lateCycles,visualTrials)),1);
        lateCycRespEaExpt{exptN,auditoryTrials} = ...
            mean(cell2mat(responseAllTrials(lateCycles,auditoryTrials)),1);
        lateSIEaExpt(exptN) = {getSelectivityIndex(...
            cell2mat(responseAllTrials(lateCycles,visualTrials)),...
            cell2mat(responseAllTrials(lateCycles,auditoryTrials)))};
        lateRespV = cell2mat(responseAllTrials(lateCycles,visualTrials));
        lateRespA = cell2mat(responseAllTrials(lateCycles,auditoryTrials));
        ttestLateAVEaExpt{exptN} = ttest2(...
            lateRespV,lateRespA,'dim',1,'alpha',0.05);

        cycMeanAllTrials = cellfun(@(x)squeeze(x)',cycMeanAllTrials,'unif',0);
        lateCycMeanEaExpt{exptN,visualTrials} = ...
            mean(cell2mat(cycMeanAllTrials(lateCycles,visualTrials)),1);
        lateCycMeanEaExpt{exptN,auditoryTrials} = ...
            mean(cell2mat(cycMeanAllTrials(lateCycles,auditoryTrials)),1);
        lateCycWinSIEaExpt{exptN} = getSelectivityIndex(...
            cell2mat(cycMeanAllTrials(lateCycles,visualTrials)),...
            cell2mat(cycMeanAllTrials(lateCycles,auditoryTrials)));
        
        lateRespV = cell2mat(cycMeanAllTrials(lateCycles,visualTrials));
        lateRespA = cell2mat(cycMeanAllTrials(lateCycles,auditoryTrials));
        ttestLateCycWinAVEaExpt{exptN} = ttest2(...
            lateRespV,lateRespA,'dim',1,'alpha',0.05);
        
        oriTuningFitEaExpt{exptN} = mouse(imouse).expt(iexp).cells.oriTuning.oriFit;
        oriTuningFitReliability{exptN} = mouse(imouse).expt(iexp).cells.oriTuning.oriFitReliability;
        
    end
end
%%
respCells = cell2mat(responsiveCellsAV');
lateStimModCells = cell2mat(lateRespCellsAV') | cell2mat(lateSuppCellsAV');
lateRespCells = cell2mat(lateRespCellsAV');

longTC = cell(1,nav);
lateCycResp = cell(1,nav);
lateCycWin = cell(1,nav);
lateWin = cell(1,nav);
for iav = 1:nav
    longTC{iav} = trialTCEaCycEaCell(nCycles,:,iav);
    lateCycResp{iav} = cell2mat(lateCycRespEaExpt(:,iav)');
    lateCycWin{iav} = cell2mat(lateCycMeanEaExpt(:,iav)');
    lateWin{iav} = cell2mat(lateWinRespAVEaExpt(:,iav)');
end
longTC = cellfun(@cell2mat,longTC,'unif',0);

siCycResp = cell2mat(lateSIEaExpt');
siCycWin = cell2mat(lateCycWinSIEaExpt');
siLateWin = cell2mat(lateWinSIEaExpt');

lateAVModCycResp = cell2mat(ttestLateAVEaExpt');
lateAVModCycWin = cell2mat(ttestLateCycWinAVEaExpt');
lateAVModLateWin = cell2mat(ttestLateWinAVEaExpt');

tcAVallCells = cell2mat(longTrialTCAV);
%%
responseLim = [-0.005 0.2];
scatLim = [-0.15 1];
scatRespLim = [-0.03 0.075];
ttLabel = 0:500:2500;
nFr = size(longTC{1,1},1);
tt = ((26:nFr)-33).*(1000/frRateHz);
ttFr = (26:nFr)-33;
ttLabelFr = ((ttLabel./1000)*frRateHz)+8;
lateWinTT = ([nFr - lateWinFrames nFr] - 33)...
    .*(1000/frRateHz);
binEdgesSI = -10:0.5:10;
siLim = [-10 10];
hmLim = [-0.1 0.1];
exCellRespLim = [-0.1 0.25];
cycTCLim = [-0.005 0.015];
cycRespLim = [0 0.02];
%%
setFigParams4Print('portrait')
figure
subplot 211
for iav = 1:nav
    y = mean(longTC{iav}(:,respCells),2);
    yerr = ste(longTC{iav}(:,respCells),2);    
    tt = ((26:length(y))-33).*(1000/frRateHz);
    h = shadedErrorBar(tt,y(26:length(y)),yerr(26:length(y)),avcolmat{iav});
    for i = 1:2
        h.edge(i).Color = 'none';
    end
    hold on
end
figXAxis([],'Time (ms)',[tt(1) tt(end)],ttLabel,ttLabel)
figYAxis([],'dF/F',responseLim)  
vline(lateWinTT,'k--')
figAxForm([],0)
title(sprintf('First Stim Responsive Cells (%s/%s)',...
    num2str(sum(respCells)),num2str(length(respCells))))

subplot 212
x = lateWin{visualTrials}(respCells);
y = lateWin{auditoryTrials}(respCells);
h = scatter(x,y,'ko');
hold on
plot(scatLim,scatLim,'k--');
[~,p] = ttest(x,y);
figXAxis([],'Visual (dF/F)',scatLim)
figYAxis([],'Auditory (dF/F)',scatLim)
figAxForm([])
title(sprintf('Responsive Cells (%s), p = %s',num2str(sum(respCells)),...
    num2str(round(p,2,'significant'))))

print([fnout '_respCellsTC&Scatter'],'-dpdf','-fillpage')

setFigParams4Print('landscape')
figure
subplot 121
si = siCycResp;
ind = lateAVModCycResp;
h = histogram(si(respCells),binEdgesSI);
h.FaceColor = [0.5 0.5 0.5];
h.EdgeColor = 'none';
h = vline(0,'k-');
h.LineWidth = 2;
h = vline(mean(si(respCells)),'r-');
h.LineWidth = 2;
hold on
h = histogram(si(ind & respCells),binEdgesSI);
h.FaceColor = 'k';
h.EdgeColor = 'none';
h = vline(mean(si(respCells & ind)),'m-');
h.LineWidth = 2;
figXAxis([],'Late Trial Selectivity Index (using resp to cycles 5-8)',siLim)
figYAxis([],'n ncells',[])
figAxForm([])
title({sprintf('Significantly Attn Modulated (%s/%s)'...
    ,num2str(sum(ind & respCells)),num2str(sum(respCells))),...
    'SI calc using cyc reposnses'})

subplot 122
x = lateCycResp{visualTrials}(respCells);
y = lateCycResp{auditoryTrials}(respCells);
h = scatter(x,y,'ko');
[~,p] = ttest(x,y);
hold on
x = lateCycResp{visualTrials}(respCells & ind & si > 0);
y = lateCycResp{auditoryTrials}(respCells & ind & si > 0);
h = scatter(x,y,'ro');
x = lateCycResp{visualTrials}(respCells & ind & si < 0);
y = lateCycResp{auditoryTrials}(respCells & ind & si < 0);
h = scatter(x,y,'bo');
plot(scatRespLim,scatRespLim,'k--');
figXAxis([],'Visual (dF/F)',scatRespLim)
figYAxis([],'Auditory (dF/F)',scatRespLim)
figAxForm([])
title(sprintf('Responsive Cells (%s), p = %s',num2str(sum(respCells)),...
    num2str(round(p,2,'significant'))))

print([fnout '_respCellsAttnSI&Scat'],'-dpdf','-fillpage')
%% example cells
if strcmp(ds,'FSAV_V1_100ms')
    tcAVAllCellsSem = cell2mat(longTrialTCAV_sem);
    exCellTC = tcAVallCells(26:end,:);
    exCellTCsem = tcAVAllCellsSem(26:end,:);

    exRespCell = 738; 
    exLateRespCell = 319; 
    exLateSuppCell = 386;
    exampleCells = [exRespCell exLateRespCell exLateSuppCell];

    setFigParams4Print('portrait')
    figure
    subplot 311
    shadedErrorBar(tt,exCellTC(:,exRespCell),exCellTCsem(:,exRespCell),'k');
    hold on
    hline(0,'k--')
    title(sprintf('respond to 1st stim, #%s',num2str(exRespCell)))
    subplot 312
    shadedErrorBar(tt,exCellTC(:,exLateRespCell),exCellTCsem(:,exLateRespCell),'k');
    hold on
    hline(0,'k--')
    title(sprintf('respond late, #%s',num2str(exLateRespCell)))
    subplot 313
    shadedErrorBar(tt,exCellTC(:,exLateSuppCell),exCellTCsem(:,exLateSuppCell),'k');
    hold on
    hline(0,'k--')
    title(sprintf('suppressed late, #%s',num2str(exLateSuppCell)))
    for i = 1:3
        subplot(3,1,i)
        figXAxis([],'time (ms)',[tt(1) tt(end)],ttLabelFr,ttLabel)
        figYAxis([],'dF/F',exCellRespLim)
        figAxForm([],0)
        hold on
        vline(0,'k--')
        vline(0,'k--')
    end

    print([fnout '_exCellTC'],'-dpdf')
else
    exampleCells = 0;
end
%% heatmaps
avTCrespCells = tcAVallCells(:,respCells);
TCrespCells = cell(1,nav);
for iav = 1:nav
    TCrespCells{iav} = longTC{iav}(:,respCells);
end

[~,sortInd] = sort(mean(avTCrespCells((end-lateWinFrames):end,:),1));
hm = cellfun(@(x) flipud(x(:,sortInd)'),TCrespCells,'unif',0);

setFigParams4Print('landscape')
figure
colormap(brewermap([],'*RdBu'));
suptitle(sprintf('Responsive Cells (%s)',num2str(sum(respCells))));

subplot 221
h = imagesc(hm{visualTrials}(:,26:end));
figXAxis([],'Time (ms)',[],ttLabelFr,ttLabel)
figYAxis([],'Cell #',[])
figAxForm([])
colorbar
caxis(hmLim)
title('Visual Trials');

subplot 222
h = imagesc(hm{auditoryTrials}(:,26:end));
figXAxis([],'Time (ms)',[],ttLabelFr,ttLabel)
figYAxis([],'Cell #',[])
figAxForm([])
colorbar
caxis(hmLim)
title('Auditory Trials');

allCells = respCells | lateStimModCells;
avTCallCells = tcAVallCells;
[~,sortInd] = sort(mean(avTCallCells((end-lateWinFrames):end,:),1));
hm = flipud(avTCallCells(:,sortInd)');
exCellInd = find(fliplr(ismember(sortInd,exampleCells)));

subplot 223
h = imagesc(hm(:,26:end));
hold on
if ~isempty(exCellInd)
    hline(exCellInd,'k-');
end
figXAxis([],'Time (ms)',[],ttLabelFr,ttLabel)
figYAxis([],'Cell #',[])
figAxForm([])
colorbar
caxis(hmLim)
title('All Trials, All Stim Modulated Cells');

print([fnout '_antiTCheatmap'],'-dpdf','-fillpage')

%% adaptation
cycTC = cell(nCycles,nav);
lateCycResp = cell(1,nav);
for iav = 1:nav
    for icyc = 1:nCycles
        cycTC{icyc,iav} = cell2mat(tcEaCellEaCycAllOutcome(icyc,:,iav));
    end
    lateCycResp{iav} = cell2mat(lateCycRespEaExpt(:,iav)');
end
cycResp = cellfun(@(x) mean(x(respwin,:),1),cycTC,'unif',0);

minRespCells = respCells & cycResp{1,visualTrials} > minRespThreshold & ...
    cycResp{1,auditoryTrials} > minRespThreshold;

earlyLateCycResp = cell(2,nav);
earlyLateCycResp(1,:) = cycResp(1,:);
earlyLateCycResp(2,:) = lateCycResp;

normCycResp = cell(nCycles,nav);
normEarlyLateCycResp = cell(2,nav);
for iav = 1:nav
    normCycResp(:,iav) = cellfun(@(x) x./cycResp{1,iav},cycResp(:,iav),'unif',0);
    normEarlyLateCycResp(:,iav) =  ...
        cellfun(@(x) x./earlyLateCycResp{1,iav},earlyLateCycResp(:,iav),'unif',0);
end
respCellsNormCycResp = cellfun(@(x) mean(x(minRespCells)),normCycResp);

%%
ttCyc = ((26:size(cycTC{1,1},1))-33).*(1000/frRateHz);
ttCycResp = ([respwin(1) respwin(end)]-33).*(1000/frRateHz);
setFigParams4Print('landscape')
figure
for icyc = 1:nCycles
    subplot(2,4,icyc)
    for iav = 1:nav
        y = mean(cycTC{icyc,iav}(26:end,minRespCells),2);
        yerr = ste(cycTC{icyc,iav}(26:end,minRespCells),2);
        h = shadedErrorBar(ttCyc,y,yerr,avcolmat{iav});
        for i = 1:2
            h.edge(i).Color = 'none';
        end
        hold on
    end
    vline(ttCycResp,'k--')
    figXAxis([],'Time (ms)',[ttCyc(1) ttCyc(end)])
    figYAxis([],'dF/F',responseLim)
    figAxForm
    title(['Stimulus #' num2str(icyc)])
end
print([fnout '_cycTC'],'-dpdf','-fillpage')

setFigParams4Print('landscape')
figure
for icyc = 1:nCycles
    subplot(2,4,icyc)
    for iav = 1:nav
        y = mean(cycTC{icyc,iav}(26:end,lateRespCells),2);
        yerr = ste(cycTC{icyc,iav}(26:end,lateRespCells),2);
        h = shadedErrorBar(ttCyc,y,yerr,avcolmat{iav});
        for i = 1:2
            h.edge(i).Color = 'none';
        end
        hold on
    end
    vline(ttCycResp,'k--')
    figXAxis([],'Time (ms)',[ttCyc(1) ttCyc(end)])
    figYAxis([],'dF/F',responseLim)
    figAxForm
    title(['Stimulus #' num2str(icyc)])
end
print([fnout '_cycTC_lateRespCells'],'-dpdf','-fillpage')

figure
suptitle(sprintf('Responsive Cells (resp must be greater than threshold (%s))',...
    num2str(sum(minRespCells))))
subplot 131
cr = cellfun(@(x) mean(x(minRespCells)),cycResp);
cr_sem = cellfun(@(x) ste(x(minRespCells),2),cycResp);
for iav = 1:nav
    h = errorbar(1:nCycles,cr(:,iav),cr_sem(:,iav),'o');
    h.Color = avcolmat{iav};
    h.MarkerFaceColor = [1 1 1];
    hold on
end
figXAxis([],'Stimulus # (ms)',[0 nCycles+1])
figYAxis([],'dF/F',cycRespLim)
figAxForm
subplot 132
cr = cellfun(@(x) mean(x(minRespCells)),normCycResp);
cr_sem = cellfun(@(x) ste(x(minRespCells),2),normCycResp);
for iav = 1:nav
    h = errorbar(1:nCycles,cr(:,iav),cr_sem(:,iav),'o');
    h.Color = avcolmat{iav};
    h.MarkerFaceColor = [1 1 1];
    hold on
end
figXAxis([],'Stimulus # (ms)',[0 nCycles+1])
figYAxis([],'Normalized dF/F',[0 1.1])
figAxForm

subplot 133
y = cellfun(@(x) mean(x(minRespCells)), normEarlyLateCycResp(2,:));
yerr = cellfun(@(x) ste(x(minRespCells),2), normEarlyLateCycResp(2,:));
for iav = 1:nav
    h = bar(iav,y(iav));
    h.BarWidth = 1;
    h.EdgeColor = 'none';
    h.FaceColor = avcolmat{iav};
    hold on
    errorbar(iav,y(iav),yerr(iav),avcolmat{iav});
end
figXAxis([],'',[0 3],1:2,{'Visual','Auditory'})
figYAxis([],'Normalized dF/F',[0 1.1])
figAxForm([],0)

print([fnout '_cycResp'],'-dpdf','-fillpage')

%% tuning
oris = 0:22.5:180;
oriLabel = (oris)+1;
oriFit = cell2mat(oriTuningFitEaExpt');
[oriFitMaxVal, oriPeak] = max(oriFit,[],1);
oriPeak = oriPeak-1;
oriBins = [0 22.5:45:157.5 180];
[~,~,oriPref] = histcounts(oriPeak,oriBins);
oriPref(oriPref == 5) = 1;
nori = 4;

isTuned = cell2mat(oriTuningFitReliability') < tuningReliabilityThresh;

oriFitLim = [-0.5 0.5];

avTCtunedCells = tcAVallCells(:,respCells & isTuned);
oriFits_tunedResponsive = oriFit(:,respCells & isTuned);

%%
oriFitNormLim = [-1 1];
figure
colormap(brewermap([],'*RdBu'));
hm = [];
n = zeros(1,nori);
for i = 1:nori
    ind = oriPref == i & isTuned & minRespCells;
    [~,sortInd] = sort(oriPeak(ind));
    fits = oriFit(:,ind);
    normFits = fits./(max(fits,[],1));
    hm = cat(1,hm,normFits(:,sortInd)');
    if i > 1
        n(i) = sum(ind)+n(i-1);
    else
        n(i) = sum(ind);
    end
end
h = imagesc(hm);
hold on
hline(n,'k-');
figXAxis([],'Orientation (deg)',[1 181],oriLabel,oris)
figYAxis([],'Cell #',[])
figAxForm([],0)
colorbar
caxis(oriFitNormLim)
title('Responsive tuned neurons sorted by pref')
print([fnout '_oriTuningFits_sorted'],'-dpdf')
%%
setFigParams4Print('landscape')
figure
suptitle(sprintf('Responsive & Tuned Cells (%s)',num2str(sum(respCells & isTuned))));

[~,sortInd] = sort(oriPeak(respCells & isTuned));
hm = oriFits_tunedResponsive(:,sortInd)';
subplot 221
colormap(brewermap([],'*RdBu'));
h = imagesc(hm);
hold on
figXAxis([],'Orientation (deg)',[1 181],oriLabel,oris)
figYAxis([],'Cell #',[])
figAxForm([])
colorbar
caxis(oriFitLim)
title('Sorted By Ori Pref')

[~,sortInd] = sort(mean(avTCtunedCells((end-lateWinFrames):end,:),1));
hm = oriFits_tunedResponsive(:,fliplr(sortInd))';
subplot 222
h = imagesc(hm);
hold on
figXAxis([],'Orientation (deg)',[1 181],oriLabel,oris)
figYAxis([],'Cell #',[])
figAxForm([])
colorbar
caxis(oriFitLim)
title('Sorted By Task Stim Resp')

subplot 223
h = polarplot(deg2rad(oriPeak(respCells & isTuned)),siCycResp(respCells & isTuned),'k.');
h.Color = [0.5 0.5 0.5];
h.MarkerSize = 12;
hold on
h = polarplot(deg2rad(oriPeak(respCells & isTuned & lateAVModCycResp)),...
    siCycResp(respCells & isTuned & lateAVModCycResp),'k.');
h.MarkerSize = 12;
title('Attn Selectivity by Ori Tuning')

subplot 224
visX = 1:3:10;
audX = visX+1;
orientations = [0 45 90 135];
leg = [];
n = cell(1,nori);
for iori = 1:nori
    ind = oriPref == iori;
    y = cellfun(@(x) mean(x(minRespCells & isTuned & ind)), normEarlyLateCycResp(2,:));
    yerr = cellfun(@(x) ste(x(minRespCells & isTuned & ind),2), normEarlyLateCycResp(2,:));
    for iav = 1:nav
        if iav == 1
            x = visX(iori);
        else
            x = audX(iori);
        end
        h = bar(x,y(iav));
        h.BarWidth = 1;
        h.EdgeColor = 'none';
        h.FaceColor = avcolmat{iav};
        hold on
        errorbar(x,y(iav),yerr(iav),avcolmat{iav});
    end
    leg(iori) = h;
    n{iori} = sprintf('%s: %s',num2str(orientations(iori)),...
        num2str(sum(minRespCells & isTuned & ind)));
end
figXAxis([],'Preferred Orientation (deg)',[0 12],visX+0.5,orientations)
figYAxis([],'Normalized dF/F',[0 1.1])
figAxForm([],0)
L = legend(leg,n,'location','northwest');
title('Responsive, tuned Cells with minimum response cutoff')
title(L,'N Cells Per Group')

print([fnout '_oriTuningFits'],'-dpdf','-fillpage')