clear all
close all
ds = 'FSAV_V1_SOM';
cellsOrDendrites = 1;
doRawData = 0;
%%
rc = behavConstsAV;
imgParams_FSAV
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
frRateHz = expt(1).frame_rate;
onTimeFr = 0.1*frRateHz;
nMonitorDelayFr = 2;
nBaselineFr = mouse(1).expt(1).pre_event_frames;
nexp = 0;
for imouse = 1:size(mouse,2)
    nexp = nexp+size(mouse(imouse).expt,2);
end
exptDates = {expt.date};
exptMice = {expt.SubNum};

basewin = 1:34;
basewin_0 = 32:34;
respwin = 36:38;
% basewinTarget = 31:33;
% respwinTarget = 35:37;
nCycles = 8;
cellGroupsAlpha = 0.01;
lateCycles = 5:nCycles;
tuningReliabilityThresh = 45;
minRespThreshold = 0.002;

%%
doExptCellPlots = false;

attnInfoExpt = struct;

exptName = cell(1,nexp);
exptSN = nan(1,nexp);
exptCycLengthFr = nan(1,nexp);
responsiveCellsAV = cell(nexp,1);
lateRespCellsAV = cell(nexp,1);
lateSuppCellsAV = cell(nexp,1);
longTrialTCAV = cell(1,nexp);
longTrialTCAV_sem = cell(1,nexp);
longTrialTCEaExpt = cell(nexp,2);
trialTCEaCycEaCell = cell(nCycles,nexp,nav);
tcEaCellEaCycAllOutcome = cell(nCycles,nexp,nav);
cycTCAllLateTrials = cell(nexp,nav);
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
        cycLengthFr = mouse(imouse).expt(iexp).info.cyc_time;
        lateWinFrames = cycLengthFr*length(lateCycles);
        exptCycLengthFr(exptN) = cycLengthFr;

%         if cycLengthFr > 11
%             respwin_cellTest = respwin(1):(respwin(1)+5);
%         else
%             respwin_cellTest = respwin;
%         end
        respwin_cellTest = respwin;
        exptName{exptN} = [mouse(imouse).expt(iexp).mouse_name '-' ...
            mouse(imouse).expt(iexp).date];
        exptSN(exptN) = str2num(mouse(imouse).expt(iexp).mouse_name);
        
        dAV_cycle1 = cat(3,...
            mouse(imouse).expt(iexp).align(pressAlign).av(visualTrials).outcome(hitTrials).cmlvCycResp{1},...
            mouse(imouse).expt(iexp).align(pressAlign).av(auditoryTrials).outcome(hitTrials).cmlvCycResp{1});
        responsiveCellsAV{exptN} = ttest(...
            squeeze(mean(dAV_cycle1(respwin_cellTest,:,:),1)),...
            squeeze(mean(dAV_cycle1(basewin_0,:,:),1)),...
            'dim',2,'tail','right','alpha',cellGroupsAlpha)' ... 
            & mean(mean(dAV_cycle1(respwin_cellTest,:,:),3),1) > minRespThreshold;
        
        if doExptCellPlots
            dAV_cycle = mean(cat(3,...
                mouse(imouse).expt(iexp).align(pressAlign).av(visualTrials).outcome(hitTrials).cmlvCycResp{4},...
                mouse(imouse).expt(iexp).align(pressAlign).av(auditoryTrials).outcome(hitTrials).cmlvCycResp{4}),3);
            dAV_cycle_err = ste(cat(3,...
                mouse(imouse).expt(iexp).align(pressAlign).av(visualTrials).outcome(hitTrials).cmlvCycResp{4},...
                mouse(imouse).expt(iexp).align(pressAlign).av(auditoryTrials).outcome(hitTrials).cmlvCycResp{4}),3);
            if exist(fullfile(rc.ashleyAnalysis,['AW' mouse(imouse).expt(iexp).mouse_name(end-1:end)]),'dir') == 7
                fnMs = fullfile(rc.ashleyAnalysis,['AW' mouse(imouse).expt(iexp).mouse_name(end-1:end)],...
                    'two-photon imaging',mouse(imouse).expt(iexp).date);
            else
                fnMs = fullfile(rc.ashleyAnalysis,mouse(imouse).expt(iexp).mouse_name,...
                    'two-photon imaging',mouse(imouse).expt(iexp).date);
            end
            plotEachCellTC(exptName{exptN},dAV_cycle2,dAV_cycle2_err,nBaselineFr,cycLengthFr,frRateHz,responsiveCellsAV{exptN},fnMs)
        end
        
        
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
            cycTC = [];
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
                if any(ismember(lateCycles,icyc))
                cycTC = cat(3,cycTC,bsxfun(@minus,tcThisCycAllOutcome(...
                    end-(cycLengthFr+nBaselineFr)+1:end,:,:),...
                    baselineAllTrials));
                end
            end
            cycTCAllLateTrials{exptN,iav} = mean(cycTC,3);
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
        
        attnInfoExpt(exptN).ms = mouse(imouse).expt(iexp).mouse_name;
        attnInfoExpt(exptN).dt = mouse(imouse).expt(iexp).date;
        attnInfoExpt(exptN).si = lateCycWinSIEaExpt{exptN};
        attnInfoExpt(exptN).avModTest = logical(ttestLateCycWinAVEaExpt{exptN});
        
        r2 =  lateCycRespEaExpt{exptN,visualTrials};
        r1 = mean(tcEaCellEaCycAllOutcome{1,exptN,visualTrials}(respwin,:),1);
        vAdapt = r2./r1;
%         vAdapt(r2 < 0) = 0;
        vAdapt(r1 < minRespThreshold) = nan;
        
        r2 = lateCycRespEaExpt{exptN,auditoryTrials};
        r1 = mean(tcEaCellEaCycAllOutcome{1,exptN,auditoryTrials}(respwin,:),1);
        aAdapt = r2./r1;
%         aAdapt(r2 < 0) = 0;
        aAdapt(r1 < minRespThreshold) = nan;
        
        attnInfoExpt(exptN).visAdapt = vAdapt;
        attnInfoExpt(exptN).audAdapt = aAdapt;
        attnInfoExpt(exptN).firstBaseRespCells = logical(responsiveCellsAV{exptN});
        attnInfoExpt(exptN).lateBaseRespCells = logical(lateRespCellsAV{exptN});
        attnInfoExpt(exptN).lateBaseSuppCells = logical(lateSuppCellsAV{exptN});
        
    end
end

% save struct with cell attention info
save([fnout,'_FSAV_attnData'],'attnInfoExpt')
%%
shortCycExptInd = exptCycLengthFr == 11;

respCells_allExpt = cell2mat(responsiveCellsAV');
lateStimModCells_allExpt = cell2mat(lateRespCellsAV') | cell2mat(lateSuppCellsAV');
lateRespCells_allExpt = cell2mat(lateRespCellsAV');


respCells_shortCycExpt = cell2mat(responsiveCellsAV(shortCycExptInd)');
lateStimModCells_shortCycExpt = cell2mat(lateRespCellsAV(shortCycExptInd)')...
    | cell2mat(lateSuppCellsAV(shortCycExptInd)');
lateRespCells_shortCycExpt = cell2mat(lateRespCellsAV(shortCycExptInd)');

longTC = cell(1,nav);
lateCycResp = cell(1,nav);
% lateCycWin = cell(1,nav);
lateWin = cell(1,nav);
for iav = 1:nav
    longTC{iav} = trialTCEaCycEaCell(nCycles,shortCycExptInd,iav);
    lateCycResp{iav} = cell2mat(lateCycRespEaExpt(:,iav)');
%     lateCycWin{iav} = cell2mat(lateCycMeanEaExpt(:,iav)');
    lateWin{iav} = cell2mat(lateWinRespAVEaExpt(:,iav)');
end
longTC = cellfun(@cell2mat,longTC,'unif',0);

siCycResp = cell2mat(lateSIEaExpt');
% siCycWin = cell2mat(lateCycWinSIEaExpt');
% siLateWin = cell2mat(lateWinSIEaExpt');

lateAVModCycResp = cell2mat(ttestLateAVEaExpt');
% lateAVModCycWin = cell2mat(ttestLateCycWinAVEaExpt');
% lateAVModLateWin = cell2mat(ttestLateWinAVEaExpt');

tcAVallCells_shortCycExpt = cell2mat(longTrialTCAV(shortCycExptInd));
tcAVAllCells_allExpt = cell2mat(cellfun(@(x) x(1:size(tcAVallCells_shortCycExpt,1),:),...
    longTrialTCAV,'unif',0));
%%
responseLim = [-0.005 0.25];
cycResponseLim = [-0.005 0.1];
scatLim = [-0.15 1];
scatRespLim = [-0.02 0.1];
ttLabel = 0:500:2500;
nFr = size(longTC{1,1},1);
tt = ((26:nFr)-33).*(1000/frRateHz);
ttFr = (26:nFr)-33;
ttLabelFr = ((ttLabel./1000)*frRateHz)+8;
lateWinTT = ([nFr - lateWinFrames nFr] - 33)...
    .*(1000/frRateHz);
binEdgesSI = -10:0.5:10;
siLim = [-15 15];
hmLim = [-0.2 0.2];
exCellRespLim = [-0.11 0.25];
cycTCLim = [-0.005 0.015];
cycRespLim = [-0.005 0.025];
normRespLim = [-6 6];
normTCLim = [-0.5 1.5];
%%
setFigParams4Print('portrait')
figure
subplot 211
for iav = 1:nav
    y = mean(longTC{iav}(:,respCells_shortCycExpt),2);
    yerr = ste(longTC{iav}(:,respCells_shortCycExpt),2);    
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
    num2str(sum(respCells_shortCycExpt)),num2str(length(respCells_shortCycExpt))))

subplot 212
x = lateWin{visualTrials}(respCells_allExpt);
y = lateWin{auditoryTrials}(respCells_allExpt);
h = scatter(x,y,'ko');
hold on
errorbar(mean(x),mean(y),ste(y,2),ste(y,2),ste(x,2),ste(x,2),'ro')
plot(scatLim,scatLim,'k--');
[~,p] = ttest(x,y);
figXAxis([],'Visual (dF/F)',scatLim)
figYAxis([],'Auditory (dF/F)',scatLim)
figAxForm([])
title(sprintf('Responsive Cells (%s), p = %s',num2str(sum(respCells_allExpt)),...
    num2str(round(p,2,'significant'))))

print([fnout '_respCellsTC&Scatter'],'-dpdf','-fillpage')

figure
subplot 211
for iav = 1:nav
    y = mean(longTC{iav}(:,respCells_shortCycExpt | lateRespCells_shortCycExpt),2);
    yerr = ste(longTC{iav}(:,respCells_shortCycExpt | lateRespCells_shortCycExpt),2);    
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
    num2str(sum(respCells_shortCycExpt | lateRespCells_shortCycExpt)),...
    num2str(length(respCells_shortCycExpt | lateRespCells_shortCycExpt))))

subplot 212
x = lateWin{visualTrials}(respCells_allExpt | lateRespCells_allExpt);
y = lateWin{auditoryTrials}(respCells_allExpt | lateRespCells_allExpt);
h = scatter(x,y,'ko');
hold on
errorbar(mean(x),mean(y),ste(y,2),ste(y,2),ste(x,2),ste(x,2),'ro')
plot(scatLim,scatLim,'k--');
[~,p] = ttest(x,y);
figXAxis([],'Visual (dF/F)',scatLim)
figYAxis([],'Auditory (dF/F)',scatLim)
figAxForm([])
title(sprintf('Responsive Cells (%s), p = %s',num2str(sum(respCells_allExpt)),...
    num2str(round(p,2,'significant'))))

print([fnout '_allRespCellsTC&Scatter'],'-dpdf','-fillpage')


setFigParams4Print('landscape')
figure
subplot 121
si = siCycResp;
ind = lateAVModCycResp;
h = histogram(si(respCells_allExpt),binEdgesSI);
h.FaceColor = [0.5 0.5 0.5];
h.EdgeColor = 'none';
h = vline(0,'k-');
h.LineWidth = 2;
h = vline(mean(si(respCells_allExpt)),'r-',...
    sprintf('all %s',num2str(mean(si(respCells_allExpt)))));
h.LineWidth = 2;
hold on
h = histogram(si(ind & respCells_allExpt),binEdgesSI);
h.FaceColor = 'k';
h.EdgeColor = 'none';
h = vline(mean(si(respCells_allExpt & ind)),'m-',...
    sprintf('sn %s',num2str(mean(si(respCells_allExpt & ind)))));
h.LineWidth = 2;
figXAxis([],'Late Trial Selectivity Index (using resp to cycles 5-8)',siLim)
figYAxis([],'n ncells',[])
figAxForm([])
title({sprintf('Significantly Attn Modulated (%s/%s/%s)'...
    ,num2str(sum(ind & respCells_allExpt)),num2str(sum(respCells_allExpt))...
    ,num2str(length(respCells_allExpt))),...
    'SI calc using cyc reposnses'})

subplot 122
x = lateCycResp{visualTrials}(respCells_allExpt);
y = lateCycResp{auditoryTrials}(respCells_allExpt);
h = scatter(x,y,'ko');
[~,p] = ttest(x,y);
hold on
x = lateCycResp{visualTrials}(respCells_allExpt & ind & si > 0);
y = lateCycResp{auditoryTrials}(respCells_allExpt & ind & si > 0);
h = scatter(x,y,'ro');
x = lateCycResp{visualTrials}(respCells_allExpt & ind & si < 0);
y = lateCycResp{auditoryTrials}(respCells_allExpt & ind & si < 0);
h = scatter(x,y,'bo');
plot(scatRespLim,scatRespLim,'k--');
figXAxis([],'Visual (dF/F)',scatRespLim)
figYAxis([],'Auditory (dF/F)',scatRespLim)
figAxForm([])
title(sprintf('Responsive Cells (%s/%s), p = %s',num2str(sum(respCells_allExpt)),...
    num2str(length(respCells_allExpt)), num2str(round(p,2,'significant'))))

print([fnout '_respCellsAttnSI&Scat'],'-dpdf','-fillpage')
%% example cells
if strcmp(ds,'FSAV_attentionV1')
%     tcAVAllCellsSem = cell2mat(longTrialTCAV_sem(shortCycExptInd));
tcAVAllCellsSem_allExpt = cell2mat(cellfun(@(x) x(1:size(tcAVallCells_shortCycExpt,1),:),...
    longTrialTCAV_sem,'unif',0));
    exCellTC = tcAVAllCells_allExpt(26:end,:);
    exCellTCsem = tcAVAllCellsSem_allExpt(26:end,:);

    exRespCell = exampleCell_1; 
    exLateRespCell = exampleCell_2; 
    exLateSuppCell = exampleCell_3;
    exampleCells = [exRespCell exLateRespCell exLateSuppCell];

    setFigParams4Print('landscape')
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
        figXAxis([],'time (ms)',[tt(1) tt(end)],ttLabel,ttLabel)
        figYAxis([],'dF/F',exCellRespLim)
        figAxForm([],0)
        hold on
        vline(0,'k--')
        vline(0:350:tt(end),'k:')
    end

    print([fnout '_exCellTC'],'-dpdf','-fillpage')
else
    exampleCells = 0;
end
%% heatmaps
avTCrespCells = tcAVallCells_shortCycExpt(:,respCells_shortCycExpt);
TCrespCells = cell(1,nav);
for iav = 1:nav
    TCrespCells{iav} = longTC{iav}(:,respCells_shortCycExpt);
end

[~,sortInd] = sort(mean(avTCrespCells((end-lateWinFrames):end,:),1));
hm = cellfun(@(x) flipud(x(:,sortInd)'),TCrespCells,'unif',0);

setFigParams4Print('landscape')
figure
colormap(brewermap([],'*RdBu'));
suptitle(sprintf('Responsive Cells (%s)',num2str(sum(respCells_shortCycExpt))));

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

allCells = respCells_shortCycExpt | lateStimModCells_shortCycExpt;
avTCallCells = tcAVAllCells_allExpt;
[~,sortInd] = sort(mean(avTCallCells((end-lateWinFrames):end,:),1));
hm = flipud(avTCallCells(:,sortInd)');
if strcmp(ds,'FSAV_attentionV1')
    exCellInd = find(fliplr(ismember(sortInd,exampleCells)));
else
    exCellInd = [];
end

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
%         cycTC{icyc,iav} = cell2mat(cellfun(@(x) x(1:41,:),...
%             tcEaCellEaCycAllOutcome(icyc,shortCycExptInd,iav)));
        
        cycTC{icyc,iav} = cell2mat(cellfun(@(x) x(1:41,:),...
            tcEaCellEaCycAllOutcome(icyc,:,iav),'unif',0));
    end
%     lateCycResp{iav} = cell2mat(lateCycRespEaExpt(shortCycExptInd,iav)');
    lateCycResp{iav} = cell2mat(lateCycRespEaExpt(:,iav)');
end
cycResp = cellfun(@(x) mean(x(respwin,:),1),cycTC,'unif',0);

% minRespCells = respCells_shortCycExpt & cycResp{1,visualTrials} > minRespThreshold & ...
%     cycResp{1,auditoryTrials} > minRespThreshold;
minRespCells = respCells_allExpt & cycResp{1,visualTrials} > minRespThreshold & ...
    cycResp{1,auditoryTrials} > minRespThreshold;

earlyLateCycResp = cell(2,nav);
earlyLateCycResp(1,:) = cycResp(1,:);
earlyLateCycResp(2,:) = lateCycResp;

normCycResp = cell(nCycles,nav);
normEarlyLateCycResp = cell(2,nav);
normCycTCLateCyc = cell(2,nav);
for iav = 1:nav
    normCycResp(:,iav) = cellfun(@(x) x./cycResp{1,iav},cycResp(:,iav),'unif',0);
    normEarlyLateCycResp(:,iav) =  ...
        cellfun(@(x) x./earlyLateCycResp{1,iav},earlyLateCycResp(:,iav),'unif',0);
    tc = cycTC{1,iav};
    normCycTCLateCyc{1,iav} = tc./earlyLateCycResp{1,iav};
    tc = cell2mat(cellfun(@(x) x(1:41,:),cycTCAllLateTrials(:,iav),'unif',0)');
    normCycTCLateCyc{2,iav} = tc./earlyLateCycResp{1,iav};
end
respCellsNormCycResp = cellfun(@(x) mean(x(minRespCells)),normCycResp);

%% plot adaptation
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
    figYAxis([],'dF/F',cycResponseLim)
    figAxForm
    title(['Stimulus #' num2str(icyc)])
end
print([fnout '_cycTC'],'-dpdf','-fillpage')

setFigParams4Print('landscape')
figure
for icyc = 1:nCycles
    subplot(2,4,icyc)
    for iav = 1:nav
        y = mean(cycTC{icyc,iav}(26:end,lateRespCells_shortCycExpt),2);
        yerr = ste(cycTC{icyc,iav}(26:end,lateRespCells_shortCycExpt),2);
        h = shadedErrorBar(ttCyc,y,yerr,avcolmat{iav});
        for i = 1:2
            h.edge(i).Color = 'none';
        end
        hold on
    end
    vline(ttCycResp,'k--')
    figXAxis([],'Time (ms)',[ttCyc(1) ttCyc(end)])
    figYAxis([],'dF/F',cycResponseLim)
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
%% plot tc of first and late cycles
figure
    earlyLateTC = cell(2,2);
    subplot 121
    for iav = 1:nav
        earlyLateTC{1,iav} = cycTC{1,iav};
        y = mean(cycTC{1,iav}(26:end,minRespCells),2);
        yerr = ste(cycTC{1,iav}(26:end,minRespCells),2);
        h = shadedErrorBar(ttCyc,y,yerr,avcolmat{iav});
        for i = 1:2
            h.edge(i).Color = 'none';
        end
        hold on
    end
    vline(ttCycResp,'k--')
    figXAxis([],'Time (ms)',[ttCyc(1) ttCyc(end)])
    figYAxis([],'dF/F',cycRespLim)
    figAxForm
    title('First Stim')
    
    subplot 122
    for iav = 1:nav
%         tc = cell2mat(cycTCAllLateTrials(shortCycExptInd,iav)');
        tc = cell2mat(cellfun(@(x) x(1:41,:),cycTCAllLateTrials(:,iav),'unif',0)');
        earlyLateTC{2,iav} = tc;
        y = mean(tc(26:end,minRespCells),2);
        yerr = ste(tc(26:end,minRespCells),2);
        h = shadedErrorBar(ttCyc,y,yerr,avcolmat{iav});
        for i = 1:2
            h.edge(i).Color = 'none';
        end
        hold on
    end
    vline(ttCycResp,'k--')
    figXAxis([],'Time (ms)',[ttCyc(1) ttCyc(end)])
    figYAxis([],'dF/F',cycRespLim)
    figAxForm
    title('Late Stim')

    earlyLateNormResp = cellfun(@(x,y) ...
        mean(x(respwin,:),1)./mean(y(respwin,:),1),earlyLateTC(2,:),...
        earlyLateTC(1,:),'unif',0);
    
    figure
    subplot 221
    for iav = 1:nav
        y = mean(normCycTCLateCyc{1,iav}(26:end,minRespCells),2);
        yerr = ste(normCycTCLateCyc{1,iav}(26:end,minRespCells),2);
        hold on
        h = shadedErrorBar(ttCyc,y,yerr);
        h.mainLine.Color = cueColor{iav};
        h.patch.FaceColor = errorPatchColor(h.mainLine.Color);
        for i = 1:2
            h.edge(i).Color = 'none';
        end
    end
    vline(ttCycResp,'k--')
    figXAxis([],'Time (ms)',[ttCyc(1) ttCyc(end)])
    figYAxis([],'dF/F',normTCLim)
    figAxForm
    title('First Stim, Normalized')
    
    subplot 222
    for iav = 1:nav
        y = mean(normCycTCLateCyc{2,iav}(26:end,minRespCells),2);
        yerr = ste(normCycTCLateCyc{2,iav}(26:end,minRespCells),2);
        hold on
        h = shadedErrorBar(ttCyc,y,yerr);
        h.mainLine.Color = cueColor{iav};
        h.patch.FaceColor = errorPatchColor(h.mainLine.Color);
        for i = 1:2
            h.edge(i).Color = 'none';
        end
    end
    vline(ttCycResp,'k--')
    figXAxis([],'Time (ms)',[ttCyc(1) ttCyc(end)])
    figYAxis([],'dF/F',normTCLim)
    figAxForm
    title('Late Stim, Normalized')
    
    subplot 223
    x = earlyLateNormResp{visualTrials}(minRespCells);
    y = earlyLateNormResp{auditoryTrials}(minRespCells);
    plot(x,y,'ko')
    hold on
    errorbar(mean(x),mean(y),ste(y,2),ste(y,2),ste(x,2),ste(x,2),'ro')
    plot(normRespLim,normRespLim,'b--')
    hline(0,'b')
    vline(0,'b')
    hline(1,'b:')
    vline(1,'b:')
    figXAxis([],'Visual',normRespLim)
    figYAxis([],'Auditory',normRespLim)
    figAxForm
    title('Late Cycle Response Normalized to First Stim')    

    subplot 224
    y = cellfun(@(x) mean(x(minRespCells)), earlyLateNormResp);
    yerr = cellfun(@(x) ste(x(minRespCells),2), earlyLateNormResp);
    d = cellfun(@(x) x(minRespCells), earlyLateNormResp,'unif',0);
    [~,p] = ttest(d{visualTrials},d{auditoryTrials},'alpha',0.05);
    h = errorbar(1:2,y,yerr);
    h.Color = 'k';
    h.LineWidth = 2;
    figXAxis([],'',[0 3],1:2,{'Visual','Auditory'})
    figYAxis([],'Normalized dF/F',[0 1.1])
    figAxForm([])
    title(sprintf('p=%s',num2str(p)))
    
    print([fnout '_adaptCycResp'],'-dpdf','-fillpage')


figure
    earlyLateTC = cell(2,2);
    subplot 121
    for iav = 1:nav
        earlyLateTC{1,iav} = cycTC{1,iav};
        y = mean(cycTC{1,iav}(26:end,minRespCells),2);
        yerr = ste(cycTC{1,iav}(26:end,minRespCells),2);
        h = shadedErrorBar(ttCyc,y,yerr,avcolmat{iav});
        for i = 1:2
            h.edge(i).Color = 'none';
        end
        hold on
    end
    vline(ttCycResp,'k--')
    figXAxis([],'Time (ms)',[ttCyc(1) ttCyc(end)])
    figYAxis([],'dF/F',cycRespLim)
    figAxForm
    title('First Stim')
    
    subplot 122
    for iav = 1:nav
%         tc = cell2mat(cycTCAllLateTrials(shortCycExptInd,iav)');
        tc = cell2mat(cellfun(@(x) x(1:41,:),cycTCAllLateTrials(:,iav),'unif',0)');
        earlyLateTC{2,iav} = tc;
        y = mean(tc(26:end,minRespCells),2);
        yerr = ste(tc(26:end,minRespCells),2);
        h = shadedErrorBar(ttCyc,y,yerr,avcolmat{iav});
        for i = 1:2
            h.edge(i).Color = 'none';
        end
        hold on
    end
    vline(ttCycResp,'k--')
    figXAxis([],'Time (ms)',[ttCyc(1) ttCyc(end)])
    figYAxis([],'dF/F',cycRespLim)
    figAxForm
    title('Late Stim')

    earlyLateNormResp = cellfun(@(x,y) ...
        mean(x(respwin,:),1)./mean(y(respwin,:),1),earlyLateTC(2,:),...
        earlyLateTC(1,:),'unif',0);
    
    earlyLateTC_rect = cell(2,2);
    for i = 1:2
       for iav = 1:2
            d = earlyLateTC{i,iav};
            d(d<0) = 0;
            earlyLateTC_rect{i,iav} = d;
        end
    end
    
    earlyLateNormResp_rect = cellfun(@(x,y) ...
        mean(x(respwin,:),1)./mean(y(respwin,:),1),earlyLateTC_rect(2,:),...
        earlyLateTC_rect(1,:),'unif',0);
    %%
setFigParams4Print('landscape')
    figure
    ind = minRespCells & lateAVModCycResp;
    suptitle(sprintf('Minimally responsive and AV mod cells (%s/%s)',...
        num2str(sum(ind)),num2str(length(ind))))
    subplot 221
    for iav = 1:nav
        y = mean(normCycTCLateCyc{1,iav}(26:end,ind),2);
        yerr = ste(normCycTCLateCyc{1,iav}(26:end,ind),2);
        hold on
        h = shadedErrorBar(ttCyc,y,yerr);
        h.mainLine.Color = cueColor{iav};
        h.patch.FaceColor = errorPatchColor(h.mainLine.Color);
        for i = 1:2
            h.edge(i).Color = 'none';
        end
    end
    vline(ttCycResp,'k--')
    figXAxis([],'Time (ms)',[ttCyc(1) ttCyc(end)])
    figYAxis([],'dF/F',normTCLim)
    figAxForm
    title('First Stim, Normalized')
    
    subplot 222
    for iav = 1:nav
        y = mean(normCycTCLateCyc{2,iav}(26:end,ind),2);
        yerr = ste(normCycTCLateCyc{2,iav}(26:end,ind),2);
        hold on
        h = shadedErrorBar(ttCyc,y,yerr);
        h.mainLine.Color = cueColor{iav};
        h.patch.FaceColor = errorPatchColor(h.mainLine.Color);
        for i = 1:2
            h.edge(i).Color = 'none';
        end
    end
    vline(ttCycResp,'k--')
    figXAxis([],'Time (ms)',[ttCyc(1) ttCyc(end)])
    figYAxis([],'dF/F',normTCLim)
    figAxForm
    title('Late Stim, Normalized')
    
    subplot 223
    x = earlyLateNormResp{visualTrials}(ind);
    y = earlyLateNormResp{auditoryTrials}(ind);
    plot(x,y,'ko')
    hold on
    errorbar(mean(x),mean(y),ste(y,2),ste(y,2),ste(x,2),ste(x,2),'ro')
    plot(normRespLim,normRespLim,'b--')
    hline(0,'b')
    vline(0,'b')
    hline(1,'b:')
    vline(1,'b:')
    figXAxis([],'Visual',normRespLim)
    figYAxis([],'Auditory',normRespLim)
    figAxForm
    title('Late Cycle Response Normalized to First Stim')    

    subplot 224
    y = cellfun(@(x) mean(x(ind)), earlyLateNormResp);
    yerr = cellfun(@(x) ste(x(ind),2), earlyLateNormResp);
    d = cellfun(@(x) x(ind),earlyLateNormResp,'unif',0);
    [~,p] = ttest(d{visualTrials},d{auditoryTrials},'alpha',0.05);
    h = errorbar(1:2,y,yerr);
    h.Color = 'k';
    h.LineWidth = 2;
    figXAxis([],'',[0 3],1:2,{'Visual','Auditory'})
    figYAxis([],'Normalized dF/F',[0 1.1])
    figAxForm([])
    title(sprintf('p=%s',num2str(p)))
    
print([fnout '_adaptCycResp_modCells'],'-dpdf','-fillpage')
%% tuning
oris = 0:22.5:180;
oriLabel = (oris)+1;
oriFit_allExpt = cell2mat(oriTuningFitEaExpt');
[~, oriPeak_allExpt] = max(oriFit_allExpt,[],1);
oriPeak_allExpt = oriPeak_allExpt-1;
oriBins = [0 22.5:45:157.5 180];
[~,~,oriPref_allExpt] = histcounts(oriPeak_allExpt,oriBins);
oriPref_allExpt(oriPref_allExpt == 5) = 1;
nori = 4;

oriFit_shortCycExpt = cell2mat(oriTuningFitEaExpt(shortCycExptInd)');
[~, oriPeak_shortCycExpt] = max(oriFit_shortCycExpt,[],1);
oriPeak_shortCycExpt = oriPeak_shortCycExpt-1;
[~,~,oriPref_shortCycExpt] = histcounts(oriPeak_shortCycExpt,oriBins);
oriPref_shortCycExpt(oriPref_shortCycExpt == 5) = 1;

isTuned_allExpt = cell2mat(oriTuningFitReliability') < tuningReliabilityThresh;
isTuned_shortCycExpt = cell2mat(oriTuningFitReliability(shortCycExptInd)') < tuningReliabilityThresh;

oriFitLim = [-0.5 0.5];

% avTCtunedCells = tcAVallCells_shortCycExpt(:,respCells_shortCycExpt & isTuned_shortCycExpt);
oriFits_tunedResponsive = oriFit_allExpt(:,respCells_allExpt & isTuned_allExpt);

%%
% oriFitNormLim = [-1 1];
% figure
% colormap(brewermap([],'*RdBu'));
% hm = [];
% n = zeros(1,nori);
% for i = 1:nori
%     ind = oriPref_shortCycExpt == i & isTuned_shortCycExpt & minRespCells;
%     [~,sortInd] = sort(oriPeak_allExpt(ind));
%     fits = oriFit_allExpt(:,ind);
%     normFits = fits./(max(fits,[],1));
%     hm = cat(1,hm,normFits(:,sortInd)');
%     if i > 1
%         n(i) = sum(ind)+n(i-1);
%     else
%         n(i) = sum(ind);
%     end
% end
% h = imagesc(hm);
% hold on
% hline(n,'k-');
% figXAxis([],'Orientation (deg)',[1 181],oriLabel,oris)
% figYAxis([],'Cell #',[])
% figAxForm([],0)
% colorbar
% caxis(oriFitNormLim)
% title('Responsive tuned neurons sorted by pref')
% print([fnout '_oriTuningFits_sorted'],'-dpdf')
%%
setFigParams4Print('landscape')
figure
suptitle(sprintf('Responsive & Tuned Cells (%s)',num2str(sum(respCells_allExpt & isTuned_allExpt))));

subplot 221
h = polarplot(deg2rad(oriPeak_allExpt(respCells_allExpt & isTuned_allExpt)),siCycResp(respCells_allExpt & isTuned_allExpt),'k.');
h.Color = [0.5 0.5 0.5];
h.MarkerSize = 12;
hold on
h = polarplot(deg2rad(oriPeak_allExpt(respCells_allExpt & isTuned_allExpt & lateAVModCycResp)),...
    siCycResp(respCells_allExpt & isTuned_allExpt & lateAVModCycResp),'k.');
h.MarkerSize = 12;
title('Attn Selectivity by Ori Tuning')

subplot 222
visX = 1:3:10;
audX = visX+1;
orientations = [0 45 90 135];
leg = [];
n = cell(1,nori);
for iori = 1:nori
    ind = oriPref_allExpt == iori;
    y = cellfun(@(x) mean(x(minRespCells & isTuned_allExpt & ind)),...
        earlyLateNormResp);
    yerr = cellfun(@(x) ste(x(minRespCells & isTuned_allExpt & ind),2),...
        earlyLateNormResp);
    d = cellfun(@(x) x(minRespCells & isTuned_allExpt & ind),...
        earlyLateNormResp,'unif',0);
    [~,p] = ttest(d{visualTrials},d{auditoryTrials},'alpha',0.05/nori);
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
    n{iori} = sprintf('%s: %s, p = %s',num2str(orientations(iori)),...
        num2str(sum(minRespCells & isTuned_allExpt & ind)),...
        num2str(round(p,2,'significant')));
end
figXAxis([],'Preferred Orientation (deg)',[0 12],visX+0.5,orientations)
figYAxis([],'Normalized dF/F',[0 1.1])
figAxForm([],0)
L = legend(leg,n,'location','northwest');
title('Responsive, tuned Cells with minimum response cutoff')
title(L,'N Cells Per Group')

subplot 224
leg = [];
n = cell(1,nori);
for iori = 1:nori
    ind = oriPref_allExpt == iori & lateAVModCycResp;
    y = cellfun(@(x) mean(x(minRespCells & isTuned_allExpt & ind)),...
        earlyLateNormResp);
    yerr = cellfun(@(x) ste(x(minRespCells & isTuned_allExpt & ind),2),...
        earlyLateNormResp);
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
        num2str(sum(minRespCells & isTuned_allExpt & ind)));
end
figXAxis([],'Preferred Orientation (deg)',[0 12],visX+0.5,orientations)
figYAxis([],'Normalized dF/F',[0 1.1])
figAxForm([],0)
L = legend(leg,n,'location','northwest');
title('Responsive, tuned, modulated cells')
title(L,'N Cells Per Group')

print([fnout '_oriTuningSIandAdapt'],'-dpdf','-fillpage')

%% each cell adaptation
figure
orientations = [0 45 90 135];
leg = [];
n = cell(1,nori);
for iori = 1:nori
    subplot(2,2,iori)
    ind = oriPref_allExpt == iori & minRespCells & isTuned_allExpt;
    y_all = cell2mat(cellfun(@(x) x(ind),...
        earlyLateNormResp,'unif',0)')';
    y = cellfun(@(x) mean(x(ind)),...
        earlyLateNormResp);
    yerr = cellfun(@(x) ste(x(ind),2),...
        earlyLateNormResp);
    d = cellfun(@(x) x(ind),...
        earlyLateNormResp,'unif',0);
    [~,p] = ttest(d{visualTrials},d{auditoryTrials},'alpha',0.05/nori);
    nc = size(y_all,1);
    rc = lateAVModCycResp(ind);
    for i = 1:nc
        h = plot(1:2,y_all(i,:),'-');
        if rc(i)
            h.Color = [0.5 0 0];
        else
            h.Color = [0.5 0.5 0.5];
        end
        hold on
    end
    h = errorbar(1:2,y,yerr,'k-');
    h.LineWidth = 2;
    leg(iori) = h;
    n{iori} = sprintf('%s: %s, p = %s',num2str(orientations(iori)),...
        num2str(sum(minRespCells & isTuned_allExpt & ind)),...
        num2str(round(p,2,'significant')));
    figXAxis([],'Trial Type',[0 3],1:2,{'Vis';'Aud'})
    figYAxis([],'Normalized dF/F',[-1 2])
    figAxForm
    title(sprintf('Pref %s',num2str(orientations(iori))))
end
% L = legend(leg,n,'location','northwest');
% title('Responsive, tuned Cells with minimum response cutoff')
% title(L,'N Cells Per Group')
%% first and late stim dF/F, tuning groups and all responsive cells, plotted with lines and errorbars
tcLim = [-0.003 0.1];
respLim = [-0.003 0.05];
respLim2 = [0 0.01];
nRows = nori;
nCols = 4;
setFigParams4Print('landscape')
figure
for iori = 1:nori
    iplot = 1+(4*(iori-1));
    subplot (nRows,nCols,iplot)
    ind = oriPref_shortCycExpt == iori & isTuned_shortCycExpt & respCells_shortCycExpt;
    for iav = 1:nav
        y = mean(longTC{iav}(:,ind),2);
        yerr = ste(longTC{iav}(:,ind),2);
        tt = ((26:length(y))-33).*(1000/frRateHz);
        hold on
        h = shadedErrorBar_chooseColor(tt,y(26:length(y)),yerr(26:length(y)),cueColor{iav});
    end
    figXAxis([],'Time (ms)',[tt(1) tt(end)],ttLabel,ttLabel)
    figYAxis([],'dF/F',tcLim)  
    vline(lateWinTT,'k--')
    figAxForm([],0)
    title(sprintf('%sdeg-Tuned Cells (%s/%s)',...
        num2str(orientations(iori)),...
        num2str(sum(ind)),num2str(length(respCells_shortCycExpt))))
    
    iplot = 2+(4*(iori-1));
    subplot (nRows,nCols,iplot)
    ind = oriPref_allExpt == iori & isTuned_allExpt & respCells_allExpt;
    for iav = 1:nav
        y = mean(earlyLateTC{1,iav}(26:end,ind),2);
        yerr = ste(earlyLateTC{1,iav}(26:end,ind),2);
        h = shadedErrorBar_chooseColor(ttCyc,y,yerr,cueColor{iav});
%         for i = 1:2
%             h.edge(i).Color = 'none';
%         end
        hold on
    end
    vline(ttCycResp,'k--')
    figXAxis([],'Time (ms)',[ttCyc(1) ttCyc(end)])
    figYAxis([],'dF/F',respLim)
    figAxForm
    title('First Stim')
    
    iplot = 3+(4*(iori-1));
    subplot (nRows,nCols,iplot)
    ind = oriPref_allExpt == iori & isTuned_allExpt & respCells_allExpt;
    for iav = 1:nav
        y = mean(earlyLateTC{2,iav}(26:end,ind),2);
        yerr = ste(earlyLateTC{2,iav}(26:end,ind),2);
        h = shadedErrorBar_chooseColor(ttCyc,y,yerr,cueColor{iav});
%         for i = 1:2
%             h.edge(i).Color = 'none';
%         end
        hold on
    end
    vline(ttCycResp,'k--')
    figXAxis([],'Time (ms)',[ttCyc(1) ttCyc(end)])
    figYAxis([],'dF/F',respLim)
    figAxForm
    title('Late Stim')
    
    iplot = 4+(4*(iori-1));
    subplot (nRows,nCols,iplot)
    y = cellfun(@(x) mean(x(ind)),earlyLateCycResp(2,:));
    yerr = cellfun(@(x) ste(x(ind),2),earlyLateCycResp(2,:));
    d = cellfun(@(x) x(ind),earlyLateCycResp(2,:),'unif',0);
    [~,p] = ttest(d{visualTrials},d{auditoryTrials},'alpha',0.05/nori);
    h = errorbar(1:2,y,yerr);
    h.Color = 'k';
    h.LineWidth = 2;
    figXAxis([],'',[0 3],1:2,{'V';'A'})
    figYAxis([],'dF/F',respLim2)
    figAxForm
    title(sprintf('p=%s',num2str(p)))
end

print([fnout '_tuningFirstAndLate'],'-dpdf','-fillpage')
%% first and late stim dF/F, tuning groups

setFigParams4Print('landscape')
figure
suptitle('Responsive, tuned Cells with minimum response cutoff')
for istim = 1:2
    subplot(1,2,istim)
    visX = 1:3:10;
    audX = visX+1;
    orientations = [0 45 90 135];
    leg = [];
    n = cell(1,nori);
    for iori = 1:nori
        ind = oriPref_allExpt == iori;
        y = cellfun(@(x) mean(x(minRespCells & isTuned_allExpt & ind)),...
            earlyLateCycResp(istim,:));
        yerr = cellfun(@(x) ste(x(minRespCells & isTuned_allExpt & ind),2),...
            earlyLateCycResp(istim,:));
        d = cellfun(@(x) x(minRespCells & isTuned_allExpt & ind),...
            earlyLateCycResp(istim,:),'unif',0);
        [~,p] = ttest(d{visualTrials},d{auditoryTrials},'alpha',0.05/nori);
        for iav = 1:nav
            if iav == 1
                x = visX(iori);
            else
                x = audX(iori);
            end
            h = bar(x,y(iav));
            h.BarWidth = 1;
            h.EdgeColor = 'none';
            h.FaceColor = cueColor{iav};
            hold on
            h = errorbar(x,y(iav),yerr(iav));
            h.Color = cueColor{iav};
        end
        leg(iori) = h;
        n{iori} = sprintf('%s: %s, p = %s',num2str(orientations(iori)),...
            num2str(sum(minRespCells & isTuned_allExpt & ind)),...
            num2str(round(p,2,'significant')));
    end
    figXAxis([],'Preferred Orientation (deg)',[0 12],visX+0.5,orientations)
    figYAxis([],'dF/F',[-0.005 0.04])
    figAxForm([],0)
    L = legend(leg,n,'location','northwest');
    if istim == 1
        title('First Stim')
    else
        title('Late Stim')
    end
    title(L,'N Cells Per Group')
end
print([fnout '_oriTuningCycResp'],'-dpdf','-fillpage')

%% each stim adaptation each tuning group

figure
suptitle(sprintf('Tuned & Min Resp Cells (%s/%s)',num2str(sum(minRespCells & isTuned_allExpt)),num2str(length(minRespCells))))
orientations = [0 45 90 135];
leg = [];
n = cell(1,nori);
for iori = 1:nori
    subplot(2,2,iori)
    ind = oriPref_allExpt == iori;
    for icyc = 1:nCycles
        y = cellfun(@(x) mean(x(minRespCells & isTuned_allExpt & ind)),...
            normCycResp(icyc,:));
        yerr = cellfun(@(x) ste(x(minRespCells & isTuned_allExpt & ind),2),...
            normCycResp(icyc,:));
%         d = cellfun(@(x) x(minRespCells & isTuned_allExpt & ind),...
%             normCycResp(2,:),'unif',0);
%         [~,p] = ttest(d{visualTrials},d{auditoryTrials},'alpha',0.05/nori);
        for iav = 1:nav
            h = errorbar(icyc,y(iav),yerr(iav),'o');
            h.Color = cueColor{iav};
            h.MarkerFaceColor = [1 1 1];
            hold on
        end
    end
    leg(iori) = h;
    n{iori} = sprintf('%s: %s, p = %s',num2str(orientations(iori)),...
        num2str(sum(minRespCells & isTuned_allExpt & ind)),...
        num2str(round(p,2,'significant')));
    figXAxis([],'Stimulus Number',[0 nCycles+1],1:nCycles,1:nCycles)
    figYAxis([],'Norm. dF/F',[0 1.1])
    figAxForm
    title(sprintf('Pref %s (n=%s)',num2str(orientations(iori)),num2str(sum(ind & minRespCells & isTuned_allExpt))))
end
% L = legend(leg,n,'location','northwest');
% title(L,'N Cells Per Group')
print([fnout '_oriTuningCycResp'],'-dpdf','-fillpage')

figure
colors = brewermap(nori+1,'PuRd');
colors = colors(2:end,:);
leg = [];
for iori = 1:nori
    ind = oriPref_allExpt == iori;
    for icyc = 1:nCycles
        y = cellfun(@(x) mean(x(minRespCells & isTuned_allExpt & ind)),...
            normCycResp(icyc,:));
        yerr = cellfun(@(x) ste(x(minRespCells & isTuned_allExpt & ind),2),...
            normCycResp(icyc,:));
        h = errorbar(icyc,y(visualTrials),yerr(visualTrials),'o');
        h.Color = colors(iori,:);
        h.MarkerFaceColor = [1 1 1];
        hold on
        if icyc == 1
            leg(iori) = h;
        end
    end
end
figXAxis([],'Stimulus Number',[0 nCycles+1],1:nCycles,1:nCycles)
figYAxis([],'Norm. dF/F',[0 1.1])
figAxForm
legend(leg,{'0','45','90','135'})
title(sprintf('Visual Trials, Ori Groups (n=%s)',num2str(sum(minRespCells & isTuned_allExpt))))
print([fnout '_oriTuningCycResp_visOnly'],'-dpdf')

%%
figure
suptitle('Responsive, tuned Cells with minimum response cutoff')
orientations = [0 45 90 135];
leg = [];
n = cell(1,nori);
for iori = 1:nori
    ind = oriPref_allExpt == iori;
    d = siCycResp(lateRespCells_allExpt & isTuned_allExpt & ind);
    y = mean(d,2);
    yerr = ste(d,2);
    [~,p] = ttest(d);
        h = bar(iori,y);
        h.BarWidth = .5;
        h.EdgeColor = 'none';
        h.FaceColor = cueColor{1};
        hold on
        h = errorbar(iori,y,yerr);
        h.Color = cueColor{1};
    leg(iori) = h;
    n{iori} = sprintf('%s: %s, p = %s',num2str(orientations(iori)),...
        num2str(sum(lateRespCells_allExpt & isTuned_allExpt & ind)),...
        num2str(round(p,2,'significant')));
end
L = legend(leg,n,'location','northwest');
figXAxis([],'Preferred Orientation (deg)',[0 nori+1],1:nori,orientations)
figYAxis([],'V-A Selectivity',[-2 2],-2:2)
figAxForm
print([fnout '_respTunedCellsSelectivity'],'-dpdf','-fillpage')

%% time-course of tuned cells
%time-courses
setFigParams4Print('portrait')
figure
for iori = 1:nori
subplot (nori,1,iori)
ind = oriPref_shortCycExpt == iori & isTuned_shortCycExpt & respCells_shortCycExpt;
for iav = 1:nav
    y = mean(longTC{iav}(:,ind),2);
    yerr = ste(longTC{iav}(:,ind),2);
    tt = ((26:length(y))-33).*(1000/frRateHz);
    hold on
    h = shadedErrorBar_chooseColor(tt,y(26:length(y)),yerr(26:length(y)),cueColor{iav});
end
figXAxis([],'Time (ms)',[tt(1) tt(end)],ttLabel,ttLabel)
figYAxis([],'dF/F',[-0.005 0.1])  
vline(lateWinTT,'k--')
figAxForm([],0)
title(sprintf('%sdeg-Tuned Cells (%s/%s)',...
    num2str(orientations(iori)),...
    num2str(sum(ind)),num2str(length(respCells_shortCycExpt))))
end
print([fnout '_tunedCellsTCs'],'-dpdf','-fillpage')

% heatmaps
avTCrespCells = tcAVallCells_shortCycExpt(:,respCells_shortCycExpt...
    & isTuned_shortCycExpt);
TCrespCells = cell(1,nav);
for iav = 1:nav
    TCrespCells{iav} = longTC{iav}(:,respCells_shortCycExpt...
    & isTuned_shortCycExpt);
end
[~,sortInd] = sort(mean(avTCrespCells((end-lateWinFrames):end,:),1));
hm = cellfun(@(x) flipud(x(:,sortInd)'),TCrespCells,'unif',0);
oriPref_hmCells = oriPref_shortCycExpt(isTuned_shortCycExpt & respCells_shortCycExpt);

setFigParams4Print('landscape')
figure
colormap(brewermap([],'*RdBu'));
suptitle(sprintf('Visual Trials,Responsive, Tuned Cells (%s)',num2str(length(oriPref_hmCells))));

for iori = 1:nori
    subplot(2,2,iori)
    h = imagesc(hm{visualTrials}(oriPref_hmCells == iori,26:end));
    figXAxis([],'Time (ms)',[],ttLabelFr,ttLabel)
    figYAxis([],'Cell #',[])
    figAxForm
    colorbar
    caxis([-0.2 0.2])
    title(sprintf('%sdeg-Tuned',num2str(orientations(iori))));
end
print([fnout '_tunedCellsHeatmap_vis'],'-dpdf','-fillpage')

setFigParams4Print('landscape')
figure
colormap(brewermap([],'*RdBu'));
suptitle(sprintf('Auditory Trials,Responsive, Tuned Cells (%s)',num2str(length(oriPref_hmCells))));

for iori = 1:nori
    subplot(2,2,iori)
    h = imagesc(hm{auditoryTrials}(oriPref_hmCells == iori,26:end));
    figXAxis([],'Time (ms)',[],ttLabelFr,ttLabel)
    figYAxis([],'Cell #',[])
    figAxForm
    colorbar
    caxis([-0.2 0.2])
    title(sprintf('%sdeg-Tuned',num2str(orientations(iori))));
end
print([fnout '_tunedCellsHeatmap_aud'],'-dpdf','-fillpage')


% heatmaps all tuned cells
avTCrespCells = tcAVallCells_shortCycExpt(:,isTuned_shortCycExpt);
TCrespCells = cell(1,nav);
for iav = 1:nav
    TCrespCells{iav} = longTC{iav}(:,isTuned_shortCycExpt);
end
[~,sortInd] = sort(mean(avTCrespCells((end-lateWinFrames):end,:),1));
hm = cellfun(@(x) flipud(x(:,sortInd)'),TCrespCells,'unif',0);
oriPref_hmCells = oriPref_shortCycExpt(isTuned_shortCycExpt);

setFigParams4Print('landscape')
figure
colormap(brewermap([],'*RdBu'));
suptitle(sprintf('Visual Trials, Tuned Cells (%s)',num2str(length(oriPref_hmCells))));

for iori = 1:nori
    subplot(2,2,iori)
    h = imagesc(hm{visualTrials}(oriPref_hmCells == iori,26:end));
    figXAxis([],'Time (ms)',[],ttLabelFr,ttLabel)
    figYAxis([],'Cell #',[])
    figAxForm
    colorbar
    caxis([-0.2 0.2])
    title(sprintf('%sdeg-Tuned',num2str(orientations(iori))));
end
print([fnout '_allTunedCellsHeatmap_vis'],'-dpdf','-fillpage')

setFigParams4Print('landscape')
figure
colormap(brewermap([],'*RdBu'));
suptitle(sprintf('Auditory Trials,Responsive, Tuned Cells (%s)',num2str(length(oriPref_hmCells))));

for iori = 1:nori
    subplot(2,2,iori)
    h = imagesc(hm{auditoryTrials}(oriPref_hmCells == iori,26:end));
    figXAxis([],'Time (ms)',[],ttLabelFr,ttLabel)
    figYAxis([],'Cell #',[])
    figAxForm
    colorbar
    caxis([-0.2 0.2])
    title(sprintf('%sdeg-Tuned',num2str(orientations(iori))));
end
print([fnout '_allTunedCellsHeatmap_aud'],'-dpdf','-fillpage')
