clear all
close all
ds = 'FSAV_attentionV1';
cellsOrDendrites = 1;
%%
rc = behavConstsAV;
imgParams_FSAV

eval(ds)
titleStr = ds(6:end);
mice = unique({expt.SubNum});
mouse_str = ['i' strjoin(mice,'_i')];

if cellsOrDendrites == 1
    load(fullfile(rc.caOutputDir,ds,...
        [mouse_str '_trOutcomeStruct_cells' ds(5:end) '.mat']));
    fnout = fullfile(rc.ashley, 'Manuscripts','Attention V1','Matlab Figs',...
        [titleStr '_anticipation_']); 
elseif cellsOrDendrites == 2
    load(fullfile(rc.caOutputDir,ds,...
        [mouse_str '_trOutcomeStruct_cells' ds(5:end) '.mat']));
    fnout = fullfile(rc.ashley, 'Manuscripts','Attention V1','Matlab Figs',...
        [titleStr '_den_anticipation_']); 
end
% 
% if cellsOrDendrites == 1
%     load(fullfile(rc.caOutputDir,ds,'FSAV Choice','FSAV_decodeData.mat'))
% end

%%
% outcomecolmat = {'k';'r'};
% avcolmat = {'k','c'};
avmat = {'visual';'auditory'};
outcomemat = {'hit';'miss'};
% pressAlign = 1;
% targetAlign = 2;
visualTrials = 1;
auditoryTrials = 2;
% prevTrialVisual = 0;
% prevTrialAuditory = 1;
nav = 2;
% hitTrials = 1;
% missTrials = 2;
% nout = 2;
% frameRateHz = expt(1).frame_rate;
% onTimeFr = 0.1*frameRateHz;
% nMonitorDelayFr = 2;
nBaselineFr = mouse(1).expt(1).info.preAlignFrames;
nFrames1s = frameRateHz;
nexp = size(expt,2);
% for imouse = 1:size(mouse,2)
%     nexp = nexp+size(mouse(imouse).expt,2);
% end
% exptDates = {expt.date};
% exptMice = {expt.SubNum};

% basewinTarget = 31:33;
% respwinTarget = 35:37;
nCycles = 8;
lateCycles = 5:nCycles;
lateWinFr = (45:88)+nBaselineFr;

orientations = [0 45 90 135];
%% pool experiment data
antiDataExpt = struct;
oriTuningExpt = struct;
for imouse = 1:size(mouse,2)
    for iexp = 1:size(mouse(imouse).expt,2)
        if imouse == 1 && iexp == 1
            exptN = 1;
        else
            exptN = exptN+1;
        end
        
        d = mouse(imouse).expt(iexp);
        
        cycLengthFr = d.info.cycTimeFrames;
        nCycLongTC = ceil(longTrialLengthFr./cycLengthFr);
        
        tc_AV = [];
        cycTC = cell(2,nCycles);
        longTC = cell(2,1);
        for iav = 1:2
            dd = d.av(iav).align(alignStart);
            if strcmp(ds,'FSAV_attentionV1')
                hits = strcmp(dd.outcome,'success');
            elseif strcmp(ds,'FSAV_V1_100ms_naive')
                hits = true(1,length(dd.outcome));
            end
            misses = strcmp(dd.outcome,'ignore');
            tc_AV = cat(3,tc_AV,dd.respTC(:,:,hits)); % make hits or misses
            for icyc = 1:nCycles
                tc = dd.respTC(:,:,dd.nCycles >= icyc & (hits | misses));
                cycStartOffset = ((icyc-1).*cycLengthFr)+nBaselineFr;
                cycTC{iav,icyc} = tc(...
                    (cycStartOffset-nBaselineFr+1):(cycStartOffset+nFrames1s),:,:);
            end
            longTC{iav} = dd.respTC(1:(longTrialLengthFr+nBaselineFr),:,...
                dd.nCycles >= nCycLongTC & (hits | misses));
        end
        
        antiDataExpt(exptN).exptName = [mouse(imouse).mouse_name '-' d.date];
        antiDataExpt(exptN).exptCycLengthFr = cycLengthFr;
        antiDataExpt(exptN).longTC = longTC;
        antiDataExpt(exptN).cycTC = cycTC;
        
        do = d.oriTuning;
        if size(do.oriResp,2) == 8
            oriTuningExpt(exptN).oriResp = do.oriResp(:,1:2:8);
            oriTuningExpt(exptN).oriRespErr = do.oriRespSem(:,1:2:8);
        elseif size(do.oriResp,2) ~= 4
            error('error in orientations used for passivie tuning')
        else
            oriTuningExpt(exptN).oriResp = do.oriResp;
            oriTuningExpt(exptN).oriRespErr = do.oriRespSem;
        end
        oriTuningExpt(exptN).fit = do.oriFit;
        oriTuningExpt(exptN).isTuned = do.oriFitReliability < tuningReliabilityThresh;
    end
end

shortCycExptInd = cell2mat({antiDataExpt.exptCycLengthFr}) == 11;
isShortCycExpt = [];
respCellsExpt = struct;
for iexp = 1:nexp
    firstTCAV = cat(3,antiDataExpt(iexp).cycTC{visualTrials,1},...
        antiDataExpt(iexp).cycTC{auditoryTrials,1});
    longTCAV = cat(3,antiDataExpt(iexp).longTC{visualTrials,1},...
        antiDataExpt(iexp).longTC{auditoryTrials,1});
    
    firstRespCells = ttest(...
        squeeze(mean(firstTCAV(respwin,:,:),1)),...
        squeeze(mean(firstTCAV(basewin_0,:,:),1)),...
        'dim',2,'tail','right','alpha',cellGroupsAlpha);
    lateRespCells = ttest(...
        squeeze(mean(longTCAV(lateWinFr,:,:),1)),...
        squeeze(mean(longTCAV(basewin,:,:),1)),...
        'dim',2,'tail','right','alpha',cellGroupsAlpha);
    lateSuppCells = ttest(...
        squeeze(mean(longTCAV(lateWinFr,:,:),1)),...
        squeeze(mean(longTCAV(basewin,:,:),1)),...
        'dim',2,'tail','left','alpha',cellGroupsAlpha);
    
    respCellsExpt(iexp).exptName = antiDataExpt(iexp).exptName;
    respCellsExpt(iexp).firstRespCells = firstRespCells;
    respCellsExpt(iexp).lateRespCells = lateRespCells;
    respCellsExpt(iexp).lateSuppCells = lateSuppCells;
    
    if shortCycExptInd(iexp)
        isShortCycExpt = cat(1,isShortCycExpt,true(length(firstRespCells),1));
    else
        isShortCycExpt = cat(1,isShortCycExpt,false(length(firstRespCells),1));
    end
        
end

cellInfo = struct;
cellInfo.firstRespCells = logical(cell2mat({respCellsExpt.firstRespCells}'));
cellInfo.lateRespCells = logical(cell2mat({respCellsExpt.lateRespCells}'));
cellInfo.lateSuppCells = logical(cell2mat({respCellsExpt.lateSuppCells}'));
cellInfo.isShortCycExpt = isShortCycExpt;
cellInfo.isTuned = logical(cell2mat({oriTuningExpt.isTuned}))';
cellInfo.oriResp = cell2mat({oriTuningExpt.oriResp}');
cellInfo.oriRespErr = cell2mat({oriTuningExpt.oriRespErr}');
cellInfo.oriFit = cell2mat({oriTuningExpt.fit})';

antiAnalysis = struct;
antiAnalysis.longTC = cell(1,3);
antiAnalysis.longTCErr = cell(1,3);
antiAnalysis.lateCycTC = cell(1,2);
antiAnalysis.lateCycSI = [];
for iexp = 1:nexp
    longTC_vis = antiDataExpt(iexp).longTC{visualTrials};
    longTC_aud = antiDataExpt(iexp).longTC{auditoryTrials};
    
    firstCycResp_vis = antiDataExpt(iexp).cycTC{visualTrials,1};
    firstCycResp_aud = antiDataExpt(iexp).cycTC{auditoryTrials,1};
    
    cycTC_vis = antiDataExpt(iexp).cycTC(visualTrials,:);
    cycTC_aud = antiDataExpt(iexp).cycTC(auditoryTrials,:);
    
    lateCycTC_vis = [];
    lateCycTC_aud = [];
    for icyc = 1:length(lateCycles)
        lateCycTC_vis = cat(3,lateCycTC_vis,...
            cycTC_vis{lateCycles(icyc)} - mean(cycTC_vis{lateCycles(icyc)}(basewin_0,:,:),1));
        lateCycTC_aud = cat(3,lateCycTC_aud,...
            cycTC_aud{lateCycles(icyc)} - mean(cycTC_aud{lateCycles(icyc)}(basewin_0,:,:),1));
    end
    
    lateCycSI = getSelectivityIndex(squeeze(mean(lateCycTC_vis(respwin,:,:),1))',...
        squeeze(mean(lateCycTC_aud(respwin,:,:),1))');
    
    antiAnalysis.longTC{visualTrials} = cat(2,antiAnalysis.longTC{visualTrials},...
        mean(longTC_vis,3));
    antiAnalysis.longTC{auditoryTrials} = cat(2,antiAnalysis.longTC{auditoryTrials},...
        mean(longTC_aud,3));
    antiAnalysis.longTC{allTrialsInd} = cat(2,antiAnalysis.longTC{allTrialsInd},...
        mean(cat(3,longTC_vis,longTC_aud),3));
    antiAnalysis.longTCErr{allTrialsInd} = cat(2,antiAnalysis.longTCErr{allTrialsInd},...
        ste(cat(3,longTC_vis,longTC_aud),3));
    antiAnalysis.lateCycTC{visualTrials} = cat(2,antiAnalysis.lateCycTC{visualTrials},...
        mean(lateCycTC_vis,3));
    antiAnalysis.lateCycTC{auditoryTrials} = cat(2,antiAnalysis.lateCycTC{auditoryTrials},...
        mean(lateCycTC_aud,3));
    
    antiAnalysis.lateCycSI = cat(2,antiAnalysis.lateCycSI,lateCycSI);
end

%% plotting params
respTCLim = [-0.005 0.05];
cycTCLim = [-0.005 0.01];
scatLim_win = [-0.2 0.6];
scatLim_cyc = [-0.035 0.085];
hmLim = [-0.1 0.1];
exCellTCLim = [-0.02 0.15];
oriRespLim = [-0.05 0.15];
siLim = [-10 10];

tcStartFrame = 26;
cycTCEndTimeMs = 350;
ttLabel_long = 0:500:2500;
ttLabel_cyc = -200:100:cycTCEndTimeMs;
nFr_long = size(antiAnalysis.longTC{1,1},1);
tt_longTC = ((tcStartFrame:nFr_long)-(nBaselineFr+nVisDelayFr)).*(1000/frameRateHz);
% ttFr_long = (tcStartFrame:nFr_long)-(nBaselineFr+nVisDelayFr);
ttLabelFr_long = ((ttLabel_long./1000)*frameRateHz)+...
    ((nBaselineFr+nVisDelayFr)-tcStartFrame+1);

nFr_cyc = size(antiAnalysis.lateCycTC{1,1},1);
tt_cycTC = ((tcStartFrame:nFr_cyc)-(nBaselineFr+nVisDelayFr)).*(1000/frameRateHz);

lateWinTT = ([lateWinFr(1) lateWinFr(end)] - (nBaselineFr+nVisDelayFr))...
    .*(1000/frameRateHz);
respWinTT = ([respwin(1) respwin(end)] - (nBaselineFr+nVisDelayFr))...
    .*(1000/frameRateHz);

%% plot anticipation analysis (Figure 2)

setFigParams4Print('landscape')

% heatmaps
lateWinRespAll = mean(antiAnalysis.longTC{allTrialsInd}(lateWinFr,:),1);
[~,lateWinSortInd] = sort(lateWinRespAll);
hm = flipud(antiAnalysis.longTC{allTrialsInd}(:,lateWinSortInd)');

figure
colormap(brewermap([],'*RdBu'));
% subplot 121
imagesc(hm(:,tcStartFrame:end))
hold on
if strcmp(ds,'FSAV_attentionV1')
    exCellInd = [exampleCell_1,exampleCell_2];
    exCellMat = zeros(1,length(cellInfo.firstRespCells));
    exCellMat(exCellInd) = 1;
    exCellSortInd = find(flip(exCellMat(lateWinSortInd)));
    hline(exCellSortInd,'k-')
end
figXAxis([],'Time from Start (ms)',[],ttLabelFr_long,ttLabelFr_long)
figYAxis([],'Cell #',[])
figAxForm
colorbar
caxis(hmLim)
title('All Trials, All Cells')
print([fnout 'heatmapAllTrialsAllCells'],'-dpdf','-fillpage')

figure
colormap gray
% subplot 122
ind = cellInfo.firstRespCells | cellInfo.lateRespCells;
indBW = ind; indBW(ind) = 0; indBW(~ind) = 1;
imagesc(flipud(indBW(lateWinSortInd)))
figYAxis([],'Cell #',[])
figAxForm
colorbar
title('Anti. Resp. Cells')
print([fnout 'cellIDforHeatmap'],'-dpdf','-fillpage')

% time-courses and quantification scatters
setFigParams4Print('portrait')
figure
subplot 311
for iav = 1:2
    y = mean(antiAnalysis.longTC{iav}...
        ((tcStartFrame:end),(cellInfo.firstRespCells | cellInfo.lateRespCells) & ...
        cellInfo.isShortCycExpt),2);
    yerr = ste(antiAnalysis.longTC{iav}...
        ((tcStartFrame:end),(cellInfo.firstRespCells | cellInfo.lateRespCells) & cellInfo.isShortCycExpt),2);
    hold on
    shadedErrorBar_chooseColor(tt_longTC,y,yerr,cueColor{iav});
end
figXAxis([],'Time (ms)',[tt_longTC(1) tt_longTC(end)],ttLabel_long,ttLabel_long)
figYAxis([],'dF/F',respTCLim)  
vline(lateWinTT,'k--')
hline(0,'k:')
figAxForm([],0)
title(sprintf('First Stim Responsive Cells (%s/%s)',...
    num2str(sum((cellInfo.firstRespCells | cellInfo.lateRespCells) & cellInfo.isShortCycExpt)),...
    num2str(sum(cellInfo.isShortCycExpt))))

subplot 323
x = mean(antiAnalysis.longTC{visualTrials}(lateWinFr,...
    (cellInfo.firstRespCells | cellInfo.lateRespCells)),1);
xerr = ste(x,2);
y = mean(antiAnalysis.longTC{auditoryTrials}(lateWinFr,...
    (cellInfo.firstRespCells | cellInfo.lateRespCells)),1);
yerr = ste(y,2);
plot(x,y,'.')
hold on
errorbar(mean(x),mean(y),yerr,yerr,xerr,xerr,'.')
plot(scatLim_win,scatLim_win,'k--')
plot(scatLim_win,[0 0],'k:')
plot([0 0],scatLim_win,'k:')
[~,p] = ttest(x,y);
figXAxis([],'Visual (dF/F)',scatLim_win)
figYAxis([],'Auditory (dF/F)',scatLim_win)
figAxForm
title(sprintf('Late Window, Resp. Cells (%s/%s), p = %s',...
    num2str(sum(cellInfo.firstRespCells | cellInfo.lateRespCells)),...
    num2str(length(cellInfo.firstRespCells)),...
    num2str(round(p,2,'significant'))))


subplot 325
for iav = 1:2
    y = mean(antiAnalysis.lateCycTC{iav}...
        ((tcStartFrame:end),(cellInfo.firstRespCells | cellInfo.lateRespCells)),2);
    yerr = ste(antiAnalysis.lateCycTC{iav}...
        ((tcStartFrame:end),(cellInfo.firstRespCells | cellInfo.lateRespCells)),2);
    hold on
    shadedErrorBar_chooseColor(tt_cycTC,y,yerr,cueColor{iav});
end
figXAxis([],'Time (ms)',[tt_cycTC(1) cycTCEndTimeMs],ttLabel_cyc,ttLabel_cyc)
figYAxis([],'dF/F',cycTCLim)  
vline(lateWinTT,'k--')
hline(0,'k:')
vline(respWinTT,'k--')
figAxForm
title(sprintf('First Stim Responsive Cells (%s/%s)',...
    num2str(sum(cellInfo.firstRespCells | cellInfo.lateRespCells)),...
    num2str(length(cellInfo.firstRespCells))))

subplot 326
x = mean(antiAnalysis.lateCycTC{visualTrials}(respwin,...
    (cellInfo.firstRespCells | cellInfo.lateRespCells)),1);
xerr = ste(x,2);
y = mean(antiAnalysis.lateCycTC{auditoryTrials}(respwin,...
    (cellInfo.firstRespCells | cellInfo.lateRespCells)),1);
yerr = ste(y,2);
plot(x,y,'.')
hold on
errorbar(mean(x),mean(y),yerr,yerr,xerr,xerr,'.')
plot(scatLim_cyc,scatLim_cyc,'k--')
plot(scatLim_cyc,[0 0],'k:')
plot([0 0],scatLim_cyc,'k:')
[~,p] = ttest(x,y);
figXAxis([],'Visual (dF/F)',scatLim_cyc)
figYAxis([],'Auditory (dF/F)',scatLim_cyc)
figAxForm
title(sprintf('Late Cycles, Resp. Cells (%s/%s), p = %s',...
    num2str(sum(cellInfo.firstRespCells | cellInfo.lateRespCells)),...
    num2str(length(cellInfo.firstRespCells)),...
    num2str(round(p,2,'significant'))))

print([fnout 'tcLongTrialsAndLateCycWithQuant_RespCells'],'-dpdf','-fillpage')

% example cells time-courses
if strcmp(ds,'FSAV_attentionV1')
    setFigParams4Print('landscape')

    figure
    subplot 221
    y = antiAnalysis.longTC{allTrialsInd}(tcStartFrame:end,exampleCell_1);
    yerr = antiAnalysis.longTCErr{allTrialsInd}(tcStartFrame:end,exampleCell_1);
    shadedErrorBar_chooseColor(tt_longTC,y,yerr,[0 0 0]);
    hold on
    figXAxis([],'Time (ms)',[tt_longTC(1) tt_longTC(end)],ttLabel_long,ttLabel_long)
    figYAxis([],'dF/F',exCellTCLim)  
    vline(lateWinTT,'k--')
    hline(0,'k:')
    figAxForm([],0)
    title(sprintf('Example Cell #%s',num2str(exampleCell_1)))

    subplot 222
    y = cellInfo.oriResp(exampleCell_1,:);
    yerr = cellInfo.oriRespErr(exampleCell_1,:);
    errorbar(orientations,y,yerr,'.')
    hold on
    x = 0:180;
    y = cellInfo.oriFit(exampleCell_1,:);
    plot(x,y,'-')
    figXAxis([],'Orienation (deg)',[-10 190])
    figYAxis([],'dF/F',oriRespLim)
    figAxForm
    title(sprintf('Passive Ori. Tuning, Cell #%s',num2str(exampleCell_1)));

    subplot 223
    y = antiAnalysis.longTC{allTrialsInd}(tcStartFrame:end,exampleCell_2);
    yerr = antiAnalysis.longTCErr{allTrialsInd}(tcStartFrame:end,exampleCell_2);
    shadedErrorBar_chooseColor(tt_longTC,y,yerr,[0 0 0]);
    hold on
    figXAxis([],'Time (ms)',[tt_longTC(1) tt_longTC(end)],ttLabel_long,ttLabel_long)
    figYAxis([],'dF/F',exCellTCLim)  
    vline(lateWinTT,'k--')
    hline(0,'k:')
    figAxForm([],0)
    title(sprintf('Example Cell #%s',num2str(exampleCell_2)))

    subplot 224
    y = cellInfo.oriResp(exampleCell_2,:);
    yerr = cellInfo.oriRespErr(exampleCell_2,:);
    errorbar(orientations,y,yerr,'.')
    hold on
    x = 0:180;
    y = cellInfo.oriFit(exampleCell_2,:);
    plot(x,y,'-')
    figXAxis([],'Orienation (deg)',[-10 190])
    figYAxis([],'dF/F',oriRespLim)
    figAxForm
    title(sprintf('Passive Ori. Tuning, Cell #%s',num2str(exampleCell_2)));

    print([fnout 'exampleCellsTCWithTuning'],'-dpdf','-fillpage')
end

% selectivity
figure
y = antiAnalysis.lateCycSI(cellInfo.firstRespCells | cellInfo.lateRespCells);
h = cdfplot(y);
hold on;
vline(mean(y),'k-')
figXAxis([],'Selectivity Index',siLim)
figYAxis([],'Fraction of Cells',[0 1])
figAxForm
title(sprintf('Resp. Cells, mean = %s, ste = %s',...
    num2str(round(mean(y),2,'significant')),...
    num2str(round(ste(y,2),2,'significant'))))

print([fnout 'selectivityCDF'],'-dpdf')