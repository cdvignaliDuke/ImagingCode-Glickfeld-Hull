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

if strcmp(ds,'FSAV_attentionV1')
    load(fullfile(rc.ashley, 'Manuscripts','Attention V1','Matlab Figs',...
        'V1_100ms_naive_anticipation_imgStats'))
    imgStats_naive = imgStats;
    clear imgStats
end

%%

nav = 2;
nBaselineFr = mouse(1).expt(1).info.preAlignFrames;
nFrames1s = frameRateHz;
nexp = size(expt,2);
nCycles = 8;
lateCycles = 5:nCycles;
lateWinFr = (45:88)+nBaselineFr;

oriBinSize = 45;
orientations = 0:oriBinSize:(180-oriBinSize);
oriBinEdges = [0, (oriBinSize/2):oriBinSize:(180-(oriBinSize/2)), 180];
nOri = length(orientations);
%% pool experiment data
antiDataExpt = struct;
oriTuningExpt = struct;
targetDataExpt = struct;
nTargets = 2; %sum(unique(targetInd) > 1);
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
            
            dd = d.av(iav).align(alignTarget);
            if iav == 1
                binEdges = oriBins;
                targets = dd.ori;
            elseif iav == 2
                binEdges = ampBins;
                targets = dd.amp;
            end
            tc = dd.respTC(1:(nBaselineFr*2),:,:);
            targetInd = discretize(targets,binEdges);
            targetDataExpt(exptN).av(iav).tc = cell(1,2);
            for itar = 1:nTargets
                ind = targetInd == itar+1;
                targetDataExpt(exptN).av(iav).tc{itar} = tc(:,:,ind);
                targetDataExpt(exptN).av(iav).targets{itar} = targets(ind);
            end
        end
        antiDataExpt(exptN).exptName = [mouse(imouse).mouse_name '-' d.date];
        antiDataExpt(exptN).exptCycLengthFr = cycLengthFr;
        antiDataExpt(exptN).longTC = longTC;
        antiDataExpt(exptN).cycTC = cycTC;
        targetDataExpt(exptN).exptName = [mouse(imouse).mouse_name '-' d.date];
        
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
        [~,oriTuningExpt(exptN).fitPeak] = max(do.oriFit,[],1);
        [~,oriPref] = histc(oriTuningExpt(exptN).fitPeak,oriBinEdges);
        oriPref(oriPref == length(orientations)+1 | oriPref == length(oriBinEdges) == 1) = 1;
        oriTuningExpt(exptN).oriPref = oriPref;
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
    
    eaTarRespCells = sum(cell2mat(cellfun(@(x) ...
        ttest(squeeze(mean(x(respwin_target,:,:),1)),...
        squeeze(mean(x(basewin_0_target,:,:),1)),...
        'dim',2,'tail','right','alpha',cellGroupsAlpha),...
        targetDataExpt(iexp).av(visualTrials).tc,'unif',0)),2) > 0;
    
    allTargetTC = [];
    for itar = 1:length(targetDataExpt(iexp).av(visualTrials).tc)
        allTargetTC = cat(3,allTargetTC,...
            targetDataExpt(iexp).av(visualTrials).tc{itar});
    end
    allTarRespCells = ttest(squeeze(mean(allTargetTC(respwin_target,:,:),1)),...
        squeeze(mean(allTargetTC(basewin_0_target,:,:),1)),...
        'dim',2,'tail','right','alpha',cellGroupsAlpha);
    
    respCellsExpt(iexp).exptName = antiDataExpt(iexp).exptName;
    respCellsExpt(iexp).firstRespCells = firstRespCells;
    respCellsExpt(iexp).lateRespCells = lateRespCells;
    respCellsExpt(iexp).lateSuppCells = lateSuppCells;
    respCellsExpt(iexp).targetRespCells = eaTarRespCells | allTarRespCells;
    
    if shortCycExptInd(iexp)
        isShortCycExpt = cat(1,isShortCycExpt,true(length(firstRespCells),1));
    else
        isShortCycExpt = cat(1,isShortCycExpt,false(length(firstRespCells),1));
    end
        
end

antiAnalysis = struct;
antiAnalysis.longTC = cell(1,3);
antiAnalysis.longTCErr = cell(1,3);
antiAnalysis.lateCycTC = cell(1,3);
antiAnalysis.firstCycTC = cell(1,3);
antiAnalysis.lateCycSI = [];
targetAnalysis = struct;
targetAnalysis.tc = cell(3,2);
targetAnalysis.targets = cell(1,2);
targetAnalysis.targetsAuROC = [];
targetAnalysis.targetsAuROCTest = [];
firstStimAuROC = [];
lateStimAuROC = [];
firstStimAuROCTest = [];
lateStimAuROCTest = [];
taskTuningTest = [];
taskTuningPref = [];
for iexp = 1:nexp
    longTC_vis = antiDataExpt(iexp).longTC{visualTrials};
    longTC_aud = antiDataExpt(iexp).longTC{auditoryTrials};
    
    nCells = size(longTC_vis,2);
    
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
    
    firstCycTC_vis = antiDataExpt(iexp).cycTC{visualTrials,1};
    firstCycTC_aud = antiDataExpt(iexp).cycTC{auditoryTrials,1};
    
    firstVisResp = squeeze(mean(firstCycTC_vis(respwin,:,:),1));
    lateVisResp = squeeze(mean(lateCycTC_vis(respwin,:,:),1));
    
    tarVisResp = cellfun(@(x) ...
        squeeze(mean(x(respwin_target,:,:),1) - mean(x(basewin_0_target,:,:),1)),...
        cat(2,targetDataExpt(iexp).av(visualTrials).tc,...
        {cat(3,targetDataExpt(iexp).av(visualTrials).tc{1},...
        targetDataExpt(iexp).av(visualTrials).tc{2})}),'unif',0);
    auroc_first = nan(nCells,3);
    auroc_first_test = nan(nCells,3);
    auroc_late = nan(nCells,3);
    auroc_late_test = nan(nCells,3);
    auroc_target = nan(nCells,1);
    auroc_target_test = nan(nCells,1);
    for icell = 1:nCells
        auroc_first(icell,:) = cellfun(@(x) ...
            roc_gh(firstVisResp(icell,:),x(icell,:)),tarVisResp);
        auroc_first_test(icell,:) = cellfun(@(x) ...
            ranksum(firstVisResp(icell,:),x(icell,:)),tarVisResp);
        auroc_late(icell,:) = cellfun(@(x) ...
            roc_gh(lateVisResp(icell,:),x(icell,:)),tarVisResp);
        auroc_late_test(icell,:) = cellfun(@(x) ...
            ranksum(lateVisResp(icell,:),x(icell,:)),tarVisResp);
        auroc_target(icell) = roc_gh(tarVisResp{1}(icell,:),tarVisResp{2}(icell,:));
        auroc_target_test(icell) = ranksum(...
            tarVisResp{1}(icell,:),tarVisResp{2}(icell,:)) < 0.05;
    end
    
    stimResp4anova = cat(2,firstVisResp,tarVisResp{1},tarVisResp{2});
    stimID4anova = cat(2,ones(1,size(firstVisResp,2)),...
        ones(1,size(tarVisResp{1},2))*2,ones(1,size(tarVisResp{2},2))*3);
    taskStimAnova = nan(1,nCells);
    for icell = 1:nCells
        taskStimAnova(icell) = anova1(stimResp4anova(icell,:),stimID4anova,'off');
    end
    [~,taskTuningID] = max(cat(2,mean(firstVisResp,2),...
        mean(tarVisResp{1},2),mean(tarVisResp{2},2)),[],2);
    
    
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
    antiAnalysis.lateCycTC{allTrialsInd} = cat(2,antiAnalysis.lateCycTC{allTrialsInd},...
        mean(cat(3,lateCycTC_vis,lateCycTC_aud),3));
    antiAnalysis.lateCycSI = cat(2,antiAnalysis.lateCycSI,lateCycSI);
    antiAnalysis.firstCycTC{visualTrials} = cat(2,antiAnalysis.firstCycTC{visualTrials},...
        mean(firstCycTC_vis,3));
    antiAnalysis.firstCycTC{auditoryTrials} = cat(2,antiAnalysis.firstCycTC{auditoryTrials},...
        mean(firstCycTC_aud,3));
    antiAnalysis.firstCycTC{allTrialsInd} = cat(2,antiAnalysis.firstCycTC{allTrialsInd},...
        mean(cat(3,firstCycTC_vis,firstCycTC_aud),3));
    
    targetAnalysis.tc{1,visualTrials} = cat(2,targetAnalysis.tc{1,visualTrials},...
        mean(targetDataExpt(iexp).av(visualTrials).tc{1},3));
    targetAnalysis.tc{2,visualTrials} = cat(2,targetAnalysis.tc{2,visualTrials},...
        mean(targetDataExpt(iexp).av(visualTrials).tc{2},3));
    targetAnalysis.tc{3,visualTrials} = cat(2,targetAnalysis.tc{3,visualTrials},...
        mean(cat(3,targetDataExpt(iexp).av(visualTrials).tc{1},...
        targetDataExpt(iexp).av(visualTrials).tc{2}),3));
    targetAnalysis.tc{1,auditoryTrials} = cat(2,targetAnalysis.tc{1,auditoryTrials},...
        mean(targetDataExpt(iexp).av(auditoryTrials).tc{1},3));
    targetAnalysis.tc{2,auditoryTrials} = cat(2,targetAnalysis.tc{2,auditoryTrials},...
        mean(targetDataExpt(iexp).av(auditoryTrials).tc{2},3));
    targetAnalysis.tc{3,auditoryTrials} = cat(2,targetAnalysis.tc{3,auditoryTrials},...
        mean(cat(3,targetDataExpt(iexp).av(auditoryTrials).tc{1},...
        targetDataExpt(iexp).av(auditoryTrials).tc{2}),3));
    
    targetAnalysis.targets{visualTrials} = cat(2,targetAnalysis.targets{visualTrials},...
        cat(1,ones(1,nCells).*mean(targetDataExpt(iexp).av(visualTrials).targets{1}),...
        ones(1,nCells).*mean(targetDataExpt(iexp).av(visualTrials).targets{2})));
    targetAnalysis.targetsAuROC = cat(1,targetAnalysis.targetsAuROC,auroc_target);
    targetAnalysis.targetsAuROCTest = cat(1,targetAnalysis.targetsAuROCTest,auroc_target_test);
    firstStimAuROC = cat(1,firstStimAuROC,auroc_first);
    lateStimAuROC = cat(1,lateStimAuROC,auroc_late);
    firstStimAuROCTest = cat(1,firstStimAuROCTest,auroc_first_test);
    lateStimAuROCTest = cat(1,lateStimAuROCTest,auroc_late_test);
    taskTuningTest = cat(1,taskTuningTest,taskStimAnova' < 0.05);
    taskTuningPref = cat(1,taskTuningPref,taskTuningID);
end
% antiAnalysis.adapt = cellfun(@(x,y) ...
%     mean(x(respwin,:),1)./mean(y(respwin,:),1),...
%     antiAnalysis.lateCycTC,antiAnalysis.firstCycTC,'unif',0);

cellInfo = struct;
cellInfo.firstRespCells = logical(cell2mat({respCellsExpt.firstRespCells}'));
cellInfo.lateRespCells = logical(cell2mat({respCellsExpt.lateRespCells}'));
cellInfo.lateSuppCells = logical(cell2mat({respCellsExpt.lateSuppCells}'));
cellInfo.minRespCells = (mean(antiAnalysis.firstCycTC{visualTrials}(respwin,:),1) > ...
    minRespThreshold & mean(antiAnalysis.firstCycTC{auditoryTrials}(respwin,:),1) > ...
    minRespThreshold)';
cellInfo.targetRespCells = logical(cell2mat({respCellsExpt.targetRespCells}'));
cellInfo.isShortCycExpt = isShortCycExpt;
cellInfo.isTuned = logical(cell2mat({oriTuningExpt.isTuned}))';
cellInfo.oriResp = cell2mat({oriTuningExpt.oriResp}');
cellInfo.oriRespErr = cell2mat({oriTuningExpt.oriRespErr}');
cellInfo.oriFit = cell2mat({oriTuningExpt.fit})';
cellInfo.oriPref = cell2mat({oriTuningExpt.oriPref})';
cellInfo.taskTuningPref = taskTuningPref;
cellInfo.taskTuningTest = taskTuningTest;

auroc_first = nan(length(cellInfo.firstRespCells),1);
auroc_late = nan(length(cellInfo.firstRespCells),1);
auroc_first_test = nan(length(cellInfo.firstRespCells),1);
auroc_late_test = nan(length(cellInfo.firstRespCells),1);
for icell = 1:length(cellInfo.firstRespCells)
    if ~cellInfo.targetRespCells(icell)
        auroc_first(icell) = firstStimAuROC(icell,3);
        auroc_late(icell) = lateStimAuROC(icell,3);
        auroc_first_test(icell) = firstStimAuROCTest(icell,3) < 0.05;
        auroc_late_test(icell) = lateStimAuROCTest(icell,3) < 0.05;
    elseif ~targetAnalysis.targetsAuROCTest(icell)
        auroc_first(icell) = firstStimAuROC(icell,3);
        auroc_late(icell) = lateStimAuROC(icell,3);
        auroc_first_test(icell) = firstStimAuROCTest(icell,3) < 0.05;
        auroc_late_test(icell) = lateStimAuROCTest(icell,3) < 0.05;
    elseif targetAnalysis.targetsAuROCTest(icell) && targetAnalysis.targetsAuROC(icell) < 0.5
        auroc_first(icell) = firstStimAuROC(icell,1);
        auroc_late(icell) = lateStimAuROC(icell,1);
        auroc_first_test(icell) = firstStimAuROCTest(icell,1) < 0.05;
        auroc_late_test(icell) = lateStimAuROCTest(icell,1) < 0.05;
    elseif targetAnalysis.targetsAuROCTest(icell) && targetAnalysis.targetsAuROC(icell) >= 0.5
        auroc_first(icell) = firstStimAuROC(icell,2);
        auroc_late(icell) = lateStimAuROC(icell,2);
        auroc_first_test(icell) = firstStimAuROCTest(icell,2) < 0.05;
        auroc_late_test(icell) = lateStimAuROCTest(icell,2) < 0.05;
    end
end
cellInfo.firstStimAuROC = auroc_first;
cellInfo.lateStimAuROC = auroc_late;
cellInfo.firstStimAuROCTest = auroc_first_test;
cellInfo.lateStimAuROCTest = auroc_late_test;
%% plotting params
respTCLim = [-0.005 0.05];
cycTCLim = [-0.005 0.01];
cycTCLim_minRespCells = [-0.005 0.025];
scatLim_win = [-0.2 0.6];
scatLim_cyc = [-0.035 0.085];
hmLim = [-0.1 0.1];
exCellTCLim = [-0.02 0.15];
oriRespLim = [-0.05 0.15];
siLim = [-10 10];
siOriLim = [-3 3];
oriBarLim_win = [0 0.08];
oriBarLim_resp = [0 0.04];
oriLim_taskResp = [-0.005 0.035];
oriNLim = [0 120];
oriTCLim = [-0.005 0.02];
targetTCLim = [-0.015 0.08];
% adaptLim = [0 1];

tcStartFrame = 26;
cycTCEndTimeMs = 350;
ttLabel_long = 0:500:2500;
ttLabel_cyc = -200:100:cycTCEndTimeMs;
ttLabel_target = -1000:250:900;
preTargetStimLabel = -700:350:0;
nFr_long = size(antiAnalysis.longTC{1,1},1);
tt_longTC = ((tcStartFrame:nFr_long)-(nBaselineFr+nVisDelayFr)).*(1000/frameRateHz);
ttLabelFr_long = ((ttLabel_long./1000)*frameRateHz)+...
    ((nBaselineFr+nVisDelayFr)-tcStartFrame+1);
ttLabelFr_target = ((ttLabel_target./1000)*frameRateHz)+...
    ((nBaselineFr+nVisDelayFr_target)+1);

nFr_cyc = size(antiAnalysis.lateCycTC{1,1},1);
tt_cycTC = ((tcStartFrame:nFr_cyc)-(nBaselineFr+nVisDelayFr)).*(1000/frameRateHz);
tt_targetTC = ((1:nFr_cyc)-(nBaselineFr+nVisDelayFr_target)).*(1000/frameRateHz);

lateWinTT = ([lateWinFr(1) lateWinFr(end)] - (nBaselineFr+nVisDelayFr))...
    .*(1000/frameRateHz);
respWinTT = ([respwin(1) respwin(end)] - (nBaselineFr+nVisDelayFr))...
    .*(1000/frameRateHz);
respWinTT_target = (...
    [respwin_target(1) respwin_target(end)] - (nBaselineFr+nVisDelayFr_target))...
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

%% create a structure with stats from analyses
imgStats = struct;
imgStats.lateWinTimeMs = round(lateWinTT,4,'significant');

imgStats.nCells.totalCells = length(cellInfo.firstRespCells);
imgStats.nCells.firstResp = sum(cellInfo.firstRespCells);
imgStats.nCells.lateOnlyResp = sum(...
    ~cellInfo.firstRespCells & cellInfo.lateRespCells);
imgStats.nCells.lateSupp = sum(cellInfo.lateSuppCells);
imgStats.nCells

x = mean(antiAnalysis.longTC{visualTrials}(lateWinFr,...
    (cellInfo.firstRespCells | cellInfo.lateRespCells)),1);
y = mean(antiAnalysis.longTC{auditoryTrials}(lateWinFr,...
    (cellInfo.firstRespCells | cellInfo.lateRespCells)),1);
imgStats.av(visualTrials).lateWinResp = mean(x);
imgStats.av(visualTrials).lateWinRespErr = ste(x,2);
imgStats.av(auditoryTrials).lateWinResp = mean(y);
imgStats.av(auditoryTrials).lateWinRespErr = ste(y,2);
[~,imgStats.lateWinTest] = ttest(x,y);

fprintf('Visual - Mean/Err Late Win: %s/%s\n',...
    num2str(round(imgStats.av(visualTrials).lateWinResp,2,'significant')),...
    num2str(round(imgStats.av(visualTrials).lateWinRespErr,2,'significant')))
fprintf('Auditory - Mean/Err Late Win: %s/%s\n',...
    num2str(round(imgStats.av(auditoryTrials).lateWinResp,2,'significant')),...
    num2str(round(imgStats.av(auditoryTrials).lateWinRespErr,2,'significant')))

x = mean(antiAnalysis.lateCycTC{visualTrials}(respwin,...
    (cellInfo.firstRespCells | cellInfo.lateRespCells)),1);
y = mean(antiAnalysis.lateCycTC{auditoryTrials}(respwin,...
    (cellInfo.firstRespCells | cellInfo.lateRespCells)),1);
imgStats.av(visualTrials).lateCycResp = mean(x);
imgStats.av(visualTrials).lateCycRespErr = ste(x,2);
imgStats.av(auditoryTrials).lateCycResp = mean(y);
imgStats.av(auditoryTrials).lateCycRespErr = ste(y,2);
[~,imgStats.lateCycTest] = ttest(x,y);

fprintf('Visual - Mean/Err Late Cyc: %s/%s\n',...
    num2str(round(imgStats.av(visualTrials).lateCycResp,2,'significant')),...
    num2str(round(imgStats.av(visualTrials).lateCycRespErr,2,'significant')))
fprintf('Auditory - Mean/Err Late Cyc: %s/%s\n',...
    num2str(round(imgStats.av(auditoryTrials).lateCycResp,2,'significant')),...
    num2str(round(imgStats.av(auditoryTrials).lateCycRespErr,2,'significant')))

imgStats.allRespCellsSI = antiAnalysis.lateCycSI(...
    cellInfo.firstRespCells | cellInfo.lateRespCells);

fprintf('Mean/Err SI: %s/%s\n',...
    num2str(round(mean(imgStats.allRespCellsSI,2),2,'significant')),...
    num2str(round(ste(imgStats.allRespCellsSI,2),2,'significant')))

if strcmp(ds,'FSAV_attentionV1')
    [~,imgStats.respCellsSITestVsNaive] = kstest2(imgStats.allRespCellsSI, ...
        imgStats_naive.allRespCellsSI,'tail','smaller');
    fprintf('Kolmogorov-Smirnov Test SI, behavior vs. naive, p = %s\n',...
        num2str(round(imgStats.respCellsSITestVsNaive,2,'significant')))
end

save([fnout 'imgStats'],'imgStats')

%% plot orientation analysis (Figure 3)
% adaptAnalysis = struct;
% ind = cellInfo.firstRespCells & cellInfo.minRespCells;
% adaptAnalysis.firstRespCells = cellfun(@(x) mean(x(ind)),antiAnalysis.adapt);
% adaptAnalysis.firstRespCellsErr = cellfun(@(x) ste(x(ind),2),antiAnalysis.adapt);
% [~,adaptAnalysis.testAV] = ttest(antiAnalysis.adapt{visualTrials}(ind),...
%     antiAnalysis.adapt{auditoryTrials}(ind),'tail','right');
% ind2 = ind & cellInfo.isTuned;
% adaptAnalysis.oriGroupsTest = anova1(antiAnalysis.adapt{3}(ind2),...
%     cellInfo.oriPref(ind2),'off');

% get first, late win, and late cycle response for passive ori tuning groups
oriGroups = struct;
oriGroups.n = nan(1,nOri);
oriGroups.nShortCyc = nan(1,nOri);
oriGroups.nTarOrDist = nan(1,nOri);
oriGroups.nTar = nan(1,nOri);
oriGroups.lateCycSI = nan(1,nOri);
oriGroups.lateCycSIErr = nan(1,nOri);
oriGroups.firstResp = nan(2,nOri);
oriGroups.firstRespErr = nan(2,nOri);
oriGroups.lateWin = nan(2,nOri);
oriGroups.lateWinErr = nan(2,nOri);
oriGroups.lateCycResp = nan(2,nOri);
oriGroups.lateCycRespErr = nan(2,nOri);
oriGroups.longTC = cell(2,nOri);
oriGroups.longTCErr = cell(2,nOri);
oriGroups.cycTC = cell(2,nOri);
oriGroups.cycTCErr = cell(2,nOri);
oriGroups.lateCycTest = nan(1,nOri);
oriGroups.lateWinTest = nan(1,nOri);
oriGroups.targetTC = cell(1,nOri);
oriGroups.targetTCErr = cell(1,nOri);
oriGroups.targetTuningResp = cell(1,nOri);
oriGroups.targetTuningRespErr = cell(1,nOri);
oriGroups.targetTuningStim = cell(1,nOri);
oriGroups.targetTuningStimErr = cell(1,nOri);
oriGroups.firstStimRespForTargetAnalysis = nan(1,nOri);
oriGroups.firstStimRespErrForTargetAnalysis = nan(1,nOri);
% adaptAnalysis.oriGroupsN = nan(1,nOri);
% adaptAnalysis.oriGroups = nan(3,nOri);
% adaptAnalysis.oriGroupsErr = nan(3,nOri);
for iori = 1:nOri
    ind = cellInfo.isTuned & cellInfo.oriPref == iori &...
        (cellInfo.firstRespCells | cellInfo.lateRespCells);
%     ind2 = cellInfo.isTuned & cellInfo.oriPref == iori & ...
%         cellInfo.firstRespCells & cellInfo.minRespCells;
    ind_tarAndDist = cellInfo.isTuned & cellInfo.oriPref == iori &...
        (cellInfo.firstRespCells | cellInfo.lateRespCells | cellInfo.targetRespCells);
    ind_tar = cellInfo.isTuned & cellInfo.oriPref == iori & cellInfo.targetRespCells;
    oriGroups.n(iori) = sum(ind);
    oriGroups.nShortCyc(iori) = sum(ind & cellInfo.isShortCycExpt);
    oriGroups.nTarOrDist(iori) = sum(ind_tarAndDist);
    oriGroups.nTar(iori) = sum(ind_tar);
%     adaptAnalysis.oriGroupsN(iori) = sum(ind2);
    oriGroups.lateCycSI(iori) = mean(antiAnalysis.lateCycSI(ind));
    oriGroups.lateCycSIErr(iori) = ste(antiAnalysis.lateCycSI(ind),2);
    [~,oriGroups.lateWinTest(iori)] = ttest(mean(antiAnalysis.longTC{visualTrials}...
        (lateWinFr,ind & cellInfo.isShortCycExpt),1),...
        mean(antiAnalysis.longTC{auditoryTrials}...
        (lateWinFr,ind & cellInfo.isShortCycExpt),1));
    [~,oriGroups.lateCycTest(iori)] = ttest(mean(antiAnalysis.lateCycTC{visualTrials}...
        (respwin,ind),1),mean(antiAnalysis.lateCycTC{auditoryTrials}...
        (respwin,ind),1));
    
    oriGroups.targetTC{iori} = cell2mat(cellfun(@(x) mean(x(:,ind_tar),2),...
        targetAnalysis.tc(1:2,visualTrials),'unif',0)');
    oriGroups.targetTCErr{iori} = cell2mat(cellfun(@(x) ste(x(:,ind_tar),2),...
        targetAnalysis.tc(1:2,visualTrials),'unif',0)');
    
    oriGroups.targetTuningResp{iori} = cat(1,mean(mean(...
            antiAnalysis.lateCycTC{visualTrials}(respwin_target,ind_tarAndDist),1)),cellfun(@(x) ...
        mean(mean(x(respwin_target,ind_tarAndDist),1) - mean(x(basewin_0,ind_tarAndDist),1)),...
        targetAnalysis.tc(1:2,visualTrials)));
    oriGroups.targetTuningRespErr{iori} = cat(1,ste(mean(...
            antiAnalysis.lateCycTC{visualTrials}(respwin_target,ind_tarAndDist),1),2),cellfun(@(x) ...
        ste(mean(x(respwin_target,ind_tarAndDist),1) - mean(x(basewin_0,ind_tarAndDist),1),2),...
        targetAnalysis.tc(1:2,visualTrials)));
    oriGroups.firstStimRespForTargetAnalysis(iori) = mean(mean(...
            antiAnalysis.firstCycTC{visualTrials}(respwin,ind_tarAndDist),1));
    oriGroups.firstStimRespErrForTargetAnalysis(iori) = ste(mean(...
            antiAnalysis.firstCycTC{visualTrials}(respwin,ind_tarAndDist),1),2);
    oriGroups.targetTuningStim{iori} = cat(1,0,...
        mean(targetAnalysis.targets{visualTrials}(:,ind_tarAndDist),2));
    oriGroups.targetTuningStimErr{iori} = cat(1,0,...
        ste(targetAnalysis.targets{visualTrials}(:,ind_tarAndDist),2));
        
    for iav = 1:2
        oriGroups.firstResp(iav,iori) = mean(mean(...
            antiAnalysis.firstCycTC{iav}(respwin,ind),1));
        oriGroups.firstRespErr(iav,iori) = ste(mean(...
            antiAnalysis.firstCycTC{iav}(respwin,ind),1),2);
        
        oriGroups.lateWin(iav,iori) = mean(mean(...
            antiAnalysis.longTC{iav}(lateWinFr,ind),1));
        oriGroups.lateWinErr(iav,iori) = ste(mean(...
            antiAnalysis.longTC{iav}(lateWinFr,ind),1),2);
        
        oriGroups.lateCycResp(iav,iori) = mean(mean(...
            antiAnalysis.lateCycTC{iav}(respwin,ind),1));
        oriGroups.lateCycRespErr(iav,iori) = ste(mean(...
            antiAnalysis.lateCycTC{iav}(respwin,ind),1),2);
        
        oriGroups.longTC{iav,iori} = mean(...
            antiAnalysis.longTC{iav}(:,ind & cellInfo.isShortCycExpt),2);
        oriGroups.longTCErr{iav,iori} = ste(...
            antiAnalysis.longTC{iav}(:,ind & cellInfo.isShortCycExpt),2);
        
        oriGroups.cycTC{iav,iori} = mean(...
            antiAnalysis.lateCycTC{iav}(:,ind),2);
        oriGroups.cycTCErr{iav,iori} = ste(...
            antiAnalysis.lateCycTC{iav}(:,ind),2);
%         adaptAnalysis.oriGroups(iav,iori) = mean(...
%             antiAnalysis.adapt{iav}(ind2),2);
%         adaptAnalysis.oriGroupsErr(iav,iori) = ste(...
%             antiAnalysis.adapt{iav}(ind2),2);
%         if iav == 2
%             adaptAnalysis.oriGroups(3,iori) = mean(...
%                 antiAnalysis.adapt{3}(ind2),2);
%             adaptAnalysis.oriGroupsErr(3,iori) = ste(...
%                 antiAnalysis.adapt{3}(ind2),2);
%         end
        
        
    end
end
ind = cellInfo.isTuned &...
    (cellInfo.firstRespCells | cellInfo.lateRespCells);
lateCycRespDiffData = ...
    mean(antiAnalysis.lateCycTC{visualTrials}(respwin,ind),1) - ...
    mean(antiAnalysis.lateCycTC{auditoryTrials}(respwin,ind),1);
oriGroups.lateCycRespTest = anova1(lateCycRespDiffData,cellInfo.oriPref(ind),'off');
oriGroups.lateCycSITest = anova1(antiAnalysis.lateCycSI(ind),cellInfo.oriPref(ind),'off');


setFigParams4Print('portrait')
figure;
suptitle('Tuned, First or Late Distractor Responsive Neurons, Short Cyc Expts')
for iori = 1:nOri
%     subplot(4,2,iori)
%     for iav = 1:2
%         y = oriGroups.longTC{iav,iori}(tcStartFrame:end);
%         yerr = oriGroups.longTCErr{iav,iori}(tcStartFrame:end);
%         hold on
%         shadedErrorBar_chooseColor(tt_longTC,y,yerr,cueColor{iav});
%     end
%     figXAxis([],'Time (ms)',[tt_longTC(1) tt_longTC(end)],ttLabel_long,ttLabel_long)
%     figYAxis([],'dF/F',oriTCLim)  
%     vline(lateWinTT,'k--')
%     hline(0,'k:')
%     figAxForm([],0)
%     title(sprintf('Pref %s, n=%s, p=%s',num2str(orientations(iori)),...
%         num2str(oriGroups.nShortCyc(iori)),num2str(round(oriGroups.lateWinTest(iori),2,'significant'))))
    
    subplot(4,2,iori)
    for iav = 1:2
        y = oriGroups.cycTC{iav,iori}(tcStartFrame:end);
        yerr = oriGroups.cycTCErr{iav,iori}(tcStartFrame:end);
        hold on
        shadedErrorBar_chooseColor(tt_cycTC,y,yerr,cueColor{iav});
    end
    figXAxis([],'Time (ms)',[tt_cycTC(1) cycTCEndTimeMs],ttLabel_cyc,ttLabel_cyc)
    figYAxis([],'dF/F',oriTCLim)  
    vline(respWinTT,'k--')
    hline(0,'k:')
    figAxForm
    title(sprintf('Pref %s, n=%s, p=%s',num2str(orientations(iori)),...
        num2str(oriGroups.n(iori)),num2str(round(oriGroups.lateCycTest(iori),2,'significant'))))
    
    subplot(4,2,iori+4)
    for itar = 1:2
        y = oriGroups.targetTC{iori}(tcStartFrame:end,itar);
        yerr = oriGroups.targetTCErr{iori}(tcStartFrame:end,itar);
        hold on
        shadedErrorBar_chooseColor(tt_cycTC,y,yerr,hiLoColor{itar});
    end
    figXAxis([],'Time (ms)',[tt_cycTC(1) cycTCEndTimeMs],ttLabel_cyc,ttLabel_cyc)
    figYAxis([],'dF/F',targetTCLim)  
    hline(0,'k:')
    vline(respWinTT,'k--')
    figAxForm
    title(sprintf('Pref %s, n=%s',num2str(orientations(iori)),...
        num2str(oriGroups.nTar(iori))))
    
end

print([fnout 'tuningAnalysis_oriGroupsTC'],'-dpdf','-fillpage')    

setFigParams4Print('portrait')
figure;
suptitle(sprintf('Tuned, Distractor-Responsive Neurons, %s/%s',...
    num2str(sum(cellInfo.isTuned & (cellInfo.firstRespCells | cellInfo.lateRespCells))),...
    num2str(length(cellInfo.isTuned))))
xsub = [-0.25 +0.25];
subplot 421
h = bar(oriGroups.firstResp','group','EdgeColor','none','BarWidth',1);
legend({'Visual','Auditory'})
figXAxis([],'Pref. Ori. (deg)',[0 nOri+1],[1:nOri],orientations)
figYAxis([],'dF/F',oriBarLim_resp)
figAxForm
title('First Stim')
subplot 422
for iav = 1:2
    hold on
    errorbar((1:nOri)+xsub(iav),oriGroups.firstResp(iav,:),...
        oriGroups.firstRespErr(iav,:),'.')
end
figXAxis([],'Pref. Ori. (deg)',[0 nOri+1],[1:nOri],orientations)
figYAxis([],'dF/F',oriBarLim_resp)
figAxForm
title('First Stim')

subplot 423
h = bar(oriGroups.lateWin','group','EdgeColor','none','BarWidth',1);
legend({'Visual','Auditory'})
figXAxis([],'Pref. Ori. (deg)',[0 nOri+1],[1:nOri],orientations)
figYAxis([],'dF/F',oriBarLim_win)
figAxForm
title('Late Window')
subplot 424
for iav = 1:2
    hold on
    errorbar((1:nOri)+xsub(iav),oriGroups.lateWin(iav,:),...
        oriGroups.lateWinErr(iav,:),'.')
end
figXAxis([],'Pref. Ori. (deg)',[0 nOri+1],[1:nOri],orientations)
figYAxis([],'dF/F',oriBarLim_win)
figAxForm
title('Late Window')

subplot 425
h = bar(oriGroups.lateCycResp','group','EdgeColor','none','BarWidth',1);
legend({'Visual','Auditory'})
figXAxis([],'Pref. Ori. (deg)',[0 nOri+1],[1:nOri],orientations)
figYAxis([],'dF/F',oriBarLim_resp)
figAxForm
title('Late Stim Resp')
subplot 426
for iav = 1:2
    hold on
    errorbar((1:nOri)+xsub(iav),oriGroups.lateCycResp(iav,:),...
        oriGroups.lateCycRespErr(iav,:),'.')
end
figXAxis([],'Pref. Ori. (deg)',[0 nOri+1],[1:nOri],orientations)
figYAxis([],'dF/F',oriBarLim_resp)
figAxForm
title(sprintf('Late Stim Resp, One-Way ANOVA p=%s',num2str(round(...
    oriGroups.lateCycRespTest,2,'significant'))))

subplot 427
h = bar(oriGroups.n,'EdgeColor','none','BarWidth',0.5);
hold on
for iori = 1:nOri
    text(iori,oriGroups.n(iori)+1,num2str(oriGroups.n(iori)))
end
figXAxis([],'Pref. Ori. (deg)',[0 nOri+1],[1:nOri],orientations)
figYAxis([],'dF/F',oriNLim)
figAxForm

subplot 427
h = bar(oriGroups.n,'EdgeColor','none','BarWidth',0.5);
hold on
for iori = 1:nOri
    text(iori,oriGroups.n(iori)+1,num2str(oriGroups.n(iori)))
end
figXAxis([],'Pref. Ori. (deg)',[0 nOri+1],[1:nOri],orientations)
figYAxis([],'N Cells',oriNLim)
figAxForm

subplot 428
h = bar(oriGroups.lateCycSI,'EdgeColor','none','BarWidth',0.5);
hold on
errorbar(1:nOri,oriGroups.lateCycSI,oriGroups.lateCycSIErr,'.')
figXAxis([],'Pref. Ori. (deg)',[0 nOri+1],[1:nOri],orientations)
figYAxis([],'Selectivity',siOriLim)
figAxForm
title(sprintf('Late Stim Resp, One-Way ANOVA p=%s',num2str(round(...
    oriGroups.lateCycSITest,2,'significant'))))

print([fnout 'tuningAnalysis_distRespCells'],'-dpdf','-fillpage')

% task tuning
setFigParams4Print('landscape')    
figure
suptitle(sprintf('Tuned and Distractor or Target Resp. Neurons %s/%s',...
    num2str(sum(cellInfo.isTuned & ...
    (cellInfo.firstRespCells | cellInfo.lateRespCells | cellInfo.targetRespCells))),...
    num2str(length(cellInfo.isTuned))))
for iori = 1:nOri
    subplot(2,2,iori)
    x = oriGroups.targetTuningStim{iori};
    xerr = oriGroups.targetTuningStimErr{iori};
    y = oriGroups.targetTuningResp{iori};
    yerr = oriGroups.targetTuningRespErr{iori};
    errorbar(0,oriGroups.firstStimRespForTargetAnalysis(iori),...
        oriGroups.firstStimRespErrForTargetAnalysis(iori),'.')
    hold on
    errorbar(x,y,yerr,yerr,xerr,xerr,'.')
    figXAxis([],'Task Ori. (deg)',[-10 100],x,oriBins([1,3,4]))
    figYAxis([],'dF/F',oriLim_taskResp)  
    vline(lateWinTT,'k--')
    hline(0,'k:')
    figAxForm([],0)
    title(sprintf('Pref %s, n=%s',num2str(orientations(iori)),...
        num2str(oriGroups.nTarOrDist(iori))))

end

print([fnout 'tuningAnalysis_oriGroupsTaskTuning'],'-dpdf','-fillpage')   

taskGroups = struct;
taskGroups.name = {'Dist.';'Tar';'n.d.'};
taskGroups.cmp(1).name = 'first:target';
taskGroups.cmp(2).name = 'late:target';
for icmp = 1:2
    taskGroups.cmp(icmp).n = nan(1,3);
    taskGroups.cmp(icmp).n4SI = nan(1,3);
    taskGroups.cmp(icmp).targetTC = cell(1,3);
    taskGroups.cmp(icmp).targetTCErr = cell(1,3);    
    taskGroups.cmp(icmp).firstTC = cell(1,3);
    taskGroups.cmp(icmp).firstTCErr = cell(1,3);
    taskGroups.cmp(icmp).lateCycSI = nan(1,3);
    taskGroups.cmp(icmp).lateCycSIErr = nan(1,3);
end
for igrp = 1:3
    if igrp == 3
        ind_first = ~cellInfo.firstStimAuROCTest & ...
            (cellInfo.firstRespCells | cellInfo.lateRespCells | cellInfo.targetRespCells);
        ind_late = ~cellInfo.lateStimAuROCTest & ...
            (cellInfo.firstRespCells | cellInfo.lateRespCells | cellInfo.targetRespCells);
    elseif igrp == 1
        ind_first = cellInfo.firstStimAuROCTest & cellInfo.firstStimAuROC < 0.5 & ...
            (cellInfo.firstRespCells | cellInfo.lateRespCells | cellInfo.targetRespCells);
        ind_late = cellInfo.lateStimAuROCTest & cellInfo.lateStimAuROC < 0.5 & ...
            (cellInfo.firstRespCells | cellInfo.lateRespCells | cellInfo.targetRespCells);
    elseif igrp == 2
        ind_first = cellInfo.firstStimAuROCTest & cellInfo.firstStimAuROC > 0.5 & ...
            (cellInfo.firstRespCells | cellInfo.lateRespCells | cellInfo.targetRespCells);
        ind_late = cellInfo.lateStimAuROCTest & cellInfo.lateStimAuROC > 0.5 & ...
            (cellInfo.firstRespCells | cellInfo.lateRespCells | cellInfo.targetRespCells);
    end
    taskGroups.cmp(1).n(igrp) = sum(ind_first);
    taskGroups.cmp(2).n(igrp) = sum(ind_late);
    taskGroups.cmp(1).n4SI(igrp) = sum(ind_first &...
        (cellInfo.firstRespCells | cellInfo.lateRespCells));
    taskGroups.cmp(2).n4SI(igrp) = sum(ind_late &...
        (cellInfo.firstRespCells | cellInfo.lateRespCells));
    taskGroups.cmp(1).targetTC{igrp} = mean(targetAnalysis.tc{3,visualTrials}...
        (:,ind_first),2);
    taskGroups.cmp(1).targetTCErr{igrp} = ste(targetAnalysis.tc{3,visualTrials}...
        (:,ind_first),2);
    taskGroups.cmp(2).targetTC{igrp} = mean(targetAnalysis.tc{3,visualTrials}...
        (:,ind_late),2);
    taskGroups.cmp(2).targetTCErr{igrp} = ste(targetAnalysis.tc{3,visualTrials}...
        (:,ind_late),2);
    taskGroups.cmp(1).firstTC{igrp} = mean(antiAnalysis.firstCycTC{visualTrials}...
        (:,ind_first),2);
    taskGroups.cmp(1).firstTCErr{igrp} = ste(antiAnalysis.firstCycTC{visualTrials}...
        (:,ind_first),2);
    taskGroups.cmp(2).firstTC{igrp} = mean(antiAnalysis.firstCycTC{visualTrials}...
        (:,ind_late),2);
    taskGroups.cmp(2).firstTCErr{igrp} = ste(antiAnalysis.firstCycTC{visualTrials}...
        (:,ind_late),2);
    taskGroups.cmp(1).lateCycSI(igrp) = mean(antiAnalysis.lateCycSI(ind_first &...
        (cellInfo.firstRespCells | cellInfo.lateRespCells)),2);
    taskGroups.cmp(1).lateCycSIErr(igrp) = ste(antiAnalysis.lateCycSI(ind_first &...
        (cellInfo.firstRespCells | cellInfo.lateRespCells)),2);
    taskGroups.cmp(2).lateCycSI(igrp) = mean(antiAnalysis.lateCycSI(ind_late &...
        (cellInfo.firstRespCells | cellInfo.lateRespCells)),2);
    taskGroups.cmp(2).lateCycSIErr(igrp) = ste(antiAnalysis.lateCycSI(ind_late &...
        (cellInfo.firstRespCells | cellInfo.lateRespCells)),2);
end
ind = (cellInfo.firstRespCells | cellInfo.lateRespCells | cellInfo.targetRespCells);
[~,sortInd] = sort(cellInfo.firstStimAuROC(ind));
tc = targetAnalysis.tc{3,visualTrials}(:,ind)';
taskGroups.cmp(1).hm = flipud(tc(sortInd,:));
[~,sortInd] = sort(cellInfo.lateStimAuROC(ind));
tc = targetAnalysis.tc{3,visualTrials}(:,ind)';
taskGroups.cmp(2).hm = flipud(tc(sortInd,:));

setFigParams4Print('portrait')
figure
suptitle('Target Aligned Resp., All Task Resp. Cells')
colormap(brewermap([],'*RdBu'))
for icmp = 1:2
    subplot(3,2,icmp)
    imagesc(taskGroups.cmp(1).hm)
    caxis(hmLim)
    colorbar
    figXAxis([],'Time (ms)',[],ttLabelFr_target,ttLabel_target)
    figYAxis([],'Cell # (auROC Sorted)',[])
    figAxForm
    title(taskGroups.cmp(icmp).name)
    subplot(3,2,icmp+2)
    L = [];
    for igrp = 1:3
        y = taskGroups.cmp(icmp).targetTC{igrp};
        yerr = taskGroups.cmp(icmp).targetTCErr{igrp};
        hold on
        h = shadedErrorBar_chooseColor(tt_targetTC,y,yerr,aurocColor{igrp});
        L(igrp) = h.mainLine;
    end
    figXAxis([],'Time (ms)',[tt_targetTC(1) tt_targetTC(end)],...
        ttLabel_target,ttLabel_target)
    figYAxis([],'dF/F',targetTCLim)  
    hline(0,'k:')
    vline(respWinTT,'k--')
    vline(preTargetStimLabel,'k:')
    figAxForm
    title(sprintf('%s, n=%s',taskGroups.cmp(icmp).name,...
        num2str(sum(taskGroups.cmp(icmp).n))))
    legend(L,strcat(taskGroups.name,repmat(' n=',[3,1]),...
        cellfun(@num2str,num2cell(taskGroups.cmp(icmp).n)','unif',0)))

    subplot(3,2,icmp+4)
    for igrp = 1:3
        y = taskGroups.cmp(icmp).firstTC{igrp}(tcStartFrame:end);
        yerr = taskGroups.cmp(icmp).firstTCErr{igrp}(tcStartFrame:end);
        hold on
        h = shadedErrorBar_chooseColor(tt_cycTC,y,yerr,aurocColor{igrp});
    end
    figXAxis([],'Time (ms)',[tt_cycTC(1) cycTCEndTimeMs],ttLabel_cyc,ttLabel_cyc)
    figYAxis([],'dF/F',targetTCLim)  
    hline(0,'k:')
    vline(respWinTT,'k--')
    figAxForm
    title(sprintf('%s, First Stim Resp.',taskGroups.cmp(icmp).name))
end
print([fnout 'tuningAnalysis_taskGroupsTC'],'-dpdf','-fillpage')  

setFigParams4Print('landscape')
figure
suptitle('Distractor Resp. Cells')
for icmp = 1:2
    subplot(1,2,icmp)
    bar(taskGroups.cmp(icmp).lateCycSI,'EdgeColor','none','BarWidth',0.5);
    hold on
    errorbar(1:3,taskGroups.cmp(icmp).lateCycSI,taskGroups.cmp(icmp).lateCycSIErr,'.')
    figXAxis([],'Pref. Task Stim.',[0 4],1:3,taskGroups.name)
    figYAxis([],'Selectivity',siOriLim)
    figAxForm
    for igrp = 1:3
        text(igrp,taskGroups.cmp(icmp).lateCycSI(igrp)+.2,num2str(taskGroups.cmp(icmp).n4SI(igrp)));
    end
    title(taskGroups.cmp(icmp).name)
end
print([fnout 'tuningAnalysis_taskGroupsSI'],'-dpdf')

