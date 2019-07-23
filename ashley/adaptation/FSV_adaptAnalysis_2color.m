clear all
close all
ds = 'FSAV_V1_SOM_temp';
cellsOrDendrites = 1;
doLoadPreviousAnalysis = false;
%%
rc = behavConstsAV;
imgParams_FSAV
bxParams_FSAV_attnV1ms

eval(ds)
titleStr = ds(6:end);
mice = unique({expt.SubNum});
mouse_str = ['i' strjoin(mice,'_i')];

load(fullfile(rc.caOutputDir,ds,...
    [mouse_str '_trOutcomeStruct_cells' ds(5:end) '.mat']));
if strcmp(ds,'FSAV_attentionV1') | strcmp(ds,'FSAV_V1_SOM_temp')
    fnout = fullfile(rc.ashleyAnalysis, 'Expt summaries','adaptation','behavior',...
        [titleStr '_']); 
    bxExpt = true;
elseif strcmp(ds,'FSAV_V1_SOM_naive_temp')
    fnout = fullfile(rc.ashleyAnalysis, 'Expt summaries','adaptation','naive',...
        [titleStr '_']);     
    bxExpt = false;
end


%%
nav = 2;
nBaselineFr = mouse(1).expt(1).info.preAlignFrames;
nFrames1s = frameRateHz;
nexp = size(expt,2);
nCycles = 7;
% lateCycles = 5:nCycles;
lateWinFr = (45:88)+nBaselineFr;
firstWinFr = (3:44)+nBaselineFr;
minTargetRT = (nVisDelayFr_target+respwin_target(1)-nBaselineFr)./frameRateHz*1000;

oriBinSize = 45;
orientations = 0:oriBinSize:(180-oriBinSize);
oriBinEdges = [0, (oriBinSize/2):oriBinSize:(180-(oriBinSize/2)), 180];
nOri = length(orientations);

nMovWin = 15;
movWinLabelFr = 30:(30+nMovWin-1);
movWinLabelFr_target = 30:(30+nMovWin-1);
% movWinLabelMs = 

% minCellN_SIFRmatch = 36;

trOutType = {'h';'m';'fa';'cr'};
trOutTypeName = {'H-All';'H-HT';'H-ET';'M-All';'M-HT';'M-ET';'FA';'CR'};
%% pool experiment data

    antiDataExpt = struct;
    antiDataExpt_passive = struct;
    tarDataExpt = struct;
    tarDataExpt_passive = struct;
%     oriTuningExpt = struct;
    for imouse = 1:size(mouse,2)
        for iexp = 1:size(mouse(imouse).expt,2)
            if imouse == 1 && iexp == 1
                exptN = 1;
            else
                exptN = exptN+1;
            end

            d = mouse(imouse).expt(iexp);
            exptID = strcmp({expt.SubNum},mouse(imouse).mouse_name) & strcmp({expt.date},d.date);

            cycLengthFr = d.info.cycTimeFrames;
            nCycLongTC = ceil(longTrialLengthFr./cycLengthFr);

            maxCycles = max(cat(2,d.tag(1).av(visualTrials).align(alignStart).nCycles,...
                d.tag(1).av(auditoryTrials).align(alignStart).nCycles));

            antiDataExpt(exptN).exptName = [mouse(imouse).mouse_name '-' d.date];
            antiDataExpt(exptN).exptCycLengthFr = cycLengthFr;
            tarDataExpt(exptN).exptName = [mouse(imouse).mouse_name '-' d.date];
            tarDataExpt(exptN).exptCycLengthFr = cycLengthFr;
            for itag = 1:2
                cycTC = cell(1,maxCycles);
    %             longTC = [];
                dd = d.tag(itag).av(visualTrials).align(alignStart);
                if bxExpt
                    hits = strcmp(dd.outcome,'success');
                    misses = strcmp(dd.outcome,'ignore');
                else
                    hits = true(1,length(dd.outcome));
                    misses = false(1,length(dd.outcome));
                end
                for icyc = 1:maxCycles
                    tc = dd.respTC(:,:,dd.nCycles >= icyc & (hits | misses));
                    cycStartOffset = ((icyc-1).*cycLengthFr)+nBaselineFr;
                    cycTC{icyc} = tc(...
                        (cycStartOffset-nBaselineFr+1):(cycStartOffset+nFrames1s),:,:);
                end
                longTC = dd.respTC(1:(longTrialLengthFr+nBaselineFr),:,...
                    dd.nCycles >= nCycLongTC & (hits | misses));

                antiDataExpt(exptN).tag(itag).name = d.tag(itag).name;
                antiDataExpt(exptN).tag(itag).longTC = longTC;
                antiDataExpt(exptN).tag(itag).cycTC = cycTC;
                
                dd = d.tag(itag).av(visualTrials).align(alignTarget);
                easyTarTC = dd.respTC(1:(nBaselineFr+nFrames1s),:,dd.ori > 60);
                tarDataExpt(exptN).tag(itag).easyTarTC = easyTarTC;
            end
            if bxExpt
                d = mousePass(imouse).expt(iexp);
                antiDataExpt_passive(exptN).exptName = [mouse(imouse).mouse_name '-' d.date];
                antiDataExpt_passive(exptN).exptCycLengthFr = cycLengthFr;
                tarDataExpt_passive(exptN).exptName = [mouse(imouse).mouse_name '-' d.date];
                tarDataExpt_passive(exptN).exptCycLengthFr = cycLengthFr;
                for itag = 1:2
                    cycTC = cell(1,maxCycles);
        %             longTC = [];
                    dd = d.tag(itag).av(visualTrials).align(alignStart);
                    hits = true(1,length(dd.outcome));
                    misses = false(1,length(dd.outcome));
                    for icyc = 1:maxCycles
                        tc = dd.respTC(:,:,dd.nCycles >= icyc & (hits | misses));
                        cycStartOffset = ((icyc-1).*cycLengthFr)+nBaselineFr;
                        cycTC{icyc} = tc(...
                            (cycStartOffset-nBaselineFr+1):(cycStartOffset+nFrames1s),:,:);
                    end
                    longTC = dd.respTC(1:(longTrialLengthFr+nBaselineFr),:,...
                        dd.nCycles >= nCycLongTC & (hits | misses));

                    antiDataExpt_passive(exptN).tag(itag).name = d.tag(itag).name;
                    antiDataExpt_passive(exptN).tag(itag).longTC = longTC;
                    antiDataExpt_passive(exptN).tag(itag).cycTC = cycTC;

                    dd = d.tag(itag).av(visualTrials).align(alignTarget);
                    easyTarTC = dd.respTC(1:(nBaselineFr+nFrames1s),:,dd.ori > 60);
                    tarDataExpt_passive(exptN).tag(itag).easyTarTC = easyTarTC;
                end
            end
        end
    end

    shortCycExptInd = cell2mat({antiDataExpt.exptCycLengthFr}) == 11;
%     isShortCycExpt = [];
    respCellsExpt = struct;
    antiAnalysis = struct;
    antiAnalysis_passive = struct;
    tarAnalysis = struct;
    tarAnalysis_passive = struct;
    cellInfo = struct;
    for itag = 1:2
        for iexp = 1:nexp
            firstTC = antiDataExpt(iexp).tag(itag).cycTC{1};
            longTC = antiDataExpt(iexp).tag(itag).longTC;
            tarTC = tarDataExpt(iexp).tag(itag).easyTarTC;

            firstRespCells = ttest(...
                squeeze(mean(firstTC(respwin,:,:),1)),...
                squeeze(mean(firstTC(basewin_0,:,:),1)),...
                'dim',2,'tail','right','alpha',cellGroupsAlpha);
            lateWinRespCells = ttest(...
                squeeze(mean(longTC(lateWinFr,:,:),1)),...
                squeeze(mean(longTC(basewin_0,:,:),1)),...
                'dim',2,'tail','right','alpha',cellGroupsAlpha);
            lateWinSuppCells = ttest(...
                squeeze(mean(longTC(lateWinFr,:,:),1)),...
                squeeze(mean(longTC(basewin_0,:,:),1)),...
                'dim',2,'tail','left','alpha',cellGroupsAlpha);
            tarRespCells = ttest(...
                squeeze(mean(tarTC(respwin,:,:),1)),...
                squeeze(mean(tarTC(basewin_0,:,:),1)),...
                'dim',2,'tail','right','alpha',cellGroupsAlpha);

            respCellsExpt(iexp).exptName = antiDataExpt(iexp).exptName;
            respCellsExpt(iexp).tag(itag).firstRespCells = firstRespCells;
            respCellsExpt(iexp).tag(itag).lateWinRespCells = lateWinRespCells;
            respCellsExpt(iexp).tag(itag).lateWinSuppCells = lateWinSuppCells;
            respCellsExpt(iexp).tag(itag).targetRespCells = tarRespCells;
        end
        
        antiAnalysis.tag(itag).name = antiDataExpt(iexp).tag(itag).name;
        antiAnalysis.tag(itag).longTC = [];
        antiAnalysis.tag(itag).longTCErr = [];
        antiAnalysis.tag(itag).cycTC = cell(1,nCycles);
        antiAnalysis.tag(itag).cycTCErr = cell(1,nCycles);
        antiAnalysis.tag(itag).lateCycTC = [];
        antiAnalysis.tag(itag).lateCycTCErr = [];
        for iexp = 1:nexp
            longTC = antiDataExpt(iexp).tag(itag).longTC;
            cycTC = antiDataExpt(iexp).tag(itag).cycTC(1:nCycles);
            lateCycTC = [];
            for icyc = 5:length(cycTC)
                lateCycTC = cat(3,lateCycTC,cycTC{icyc} - mean(cycTC{icyc}(basewin_0,:,:),1));
            end

            antiAnalysis.tag(itag).longTC = cat(2,antiAnalysis.tag(itag).longTC,...
                mean(longTC,3));
            antiAnalysis.tag(itag).longTCErr = cat(2,antiAnalysis.tag(itag).longTCErr,...
                ste(longTC,3));
            antiAnalysis.tag(itag).cycTC = cellfun(@(x,y) cat(2,x,mean(y,3)),antiAnalysis.tag(itag).cycTC,...
                cycTC,'unif',0);
            antiAnalysis.tag(itag).cycTCErr = cellfun(@(x,y) cat(2,x,ste(y,3)),antiAnalysis.tag(itag).cycTCErr,...
                cycTC,'unif',0);
            antiAnalysis.tag(itag).lateCycTC = cat(2,antiAnalysis.tag(itag).lateCycTC,...
                mean(lateCycTC,3));
            antiAnalysis.tag(itag).lateCycTCErr = cat(2,antiAnalysis.tag(itag).lateCycTCErr,...
                ste(lateCycTC,3));
            
            lateCycRespCells = ttest(...
                squeeze(mean(lateCycTC(respwin,:,:),1)),...
                squeeze(mean(lateCycTC(basewin_0,:,:),1)),...
                'dim',2,'tail','right','alpha',cellGroupsAlpha);
            
            respCellsExpt(iexp).tag(itag).lateCycRespCells = lateCycRespCells;
        end
        
        if bxExpt
            antiAnalysis_passive.tag(itag).longTC = [];
            antiAnalysis_passive.tag(itag).longTCErr = [];
            antiAnalysis_passive.tag(itag).cycTC = cell(1,nCycles);
            antiAnalysis_passive.tag(itag).cycTCErr = cell(1,nCycles);
            antiAnalysis_passive.tag(itag).lateCycTC = [];
            antiAnalysis_passive.tag(itag).lateCycTCErr = [];
            for iexp = 1:nexp
            longTC = antiDataExpt_passive(iexp).tag(itag).longTC;
            cycTC = antiDataExpt_passive(iexp).tag(itag).cycTC(1:nCycles);
            lateCycTC = [];
            for icyc = 5:length(cycTC)
                lateCycTC = cat(3,lateCycTC,cycTC{icyc} - mean(cycTC{icyc}(basewin_0,:,:),1));
            end

            antiAnalysis_passive.tag(itag).longTC = cat(2,antiAnalysis_passive.tag(itag).longTC,...
                mean(longTC,3));
            antiAnalysis_passive.tag(itag).longTCErr = cat(2,antiAnalysis_passive.tag(itag).longTCErr,...
                ste(longTC,3));
            antiAnalysis_passive.tag(itag).cycTC = cellfun(@(x,y) cat(2,x,mean(y,3)),antiAnalysis_passive.tag(itag).cycTC,...
                cycTC,'unif',0);
            antiAnalysis_passive.tag(itag).cycTCErr = cellfun(@(x,y) cat(2,x,ste(y,3)),antiAnalysis_passive.tag(itag).cycTCErr,...
                cycTC,'unif',0);
            antiAnalysis_passive.tag(itag).lateCycTC = cat(2,antiAnalysis_passive.tag(itag).lateCycTC,...
                mean(lateCycTC,3));
            antiAnalysis_passive.tag(itag).lateCycTCErr = cat(2,antiAnalysis_passive.tag(itag).lateCycTCErr,...
                ste(lateCycTC,3));
            end
        end
        
        tarAnalysis.tag(itag).easyTC = [];
        tarAnalysis_passive.tag(itag).easyTC = [];
        for iexp = 1:nexp
            tc = tarDataExpt(iexp).tag(itag).easyTarTC;
            tarAnalysis.tag(itag).easyTC = cat(2,tarAnalysis.tag(itag).easyTC,...
                mean(tc,3));
            if bxExpt
                tc = tarDataExpt_passive(iexp).tag(itag).easyTarTC;
                tarAnalysis_passive.tag(itag).easyTC = cat(2,tarAnalysis_passive.tag(itag).easyTC,...
                    mean(tc,3));
            end
        end
        
        cellInfo.tag(itag).name = antiDataExpt(1).tag(itag).name;
        cellInfo.tag(itag).firstRespCells = [];
        cellInfo.tag(itag).lateWinRespCells = [];
        cellInfo.tag(itag).lateWinSuppCells = [];
        cellInfo.tag(itag).targetRespCells = [];
        cellInfo.tag(itag).lateCycRespCells = [];
        for iexp = 1:nexp
            cellInfo.tag(itag).firstRespCells = cat(1,cellInfo.tag(itag).firstRespCells,...
                logical(cell2mat({respCellsExpt(iexp).tag(itag).firstRespCells}')));
            cellInfo.tag(itag).lateWinRespCells = cat(1,cellInfo.tag(itag).lateWinRespCells,...
                logical(cell2mat({respCellsExpt(iexp).tag(itag).lateWinRespCells}')));
            cellInfo.tag(itag).lateWinSuppCells = cat(1,cellInfo.tag(itag).lateWinSuppCells,...
                logical(cell2mat({respCellsExpt(iexp).tag(itag).lateWinSuppCells}')));
            cellInfo.tag(itag).targetRespCells = cat(1,cellInfo.tag(itag).targetRespCells,...
                logical(cell2mat({respCellsExpt(iexp).tag(itag).targetRespCells}')));
            cellInfo.tag(itag).lateCycRespCells = cat(1,cellInfo.tag(itag).lateCycRespCells,...
                logical(cell2mat({respCellsExpt(iexp).tag(itag).lateCycRespCells}')));
        end
        cellInfo.tag(itag).minRespCells = (mean(antiAnalysis.tag(itag).cycTC{1}(respwin,:),1) > ...
            minRespThreshold)';
    end
%     cellInfo.tag(itag).isTuned = logical(cell2mat({oriTuningExpt.isTuned}))';
%     cellInfo.tag(itag).oriResp = cell2mat({oriTuningExpt.oriResp}');
%     cellInfo.tag(itag).oriRespErr = cell2mat({oriTuningExpt.oriRespErr}');
%     cellInfo.tag(itag).oriFit = cell2mat({oriTuningExpt.fit})';
%     cellInfo.tag(itag).oriPref = cell2mat({oriTuningExpt.oriPref})';
%     cellInfo.tag(itag).hwhm = hwhmFromOriFit(cellInfo.oriFit(:,1:180)',1:180)';

save([fnout 'adaptAnalysis'],'antiAnalysis','tarAnalysis','antiAnalysis_passive','cellInfo','respCellsExpt')
%% plotting params
respTCLim = [-0.005 0.05];
cycTCLim = [-0.01 0.12];
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
oriTCLim = [-0.005 0.08];
targetTCLim = [-0.015 0.08];
outTCLim = [-0.005 0.04];
firstTCLim = [-0.005 0.04];
adaptLim = [0 1];
suppTCLim = [-0.05 0.005];
suppScatLim_win = [-0.2 0.1];
suppScatLim_cyc = [-0.015 0.015];

tcStartFrame = 26;
cycTCEndTimeMs = 350;
cycTCEndFr = 45;
ttLabel_long = 0:500:2500;
ttLabel_cyc = -200:100:cycTCEndTimeMs;
ttLabel_target = -1000:250:900;
preTargetStimLabel = -700:350:0;
nFr_long = size(antiAnalysis.tag(1).longTC,1);
tt_longTC = ((tcStartFrame:nFr_long)-(nBaselineFr+nVisDelayFr)).*(1000/frameRateHz);
ttLabelFr_long = ((ttLabel_long./1000)*frameRateHz)+...
    ((nBaselineFr+nVisDelayFr)-tcStartFrame+1);
ttLabelFr_cyc = ((ttLabel_cyc./1000)*frameRateHz)+...
    ((nBaselineFr+nVisDelayFr)-tcStartFrame+1);
ttLabelFr_target = ((ttLabel_target./1000)*frameRateHz)+...
    ((nBaselineFr+nVisDelayFr_target)+1);

nFr_cyc = size(antiAnalysis.tag(itag).cycTC{1,1},1);
tt_cycTC = ((tcStartFrame:nFr_cyc)-(nBaselineFr+nVisDelayFr)).*(1000/frameRateHz);
% tt_targetTC = ((1:nFr_cyc)-(nBaselineFr+nVisDelayFr_target)).*(1000/frameRateHz);

lateWinTT = ([lateWinFr(1) lateWinFr(end)] - (nBaselineFr+nVisDelayFr))...
    .*(1000/frameRateHz);
respWinTT = ([respwin(1) respwin(end)] - (nBaselineFr+nVisDelayFr))...
    .*(1000/frameRateHz);
respWinTT_target = (...
    [respwin_target(1) respwin_target(end)] - (nBaselineFr+nVisDelayFr_target))...
    .*(1000/frameRateHz);
baseWinTT = (...
    [basewin_0(1) basewin_0(end)] - (nBaselineFr+nVisDelayFr))...
    .*(1000/frameRateHz);

movWinLabelFr = 30:(30+nMovWin-1);
movWinLabelMs = (movWinLabelFr - (nBaselineFr+nVisDelayFr)).*(1000/frameRateHz);

weightLim = [-3 4.2];
binnedWeightLim = [-0.4 0.4];
weightLimSum = [-0.8 0.8];
siLimSum = [-0.5 2.5];

% lateCycRespAll = mean(antiAnalysis.lateCycTC{1}(respwin,:),1);

%%
for itag = 1:2
    ind = cellInfo.tag(itag).firstRespCells & cellInfo.tag(itag).minRespCells;
    cycResp = cellfun(@(x) mean(x(respwin,:),1) - mean(x(basewin_0,:),1),antiAnalysis.tag(itag).cycTC,'unif',0);

    adaptResp = cellfun(@(x) x./cycResp{1},cycResp,'unif',0);

    setFigParams4Print('landscape')
    figure
    if bxExpt
        suptitle(sprintf('behav, %s, Min. & First Stim Resp. Neurons',...
            cellInfo.tag(itag).name))
    else
        suptitle(sprintf('naive, %s, Min. & First Stim Resp. Neurons',...
            cellInfo.tag(itag).name))
    end
    for icyc = 1:nCycles
        subplot(2,nCycles,icyc)
        bl = mean(antiAnalysis.tag(itag).cycTC{icyc}(basewin_0,ind),1);
        y = antiAnalysis.tag(itag).cycTC{icyc}(26:end,ind) - bl;
        yerr = ste(antiAnalysis.tag(itag).cycTC{icyc}(26:end,ind),2);
        shadedErrorBar_chooseColor(tt_cycTC,mean(y,2),yerr,[0 0 0]);
        figXAxis([],'Time from Stim (ms)',[tt_cycTC(1) 350])
        figYAxis([],'dF/F',cycTCLim)
        figAxForm
        hline(0,'k:')
        hold on
        vline(respWinTT,'k--');
        if icyc == 1
            title(sprintf('Stim #%s (%s/%s)',num2str(icyc),num2str(sum(ind)),...
                num2str(length(ind))))
        else
            title(sprintf('Stim #%s',num2str(icyc)))
        end
    end

    subplot 223
    y = cellfun(@(x) mean(x(ind)),cycResp);
    yerr = cellfun(@(x) ste(x(ind),2),cycResp);
    errorbar(1:nCycles,y,yerr,'.')
    figXAxis([],'Stim #',[0 nCycles+1],1:nCycles,1:nCycles)
    figYAxis([],'dF/F',cycTCLim)
    figAxForm
    title(sprintf('First Stim Resp. Cells (%s/%s)',num2str(sum(ind)),num2str(length(ind))))
    subplot 224
    y = cellfun(@(x) mean(x(ind)),adaptResp);
    yerr = cellfun(@(x) ste(x(ind),2),adaptResp);
    errorbar(1:nCycles,y,yerr,'.')
    figXAxis([],'Stim #',[0 nCycles+1],1:nCycles,1:nCycles)
    figYAxis([],'Norm dF/F',[0 1.5])
    figAxForm
    title(sprintf('All Resp. Cells, n=%s', num2str(sum(ind))))

    print([fnout 'adaptation_' cellInfo.tag(itag).name],'-dpdf','-fillpage')
end
%% 
if bxExpt
for itag = 1:2
    ind = cellInfo.tag(itag).firstRespCells & cellInfo.tag(itag).minRespCells;
    cycResp = cellfun(@(x) mean(x(respwin,:),1) - mean(x(basewin_0,:),1),antiAnalysis_passive.tag(itag).cycTC,'unif',0);

    adaptResp = cellfun(@(x) x./cycResp{1},cycResp,'unif',0);

    setFigParams4Print('landscape')
    figure
    suptitle(sprintf('behav - passive expt, %s, Min. & First Stim Resp. Neurons',...
        cellInfo.tag(itag).name))
    
    for icyc = 1:nCycles
        subplot(2,nCycles,icyc)
        bl = mean(antiAnalysis_passive.tag(itag).cycTC{icyc}(basewin_0,ind),1);
        y = antiAnalysis_passive.tag(itag).cycTC{icyc}(26:end,ind) - bl;
        yerr = ste(antiAnalysis_passive.tag(itag).cycTC{icyc}(26:end,ind),2);
        shadedErrorBar_chooseColor(tt_cycTC,mean(y,2),yerr,[0 0 0]);
        figXAxis([],'Time from Stim (ms)',[tt_cycTC(1) 350])
        figYAxis([],'dF/F',cycTCLim)
        figAxForm
        hline(0,'k:')
        hold on
        vline(respWinTT,'k--');
        if icyc == 1
            title(sprintf('Stim #%s (%s/%s)',num2str(icyc),num2str(sum(ind)),...
                num2str(length(ind))))
        else
            title(sprintf('Stim #%s',num2str(icyc)))
        end
    end

    subplot 223
    y = cellfun(@(x) mean(x(ind)),cycResp);
    yerr = cellfun(@(x) ste(x(ind),2),cycResp);
    errorbar(1:nCycles,y,yerr,'.')
    figXAxis([],'Stim #',[0 nCycles+1],1:nCycles,1:nCycles)
    figYAxis([],'dF/F',cycTCLim)
    figAxForm
    title(sprintf('First Stim Resp. Cells (%s/%s)',num2str(sum(ind)),num2str(length(ind))))
    subplot 224
    y = cellfun(@(x) mean(x(ind)),adaptResp);
    yerr = cellfun(@(x) ste(x(ind),2),adaptResp);
    errorbar(1:nCycles,y,yerr,'.')
    figXAxis([],'Stim #',[0 nCycles+1],1:nCycles,1:nCycles)
    figYAxis([],'Norm dF/F',[0 1.5])
    figAxForm
    title(sprintf('All Resp. Cells, n=%s', num2str(sum(ind))))

    print([fnout 'adaptation_passive_' cellInfo.tag(itag).name],'-dpdf','-fillpage')
end
end

%%

hmLim = [-0.3 0.3];
figure
suptitle('Anticipation: Late Resp. Cells')
colormap(brewermap([],'*RdBu'));
for itag = 1:2
    subplot(2,2,itag)
%     ind = logical(cellInfo.tag(itag).lateWinRespCells);
%     ind = logical(cellInfo.tag(itag).lateWinSuppCells);
    ind = cellInfo.tag(itag).lateWinRespCells | ...
        cellInfo.tag(itag).lateWinSuppCells | ...
        cellInfo.tag(itag).firstRespCells | ...
        cellInfo.tag(itag).lateCycRespCells;
    lateWinTC = antiAnalysis.tag(itag).longTC(:,ind);
    lateWinResp = mean(lateWinTC(lateWinFr,:),1);
    [~,lateWinSortInd] = sort(lateWinResp);
    hm = flipud(lateWinTC(:,lateWinSortInd)');
    imagesc(hm(:,tcStartFrame:end))
    hold on
    figXAxis([],'Time from Start (fr)',[],ttLabelFr_long,ttLabelFr_long)
    figYAxis([],'Cell #',[])
    figAxForm
    colorbar
    caxis(hmLim)
    title(sprintf('Behaving, %s Cells',antiAnalysis.tag(itag).name))
    
    if bxExpt
        subplot(2,2,itag+2)
        lateWinTC = antiAnalysis_passive.tag(itag).longTC(:,ind);
        hm = flipud(lateWinTC(:,lateWinSortInd)');
        imagesc(hm(:,tcStartFrame:end))
        hold on
        if strcmp(ds,'FSAV_attentionV1')
            exCellInd = [exampleCell_1,exampleCell_2];
            exCellMat = zeros(1,length(cellInfo.firstRespCells));
            exCellMat(exCellInd) = 1;
            exCellSortInd = find(flip(exCellMat(lateWinSortInd)));
            hline(exCellSortInd,'k-')
        end
        figXAxis([],'Time from Start (fr)',[],ttLabelFr_long,ttLabelFr_long)
        figYAxis([],'Cell #',[])
        figAxForm
        colorbar
        caxis(hmLim)
        title(sprintf('Passive, %s Cells',antiAnalysis.tag(itag).name))
    end
end
print([fnout 'heatmapsAllCells'],'-dpdf','-fillpage')

figure
suptitle('Easy Targets (>60 deg): All Cells')
colormap(brewermap([],'*RdBu'));
for itag = 1:2
    subplot(2,2,itag)
    ind = logical(cellInfo.tag(itag).targetRespCells);
%     ind = logical(cellInfo.tag(itag).lateWinSuppCells);
%     ind = cellInfo.tag(itag).lateWinRespCells | cellInfo.tag(itag).lateWinSuppCells;
    tc = tarAnalysis.tag(itag).easyTC(:,ind);
    resp = mean(tc(respwin,:),1);
    [~,lateWinSortInd] = sort(resp);
    hm = flipud(tc(:,lateWinSortInd)');
    imagesc(hm(:,tcStartFrame:end))
    hold on
    figXAxis([],'Time from Start (fr)',[],ttLabelFr_target,ttLabelFr_target)
    figYAxis([],'Cell #',[])
    figAxForm
    colorbar
    caxis(hmLim)
    title(sprintf('Behaving, %s Cells',antiAnalysis.tag(itag).name))
    
    if bxExpt
        subplot(2,2,itag+2)
        tc = tarAnalysis_passive.tag(itag).easyTC(:,ind);
        hm = flipud(tc(:,lateWinSortInd)');
        imagesc(hm(:,tcStartFrame:end))
        hold on
        figXAxis([],'Time from Start (fr)',[],ttLabelFr_target,ttLabelFr_target)
        figYAxis([],'Cell #',[])
        figAxForm
        colorbar
        caxis(hmLim)
        title(sprintf('Passive, %s Cells',antiAnalysis.tag(itag).name))
    end
end

hmLim = [-0.1 0.1];
figure
suptitle('Late Resp. Cells')
colormap(brewermap([],'*RdBu'));
for itag = 1:2
    subplot(2,2,itag)
    ind = logical(cellInfo.tag(itag).lateCycRespCells);
%     ind = logical(cellInfo.tag(itag).lateWinSuppCells);
%     ind = cellInfo.tag(itag).lateWinRespCells | cellInfo.tag(itag).lateWinSuppCells;
    tc = antiAnalysis.tag(itag).lateCycTC(:,ind);
    resp = mean(tc(respwin,:),1);
    [~,lateWinSortInd] = sort(resp);
    hm = flipud(tc(:,lateWinSortInd)');
    imagesc(hm(:,tcStartFrame:end))
    hold on
    figXAxis([],'Time from Start (fr)',[],ttLabelFr_target,ttLabelFr_target)
    figYAxis([],'Cell #',[])
    figAxForm
    colorbar
    caxis(hmLim)
    title(sprintf('Behaving, %s Cells',antiAnalysis.tag(itag).name))
    
    if bxExpt
        subplot(2,2,itag+2)
        tc = antiAnalysis_passive.tag(itag).lateCycTC(:,ind);
        hm = flipud(tc(:,lateWinSortInd)');
        imagesc(hm(:,tcStartFrame:end))
        hold on
        figXAxis([],'Time from Start (fr)',[],ttLabelFr_target,ttLabelFr_target)
        figYAxis([],'Cell #',[])
        figAxForm
        colorbar
        caxis(hmLim)
        title(sprintf('Passive, %s Cells',antiAnalysis.tag(itag).name))
    end
end

%% cluster analysis
allCells = [];
allCells_pass = [];
cellTag = [];
for itag = 1:2
    ind = cellInfo.tag(itag).lateWinRespCells | ...
        cellInfo.tag(itag).lateWinSuppCells | ...
        cellInfo.tag(itag).firstRespCells | ...
        cellInfo.tag(itag).lateCycRespCells;
    lateWinTC = antiAnalysis.tag(itag).longTC(:,ind);
    lateWinResp = mean(lateWinTC(lateWinFr,:),1);
    [~,lateWinSortInd] = sort(lateWinResp);
    lateWinTC_pass = antiAnalysis_passive.tag(itag).longTC(:,ind);
    allCells = cat(2,allCells,lateWinTC(:,lateWinSortInd));
    allCells_pass = cat(2,allCells_pass,lateWinTC_pass(:,lateWinSortInd));
    cellTag = cat(2,cellTag,repmat({cellInfo.tag(itag).name},1,sum(ind)));
end

nc = size(allCells,2);
allCells_norm = allCells./max(allCells,[],1);
allCells_pass_norm = allCells_pass./max(allCells,[],1);

% allCells_pca = pca(allCells_norm);

nCluster = 3;
[clustID, clustCentroid] = kmeans(allCells_norm',nCluster,'MaxIter',1000000,'Replicates',1000);
[~,clustID_pass] = pdist2(clustCentroid,allCells_pass_norm','euclidean','smallest',1);

figure
subplot 221
histogram(clustID,0:(nCluster+1),'Normalization','probability')
figXAxis([],'Cluster #',[0 (nCluster+2)])
figYAxis([],'Fraction of Cells',[0 1])
figAxForm
subplot 222
histogram(clustID_pass,0:(nCluster+1),'Normalization','probability')
figXAxis([],'Cluster #',[0 (nCluster+2)])
figYAxis([],'Fraction of Cells',[0 1])
figAxForm
subplot 223
for ic = 1:nCluster
    cInd = clustID == ic;
    hold on
    plot(ones(1,sum(cInd)).*ic,find(cInd),'.','MarkerSize',10);
end
tagInd = strcmp(cellTag,'SOM+');
plot(ones(1,sum(tagInd)).*(nCluster+1),find(tagInd),'k.','MarkerSize',10)
figXAxis([],'Cluster #',[0 (nCluster+2)],1:(nCluster+1),...
    cat(2,cellfun(@num2str,num2cell(1:nCluster),'unif',0),{'SOM+'}))
figYAxis([],'Fraction of Cells',[0 nc])
figAxForm
subplot 224
for ic = 1:nCluster
    cInd = clustID_pass == ic;
    hold on
    plot(ones(1,sum(cInd)).*ic,find(cInd),'.','MarkerSize',10);
end
figXAxis([],'Cluster #',[0 (nCluster+1)])
figYAxis([],'Fraction of Cells',[0 nc])
figAxForm
print([fnout 'clusterAnalysisBxPass_clusterID'],'-dpdf','-fillpage')

clusterTC = cell(1,nCluster);
clusterTC_pass = cell(1,nCluster);
for ic = 1:nCluster
    cInd = clustID == ic;
    clusterTC{ic} = allCells(:,cInd);
    cInd = clustID_pass == ic;
    clusterTC_pass{ic} = allCells_pass(:,cInd);    
end

figure
for ic = 1:nCluster
    subplot(2,nCluster,ic)
    y = mean(clusterTC{ic}(tcStartFrame:end,:),2);
    yerr = ste(clusterTC{ic}(tcStartFrame:end,:),2);
    if all(~isnan(yerr))
        shadedErrorBar_chooseColor(tt_longTC,y,yerr,[0 0 0]);
    else
        plot(tt_longTC,y,'k')
    end
    figXAxis([],'Time from start (ms)',[tt_longTC(1) tt_longTC(end)],ttLabel_long)
    figYAxis([],'dF/F',[-0.1 0.1])
    figAxForm
    title(sprintf('Cluster %s, n=%s',num2str(ic),num2str(size(clusterTC{ic},2))))
    
    subplot(2,nCluster,ic+nCluster)
    y = mean(clusterTC_pass{ic}(tcStartFrame:end,:),2);
    yerr = ste(clusterTC_pass{ic}(tcStartFrame:end,:),2);
    if all(~isnan(yerr))
        shadedErrorBar_chooseColor(tt_longTC,y,yerr,[0 0 0]);
    else
        plot(tt_longTC,y,'k')
    end
    figXAxis([],'Time from start (ms)',[tt_longTC(1) tt_longTC(end)],ttLabel_long)
    figYAxis([],'dF/F',[-0.1 0.1])
    figAxForm
    title(sprintf('Cluster %s, n=%s',num2str(ic),num2str(size(clusterTC_pass{ic},2))))
end
print([fnout 'clusterAnalysisBxPass_clusterTC'],'-dpdf','-fillpage')

switchClustID = clustID'-clustID_pass ~= 0;

figure
subplot 121
y = sum(switchClustID & ~tagInd)./sum(~tagInd);
plot(1,y,'.','MarkerSize',20)
hold on
y = sum(switchClustID & tagInd)./sum(tagInd);
plot(2,y,'.','MarkerSize',20);
figXAxis([],'',[0 3],1:2,{'SOM-','SOM+'})
figYAxis([],'Fraction Switch Cluster',[0 0.2])
figAxForm
subplot 122
hold on
for ic = 1:nCluster
    cInd = clustID' == ic & switchClustID;
    n = histcounts(clustID_pass(cInd),1:nCluster+1);
    scatter(ic.*(ones(1,nCluster)),1:nCluster,(n*5)+1)
end
figXAxis([],'Behavior Cluster #',[0 nCluster+1])
figYAxis([],'Passive Cluster #',[0 nCluster+1])
figAxForm
title(sprintf('Cluster Switch Cells, n=%s',num2str(sum(switchClustID))))
print([fnout 'clusterAnalysisBxPass_nClusterSwitch'],'-dpdf','-fillpage')

%%
    figure
    suptitle('Target Resp.')
for itag = 1:2
    ind = logical(cellInfo.tag(itag).targetRespCells);
    tarResp_bx = mean(tarAnalysis.tag(itag).easyTC(respwin,:),1) - mean(tarAnalysis.tag(itag).easyTC(basewin_0,:),1);
    tarResp_pass = mean(tarAnalysis_passive.tag(itag).easyTC(respwin,:),1) - mean(tarAnalysis_passive.tag(itag).easyTC(basewin_0,:),1);

    tarRespLim = [-0.01 0.12];
%     for icyc = 1:nCycles
        subplot(1,2,itag)
        x = tarResp_bx(ind);
        y = tarResp_pass(ind);
        plot(x,y,'.','MarkerSize',10)
        hold on
        plot(tarRespLim,tarRespLim,'k--')
        figXAxis([],'Behavior',tarRespLim)
        figYAxis([],'Passive',tarRespLim)
        figAxForm
        title(sprintf('%s Cells',antiAnalysis.tag(itag).name))
%     end
end

figure
    suptitle('Late Resp.')
for itag = 1:2
    ind = logical(cellInfo.tag(itag).lateCycRespCells);
    lateResp_bx = mean(antiAnalysis.tag(itag).lateCycTC(respwin,:),1) - mean(antiAnalysis.tag(itag).lateCycTC(basewin_0,:),1);
    lateResp_pass = mean(antiAnalysis_passive.tag(itag).lateCycTC(respwin,:),1) - mean(antiAnalysis_passive.tag(itag).lateCycTC(basewin_0,:),1);

    lateRespLim = [-0.01 0.1];
%     for icyc = 1:nCycles
        subplot(1,2,itag)
        x = lateResp_bx(ind);
        y = lateResp_pass(ind);
        plot(x,y,'.','MarkerSize',10)
        hold on
        plot(lateRespLim,lateRespLim,'k--')
        figXAxis([],'Behavior',lateRespLim)
        figYAxis([],'Passive',lateRespLim)
        figAxForm
        title(sprintf('%s Cells',antiAnalysis.tag(itag).name))
%     end
end