clear all
close all
ds = '';
%%
rc = behavConstsAV;
if strcmp(rc.name,'ashle')
    dataGroup = ['awFSAVdatasets' ds];
else
    dataGroup = [];
end
eval(dataGroup)
titleStr = ds;
if strcmp(titleStr, '')
    titleStr = 'V1_100ms';
else
    titleStr = titleStr(2:end);
end
str = unique({expt.SubNum});
mouse_str = ['i' strjoin(str,'_i')];

load(fullfile(rc.caOutputDir,dataGroup,[mouse_str '_CaSummary' ds '.mat']));
% load(fullfile(rc.caOutputDir,dataGroup,[titleStr '_' mouse_str '_modCells.mat']));
fnout = fullfile(rc.caOutputDir, dataGroup,[titleStr '_startAlign']); 

%%
outcomecolmat = {'k';'r'};
avcolmat = {'k','c'};
avmat = {'visual';'auditory'};
outcomemat = {'hit';'miss'};
pressAlign = 1;
visualTrials = 1;
auditoryTrials = 2;
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

%% response to each of first 5 cyc
basewin = 32:34;
respwin = 36:38;
nCycles = 8;
tcEaCellEaCycAllOutcome = cell(nCycles,nexp,nav);
trialTCEaCycEaCell = cell(nCycles,nexp,nav);
nTrialsEaCycAllOutcome = zeros(nCycles,nexp,nav);
tcEaCellEaCycByOutcome = cell(nCycles,nexp,nav,nout);
nTrialsEaCycByOutcome = zeros(nCycles,nexp,nav,nout);
exptName = cell(1,nexp);
h = cell(1,nexp);
p = cell(nexp,nCycles);
responsiveCellsTtestEaCyc = cell(nexp,nCycles);
responsive2itiCellsTtestEaCyc = cell(nexp,nCycles);
suppressed2itiCellsTtestEaCyc = cell(nexp,nCycles);
for imouse = 1:size(mouse,2)
    for iexp = 1:size(mouse(imouse).expt,2)
            if imouse == 1 && iexp == 1
                exptN = 1;
            else
                exptN = exptN+1;
            end
        exptName{exptN} = [mouse(imouse).expt(iexp).mouse_name '-' ...
            mouse(imouse).expt(iexp).date];
        for iav = 1:nav
            d = mouse(imouse).expt(iexp).align(pressAlign).av(iav);        
            for icyc = 1:nCycles
                tcThisCycAllOutcome = cat(3,d.outcome(hitTrials).cmlvCycResp{icyc},...
                    d.outcome(missTrials).cmlvCycResp{icyc});
                extraFrames = cycLengthFr*(icyc-1);
                trialTCEaCycEaCell{icyc,exptN,iav} = mean(tcThisCycAllOutcome,3);
                for ioutcome = 1:nout
                    tcThisOutcome = d.outcome(ioutcome).cmlvCycResp{icyc};
                    nTrialsEaCycByOutcome(icyc,exptN,iav,ioutcome) = size(tcThisOutcome,3);
                    baselineByOutcome = mean(tcThisOutcome(basewin+extraFrames,:,:),1);
                    tcEaCellEaCycByOutcome{icyc,exptN,iav,ioutcome}...
                        = mean(bsxfun(@minus,tcThisOutcome(...
                        end-(cycLengthFr+nBaselineFr)+1:end,:,:),baselineByOutcome),3);
                end
                if iav == 1
                    if icyc == 1
                        base_resp_thisCyc = mean(tcThisCycAllOutcome(basewin,:,:),1);
                        base_resp_iti = base_resp_thisCyc;
                    else
                        base_resp_thisCyc = mean(tcThisCycAllOutcome(basewin+extraFrames,:,:),1);
                        base_resp_iti = mean(tcThisCycAllOutcome(basewin,:,:),1);
                    end
                    responsiveCellsTtestEaCyc{exptN,icyc} = zeros(1,size(tcThisCycAllOutcome,2));
                    responsive2itiCellsTtestEaCyc{exptN,icyc} = zeros(1,size(tcThisCycAllOutcome,2));
                    suppressed2itiCellsTtestEaCyc{exptN,icyc} = zeros(1,size(tcThisCycAllOutcome,2));
                    hResp{exptN} = zeros(1,size(tcThisCycAllOutcome,2));
                    for iCell = 1:size(tcThisCycAllOutcome,2)
                        if icyc == 1
                            [hResp{exptN}(1,iCell), ~] = ttest(squeeze(mean(...
                                tcThisCycAllOutcome(respwin,iCell,:),1)),...
                                squeeze(base_resp_thisCyc(:,iCell,:)),'tail','right');
                            if hResp{exptN}(1,iCell) == 1
                                [max_val, max_ind{exptN}(1,iCell)] = ...
                                    max(diff(mean(...
                                    tcThisCycAllOutcome(26:end,iCell,:),3),1),[],1);
                                if max_ind{exptN}(1,iCell) < respwin(end)-25
                                    responsiveCellsTtestEaCyc{exptN,icyc}(1,iCell) = 1;
                                end
                            end
                        else
                            [responsiveCellsTtestEaCyc{exptN,icyc}(1,iCell), ~]...
                                = ttest2(squeeze(mean(...
                                tcThisCycAllOutcome(respwin,iCell,:),1)),...
                                squeeze(base_resp_thisCyc(:,iCell,:)),'tail','right');
                        end
                        [responsive2itiCellsTtestEaCyc{exptN,icyc}(1,iCell), ~]...
                            = ttest2(squeeze(mean(tcThisCycAllOutcome(respwin+extraFrames,iCell,:),1)),...
                            squeeze(base_resp_iti(:,iCell,:)),'tail','right');
                        [suppressed2itiCellsTtestEaCyc{exptN,icyc}(1,iCell), ~]...
                            = ttest2(squeeze(mean(tcThisCycAllOutcome(respwin+extraFrames,iCell,:),1)),...
                            squeeze(base_resp_iti(:,iCell,:)),'tail','left');
                        
                    end
                end
                thisExptRespEaCyc = nan(nCycles,size(tcThisCycAllOutcome,2));
                responseAllTrials = mean(tcThisCycAllOutcome(...
                    respwin+extraFrames,:,:),1);
                baselineAllTrials = mean(tcThisCycAllOutcome(basewin+extraFrames,:,:),1);
                thisExptRespEaCyc(icyc,:) = squeeze(mean(responseAllTrials - ...
                    baselineAllTrials,3));
                nTrialsEaCycAllOutcome(icyc,exptN,iav) = size(tcThisCycAllOutcome,3);
                tcEaCellEaCycAllOutcome{icyc,exptN,iav} = mean(bsxfun(@minus,tcThisCycAllOutcome(...
                    end-(cycLengthFr+nBaselineFr)+1:end,:,:),...
                    baselineAllTrials),3);

            end
        end
    end
end
%% figures

figure; 
suptitle('response to each vis stim, each experiment (n cells)')
for i = 1:nexp
    subplot(4,5,i);
    ind = find(responsiveCellsTtestEaCyc{i,1});
    nCells = length(ind);
    for icyc = 1:5
        plot(mean(tcEaCellEaCycAllOutcome{icyc,i,visualTrials}(25:end,ind),2)); 
        hold on; 
        title(num2str(length(ind)))
    end 
end

[tcEaCycResponsiveCells,avgResponseEaCycAV,semResponseEaCycAV] = ...
    getCycResponse4CellIndAV(tcEaCellEaCycAllOutcome,...
    responsiveCellsTtestEaCyc,visualTrials,respwin,1);
nCellsResp1Stim = size(tcEaCycResponsiveCells{1,1},2);

lateResponseCellInd = combineCellsWithLogicals(...
    responsive2itiCellsTtestEaCyc(:,6:nCycles));
[tcEaCycLateRespCells,avgLateRespEaCycAV,semLateRespEaCycAV] = ...
    getCycResponse4CellIndAV(tcEaCellEaCycAllOutcome,...
    lateResponseCellInd,visualTrials,respwin,0);
nCellsResp6Stim = size(tcEaCycLateRespCells{1,1},2);

lateSuppCellInd = combineCellsWithLogicals(...
    suppressed2itiCellsTtestEaCyc(:,6:nCycles));
[tcEaCycLateSuppCells,avgLateSuppEaCycAV,semLateSuppEaCycAV] = ...
    getCycResponse4CellIndAV(tcEaCellEaCycAllOutcome,...
    lateSuppCellInd,visualTrials,respwin,0);
nCellsSupp6Stim = size(tcEaCycLateSuppCells{1,1},2);

figure
suptitle('responsive cells, adaptation analysis window, visual vs auditory trials')
subplot(2,2,1)
for icyc = 1:nCycles
    plot(((26:size(tcEaCycResponsiveCells{icyc,visualTrials},1))-33).*(1000/frRateHz),mean(tcEaCycResponsiveCells{icyc,visualTrials}(26:end,:),2))
    hold on
end
title(['n = ' num2str(nCellsResp1Stim)])
xlabel('Time (ms)')
ylabel('dF/F')
vline((respwin([1 end])-33).*(1000/frRateHz))

subplot(2,2,3)
for iav = 1:nav
    errorbar(1:nCycles, avgResponseEaCycAV(:,iav), semResponseEaCycAV(:,iav));
    hold on
end
xlim([0 nCycles+1])
xlabel('Stimulus number')
ylabel('dF/F')
legend({'visual','auditory'})

normTCEaCycResponsiveCells = getNormTCByCondition(tcEaCycResponsiveCells,1,...
    visualTrials,respwin);

responsiveCellsPctChangeEaCyc = getPctChangeFromCyc1AV(...
    avgResponseEaCycAV,visualTrials);
lateRespCellsPctChangeEaCyc = getPctChangeFromCyc1AV(...
    avgLateRespEaCycAV,visualTrials);
lateSuppCellsPctChangeEaCyc = getPctChangeFromCyc1AV(...
    1+avgLateSuppEaCycAV,visualTrials);

[avgNormTCEaCycAV,avgNormResponseEaCycAV,semNormResponseEaCycAV] = ...
    getNormCycRespAV(normTCEaCycResponsiveCells,visualTrials,respwin);
% [avgNormTCLateRespEaCycAV,avgNormLateRespEaCycAV,semNormLateRespEaCycAV] = ...
%     getNormCycRespAV(normTCEaCycLateRespCells,visualTrials,respwin);
% [avgNormTCLateSuppEaCycAV,avgNormLateSuppEaCycAV,semNormLateSuppEaCycAV] = ...
%     getNormCycRespAV(normTCEaCycLateSuppCells,visualTrials,respwin);



subplot(2,2,2)
for icyc = 1:nCycles
    plot(((26:size(avgNormTCEaCycAV,1))-33).*(1000/frRateHz),avgNormTCEaCycAV(26:end,icyc,visualTrials))
    hold on
end
title(['n = ' num2str(nCellsResp1Stim)])
xlabel('Time (ms)')
ylabel('Normalized dF/F')
vline((respwin([1 end])-33).*(1000/frRateHz))

subplot(2,2,4)
for iav = 1:nav
    errorbar(1:nCycles, avgNormResponseEaCycAV(:,iav), semNormResponseEaCycAV(:,iav));
    hold on
end
ylim([0 1.1])
xlim([0 nCycles+1])
xlabel('Stimulus number')
ylabel('Normalized dF/F')
legend({'visual','auditory'})


% print('Z:\home\lindsey\Analysis\2P\Adaptation\Adaptation_allCells_AV_AWdataset.pdf','-dpdf','-bestfit')
%% cell groups AV
setFigParams4Print('portrait')
responseLim = [-0.01 0.02];
normRespLim = [0 1.1];
cellType = {'respond to 1st stim';'responsive late in trial';'suppressed late in trial'};
figure;

subplot(3,2,1)
for icyc = 1:nCycles
    plot(((26:size(tcEaCycLateRespCells{icyc,visualTrials},1))-33).*(1000/frRateHz),...
        mean(tcEaCycLateRespCells{icyc,visualTrials}(26:end,:),2))
    hold on
end
title(sprintf('%s (%s)',cellType{2}, num2str(nCellsResp6Stim)))

subplot(3,2,2)
for icyc = 1:nCycles
    plot(((26:size(tcEaCycLateSuppCells{icyc,visualTrials},1))-33).*(1000/frRateHz),...
        mean(tcEaCycLateSuppCells{icyc,visualTrials}(26:end,:),2))
    hold on
end
title(sprintf('%s (%s)',cellType{3}, num2str(nCellsSupp6Stim)))

subplot(3,2,3)
for iav = 1:nav
    errorbar(1:nCycles, avgLateRespEaCycAV(:,iav), semLateRespEaCycAV(:,iav));
    hold on
end
subplot(3,2,4)
for iav = 1:nav
    errorbar(1:nCycles, avgLateSuppEaCycAV(:,iav), semLateSuppEaCycAV(:,iav));
    hold on
end

for itype = 1:2
    subplot(3,2,itype)
    figXAxis([],'Time (ms)',[])
    figYAxis([],'dF/F',responseLim)
    figAxForm([])
    
    subplot(3,2,itype+2)
    figXAxis([],'Stimulus number',[0 nCycles+1])
    figYAxis([],'dF/F',responseLim)
    figAxForm([])
    legend({'visual','auditory'})
    title(cellType{itype+1})
end

subplot(3,2,5)
normRespLim = [0 1.1];
for iav = 1:nav
    plot(1:nCycles, lateRespCellsPctChangeEaCyc(:,iav));
    hold on
end
figXAxis([],'Stimulus number',[0 nCycles+1])
figYAxis([],'% change',normRespLim)
figAxForm([])
legend({'visual','auditory'})
title(cellType{2})
subplot(3,2,6)
normRespLim = [0.99 1.01];
for iav = 1:nav
    plot(1:nCycles, lateSuppCellsPctChangeEaCyc(:,iav));
    hold on
end
figXAxis([],'Stimulus number',[0 nCycles+1])
figYAxis([],'% change',normRespLim)
figAxForm([])
legend({'visual','auditory'})
title(cellType{3})
print([fnout 'adaptationAVLateCellGroups'],'-dpdf','-fillpage')

%% cell group time-courses
trialTCEaCycResponsiveCells = getCycResponse4CellIndAV(trialTCEaCycEaCell,...
    responsiveCellsTtestEaCyc,visualTrials,respwin,1);
trialTCEaCycLateRespCells = getCycResponse4CellIndAV(trialTCEaCycEaCell,...
    lateResponseCellInd,visualTrials,respwin,1);
trialTCEaCycLateSuppCells = getCycResponse4CellIndAV(trialTCEaCycEaCell,...
    lateSuppCellInd,visualTrials,respwin,1);

responseLim = [-0.025 0.06];
setFigParams4Print('portrait')
figure
subplot 321
for iav = 1:nav
    y = mean(trialTCEaCycResponsiveCells{6,iav},2);
    yerr = ste(trialTCEaCycResponsiveCells{6,iav},2);
    tt = ((26:length(y))-33).*(1000/frRateHz);
    shadedErrorBar(tt,y(26:length(y)),yerr(26:length(y)),avcolmat{iav})
    hold on
end
title(cellType{1})
subplot 322
for iav = 1:nav
    y = mean(trialTCEaCycResponsiveCells{nCycles,iav},2);
    yerr = ste(trialTCEaCycResponsiveCells{nCycles,iav},2);
    tt = ((26:length(y))-33).*(1000/frRateHz);
    shadedErrorBar(tt,y(26:length(y)),yerr(26:length(y)),avcolmat{iav})
    hold on
end
title(cellType{1})

subplot 323
for iav = 1:nav
    y = mean(trialTCEaCycLateRespCells{6,iav},2);
    yerr = ste(trialTCEaCycLateRespCells{6,iav},2);
    tt = ((26:length(y))-33).*(1000/frRateHz);
    shadedErrorBar(tt,y(26:length(y)),yerr(26:length(y)),avcolmat{iav})
    hold on
end
title(cellType{2})
subplot 324
for iav = 1:nav
    y = mean(trialTCEaCycLateRespCells{nCycles,iav},2);
    yerr = ste(trialTCEaCycLateRespCells{nCycles,iav},2);
    tt = ((26:length(y))-33).*(1000/frRateHz);
    shadedErrorBar(tt,y(26:length(y)),yerr(26:length(y)),avcolmat{iav})
    hold on
end
title(cellType{2})

subplot 325
for iav = 1:nav
    y = mean(trialTCEaCycLateSuppCells{6,iav},2);
    yerr = ste(trialTCEaCycLateSuppCells{6,iav},2);
    tt = ((26:length(y))-33).*(1000/frRateHz);
    shadedErrorBar(tt,y(26:length(y)),yerr(26:length(y)),avcolmat{iav})
    hold on
end
title(cellType{3})
subplot 326
for iav = 1:nav
    y = mean(trialTCEaCycLateSuppCells{nCycles,iav},2);
    yerr = ste(trialTCEaCycLateSuppCells{nCycles,iav},2);
    tt = ((26:length(y))-33).*(1000/frRateHz);
    shadedErrorBar(tt,y(26:length(y)),yerr(26:length(y)),avcolmat{iav})
    hold on
end
title(cellType{3})

for iplot = 1:6
    subplot(3,2,iplot)
    figXAxis([],'time (ms)',[tt(1) tt(end)])
    figYAxis([],'dF/F',responseLim)    
end
print([fnout 'tcAVLateCellGroups'],'-dpdf','-fillpage')
%% baseline

[avgBaselineEaCycAV,semBaselineEaCycAV] = ...
    getCycBaseline4CellIndAV(trialTCEaCycEaCell,...
    responsiveCellsTtestEaCyc,visualTrials,basewin,respwin,cycLengthFr,1);

[avgBaselineLateRespEaCycAV,semBaselineLateRespEaCycAV] = ...
    getCycBaseline4CellIndAV(trialTCEaCycEaCell,...
    lateResponseCellInd,visualTrials,basewin,respwin,cycLengthFr,0);

[avgBaselineLateSuppEaCycAV,semBaselineLateSuppEaCycAV] = ...
    getCycBaseline4CellIndAV(trialTCEaCycEaCell,...
    lateSuppCellInd,visualTrials,basewin,respwin,cycLengthFr,0);

setFigParams4Print('landscape')
responseLim = [-0.025 0.04];
normRespLim = [0 1.1];
figure
suptitle('baseline each cycle')
for iav = 1:nav
    subplot(1,3,1)
    errorbar(1:nCycles, avgBaselineEaCycAV(:,iav), ...
        semBaselineEaCycAV(:,iav));
    hold on
    subplot(1,3,2)
    errorbar(1:nCycles, avgBaselineLateRespEaCycAV(:,iav), ...
        semBaselineLateRespEaCycAV(:,iav));
    hold on
    subplot(1,3,3)
    errorbar(1:nCycles, avgBaselineLateSuppEaCycAV(:,iav), ...
        semBaselineLateSuppEaCycAV(:,iav));
    hold on
end

for itype = 1:3
    subplot(1,3,itype)
    figXAxis([],'Stimulus number',[0 nCycles+1])
    figYAxis([],'dF/F',responseLim)
    figAxForm([])
    legend({'visual','auditory'},'location','northeastoutside')
    title(cellType{itype})
end
print([fnout 'baseAVCellGroups'],'-dpdf','-fillpage')

%% second stim
secondStimRespInd = responsiveCellsTtestEaCyc(:,2);%cellfun(@(x,y) x == 0 & y == 1,responsiveCellsTtestEaCyc(:,1),responsiveCellsTtestEaCyc(:,2),'unif',0);

[tcEaCycSecondRespCells,avgSecondRespEaCycAV,semSecondRespEaCycAV] = ...
    getCycResponse4CellIndAV(tcEaCellEaCycAllOutcome,...
    secondStimRespInd,visualTrials,respwin,1);
nCellsResp2Stim = size(tcEaCycSecondRespCells{1,1},2);

normTCEaCycSecondRespCells = getNormTCByCondition(tcEaCycSecondRespCells,1,...
    visualTrials,respwin);
[avgNormSecondRespTCEaCycAV,avgNormSecondRespEaCycAV,semNormSecondRespEaCycAV] = ...
    getNormCycRespAV(normTCEaCycSecondRespCells,visualTrials,respwin);

figure
suptitle('2nd stim responsive cells')
subplot(2,2,1)
for icyc = 1:nCycles
    plot(((26:size(tcEaCycSecondRespCells{icyc,visualTrials},1))-33).*(1000/frRateHz),...
        mean(tcEaCycSecondRespCells{icyc,visualTrials}(26:end,:),2))
    hold on
end
title(['n = ' num2str(nCellsResp2Stim)])
xlabel('Time (ms)')
ylabel('dF/F')
vline((respwin([1 end])-33).*(1000/frRateHz))

subplot(2,2,3)
for iav = 1:nav
    errorbar(1:nCycles, avgSecondRespEaCycAV(:,iav), semSecondRespEaCycAV(:,iav));
    hold on
end
xlim([0 nCycles+1])
xlabel('Stimulus number')
ylabel('dF/F')
legend({'visual','auditory'})

subplot(2,2,2)
for icyc = 1:nCycles
    plot(((26:size(avgNormSecondRespTCEaCycAV,1))-33).*(1000/frRateHz),...
        avgNormSecondRespTCEaCycAV(26:end,icyc,visualTrials))
    hold on
end
title(['n = ' num2str(nCellsResp2Stim)])
xlabel('Time (ms)')
ylabel('Normalized dF/F')
vline((respwin([1 end])-33).*(1000/frRateHz))

subplot(2,2,4)
for iav = 1:nav
    errorbar(1:nCycles, avgNormSecondRespEaCycAV(:,iav), ...
        semNormSecondRespEaCycAV(:,iav));
    hold on
end
ylim([0.8 1.1])
xlim([0 nCycles+1])
xlabel('Stimulus number')
ylabel('Normalized dF/F')
legend({'visual','auditory'})
print([fnout 'AVrespSecondStim'],'-dpdf','-fillpage')
%% hits vs miss
tcEaCycByOutcome = cell(nCycles,nav,nout);
for icyc = 1:nCycles
    for iexp = 1:nexp
        for iav = 1:nav
            for ioutcome = 1:nout
                tcEaCycByOutcome{icyc,iav,ioutcome} = [tcEaCycByOutcome{icyc,iav,ioutcome} tcEaCellEaCycByOutcome{icyc,iexp,iav,ioutcome}]; 
            end
        end
    end
end


for iexp = 1:nexp
    if iexp == 1
        offset = 0;
        cellTypeInd = cell(1,3);
    else
        offset = size(tcEaCellEaCycByOutcome{1,iexp-1,1,1},2);
    end
    cellTypeInd{1} = [cellTypeInd{1} responsiveCells1StimInd{iexp}+offset];
    cellTypeInd{2} = [cellTypeInd{2} responsiveCells6StimInd{iexp}+offset];
    cellTypeInd{3} = [cellTypeInd{3} suppressedCells6StimInd{iexp}+offset];
end

avgResponseEaCycOutcome = cell(1,3);
semResponseEaCycOutcome = cell(1,3);
for itype = 1:3
    ind = cellTypeInd{itype};
    avgResponseEaCycOutcome{itype} = cell2mat(cellfun(@(x) nanmean(nanmean(x(respwin,ind),1),2),tcEaCycByOutcome,'unif',0));
    semResponseEaCycOutcome{itype} = cell2mat(cellfun(@(x) ste(nanmean(x(respwin,ind),1),2),tcEaCycByOutcome,'unif',0));
end
avgTCEaCycByOutcome = cellfun(@(x) mean(x(:,cellTypeInd{1}),2), tcEaCycByOutcome,'unif',0);
normTCEaCycByOutcome = cellfun(@(x) x(:,cellTypeInd{1})./nanmean(tcEaCycByOutcome{1,1,1}(respwin,cellTypeInd{1}),1),tcEaCycByOutcome,'unif',0);
avgNormTCEaCycByOutcome = cellfun(@(x) nanmean(x,2),normTCEaCycByOutcome,'unif',0);
avgNormResponseEaCycByOutcome = cell2mat(cellfun(@(x) nanmean(nanmean(x(respwin,:),1),2),normTCEaCycByOutcome,'unif',0));
semNormResponseEaCycByOutcome = cell2mat(cellfun(@(x) ste(nanmean(x(respwin,:),1),2),normTCEaCycByOutcome,'unif',0));

setFigParams4Print('landscape')
figure;
suptitle('visual trials, hits')
subplot 121
for icyc = 1:nCycles
    plot(((26:size(avgNormTCEaCycAV,1))-33).*(1000/frRateHz),avgTCEaCycByOutcome{icyc,visualTrials,hitTrials}(26:end));
    hold on
end
figXAxis([],'time (ms)',[])
figYAxis([],'dF/F',[])
subplot 122
for icyc = 1:nCycles
    plot(((26:size(avgNormTCEaCycAV,1))-33).*(1000/frRateHz),avgNormTCEaCycByOutcome{icyc,visualTrials,hitTrials}(26:end));
    hold on
end
figXAxis([],'time (ms)',[])
figYAxis([],'Normalized dF/F',[])


outcomecolmat = {'k';'r'};
avcolmat = {'k','c'};
avmat = {'visual';'auditory'};
outcomemat = {'hit';'miss'};
responseLim = [-0.002 0.02];
setFigParams4Print('portrait')
figure;
suptitle('response normalized to visual hits')
for iav = 1:nav
    subplot(2,nav,iav)
    for ioutcome = 1:nout
        hold on
        x = avgResponseEaCycOutcome(:,iav,ioutcome);
        xerr = semResponseEaCycOutcome(:,iav,ioutcome);
        errorbar(1:nCycles,x,xerr,outcomecolmat{ioutcome})
    end
    title(titlemat{iav})
    figXAxis([],'Stimulus number',[0 nCycles+1],1:nCycles,1:nCycles)
    figYAxis([],'dF/F',responseLim)
    legend(outcomemat)
    subplot(2,nav,iav+2)
    for ioutcome = 1:nout
        hold on
        x = avgNormResponseEaCycByOutcome(:,iav,ioutcome);
        xerr = semNormResponseEaCycByOutcome(:,iav,ioutcome);
        errorbar(1:nCycles,x,xerr,outcomecolmat{ioutcome})
    end
    title(titlemat{iav})
    figXAxis([],'Stimulus number',[0 nCycles+1],1:nCycles,1:nCycles)
    figYAxis([],'Normalized dF/F',[-0.3 1.5])
    legend(outcomemat)
end
print([fnout 'adaptationAVByOutcome'],'-dpdf','-fillpage')

setFigParams4Print('landscape')
figure;
for itype = 1:3
    subplot(1,3,itype)
    for iav = 1:nav
        x = avgResponseEaCycOutcome{itype}(:,iav,hitTrials);
        xerr = semResponseEaCycOutcome{itype}(:,iav,hitTrials);
        errorbar(1:nCycles,x,xerr,avcolmat{iav})
        hold on
    end
    title(cellType{itype});
    figXAxis([],'Stimulus number',[0 nCycles+1],1:nCycles,1:nCycles)
%     figYAxis([],'Normalized dF/F',[0 1.1])
    figAxForm([])
    legend(avmat)
end
%%
figure; 
colmat = ['k','c'];
norm_tcs_cum = cell(1,nav);
temp_tcs_cum = cell(1,nav);
for icyc = 1:nCycles
    subplot(3,4,icyc)
    for iav = 1:nav
        shadedErrorBar(((26:size(avgNormTCEaCycAV,1))-33).*(1000/frRateHz), mean(normTCEaCycResponsiveCells{icyc,iav}(26:end,:),2), ste(normTCEaCycResponsiveCells{icyc,iav}(26:end,:),2),colmat(:,iav));
        hold on
        if icyc > 2 & icyc < 7
            norm_tcs_cum{iav} = cat(3, norm_tcs_cum{iav}, normTCEaCycResponsiveCells{icyc,iav});
            temp_tcs_cum{iav} = cat(3, temp_tcs_cum{iav}, tcEaCycResponsiveCells{icyc,iav}(:,big_ind));
        end
    end
    title(['Cycle #' num2str(icyc)])
end
subplot(3,2,5)
for iav = 1:nav
    shadedErrorBar(((26:size(avgNormTCEaCycAV,1))-33).*(1000/frRateHz), mean(mean(norm_tcs_cum{iav}(26:end,:,:),3),2), ste(mean(normTCEaCycResponsiveCells{icyc,iav}(26:end,:,:),3),2),colmat(:,iav));
    hold on
end
title(['Cycle #3-6'])

subplot(3,2,6)
scatter(mean(mean(temp_tcs_cum{1}(respwin,:,:),1),3),mean(mean(temp_tcs_cum{2}(respwin,:,:),1),3),'ok')
[h_ppr, p_ppr] = ttest(mean(mean(temp_tcs_cum{1}(respwin,:,:),1),3),mean(mean(temp_tcs_cum{2}(respwin,:,:),1),3));
xlabel('Visual')
ylabel('Auditory')
xlim([-0.01 0.05])
ylim([-0.01 0.05])
refline(1,0) 
axis square
title(['p = ' num2str(chop(p_ppr,2))])

print('Z:\home\lindsey\Analysis\2P\Adaptation\Adaptation_allCellsByCycle_AV_AWdataset.pdf','-dpdf','-bestfit')

