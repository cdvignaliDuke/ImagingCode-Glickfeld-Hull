function plotFSAVadaptation(expt,mouse,fnout)
pressAlign = 1;
visualTrials = 1;
hitTrials = 1;
missTrials = 2;
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
exptColors = brewermap(nexp+1,'YlOrRd');
exptColors = exptColors(2:nexp+1,:);
%%
% trials that are at least 1 cyc long - average for each experiment
meanTC_1cyc = cell(1,nexp);
errTC_1cyc = cell(1,nexp);
for imouse = 1:size(mouse,2)
    for iexp = 1:size(mouse(imouse).expt,2)
        if imouse == 1 && iexp == 1
            exptN = 1;
        else
            exptN = exptN+1;
        end
        responsiveCells = mouse(imouse).expt(iexp).cells(8).ind;
        d = mouse(imouse).expt(iexp).align(pressAlign).av(visualTrials);
        tc_1cyc = cat(3,d.outcome(hitTrials).cmlvCycResp{1},...
            d.outcome(missTrials).cmlvCycResp{1});
        meanTC_1cyc{exptN} = mean(tc_1cyc(:,responsiveCells,:),3);
        errTC_1cyc{exptN} = ste(tc_1cyc(:,responsiveCells,:),3);        
    end
end

meanTCacrossCells_1cyc = cell2mat(cellfun(@(x) mean(x,2),meanTC_1cyc,...
    'unif',0));
figure;
plot(timestamp_1cyc,meanTCacrossCells_1cyc)
hold on
vline([0 onTimeFr],'k--')
vline(nMonitorDelayFr,'r--')
figXAxis([],'time (frames)',plotTimeRange_1cyc)
figYAxis([],'dF/F',[])
title('response to 1st stim, all responsive cells')


% choose response window
nRespFr = 3;
nDelayFr = 4;
responseWindow = nBaselineFr+nDelayFr:nBaselineFr+nDelayFr+nRespFr;
baselineWindow = responseWindow(1)-4:responseWindow(1)-2;

% response to each of first 5 cyc
nCycles = 5;
tcEaCellEaCyc = cell(nCycles,nexp);
meanRespEaCell = cell(1,nexp);
nTrialsEaCyc = zeros(nCycles,nexp);
exptName = cell(1,nexp);
for imouse = 1:size(mouse,2)
    for iexp = 1:size(mouse(imouse).expt,2)
            if imouse == 1 && iexp == 1
                exptN = 1;
            else
                exptN = exptN+1;
            end
        exptName{exptN} = [mouse(imouse).expt(iexp).mouse_name '-' ...
            mouse(imouse).expt(iexp).date];
        responsiveCells = mouse(imouse).expt(iexp).cells(8).ind;
        d = mouse(imouse).expt(iexp).align(pressAlign).av(visualTrials);        
        for icyc = 1:nCycles
            tcThisCyc = cat(3,d.outcome(hitTrials).cmlvCycResp{icyc}...
                (:,responsiveCells,:),d.outcome(missTrials).cmlvCycResp...
                {icyc}(:,responsiveCells,:));
            if icyc == 1
                thisExptRespEaCyc = nan(nCycles,size(tcThisCyc,2));
            end
            extraFrames = cycLengthFr*(icyc-1);
            responseAllTrials = mean(tcThisCyc(...
                responseWindow+extraFrames,:,:),1);
            baselineAllTrials = mean(tcThisCyc(...
                baselineWindow+extraFrames,:,:),1);
            thisExptRespEaCyc(icyc,:) = squeeze(mean(responseAllTrials - ...
                baselineAllTrials,3));
            nTrialsEaCyc(icyc,exptN) = size(responseAllTrials,3);
            tcEaCellEaCyc{icyc,exptN} = mean(bsxfun(@minus,tcThisCyc(...
                end-(cycLengthFr+nBaselineFr)+1:end,:,:),...
                baselineAllTrials),3);
                
        end
        meanRespEaCell{exptN} = thisExptRespEaCyc;
    end
end
meanRespEaCellNorm = cellfun(@(x) bsxfun(@rdivide,x,x(1,:)),...
    meanRespEaCell,'unif',0);
tcEaCellEaCycNorm = cellfun(@(x,y) x./y(1,:),tcEaCellEaCyc,...
    repmat(meanRespEaCell,nCycles,1),'unif',0);
meanRespEaCycEaExpt = cell2mat(cellfun(@(x) mean(x,2),meanRespEaCellNorm,'unif',0));
errRespEaCycEaExpt = cell2mat(cellfun(@(x) ste(x,2),meanRespEaCellNorm,'unif',0));
meanTcEaCycEaExpt = cellfun(@(x) mean(x,2),tcEaCellEaCycNorm,'unif',0);
errTcEaCycEaExpt = cellfun(@(x) ste(x,2),tcEaCellEaCycNorm,'unif',0);
nCellsEaExpt = cellfun(@(x) size(x,2),meanRespEaCell(1,:),'unif',0);
exptNameNCells = cellfun(@(x,y) sprintf('%s (%s)',x,num2str(y)),exptName,...
    nCellsEaExpt,'unif',0);
nCellsCutoff = 3;
expt2Use4Analysis = cell2mat(nCellsEaExpt) > nCellsCutoff;
meanRespEaCyc = nanmean(meanRespEaCycEaExpt(:,expt2Use4Analysis),2);
errRespEaCyc = ste(meanRespEaCycEaExpt(:,expt2Use4Analysis),2);

meanTcEaCyc = nan(nBaselineFr+cycLengthFr,nCycles);
errTcEaCyc = nan(nBaselineFr+cycLengthFr,nCycles);
for icyc = 1:nCycles
    meanTcEaCyc(:,icyc) = mean(cell2mat(meanTcEaCycEaExpt...
        (icyc,expt2Use4Analysis)),2);
    errTcEaCyc(:,icyc) = ste(cell2mat(meanTcEaCycEaExpt...
        (icyc,expt2Use4Analysis)),2);
end
%%
resp_lim = [0 1.1];
setFigParams4Print('landscape')
FSAdaptFig = figure;
subplot 224
h = errorbar(1:nCycles,meanRespEaCyc,errRespEaCyc,'ko');
h.MarkerFaceColor = [1 1 1];
figXAxis([],'cycle number',[0 nCycles+1],1:nCycles,1:nCycles)
figYAxis([],'norm dF/F',resp_lim);
figAxForm([])
title({'response to each stim';'cells responsive to 1st baseline stim';'visual trials only'})



figure(FSAdaptFig);
subplot 222
for iexp = 1:nexp
    if ~expt2Use4Analysis(iexp)
        continue
    end
    hold on
    h = errorbar(1:nCycles,meanRespEaCycEaExpt(:,iexp),errRespEaCycEaExpt(:,iexp),'o-');
    h.Color = exptColors(iexp,:);
    h.MarkerFaceColor = exptColors(iexp,:);
    h.MarkerEdgeColor = [1 1 1];
end
resp_lim = [-6 6];
figXAxis([],'cycle number',[0 nCycles+1],1:nCycles,1:nCycles)
figYAxis([],'dF/F',resp_lim);
figAxForm([])
title({'cycle response for each expt'})

subplot 221
resp_lim = [-0.03 0.12];
[respPatchCoordinates(:,1),respPatchCoordinates(:,2)] = respPatchCoord(...
    responseWindow-nBaselineFr,resp_lim);
p = patch(respPatchCoordinates(:,1), respPatchCoordinates(:,2), [0.75 0.75 0.75]);
p.EdgeColor = [1 1 1];
hold on
[baselinePatchCoordinates(:,1),baselinePatchCoordinates(:,2)] = ...
    respPatchCoord(baselineWindow-nBaselineFr,resp_lim);
p = patch(baselinePatchCoordinates(:,1), baselinePatchCoordinates(:,2), [0.9 0.9 0.9]);
p.EdgeColor = [1 1 1];
hold on
clear linesForLegend
for iexp = 1:nexp
    if ~expt2Use4Analysis(iexp)
        linesForLegend(iexp) = hline(0,'k-');
        continue
    end
    h = plot(timestamp_1cyc,meanTCacrossCells_1cyc(:,iexp));
    h.Color = exptColors(iexp,:);
    linesForLegend(iexp) = h;
    hold on
end
vline([0 onTimeFr],'k--')
figXAxis([],'time (frames)',plotTimeRange_1cyc)
figYAxis([],'dF/F',resp_lim)
figAxForm([])
title({'response to 1st stim, all responsive cells';'dotted line vis stim on/off'})
legend(linesForLegend(expt2Use4Analysis),exptNameNCells(expt2Use4Analysis),'location','northeastoutside')

subplot 223
cycColors = brewermap(nCycles+3,'Blues');
cycColors = cycColors(4:nCycles+3,:);
resp_lim = [-0.6 2.5];
[respPatchCoordinates(:,1),respPatchCoordinates(:,2)] = respPatchCoord(...
    responseWindow-nBaselineFr,resp_lim);
p = patch(respPatchCoordinates(:,1), respPatchCoordinates(:,2), [0.75 0.75 0.75]);
p.EdgeColor = [1 1 1];
hold on
[baselinePatchCoordinates(:,1),baselinePatchCoordinates(:,2)] = ...
    respPatchCoord(baselineWindow-nBaselineFr,resp_lim);
p = patch(baselinePatchCoordinates(:,1), baselinePatchCoordinates(:,2), [0.9 0.9 0.9]);
p.EdgeColor = [1 1 1];
hold on
clear linesForLegend
for icyc = 1:nCycles
    hold on
    h = plot(timestamp_1cyc,meanTcEaCyc(:,icyc));    
    h.Color = cycColors(icyc,:);    
    linesForLegend(icyc) = h;
end
legend(linesForLegend,strsplit(num2str(1:5)),'location','northeastoutside')
vline([0 onTimeFr],'k--')
figXAxis([],'time (frames)',plotTimeRange_1cyc)
figYAxis([],'dF/F',resp_lim)
figAxForm([])
title('response aligned to each cycle')
print([fnout 'eaCycleAdaptation'],'-dpdf','-fillpage');
end
    