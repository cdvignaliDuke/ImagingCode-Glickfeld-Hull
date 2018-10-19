clear all
close all
ds = 'FSAV_attentionV1';
cellsOrDendrites = 1;
%%
rc = behavConstsAV;
imgParams_FSAV
dataGroup = ds;
eval(dataGroup)

titleStr = ds;
str = unique({expt.SubNum});
mouse_str = ['i' strjoin(str,'_i')];

if cellsOrDendrites == 1
    load(fullfile(rc.caOutputDir,dataGroup,[mouse_str '_trOutcomeStruct_cells' ds(5:end) '.mat']));    
elseif cellsOrDendrites == 2
    load(fullfile(rc.caOutputDir,dataGroup,[mouse_str '_trOutcomeStruct_dendrites' ds(5:end) '.mat']));    
end        
if cellsOrDendrites == 1
    fnout = fullfile(rc.caOutputDir, dataGroup,[titleStr '_trOut_']); 
elseif cellsOrDendrites == 2
    fnout = fullfile(rc.caOutputDir, dataGroup,[titleStr 'trOut_den_']); 
end

%%
nmice = size(mouse,2);
minTrN = 5;
fitReliabilityCutoff = 30;
if strcmp(ds,'FSAV_V1_100ms_naive')
    visRT = [0 570];
    audRT = [0 570];
else
    visRT = [250 570];
    audRT = [150 450];
end
nFrMonitorDelay = 3; %using cLeverDown
nFrMonitorDelay_target = 2; %using cTargetOn
nFrVisDelay = 3; 
nFrRespwin = 3;
respCutoffDff = 0.002;
respAlpha = 0.01;
frameRateHz = 30;
%%
colors_FACR = {'m','b'};
colors_HM = {'k','r'};
lines_VA = {'-',':'};
titles_VA = {'Visual','Auditory'};
visualTrials = 1;
auditoryTrials = 2;
hitTrials = 1;
missTrials = 2;
nav = 2;
firstAlign = 1;
faAlign = 2;
crAlign = 3;
targetAlign = 4;
nalign = 4;
faInd = 1;
crInd = 2;
hInd = 1;
mInd = 2;

nexp = 0;
for imouse = 1:size(mouse,2)
    nexp = nexp+size(mouse(imouse).expt,2);
end
%%
nBaselineFr = mouse(1).expt(1).info.preAlignFrames;
cycTimeFr = mouse(1).expt(1).info.cycTimeFrames;

respwin = (nBaselineFr+nFrVisDelay+nFrMonitorDelay):...
    (nBaselineFr+nFrVisDelay+nFrMonitorDelay+nFrRespwin - 1);
basewin = (nBaselineFr-nFrRespwin+nFrMonitorDelay+1):...
    (nBaselineFr+nFrMonitorDelay);

respwin_target = (nBaselineFr+nFrVisDelay+nFrMonitorDelay_target):...
    (nBaselineFr+nFrVisDelay+nFrMonitorDelay_target+nFrRespwin - 1);
basewin_target = (nBaselineFr-nFrRespwin+nFrMonitorDelay_target+1):...
    (nBaselineFr+nFrMonitorDelay_target);
%% compare FA and CR
load(fullfile(rc.ashleyAnalysis,'FSAV Summaries',ds,...
    [ds(6:end) '_startAlign_anticipationFirstRespCells']))
if strcmp(ds,'FSAV_attentionV1')
respTC = cell(2,nav);
firstRespCellsInd_FACR = [];
firstRespCells_FACR = [];
firstRespCellsInd = [];
for imouse = 1:nmice
    for iexp = 1:size(mouse(imouse).expt,2)

        d = mouse(imouse).expt(iexp);

        firstRespTCs = cat(3,d.av(visualTrials).align(firstAlign).respTC,...
                            d.av(auditoryTrials).align(firstAlign).respTC);

        firstResponsiveCells = logical(...
                        ttest(squeeze(mean(firstRespTCs(respwin,:,:),1)),...
                        squeeze(mean(firstRespTCs(basewin,:,:),1)),...
                        'dim',2,'tail','right','alpha',respAlpha))';
        respCutoffMet = squeeze(mean(mean(firstRespTCs(respwin,:,:),1),3))...
                        > respCutoffDff;
                    
        for ialign = 1:2
            ntr_v = size(d.av(visualTrials).align(ialign+1).respTC,3);
            ntr_a = size(d.av(auditoryTrials).align(ialign+1).respTC,3);
            if ialign == 1
                trInd_v = d.av(visualTrials).align(ialign+1).reactTime > visRT(1) & ...
                    d.av(visualTrials).align(ialign+1).reactTime < visRT(2);
                trInd_a = d.av(auditoryTrials).align(ialign+1).reactTime > audRT(1) & ...
                    d.av(auditoryTrials).align(ialign+1).reactTime < audRT(2);
                ntr_v = sum(trInd_v);
                ntr_a = sum(trInd_a);
            end
            if ntr_v >= minTrN && ntr_a >= minTrN
                rv = mean(d.av(visualTrials).align(ialign+1).respTC(:,:,trInd_v),3) - ...
                    mean(mean(d.av(visualTrials).align(ialign+1).respTC(basewin,:,trInd_v),3),1);
                ra = mean(d.av(auditoryTrials).align(ialign+1).respTC(:,:,trInd_a),3) - ...
                    mean(mean(d.av(auditoryTrials).align(ialign+1).respTC(basewin,:,trInd_a),3),1);
                ind = firstResponsiveCells & respCutoffMet;
                if ialign == 1
                nc = size(rv,2);
                    firstRespCells_FACR = cat(2,firstRespCells_FACR,true(1,nc));
                end
            else
                rv = [];
                ra = [];
                ind = [];
                if ialign == 1
                    nc = size(d.av(visualTrials).align(1).respTC,2);
                    firstRespCells_FACR = cat(2,firstRespCells_FACR,false(1,nc));
                end
            end
            if ialign == 1
                firstRespCellsInd_FACR = cat(2,firstRespCellsInd_FACR,ind);
                firstRespCellsInd = cat(2,firstRespCellsInd,...
                    firstResponsiveCells & respCutoffMet);
            end
            respTC{ialign,visualTrials} = cat(2,respTC{ialign,visualTrials},rv);
            respTC{ialign,auditoryTrials} = cat(2,respTC{ialign,auditoryTrials},ra);
        end
    end
end
firstRespCellsInd = logical(firstRespCellsInd);
firstRespCellsInd_FACR = logical(firstRespCellsInd_FACR);
firstRespCells4TC = firstRespCells(logical(firstRespCells_FACR));
% plot
[nTCFr,nCells] = size(respTC{1,1});
tt_all = ((1:nTCFr)-(nBaselineFr+nFrMonitorDelay+1))./frameRateHz;
tt_basewin = ([basewin(1) basewin(end)]-(nBaselineFr+nFrMonitorDelay+1))...
    ./frameRateHz;
tt_respwin = ([respwin(1) respwin(end)]-(nBaselineFr+nFrMonitorDelay+1))...
    ./frameRateHz;

ttWinInd = find(tt_all < -0.5,1,'last'):find(tt_all > 0.5,1,'first');
tt = tt_all(ttWinInd);
tt_label = -0.5:0.25:0.5;
respTC_lim = [-0.004 0.06];
resp_lim = [-0.1 0.16];

setFigParams4Print('landscape')
figure
suptitle(sprintf('First Stim Responsive Cells (%s)',num2str(sum(firstRespCells4TC))))
for iav = 1:nav
    leg = [];
    subplot(2,2,iav)
    x = mean(respTC{faInd,iav}(ttWinInd,firstRespCells4TC),2);
    xerr = ste(respTC{faInd,iav}(ttWinInd,firstRespCells4TC),2);
    h = shadedErrorBar(tt,x,xerr,colors_FACR{faInd});
    h.mainLine.LineStyle = lines_VA{iav};
    for i = 1:2
        h.edge(i).Color = h.patch.FaceColor;
    end
    leg(1) = h.mainLine;
    hold on
    x = mean(respTC{crInd,iav}(ttWinInd,firstRespCells4TC),2);
    xerr = ste(respTC{crInd,iav}(ttWinInd,firstRespCells4TC),2);
    h = shadedErrorBar(tt,x,xerr,colors_FACR{crInd});
    h.mainLine.LineStyle = lines_VA{iav};
    for i = 1:2
        h.edge(i).Color = h.patch.FaceColor;
    end
    leg(2) = h.mainLine;
    figXAxis([],'Time from Vis Stim (s)',[tt(1) tt(end)],tt_label,tt_label)
    figYAxis([],'dF/F',respTC_lim)
    figAxForm
    vline(tt_basewin,'k:')
    vline(tt_respwin,'k:')
    hline(0,'k:')
    title(sprintf('%s Trials',titles_VA{iav}))
    legend(leg,{'FA','CR'},'location','northwest')
    
    
    subplot(2,2,iav+2)
    x = nanmean(respTC{faInd,iav}(respwin,firstRespCells4TC),1);
    y = nanmean(respTC{crInd,iav}(respwin,firstRespCells4TC),1);
    h = plot(x,y,'o');
    h.Color = [0.5 0.5 0.5];
    h.MarkerFaceColor = [0.5 0.5 0.5];
    hold on
    h = errorbar(mean(x),mean(y),ste(y,2),ste(y,2),ste(x,2),ste(x,2),'ro');
    h.MarkerFaceColor = 'r';
    plot(resp_lim,resp_lim,'k--')
    [~,p] = ttest(x,y,'alpha',0.05);
    p = round(p,2,'significant');
    figXAxis([],'False Alarm',resp_lim)
    figYAxis([],'Correct Reject',resp_lim)
    figAxForm
    title(sprintf('All Target Responsive, p = %s, n = %s',...
        num2str(p),num2str(sum(firstRespCells4TC))))
end
print([fnout 'respTC_FACR'],'-dpdf','-fillpage')

%% compare RT    
vrts = cell(1,size(mouse,2));
arts = cell(1,size(mouse,2));
for imouse = 1:nmice
    for iexp = 1:size(mouse(imouse).expt,2)
        d =  mouse(imouse).expt(iexp);
        h = strcmp(d.av(visualTrials).align(targetAlign).outcome,'success');
        v = d.av(visualTrials).align(targetAlign).reactTime(h);
        h = strcmp(d.av(auditoryTrials).align(targetAlign).outcome,'success');
        a = d.av(auditoryTrials).align(targetAlign).reactTime(h);
        vrts{imouse} = cat(2,vrts{imouse},mean(v));
        arts{imouse} = cat(2,arts{imouse},mean(a));
    end
end
 
figure
x = ones(1,size(mouse,2));
y = cellfun(@mean,vrts);
h = scatter(x,y,'ko');
h.MarkerFaceColor = 'k';
hold on
h = errorbar(1,mean(y),ste(y,2),'ko');
h.MarkerFaceColor = [1 1 1];
x = ones(1,size(mouse,2))*2;
y = cellfun(@mean,arts);
h = scatter(x,y,'ko');
h.MarkerFaceColor = 'k';
h = errorbar(2,mean(y),ste(y,2),'ko');
h.MarkerFaceColor = [1 1 1];
figXAxis([],'',[0 3],[1 2],{'Visual','Auditory'})
figYAxis([],'RT (ms)',[])
figAxForm
title('Mean Reaction Time, Each Mouse')
       
print([fnout 'RTeaMouse'],'-dpdf','-fillpage')
end
%% compare H and M    
oris = [];
amps = [];
for imouse = 1:nmice
    for iexp = 1:size(mouse(imouse).expt,2)
        oris = cat(2,oris,mouse(imouse).expt(iexp).info.visTargets);
        amps = cat(2,amps,mouse(imouse).expt(iexp).info.audTargets);
    end
end
oris = unique(round(oris,2,'significant'));
oris = oris(2:end);
nori = length(oris);
amps = unique(round(amps,2,'significant'));
amps = amps(2:end);
namp = length(amps);

respTCmatched = cell(2,nav);
respTCmatchedEaMs = cell(2,nav,nmice);
targetRespCellsInd = [];
nTargetRespCellsExpt = nan(1,nexp);
nCellsExpt = nan(1,nexp);
exptName = cell(1,nexp);
targetRespCellsEaMs = cell(1,nmice);
detectSelectivity_HM = [];
HM_ttest = [];
matchTargetSzMs = zeros(nav,nmice);
for imouse = 1:nmice
    hMatchMsVis = [];
    mMatchMsVis = [];
    hMatchMsAud = [];
    mMatchMsAud = [];
    targetSz = zeros(nav,size(mouse(imouse).expt,2));
    for iexp = 1:size(mouse(imouse).expt,2)        
        if imouse == 1 && iexp == 1
            exptN = 1;
        else
            exptN = exptN+1;
        end
        exptName{exptN} = [mouse(imouse).mouse_name '-' mouse(imouse).expt(iexp).date];
        for iav = 1:nav
            d = mouse(imouse).expt(iexp).av(iav).align(targetAlign);
            rt = d.reactTime;
            if iav == 1
                h = strcmp(d.outcome,'success') & rt > visRT(1) & rt < visRT(2);
                m = strcmp(d.outcome,'ignore') | rt > visRT(2);
                tOri = round(d.ori,2,'significant');
                orientations = round(unique(tOri),2,'significant');
                hitResp = cell(1,length(orientations));
                missResp = cell(1,length(orientations));
                targetRespOri = cell(1,length(orientations));
                for i = 1:length(orientations)
                   oriInd = tOri == orientations(i);  
                   if isempty(oriInd)
                       hitResp{i} = nan(size(d.respTC,1),size(d.respTC,2),0);
                       missResp{i} = nan(size(d.respTC,1),size(d.respTC,2),0);
                   else
                       hitResp{i} = d.respTC(:,:,oriInd & h);
                       missResp{i} = d.respTC(:,:,oriInd & m);
                       if size(hitResp{i},3) + size(missResp{i},3) >= minTrN
                           targetRespOri{i} = cat(3,hitResp{i},missResp{i});
                       else
                           targetRespOri{i} = nan(size(d.respTC,1),size(d.respTC,2),0);
                       end
                   end           
                end
                targetResp = d.respTC(:,:,h | m);

                nHit = cellfun(@(x) size(x,3),hitResp);
                nMiss = cellfun(@(x) size(x,3),missResp);

                hitRespMatch = [];
                missRespMatch = [];
                vistargets = [];
                for i = 1:length(orientations)
                    if nHit(i) > nMiss(i)
                        ind = randsample(nHit(i),nMiss(i));
                        hitRespMatch = cat(3,hitResp{i}(:,:,ind),hitRespMatch);
                        missRespMatch = cat(3,missResp{i},missRespMatch);                
                    else
                        ind = randsample(nMiss(i),nHit(i));
                        hitRespMatch = cat(3,hitResp{i},hitRespMatch);
                        missRespMatch = cat(3,missResp{i}(:,:,ind),missRespMatch); 
                    end
                    if nHit(i) > 0 && nMiss(i) > 0
                        vistargets = [vistargets orientations(i)];
                    end
                end
                respTCmatched{hitTrials,iav} = cat(2,respTCmatched{hitTrials,iav},...
                    mean(hitRespMatch,3) - ...
                    mean(mean(hitRespMatch(basewin_target,:,:),1),3));
                respTCmatched{missTrials,iav} = cat(2,respTCmatched{missTrials,iav},...
                    mean(missRespMatch,3) - ...
                    mean(mean(missRespMatch(basewin_target,:,:),1),3));
                hMatchMsVis = cat(2,hMatchMsVis,mean(hitRespMatch,3) - ...
                    mean(mean(hitRespMatch(basewin_target,:,:),1),3));
                mMatchMsVis = cat(2,mMatchMsVis,mean(missRespMatch,3) - ...
                    mean(mean(missRespMatch(basewin_target,:,:),1),3));
                targetSz(iav,iexp) = nanmean(vistargets);
                
                hresp = squeeze(mean(hitRespMatch(respwin_target,:,:)))';
                mresp = squeeze(mean(missRespMatch(respwin_target,:,:)))';
                si = getSelectivityIndex(hresp,mresp);
                hmresptest = ttest2(hresp,mresp,'dim',1,'tail','both');
                if all(isnan(hmresptest))
                    hmresptest = false(1,length(hmresptest));
                end
                detectSelectivity_HM = cat(2,detectSelectivity_HM,si);
                HM_ttest = cat(2,HM_ttest,hmresptest);
                    

                targetRespAlpha = 0.05./(length(orientations)-1);
                targetsInd = cellfun(@(x) ~isempty(x),targetRespOri);
                targetTtest = cellfun(@(x) ttest(...
                    squeeze(mean(x(respwin_target,:,:),1)),...
                    squeeze(mean(x(basewin_target,:,:),1)),...
                    'dim',2,'tail','right','alpha',targetRespAlpha),...
                    targetRespOri(targetsInd),'unif',0);
                if ~any(targetsInd)
                    targetRespCellsOri = false(size(targetResp,2),1);
                else
                    targetRespCellsOri = sum(cell2mat(targetTtest),2) > 1;
                end
                targetRespCells = ttest(...
                    squeeze(mean(targetResp(respwin_target,:,:),1)),...
                    squeeze(mean(targetResp(basewin_target,:,:),1)),...
                    'dim',2,'tail','right');
                targetRespCellsInd = cat(2,targetRespCellsInd,...
                    (targetRespCellsOri | logical(targetRespCells))');
                targetRespCellsEaMs{imouse} = cat(2,targetRespCellsEaMs{imouse},...
                    (targetRespCellsOri | logical(targetRespCells))');
                nTargetRespCellsExpt(exptN) = sum((targetRespCellsOri | logical(targetRespCells)));
                nCellsExpt(exptN) = length(targetRespCellsOri);
            elseif iav == 2
                h = strcmp(d.outcome,'success') & rt > audRT(1) & rt < audRT(2);
                m = strcmp(d.outcome,'ignore') | rt > audRT(2);
                tAmp = round(d.amp,2,'significant');
                amplitudes = round(unique(tAmp),2,'significant');
                hitResp = cell(1,length(amplitudes));
                missResp = cell(1,length(amplitudes));
                for i = 1:length(amplitudes)
                   ampInd = tAmp == amplitudes(i);  
                   if isempty(ampInd)
                       hitResp{i} = nan(size(d.respTC,1),size(d.respTC,2),0);
                       missResp{i} = nan(size(d.respTC,1),size(d.respTC,2),0);
                   else
                       hitResp{i} = d.respTC(:,:,ampInd & h);
                       missResp{i} = d.respTC(:,:,ampInd & m);
                   end           
                end
                
                nHit = cellfun(@(x) size(x,3),hitResp);
                nMiss = cellfun(@(x) size(x,3),missResp);

                hitRespMatch = [];
                missRespMatch = [];
                audTargets = [];
                for i = 1:length(amplitudes)
                    if nHit(i) > nMiss(i)
                        ind = randsample(nHit(i),nMiss(i));
                        hitRespMatch = cat(3,hitResp{i}(:,:,ind),hitRespMatch);
                        missRespMatch = cat(3,missResp{i},missRespMatch);                
                    else
                        ind = randsample(nMiss(i),nHit(i));
                        hitRespMatch = cat(3,hitResp{i},hitRespMatch);
                        missRespMatch = cat(3,missResp{i}(:,:,ind),missRespMatch); 
                    end
                    if nHit(i) > 0 && nMiss(i) > 0
                        audTargets = [audTargets amplitudes(i)];
                    end
                end
                respTCmatched{hitTrials,iav} = cat(2,respTCmatched{hitTrials,iav},...
                    mean(hitRespMatch,3) - ...
                    mean(mean(hitRespMatch(basewin_target,:,:),1),3));
                respTCmatched{missTrials,iav} = cat(2,respTCmatched{missTrials,iav},...
                    mean(missRespMatch,3) - ...
                    mean(mean(missRespMatch(basewin_target,:,:),1),3));
                hMatchMsAud = cat(2,hMatchMsAud,mean(hitRespMatch,3) - ...
                    mean(mean(hitRespMatch(basewin_target,:,:),1),3));
                mMatchMsAud = cat(2,mMatchMsAud,mean(missRespMatch,3) - ...
                    mean(mean(missRespMatch(basewin_target,:,:),1),3));
                targetSz(iav,iexp) = nanmean(audTargets);
                
            end
        end
    end
    respTCmatchedEaMs{hInd,visualTrials,imouse} = hMatchMsVis;
    respTCmatchedEaMs{mInd,visualTrials,imouse} = mMatchMsVis;
    respTCmatchedEaMs{hInd,auditoryTrials,imouse} = hMatchMsAud;
    respTCmatchedEaMs{mInd,auditoryTrials,imouse} = mMatchMsAud;
    matchTargetSzMs(:,imouse) = mean(targetSz,2);
end
targetRespCellsInd = logical(targetRespCellsInd);
HM_ttest = logical(HM_ttest);

% plot
[nTCFr,nCells] = size(respTCmatched{1,1});
tt_all = ((1:nTCFr)-(nBaselineFr+nFrMonitorDelay_target+1))./frameRateHz;
tt_basewin = ([basewin_target(1) basewin_target(end)]-(nBaselineFr+nFrMonitorDelay_target+1))...
    ./frameRateHz;
tt_respwin = ([respwin_target(1) respwin_target(end)]-(nBaselineFr+nFrMonitorDelay_target+1))...
    ./frameRateHz;

ttWinInd = find(tt_all < -0.5,1,'last'):find(tt_all > 0.5,1,'first');
tt = tt_all(ttWinInd);
tt_label = -0.5:0.25:0.5;

figure
suptitle('Matched Target Response')
for iav = 1:nav
    subplot(2,2,iav)
    leg = [];
    if iav == 1
        ind = targetRespCellsInd;
    elseif iav == 2
        ind = firstRespCellsInd;
    end
    x = nanmean(respTCmatched{hInd,iav}(ttWinInd,ind),2);
    xerr = ste(respTCmatched{hInd,iav}(ttWinInd,ind),2);
    h = shadedErrorBar(tt,x,xerr,colors_HM{hInd});
    for i = 1:2
        h.edge(i).Color = h.patch.FaceColor;
    end
    leg(1) = h.mainLine;
    hold on
    x = nanmean(respTCmatched{mInd,iav}(ttWinInd,ind),2);
    xerr = ste(respTCmatched{mInd,iav}(ttWinInd,ind),2);
    h = shadedErrorBar(tt,x,xerr,colors_HM{crInd});
    for i = 1:2
        h.edge(i).Color = h.patch.FaceColor;
    end
    leg(2) = h.mainLine;
    figXAxis([],'Time from Target Stim (s)',[tt(1) tt(end)],tt_label,tt_label)
    figYAxis([],'dF/F',respTC_lim)
    figAxForm
    vline(tt_basewin,'k:')
    vline(tt_respwin,'k:')
    if iav == 1
        leg(3) = vline((visRT(1)/1000) - (nFrMonitorDelay_target/frameRateHz),'b-');
        title({sprintf('Target Responsive Cells (%s)',...
            num2str(sum(targetRespCellsInd)));...
            sprintf('%s Trials',titles_VA{iav})});
    elseif iav == 2
        leg(3) = vline((audRT(1)/1000) - (nFrMonitorDelay_target/frameRateHz),'b-');
        title({sprintf('First Stim Responsive Cells (%s)',...
            num2str(sum(firstRespCellsInd)));...
            sprintf('%s Trials',titles_VA{iav})});
    end
    hline(0,'k:')    
    legend(leg,{'H','M','min RT'},'location','northwest')
    
    
    subplot(2,2,iav+2)
    x = nanmean(respTCmatched{hInd,iav}(respwin_target,ind),1);
    xerr = ste(respTCmatched{hInd,iav}(respwin_target,ind),1);
    y = nanmean(respTCmatched{mInd,iav}(respwin_target,ind),1);
    yerr = ste(respTCmatched{mInd,iav}(respwin_target,ind),1);
    h = plot(x,y,'o');
    h.Color = [0.5 0.5 0.5];
    h.MarkerFaceColor = [0.5 0.5 0.5];
    hold on
    h = errorbar(mean(x),mean(y),ste(y,2),ste(y,2),ste(x,2),ste(x,2),'ro');
    h.MarkerFaceColor = 'r';
    plot(resp_lim,resp_lim,'k--')
    [~,p] = ttest(x,y,'alpha',0.05);
    p = round(p,2,'significant');
    figXAxis([],'Hit',resp_lim)
    figYAxis([],'Miss',resp_lim)
    figAxForm
    title(sprintf('All Target Responsive, p = %s, n = %s',...
        num2str(p),num2str(sum(ind))))
    
end

print([fnout 'respTC_HM'],'-dpdf','-fillpage')

figure;
[nrow,ncol] = optimizeSubplotDim(nmice);
for imouse = 1:nmice
    subplot(nrow,ncol,imouse);
    y = nanmean(respTCmatchedEaMs{hInd,visualTrials,imouse}...
        (ttWinInd,logical(targetRespCellsEaMs{imouse})),2);
    yerr = ste(respTCmatchedEaMs{hInd,visualTrials,imouse}...
        (ttWinInd,logical(targetRespCellsEaMs{imouse})),2);
    h = shadedErrorBar(tt,y,yerr,colors_HM{hInd});
    for i = 1:2
        h.edge(i).Color = h.patch.FaceColor;
    end
    hold on
    y = nanmean(respTCmatchedEaMs{mInd,visualTrials,imouse}...
        (ttWinInd,logical(targetRespCellsEaMs{imouse})),2);
    yerr = ste(respTCmatchedEaMs{mInd,visualTrials,imouse}...
        (ttWinInd,logical(targetRespCellsEaMs{imouse})),2);
    h = shadedErrorBar(tt,y,yerr,colors_HM{mInd});
    for i = 1:2
        h.edge(i).Color = h.patch.FaceColor;
    end
    figXAxis([],'Time from Target Stim (s)',[tt(1) tt(end)],tt_label,tt_label)
    figYAxis([],'dF/F',respTC_lim)
    figAxForm
    vline(tt_basewin,'k:')
    vline(tt_respwin,'k:')
    title({sprintf('N = %s',...
        num2str(sum(targetRespCellsEaMs{imouse})));...
        sprintf('Mouse %s',mouse(imouse).mouse_name)});
end
figure;
[nrow,ncol] = optimizeSubplotDim(nmice);
for imouse = 1:nmice
    subplot(nrow,ncol,imouse);
    x = nanmean(respTCmatchedEaMs{hInd,visualTrials,imouse}...
        (respwin_target,logical(targetRespCellsEaMs{imouse})),1);
    y = nanmean(respTCmatchedEaMs{mInd,visualTrials,imouse}...
        (respwin_target,logical(targetRespCellsEaMs{imouse})),1);
    h = scatter(x,y,'o');
    
    h.MarkerFaceColor = 'k';
    h.MarkerEdgeColor = [1 1 1];
    hold on
    xerr = ste(x,2);
    yerr = ste(y,2);
    h = errorbar(mean(x),mean(y),yerr,yerr,xerr,xerr,'ro');
    h.MarkerFaceColor = 'r';
    plot(resp_lim,resp_lim,'k--')
    [~,p] = ttest(x,y,'alpha',0.05);
    p = round(p,2,'significant');
    figXAxis([],'Hit',resp_lim)
    figYAxis([],'Miss',resp_lim)
    figAxForm
    title(sprintf('All Target Responsive, p = %s, n = %s',...
        num2str(p),num2str(sum(ind))))
    title({sprintf('N = %s, p = %s',...
        num2str(sum(targetRespCellsEaMs{imouse})),...
        num2str(p));...
        sprintf('Mouse %s',mouse(imouse).mouse_name)});
end
%%
targetRespMatched = cellfun(@(x) mean(x(respwin_target,:)),respTCmatched,'unif',0);

%% ori tuning 
orientations = [0 45 90 135];
nOrientations = length(orientations);
oriPrefBins = [0 22.5 67.5 112.5 157.5 180];

oriFitPeak = [];
isTuned = [];
for imouse = 1:size(mouse,2)
    for iexp = 1:size(mouse(imouse).expt,2)
        d = mouse(imouse).expt(iexp).oriTuning;
        [~,fitPeak] = max(d.oriFit);
        oriFitPeak = cat(2,oriFitPeak,fitPeak-1);
        isTuned = cat(2,isTuned,d.oriFitReliability < fitReliabilityCutoff);
    end
end
isTuned = logical(isTuned);

nbins = 16;
figure
h = histogram(oriFitPeak(isTuned & firstRespCellsInd),nbins);
h.EdgeColor = 'none';
hold on
h = histogram(oriFitPeak(isTuned & targetRespCellsInd),nbins);
h.EdgeColor = 'none';
L = legend({'First Stim','Target'});
title(L,'Cell Type')
figXAxis([],'Pref Orientations (deg)',[-1 181]);
figYAxis([],'N Cells',[])
figAxForm
title('Orientation Preference of Tuned and Responsive Cells')

[~,~,oriPref] = histcounts(oriFitPeak,oriPrefBins);
oriPref(oriPref == 5) = 1;
oriPref(~isTuned) = 0;

% plot
resp_lim = [-0.06 0.28];


figure
suptitle(sprintf('Tuned, Target Responsive Neurons, Target Aligned (%s)',...
    num2str(sum(isTuned & targetRespCellsInd))))
for i = 1:nOrientations
    subplot (2,2,i)
    ind = oriPref == i & targetRespCellsInd;
    x = targetRespMatched{hitTrials,visualTrials}...
        (ind);
    xerr = ste(x,2);
    y = targetRespMatched{missTrials,visualTrials}...
        (ind);
    yerr = ste(y,2);
    h = scatter(x,y,'o');
    h.MarkerFaceColor = 'k';
    h.MarkerEdgeColor = [1 1 1];
    hold on
    h = errorbar(mean(x),mean(y),yerr,yerr,xerr,xerr,'ro');
    h.MarkerFaceColor = 'r';
    plot(resp_lim,resp_lim,'k--')
    [~,p] = ttest(x,y,'alpha',0.05/nOrientations);
    p = round(p,2,'significant');
    figXAxis([],'Hit',resp_lim)
    figYAxis([],'Miss',resp_lim)
    figAxForm
    title(sprintf('Pref %s deg, p = %s, n = %s',...
        num2str(orientations(i)),num2str(p),num2str(sum(ind))))
end

print([fnout 'HM_scat_oriPref'],'-dpdf','-fillpage')


figure
subplot 221
ind = targetRespCellsInd;
x = targetRespMatched{hInd,visualTrials}(ind);
y = targetRespMatched{mInd,visualTrials}(ind);
h = scatter(x,y,'o');
h.MarkerFaceColor = 'k';
h.MarkerEdgeColor = [1 1 1];
hold on
xerr = ste(x,2);
yerr = ste(y,2);
h = errorbar(mean(x),mean(y),yerr,yerr,xerr,xerr,'ro');
h.MarkerFaceColor = 'r';
plot(resp_lim,resp_lim,'k--')
[~,p] = ttest(x,y,'alpha',0.05);
p = round(p,2,'significant');
figXAxis([],'Hit',resp_lim)
figYAxis([],'Miss',resp_lim)
figAxForm
title(sprintf('All Target Responsive, p = %s, n = %s',...
    num2str(p),num2str(sum(ind))))

subplot 222
ind = targetRespCellsInd & isTuned;
x = targetRespMatched{hInd,visualTrials}(ind);
y = targetRespMatched{mInd,visualTrials}(ind);
h = scatter(x,y,'o');
h.MarkerFaceColor = 'k';
h.MarkerEdgeColor = [1 1 1];
hold on
xerr = ste(x,2);
yerr = ste(y,2);
h = errorbar(mean(x),mean(y),yerr,yerr,xerr,xerr,'ro');
h.MarkerFaceColor = 'r';
plot(resp_lim,resp_lim,'k--')
[~,p] = ttest(x,y,'alpha',0.05/nOrientations);
p = round(p,2,'significant');
figXAxis([],'Hit',resp_lim)
figYAxis([],'Miss',resp_lim)
figAxForm
title(sprintf('Tuned, p = %s, n = %s',...
    num2str(p),num2str(sum(ind))))

% subplot 223
% ind = isTuned & targetRespCellsInd;
% h = polarplot(deg2rad(oriFitPeak(ind)), detectSelectivity_HM(ind),'ko');
% title('Ori Tuning x Detect Selectivity')

print([fnout 'HM_scat_polar_tuned'],'-dpdf','-fillpage')

%%
targetOriBins = [0 40 60 90];%0:30:90;
nbin = 3;

binTargetTC = cell(2,nbin);
binTargetRespCells = cell(1,nbin);
for imouse = 1:size(mouse,2)
    for iexp = 1:size(mouse(imouse).expt,2)
        d = mouse(imouse).expt(iexp).av(visualTrials).align(targetAlign);
        rt = d.reactTime;
        h = strcmp(d.outcome,'success') & rt > visRT(1) & rt < visRT(2);
        m = strcmp(d.outcome,'ignore') | rt > visRT(2);
        tc = d.respTC;
        [~,~,binID] = histcounts(d.ori,targetOriBins);
        for i = 1:nbin
            if sum(binID == i & h) >= minTrN && sum(binID == i & m) >= minTrN
                hitResp = tc(:,:,binID == i & h);
                missResp = tc(:,:,binID == i & m);
                allResp = cat(3,hitResp,missResp);
            else
                hitResp = [];
                missResp = [];
                allResp = [];
            end
            binTargetTC{hitTrials,i} = cat(2,binTargetTC{hitTrials,i},...
                mean(hitResp,3));
            binTargetTC{missTrials,i} = cat(2,binTargetTC{missTrials,i},...
                mean(missResp,3));
            if ~isempty(allResp)
                targetTest = ttest(...
                    squeeze(mean(allResp(respwin_target,:,:),1)),...
                    squeeze(mean(allResp(basewin_target,:,:),1)),...
                    'dim',2,'tail','right','alpha',0.05./(nbins-1));
            else
                targetTest = [];
            end
            binTargetRespCells{i} = cat(1,binTargetRespCells{i},targetTest);
        end
    end
end
binTargetRespCells = cellfun(@logical,binTargetRespCells,'unif',0);

% plot
[nTCFr,~] = size(binTargetTC{1,1});
tt_all = ((1:nTCFr)-(nBaselineFr+nFrMonitorDelay_target+1))./frameRateHz;
tt_basewin = ([basewin_target(1) basewin_target(end)]-(nBaselineFr+nFrMonitorDelay_target+1))...
    ./frameRateHz;
tt_respwin = ([respwin_target(1) respwin_target(end)]-(nBaselineFr+nFrMonitorDelay_target+1))...
    ./frameRateHz;

ttWinInd = find(tt_all < -0.5,1,'last'):find(tt_all > 0.5,1,'first');
tt = tt_all(ttWinInd);
tt_label = -0.5:0.25:0.5;
respTC_lim = [-0.004 0.07];

figure
suptitle('Binned Target Response: cells are responsive within bin')
for ibin = 1:nbin
    subplot(1,nbin,ibin)
    leg = [];
    ind = binTargetRespCells{ibin};
    x = mean(binTargetTC{hInd,ibin}(ttWinInd,ind),2);
    xerr = ste(binTargetTC{hInd,ibin}(ttWinInd,ind),2);
    h = shadedErrorBar(tt,x,xerr,colors_HM{hInd});
    for i = 1:2
        h.edge(i).Color = h.patch.FaceColor;
    end
    leg(1) = h.mainLine;
    hold on
    x = mean(binTargetTC{mInd,ibin}(ttWinInd,ind),2);
    xerr = ste(binTargetTC{mInd,ibin}(ttWinInd,ind),2);
    h = shadedErrorBar(tt,x,xerr,colors_HM{crInd});
    for i = 1:2
        h.edge(i).Color = h.patch.FaceColor;
    end
    leg(2) = h.mainLine;
    figXAxis([],'Time from Target Stim (s)',[tt(1) tt(end)],tt_label,tt_label)
    figYAxis([],'dF/F',respTC_lim)
    figAxForm
    vline(tt_basewin,'k:')
    vline(tt_respwin,'k:')
    title(sprintf('%s-%s Targets, n = %s/%s)',...
        num2str(targetOriBins(ibin)),num2str(targetOriBins(ibin+1)),...
        num2str(sum(ind)),num2str(length(ind))));
    hline(0,'k:')  
    if ibin == 1
        legend(leg,{'H','M'},'location','northwest')
    end
end

print([fnout 'binnedTC_HM'],'-dpdf','-fillpage');

%% Press align H & M
minCycN = 6;
minTrN = 2;

startAlignTC = cell(3,nav);
cellsUsed = cell(1,nav);
for imouse = 1:size(mouse,2)
    for iexp = 1:size(mouse(imouse).expt,2)
        for iav = 1:nav
            d = mouse(imouse).expt(iexp).av(iav).align(firstAlign);
            cInd = d.nCycles >= minCycN;
            h = strcmp(d.outcome,'success');
            m = strcmp(d.outcome,'ignore');
            fa = strcmp(d.outcome,'failure');

            tc = d.respTC;
            
            if sum(h & cInd) > minTrN && sum(m & cInd) > minTrN && sum(fa & cInd) > minTrN
                startAlignTC{1,iav} = cat(2,startAlignTC{1,iav},...
                    mean(tc(:,:,h & cInd),3));
                startAlignTC{2,iav} = cat(2,startAlignTC{2,iav},...
                    mean(tc(:,:,m & cInd),3));
                startAlignTC{3,iav} = cat(2,startAlignTC{3,iav},...
                    mean(tc(:,:,fa & cInd),3));
                cellsUsed{iav} = cat(2,cellsUsed{iav},true(1,size(tc,2)));
            else
                cellsUsed{iav} = cat(2,cellsUsed{iav},false(1,size(tc,2)));
            end
        end
    end
end


[nTCFr,~] = size(startAlignTC{1,1});
tt_all = ((1:nTCFr)-(nBaselineFr+nFrMonitorDelay+1))./frameRateHz;

tt_cyc = (0:cycTimeFr:((minCycN-1)*cycTimeFr))./frameRateHz;

ttWinInd = find(tt_all < -0.5,1,'last'):find(tt_all > 1.8,1,'first');
tt = tt_all(ttWinInd);
tt_label = -0.5:0.25:2;
respTC_lim = [-0.004 0.2];

figure

for iav = 1:nav
    subplot(1,nav,iav)
    leg = [];
    ind = firstRespCellsInd(logical(cellsUsed{iav}));
    for iout = 1:3
        y = mean(startAlignTC{iout,iav}(ttWinInd,ind),2);
        yerr = ste(startAlignTC{iout,iav}(ttWinInd,ind),2);
        if iout < 3
            c = colors_HM{iout};
        else
            c = colors_FACR{1};
        end
        h = shadedErrorBar(tt,y,yerr,c);
        h.mainLine.LineStyle = lines_VA{iav};
        for i = 1:2
            h.edge(i).LineStyle = 'none';
        end
        leg(iout) = h.mainLine;
        hold on
    end
    figXAxis([],'Time from Target Stim (s)',[tt(1) tt(end)],tt_label,tt_label)
    figYAxis([],'dF/F',respTC_lim)
    figAxForm
    vline(tt_cyc,'k:')
    title({sprintf('First Stim Responsive Cells (%s)',...
            num2str(sum(ind)));...
            sprintf('%s Trials',titles_VA{iav})});
end

%% heatmaps and example cells and tuning
load(fullfile(rc.ashleyAnalysis,'FSAV Summaries',ds,...
    [ds(6:end) '_startAlign_anticipationHeatmapSorting']))
targetOriBins = [0 45 90];
nbin = 2;
targetTC = cell(1,nbin);
targetTC_sem = cell(1,nbin);
oriFits = [];
oriResp = [];
oriResp_sem = [];
for imouse = 1:size(mouse,2)
    for iexp = 1:size(mouse(imouse).expt,2)
        d = mouse(imouse).expt(iexp).av(visualTrials).align(targetAlign);
        rt = d.reactTime;
        h = strcmp(d.outcome,'success') & rt > visRT(1) & rt < visRT(2);
        m = strcmp(d.outcome,'ignore') | rt > visRT(2);
        tc = d.respTC;
        [~,~,binID] = histcounts(d.ori,targetOriBins);
        for i = 1:nbin
            if sum(binID == i & (h|m)) >= minTrN
                allResp = tc(:,:,binID == i & (h|m));
            else
                allResp = nan(size(tc,1),size(tc,2),1);
            end
            targetTC{1,i} = cat(2,targetTC{1,i},...
                mean(allResp,3));
            targetTC_sem{1,i} = cat(2,targetTC_sem{1,i},...
                ste(allResp,3));
        end
        
        d = mouse(imouse).expt(iexp).oriTuning;
        oriFits = cat(2,oriFits,d.oriFit);
        if size(d.oriResp,2) == 8
            or = d.oriResp(:,1:2:8);
            ors = d.oriRespSem(:,1:2:8);
            oriResp = cat(1,oriResp,or);
            oriResp_sem = cat(1,oriResp_sem,ors);   
        else                        
            oriResp = cat(1,oriResp,d.oriResp);
            oriResp_sem = cat(1,oriResp_sem,d.oriRespSem);     
        end
    end
end


[nTCFr,~] = size(targetTC{1,1});
tt_all_fr = (1:nTCFr)-(nBaselineFr+nFrMonitorDelay+1);
tt_all = tt_all_fr./frameRateHz;

ttWinInd = find(tt_all < -0.5,1,'last'):find(tt_all > 0.51,1,'first');
tt = tt_all(ttWinInd);
tt_label = -0.5:0.25:2;
tt_label_fr = find(ismember(tt,tt_label));
tt_label_s = tt(tt_label_fr);
hmLim =  [-0.1 0.1];


[~,sortInd] = sort(mean(targetTC{end}(respwin,:),1));
if strcmp(ds,'FSAV_attentionV1')
    exampleCells = [exampleCell_1, exampleCell_2, exampleCell_3];
    exCellInd = find(fliplr(ismember(antiSort,exampleCells)));
else
    exCellInd = [];
end

setFigParams4Print('landscape')
figure
colormap(brewermap([],'*RdBu'));
suptitle('Target Aligned Response, Sorted by Response Late Anti, All Cells');
subplot 221
hm = flipud(targetTC{1}(:,antiSort)');
imagesc(hm(:,ttWinInd));
hold on
if ~isempty(exCellInd)
    hline(exCellInd,'k-');
end
figXAxis([],'Time (s)',[],tt_label_fr,tt_label_s)
figYAxis([],'Cell #',[])
figAxForm([])
caxis(hmLim)
colorbar
title('Difficult Targets (<45 deg)')

subplot 222
hm = flipud(targetTC{2}(:,antiSort)');
imagesc(hm(:,ttWinInd));
hold on
if ~isempty(exCellInd)
    hline(exCellInd,'k-');
end
figXAxis([],'Time (s)',[],tt_label_fr,tt_label_s)
figYAxis([],'Cell #',[])
figAxForm([])
caxis(hmLim)
colorbar
title('Easy Targets (>=45 deg)')

subplot 223
hm = flipud(oriFits(:,antiSort)');
imagesc(hm)
hold on
if ~isempty(exCellInd)
    hline(exCellInd,'k-');
end
figXAxis([],'Orientation (deg)',[],45:45:180,45:45:180)
figYAxis([],'Cell #',[])
figAxForm
caxis([-0.25 0.25])
colorbar
title('Orienation Tuning')
print([fnout 'binnedTC_heatmap_allCells'],'-dpdf','-fillpage');


colors = brewermap(nbin,'Greys');
respTC_lim = [-0.11 0.25];
setFigParams4Print('landscape')
figure
exampleCellName = {'First-Stim Responsive';'Late Anti Responsive';'Late Anti Suppressed'};
suptitle('Target Aligned, Example Cells,Dark/Light: >/< 45deg')
for icell = 1:3
    subplot(3,1,icell)
    for i = 1:nbin
        y = targetTC{i}(ttWinInd,exampleCells(icell));
        yerr = targetTC_sem{i}(ttWinInd,exampleCells(icell));
        hold on
        h = shadedErrorBar(tt,y,yerr);
        h.mainLine.Color = colors(i,:);
        h.edge(1).Color = 'none';
        h.edge(2).Color = 'none';
        h.patch.FaceColor = h.mainLine.Color +((1 - h.mainLine.Color).*0.5);
    end
    figXAxis([],'Time from Target Stim (s)',[tt(1) tt(end)],tt_label,tt_label)
    figYAxis([],'dF/F',respTC_lim)
    figAxForm([],0)
    title(sprintf('%s #%s',exampleCellName{icell},num2str(exampleCells(icell))))
end
print([fnout 'binnedTC_exampleCells'],'-dpdf','-fillpage');

figure
suptitle('Tuning, Example Cells')
for icell = 1:3
    subplot(3,1,icell)
    for i = 1:nbin
        y = oriResp(exampleCells(icell),:);
        yerr = oriResp_sem(exampleCells(icell),:);
        y(y<0) = 0;
        yerr(y<=0) = 0;
        hold on
        h = errorbar([0:45:135],y,yerr,'ko');
        hold on
        h = plot(0:180,oriFits(:,exampleCells(icell)),'k-');
        
    end
    figXAxis([],'Orientation (deg)',[-45 225],45:45:180,45:45:180)
    figYAxis([],'Cell #',[-0.1 0.1])
    figAxForm([],0)
    title(sprintf('%s #%s',exampleCellName{icell},num2str(exampleCells(icell))))
end
print([fnout 'oriTuning_exampleCells'],'-dpdf','-fillpage');

