clear all
close all
ds = 'ConAV_V1_naive';
rc = behavConstsAV;
eval(ds)
iexp = 1;
%% load data
mouse = expt(iexp).mouse;
expDate = expt(iexp).date;

fn = fullfile(rc.ashleyAnalysis,mouse,'two-photon imaging',expDate,'data processing');
load(fullfile(fn,'timecourses_bx_cells'))
dataTC = data_bx_tc_subnp;
clear data_bx_tc_subnp data_bx_tc
for irun = 1:expt(iexp).nrun
    if irun == 1
        mw_temp = loadMworksFile(mouse,expDate,expt(iexp).time_mat(irun,:));
    else
        mw_temp = [mw_temp loadMworksFile(mouse,expDate,expt(iexp).time_mat(irun,:))];
    end
end
mw = concatenateDataBlocks(mw_temp);

%% image info
[nfr,ncells] = size(dataTC);
%% analysis params

frHz= 30;
nvisdelayframes = 3;
nvisdelayframes_target = 2;
nboot = 1000;

nblframes = 30;
ndistractframes = 77;
minCycN_longTC = 7;

basewin = 26:30;
respwin = 34:38;

startInd = 1;
targetInd = 2;

tuningcutoff = 0.8;

tt_long = ((1:(nblframes+ndistractframes)) - (nblframes+nvisdelayframes))...
    .*(1000./frHz);
ttLabel_long = 0:500:2500;
tt_target = ((1:(nblframes+nblframes)) - (nblframes+nvisdelayframes))...
    .*(1000./frHz);
ttLabel_target = 0:200:1000;
tt_late = ((1:(nblframes+nblframes)) - (nblframes+nvisdelayframes)-cycTime)...
    .*(1000./frHz);
ttLabel_late = -1000:200:0;
%% add run block offsets

%% task info

tStart = cell2mat(mw.cLeverDown);
tTargetOn = cell2mat(mw.cTargetOn);
tCycN = cell2mat(mw.tCyclesOn);
tOri = celleqel2mat_padded(mw.tBaseGratingDirectionDeg);
nTrials = length(tStart);
tInvTarget = celleqel2mat_padded(mw.cCatchOn);
tTargetCon = cell2mat(mw.tGratingContrast);
baseCon = unique(cell2mat(mw.tBaseGratingContrast));
tTargetConChange = tTargetCon - baseCon;
conChanges = unique(tTargetConChange);
tIsVis = tTargetConChange ~= 0;
cycTime = unique(cell2mat(mw.nFramesOn)) + unique(cell2mat(mw.nFramesOff));
%% align data to trial start 
dataStruct = struct;

ind = tCycN >= minCycN_longTC & isnan(tInvTarget);
ntr = sum(ind);
dataStruct.start.longTC = nan(nblframes+ndistractframes,ntr,ncells);
dataStruct.start.tOri_longTC = nan(1,ntr);
dataStruct.start.tIsVis_longTC = nan(1,ntr);
dataStruct.start.firstTC = nan(nblframes.*2,nTrials,ncells);
dataStruct.start.lateTC = nan(nblframes.*2,nTrials,ncells);
dataStruct.target.tc = nan(nblframes.*2,nTrials,ncells);
dataStruct.tOri = nan(1,nTrials);
dataStruct.tIsVis = nan(1,nTrials);
for itrial = 1:nTrials
    % long TC
    if itrial == 1
        trfillind = 0;
    end
    if ind(itrial)
        trfillind = trfillind+1;
        trframeind = (tStart(itrial)-nblframes+1):(tStart(itrial)+ndistractframes);
        f0 = mean(dataTC(trframeind(1:nblframes),:),1);
        dff = (dataTC(trframeind,:) - f0)./f0;
        dataStruct.start.longTC(:,trfillind,:) = dff;
        dataStruct.start.tOri_longTC(trfillind) = tOri(itrial);
        dataStruct.start.tIsVis_longTC(trfillind) = tIsVis(itrial);
    end
    
    % response to first stim
    trframeind = (tStart(itrial)-nblframes+1):(tStart(itrial)+nblframes);
    f0 = mean(dataTC(trframeind(1:nblframes),:),1);
    dff = (dataTC(trframeind,:) - f0)./f0;
    dataStruct.start.firstTC(:,itrial,:) = dff;
    dataStruct.tOri(itrial) = tOri(itrial);
    dataStruct.tIsVis(itrial) = tIsVis(itrial);   
    
    % response to target stim
    trframeind = (tTargetOn(itrial)-nblframes+1):(tTargetOn(itrial)+nblframes);
    dff = (dataTC(trframeind,:) - f0)./f0;
    dataStruct.target.tc(:,itrial,:) = dff - mean(dff(basewin,:),1);
    
    % response to late stim
    trframeind = ((tTargetOn(itrial)-nblframes+1):(tTargetOn(itrial)+nblframes))...
        - cycTime;
    dff = (dataTC(trframeind,:) - f0)./f0;
    dataStruct.start.lateTC(:,itrial,:) = dff - mean(dff(basewin,:),1);
end
dataStruct.start.firstResp = squeeze(mean(dataStruct.start.firstTC(respwin,:,:),1) - ...
    mean(dataStruct.start.firstTC(basewin,:,:),1));
dataStruct.target.resp = squeeze(mean(dataStruct.target.tc(respwin,:,:),1) - ...
    mean(dataStruct.target.tc(basewin,:,:),1));
dataStruct.start.lateResp = squeeze(mean(dataStruct.start.lateTC(respwin,:,:),1) - ...
    mean(dataStruct.target.tc(basewin,:,:),1));

%% cell Info

% ori tuning
oris = unique(dataStruct.tOri);
nori = length(oris);

oriTuning = struct;
oriTuning.start.resp = nan(nori,ncells);
oriTuning.start.respErr = nan(nori,ncells);
oriTuning.target.resp = nan(nori,ncells);
oriTuning.target.respErr = nan(nori,ncells);
for iori = 1:nori
    ind = dataStruct.tOri == oris(iori);
    oriTuning.start.resp(iori,:) = mean(dataStruct.start.firstResp(ind,:),1);
    oriTuning.start.respErr(iori,:) = ste(dataStruct.start.firstResp(ind,:),1);
    oriTuning.target.resp(iori,:) = mean(dataStruct.target.resp(ind,:),1);
    oriTuning.target.respErr(iori,:) = ste(dataStruct.target.resp(ind,:),1);
end

[oriTuning.start.fit,oriTuning.start.fitRsquared] = vonmisesFit(...
    oriTuning.start.resp',oris);
[oriTuning.target.fit,oriTuning.target.fitRsquared] = vonmisesFit(...
    oriTuning.target.resp',oris);

% responsive cells
startRespCells = ttest(...
    squeeze(mean(dataStruct.start.longTC(basewin,:,:),1)),...
    squeeze(mean(dataStruct.start.longTC(respwin,:,:),1)),'tail','left');
lateRespCells = ttest(...
    squeeze(mean(dataStruct.start.lateTC(basewin,:,:),1)),...
    squeeze(mean(dataStruct.start.lateTC(respwin,:,:),1)),'tail','left');
tarRespCells = ttest(...
    squeeze(mean(dataStruct.target.tc(basewin,:,:),1)),...
    squeeze(mean(dataStruct.target.tc(respwin,:,:),1)),'tail','left');
    
% structure
cellInfo = struct;
cellInfo.startResp = logical(startRespCells);
cellInfo.lateResp = logical(lateRespCells);
cellInfo.targetResp = logical(tarRespCells);
cellInfo.isTuned = oriTuning.start.fitRsquared > tuningcutoff;
%% decode contrast
nDCCells = 20;
cellInd = find(cellInfo.lateResp | cellInfo.targetResp);
% sampCells = randsample(cellInd,nDCCells);
easyVisTrialInd = tTargetConChange == conChanges(end);
hardVisTrialInd = tTargetConChange == conChanges(end-1);
distResp = dataStruct.start.lateResp(easyVisTrialInd,cellInd);
tarResp = dataStruct.target.resp(easyVisTrialInd,cellInd);
% distResp = dataStruct.start.lateResp(hardVisTrialInd,cellInd);
% tarResp = dataStruct.target.resp(hardVisTrialInd,cellInd);
resp = zscore(cat(1,distResp,tarResp));
targetTrialInd = cat(1,zeros(size(distResp,1),1),ones(size(tarResp,1),1));
C = eye(size(resp,2));
[~,~,targetGLM] = glmfit(resp*C,targetTrialInd,'binomial');
pctCorrHO = getPctCorr_hoData(resp,targetTrialInd,0.5);

figure
histogram(targetGLM.beta(2:end))
%% 
d = dataStruct.start;
L = [];
figure
subplot 411
ind = logical(d.tIsVis_longTC);
y = mean(mean(d.longTC(:,ind,cellInfo.startResp),2),3);
yerr = ste(mean(d.longTC(:,ind,cellInfo.startResp),2),3);
h=shadedErrorBar_chooseColor(tt_long,y,yerr,[0 0 0]);
L(1) = h.mainLine;
hold on
ind = ~d.tIsVis_longTC;
y = mean(mean(d.longTC(:,ind,cellInfo.startResp),2),3);
yerr = ste(mean(d.longTC(:,ind,cellInfo.startResp),2),3);
h=shadedErrorBar_chooseColor(tt_long,y,yerr,[0.75 0.25 0.5]);
L(2) = h.mainLine;
figXAxis([],'Time from Start (ms)',[-500 max(tt_long)],ttLabel_long,ttLabel_long)
figYAxis([],'dF/F',[])
figAxForm([],0)
hline(0,'k--')
title(sprintf('First Stim Resp Cells (n=%s)',num2str(sum(cellInfo.startResp))))
legend(L,{'Vis','Aud'})

L = [];
subplot 423
ind = logical(dataStruct.tIsVis & isnan(tInvTarget));
y = mean(mean(d.lateTC(:,ind,cellInfo.lateResp),2),3);
yerr = ste(mean(d.lateTC(:,ind,cellInfo.lateResp),2),3);
h=shadedErrorBar_chooseColor(tt_late,y,yerr,[0 0 0]);
L(1) = h.mainLine;
hold on
ind = logical(~dataStruct.tIsVis & isnan(tInvTarget));
y = mean(mean(d.lateTC(:,ind,cellInfo.startResp),2),3);
yerr = ste(mean(d.lateTC(:,ind,cellInfo.startResp),2),3);
h=shadedErrorBar_chooseColor(tt_late,y,yerr,[0.75 0.25 0.5]);
L(2) = h.mainLine;
figXAxis([],'Time from Target (ms)',[-1000 0],ttLabel_late,ttLabel_late)
figYAxis([],'dF/F',[])
figAxForm([],0)
hline(0,'k--')
title(sprintf('Late Resp Cells (n=%s)',num2str(sum(cellInfo.lateResp))))
legend(L,{'Vis','Aud'})


d = dataStruct.target;
tarcolors = {[0.75, 0.25,0.5],[0.5 0.5 0.5],[0 0 0]};
L = [];
subplot 424
for icon = 1:length(conChanges)
    ind = tTargetConChange == conChanges(icon);
    y = mean(mean(d.tc(:,ind,cellInfo.targetResp),2),3);
    yerr = ste(mean(d.tc(:,ind,cellInfo.startResp),2),3);
    h=shadedErrorBar_chooseColor(tt_target,y,yerr,tarcolors{icon});
    L(icon) = h.mainLine;
    hold on
end
figXAxis([],'Time from Target (ms)',[-200 max(tt_target)],ttLabel_target,ttLabel_target)
figYAxis([],'dF/F',[])
figAxForm([],0)
hline(0,'k--')
title(sprintf('Target Resp Cells (n=%s)',num2str(sum(cellInfo.targetResp))))
legend(L,{'Aud',num2str(conChanges(2)),num2str(conChanges(end))})

exCells = [6,18,68,70];
% exCells = [1,12,36,76];

for iplot = 1:4
    subplot(4,2,4+iplot)
    x = oris;
    y = oriTuning.start.resp(:,exCells(iplot));
    yerr = oriTuning.start.respErr(:,exCells(iplot));
    h = errorbar(x,y,yerr,'o');
    h.Color = [0.5 0.5 0.5];
    hold on
    y = oriTuning.target.resp(:,exCells(iplot));
    yerr = oriTuning.target.respErr(:,exCells(iplot));
    h = errorbar(x,y,yerr,'ko');
    x = 0:180;
    y = oriTuning.start.fit(:,exCells(iplot));
    h = plot(x,y,'-');
    h.Color = [0.5 0.5 0.5];
    y = oriTuning.target.fit(:,exCells(iplot));
    h = plot(x,y,'k-');
    figXAxis([],'Orientation',[-10 190],oris,oris)
    figYAxis([],'dF/F',[-0.01 0.09])
    figAxForm
    title(sprintf('Cell #%s',num2str(exCells(iplot))))
end

print(fullfile(fn,'conAVsummary'),'-dpdf','-fillpage')