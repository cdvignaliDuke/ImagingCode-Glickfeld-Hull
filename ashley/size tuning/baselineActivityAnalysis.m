clear all
close all
ds = 'szTuning_dreadds_PM';
rc = behavConstsAV;
eval(ds)
slct_expt = 2;
%%
iexp = slct_expt;

tcLengthMin = 5;
maxTcTimeMin = 60;
tcLengthFr = (tcLengthMin)*60*params.frameRate;

mouse = expt(iexp).mouse;
subnum = mouse;
expDate = expt(iexp).date;
runs = expt(iexp).sizeTuningFolder;
nrun = length(runs);
nItiFr = round(params.itiTimeS.*params.frameRate);
nBaselineFr = round(params.nBaselineMs.*params.frameRate./1000);

fnout = fullfile(rc.ashleyAnalysis,mouse,'two-photon imaging',expDate);

%% load data
tc = cell(1,nrun);
mw = cell(1,nrun);
for irun = 1:nrun
    dataFolder = runs{irun};
    dataTime = expt(iexp).sizeTuningTime{irun};
    load(fullfile(fnout,dataFolder,'timecourses'))
    tc{irun} = dataTC;
    mw{irun} = loadMworksFile(subnum,expDate,dataTime);
end

%%
drugTimes = cellfun(@num2str,expt(iexp).sizeTuningTimeFromDrugMin,'unif',0);
drugTimesS = cellfun(@(x) (str2double(x)).*60,drugTimes,'unif',0);
timeFromDrugS = cellfun(@(x,y) (((0:length(x)-1)./params.frameRate)+y),...
    tc,drugTimesS','unif',0);
timeFromDrugS{1} = (-length(timeFromDrugS{1}):-1)./params.frameRate;
timeFromDrugMin = cellfun(@(x) x./60,timeFromDrugS,'unif',0);


anaylsisStartTimesMin = 0:tcLengthMin:maxTcTimeMin;
nTimes = length(anaylsisStartTimesMin)+1;

%%
%% extract spontaneous activity

spontTC = [];
spontTimePointsMin = [];

for irun = 1:nrun
    mwRun = mw{irun};
    tcRun = tc{irun};
    nOn = mwRun.nScansOn;
    nOff = mwRun.nScansOff;
    [nFr,nCells] = size(tcRun);
    nTrials = nFr./(nOn+nOff);
    tcTrials = reshape(tcRun,[nOn+nOff, nTrials, nCells]);
    drugTimeTrials = reshape(timeFromDrugMin{irun},[nOn+nOff, nTrials]);
            
    trialInd = 1:nTrials;
    tSize = mwRun.tGratingDiameterDeg;
    if length(tSize) < nTrials
        trialInd = 1:length(tSize);
    end
    
    spontInd = setdiff(1:(nOff+nOn),[1:nBaselineFr nOff+1:(nOff+nOn)]);
    
    spontTCRun = tcTrials(spontInd,trialInd,:);
    spontTCRun = reshape(spontTCRun, [length(spontInd)*length(trialInd) nCells]);
    
    timepointsRun = drugTimeTrials(spontInd,trialInd);
    timepointsRun = timepointsRun(:);
    
    if length(tSize) < nTrials
        extraFrames = tcTrials(:,(trialInd(end)+1):end,:);
        spontTC = cat(1,spontTC,spontTCRun,...
            reshape(extraFrames,[(nOn+nOff)*(nTrials - length(tSize)),nCells]));
        extraTimepoints = drugTimeTrials(:,trialInd(end)+1);
        spontTimePointsMin = cat(1,spontTimePointsMin,timepointsRun,...
            extraTimepoints(:));
    else
        spontTC = cat(1,spontTC,spontTCRun);
        spontTimePointsMin = cat(1,spontTimePointsMin,timepointsRun);
    end
    
end


meanCell = mean(spontTC,1);
stdCell = std(spontTC,1);
twoStdThreshold = (2.*stdCell)+meanCell;
spontEvents = spontTC > twoStdThreshold;

timeBinEdges = [-10:10:70];
[nTimepointsInBin,~,timeBinID] = histcounts(spontTimePointsMin,timeBinEdges);
nBins = max(timeBinID);
nTimepointsInBin = nTimepointsInBin(1:nBins);
minTrN = 10;
minTimepoints = min(nTimepointsInBin);

matchedBinnedSpontEvents = cell(1,nBins);
for ibin = 1:nBins
    ind = sort(randsample(find(timeBinID == ibin),minTimepoints));
    matchedBinnedSpontEvents{ibin} = spontEvents(ind,:);    
end
spontEventsSum = cell2mat(cellfun(@sum,matchedBinnedSpontEvents,'unif',0)');

spontEventsSumNorm = spontEventsSum./spontEventsSum(1,:);
%%
drugTimes = cellfun(@num2str,num2cell(timeBinEdges(1:nBins)),'unif',0);
figure
subplot 131
eventLim = [-20 900];
colors = brewermap(nBins+1,'Blues');
colors = colors(3:end,:);
x = spontEventsSum(1,:);
for i = 1:nBins-1
    hold on
    y = spontEventsSum(i+1,:);
    h = plot(x,y,'.');
    h.MarkerSize = 15;
    h.Color = colors(i,:);
end
plot(eventLim, eventLim,'k--')
figXAxis([],'Pre-CNO Events',eventLim)
figYAxis([],'Post-CNO Events',eventLim)
figAxForm
legend(cat(1,drugTimes(2:end)',{'Unity'}),'location','northeastoutside')

subplot 132
x = timeBinEdges(1:nBins);
y = mean(spontEventsSum,2);
yerr = ste(spontEventsSum,2);
h = errorbar(x,y,yerr,'o-');
h.MarkerFaceColor = 'none';
figXAxis([],'Time From Drug Delivery',[x(1)-10 x(end)+10])
figYAxis([],'Average N Events',[])
figAxForm

subplot 133
x = timeBinEdges(1:nBins);
y = mean(spontEventsSumNorm,2);
yerr = ste(spontEventsSumNorm,2);
h = errorbar(x,y,yerr,'o-');
h.MarkerFaceColor = 'none';
figXAxis([],'Time From Drug Delivery',[x(1)-10 x(end)+10])
figYAxis([],'Normalized N Events',[])
figAxForm
%% entire tc events
filterWindowFr = 50;
movAvgFilt = ones(filterWindowFr+1,1)/(filterWindowFr+1);
eventsTC = cell(1,nrun);
eventsSum = cell(1,nrun);
for irun = 1:nrun
    tcRun = tc{irun};
    [nFr,nCells] = size(tcRun);
    meanCell = mean(tcRun,1);
    stdCell = std(tcRun,1);
    twoStdThreshold = (2.*stdCell)+meanCell;
    events = tcRun > twoStdThreshold;
    smoothEvents = nan(nFr,nCells);
    for i = 1:nCells
        smoothEvents(:,i) = conv(events(:,i),movAvgFilt,'same');
    end
    eventsTC{irun} = smoothEvents;
    eventsSum{irun} = sum(events,1);
end

%% extract iti activity
itiEventsSum = cell(1,nrun);
itiTC = cell(1,nrun);
itiTimeFromDrugMin = cell(1,nrun);

for irun = 1:nrun
    mwRun = mw{irun};
    tcRun = tc{irun};
    nOn = mwRun.nScansOn;
    nOff = mwRun.nScansOff;
    [nFr,nCells] = size(tcRun);
    nTrials = nFr./(nOn+nOff);
    tcTrials = reshape(tcRun,[nOn+nOff, nTrials, nCells]);
    itiFrameInd = (nOff-nItiFr+1):nOff;
%     itiTC{irun} = tcTrials(itiFrameInd,:,:);
    itiTC{irun} = reshape(tcTrials(itiFrameInd,:,:),[nTrials*nItiFr, nCells]);
    
    meanCell = mean(itiTC{irun},1);
    stdCell = std(itiTC{irun},1);
    twoStdThreshold = (2.*stdCell)+meanCell;
    itiEvents = itiTC{irun} > twoStdThreshold;
    itiEventsSum{irun} = sum(itiEvents,1);
    
    drugTimeTrials = reshape(timeFromDrugMin{irun},[nOn+nOff, nTrials]);
    itiTimeFromDrugMin{irun} = reshape(drugTimeTrials(itiFrameInd,:),[nItiFr*nTrials,1]);
%     itiTimeFromDrugMin{irun} = drugTimeTrials(itiFrameInd,:);
end

%%

%% stimulus-driven activity
sizeTC = cell(1,nrun);
tcTimePoints = [];
respTimePoints = [];
drugTimePoints = [];
respStimType = [];
for irun = 1:nrun
    mwRun = mw{irun};
    tcRun = tc{irun};
    timesRun = timeFromDrugMin{irun};
    nOn = mwRun.nScansOn;
    nOff = mwRun.nScansOff;
    [nFr,nCells] = size(tcRun);
    nTrials = nFr./(nOn+nOff);    
    tcTrials = reshape(tcRun,[nOn+nOff, nTrials, nCells]);
    trialTimes = reshape(timesRun,[nOn+nOff,nTrials]);
    trialTimes = trialTimes(nOff+1,:);
    f = mean(tcTrials((nOff-nBaselineFr+1):nOff,:,:),1);
    dff = (tcTrials-f)./f;
    tSize = round(celleqel2mat_padded(mwRun.tGratingDiameterDeg),2,'significant');
    sizes = unique(tSize);
    nSize = length(sizes);
    
    if nTrials ~= length(tSize)
        n = length(tSize);
        dff = dff(:,1:n,:);
        trialTimes = trialTimes(1:n);
    end
    tcTimePoints = cat(2,tcTimePoints,dff);
    respwin = (nBaselineFr+1):(nOn+nBaselineFr);
    respTrials = squeeze(mean(dff(respwin,:,:),1));
    respTimePoints = cat(1,respTimePoints,respTrials);
    drugTimePoints = cat(2,drugTimePoints,trialTimes);
    respStimType = cat(2,respStimType,tSize);
    
    stimFrInd = (nOff-nBaselineFr+1):(nOff+nOn);
    sizeTCRun = nan(nSize,nOn+nBaselineFr,nCells);
    if irun == 1
        responsiveCellsSizes = nan(nSize,nCells);
        alpha = 0.05./nSize;
    end
    for i = 1:nSize
        ind = tSize == sizes(i);
        sizeTCRun(i,:,:) = squeeze(mean(dff(stimFrInd,ind,:),2));
        if irun == 1
            basewin = 1:nBaselineFr;
            baseResp = squeeze(mean(dff(basewin,ind,:),1));
            respResp = squeeze(mean(dff(respwin,ind,:),1));
            responsiveCellsSizes(i,:) = logical(ttest(...
                baseResp,respResp,'tail','left','alpha',alpha));
        end
    end
    sizeTC{irun} = sizeTCRun;
    
    
    
end

responsiveCells = find(any(responsiveCellsSizes,1));

%%

exptTimeLim = [-10 70];
nExampleCells = 5;
if length(responsiveCells) < nExampleCells
    nRandCells = nExampleCells - length(responsiveCells);
    exampleCells = cat(2,responsiveCells,randsample(setdiff(1:nCells,responsiveCells),nRandCells));
else
    exampleCells = randsample(responsiveCells,nExampleCells);
end

timeBinEdges = [-10:10:70];
[~,~,timeBinID] = histcounts(drugTimePoints,timeBinEdges);
minTrN = 10;
sizes = unique(respStimType);
nSize = length(sizes);
tcSizeTimePoints = cell(max(timeBinID),nSize);
respCells = cell(1,nSize);

for i = 1:nSize
    ind = respStimType == sizes(i);
    tcSize = tcTimePoints(:,ind,:);
    respSize = squeeze(mean(tcSize(respwin,:,:),1));
    baseSize = squeeze(mean(tcSize(basewin,:,:),1));
    respCells{i} = logical(ttest(baseSize',respSize','tail','left'));
    for ibin = 1:max(timeBinID)
        ind = respStimType == sizes(i) & timeBinID == ibin;
        if sum(ind) > minTrN
            tcSizeTimePoints{ibin,i} = squeeze(mean(tcTimePoints(:,ind,:),2));
        end            
    end
end

%%
tt = double(-nBaselineFr+1:nOn)./params.frameRate.*1000;
ttLabel = -500:250:max(tt);
colors = brewermap(max(timeBinID)+1,'Blues');
colors = colors(2:end,:);
figure
for icell = 1:5
    subplot(2,3,icell)
    xc = exampleCells(icell);
    for ibin = 1:max(timeBinID)
        if ~isempty(tcSizeTimePoints{ibin,2})
            y = tcSizeTimePoints{ibin,2}(:,xc);
            x = double((-nBaselineFr+1):(length(y)-nBaselineFr))./params.frameRate.*1000;
            hold on
            h = plot(x,y,'-');
            if ibin == 1
                h.Color = 'k';
            else
                h.Color = colors(ibin-1,:);
            end        
        end
    end
    figXAxis([],'Time From Vis Stim (ms)',[min(x) max(x)],ttLabel,ttLabel)
    figYAxis([],'dF/F',[])
    figAxForm
    title(sprintf('Cell %s',num2str(xc)))
end

%%
exptTimeLim = [-10 70];
nExampleCells = 5;
if length(responsiveCells) < nExampleCells
    nRandCells = nExampleCells - length(responsiveCells);
    exampleCells = cat(2,responsiveCells,randsample(setdiff(1:nCells,responsiveCells),nRandCells));
else
    exampleCells = randsample(responsiveCells,nExampleCells);
end

colors = brewermap(nSize+1,'Reds');
colors = colors(2:end,:);
figure
for icell = 1:nExampleCells
    subplot(nExampleCells,1,icell)
    xc = exampleCells(icell);
    for i = nSize
        ind = respStimType == sizes(i);
        x = drugTimePoints(ind);
        y = respTimePoints(ind,xc);
        hold on
        h = plot(x,y,'o');
        h.Color = colors(i,:);
        h.MarkerFaceColor = colors(i,:);
    end    
    title(sprintf('Cell %s',num2str(xc)))
    xLabel = sprintf('Time From %s Delivery (min)',expt(slct_expt).drug);
    figXAxis([],xLabel,exptTimeLim)
    figYAxis([],'dF/F',[])
    figAxForm([],0)
    hline(0,'k:')
end

%% plot iti activity


setFigParams4Print('landscape')
figure
suptitle({sprintf('%s-%s-%s',mouse,expDate,expt(iexp).img_loc{1});...
    'Event = F > 2std of mean';['ITI Only' num2str(nCells) ' Cells']})

subplot 121
eventLim = [0 300];
colors = brewermap(nrun+1,'Blues');
colors = colors(3:end,:);
x = itiEventsSum{1};
for i = 1:(nrun-1)
    hold on
    y = itiEventsSum{i+1};
    h = plot(x,y,'.');
    h.Color = colors(i,:);
end
plot(eventLim, eventLim,'k--')
figXAxis([],'Pre-CNO Events',eventLim)
figYAxis([],'Post-CNO Events',eventLim)
figAxForm
legend(cat(1,drugTimes(2:end),{'Unity'}),'location','northeastoutside')

subplot 122
x = 1:nrun;
y = cellfun(@mean,itiEventsSum);
yerr = cellfun(@(x) ste(x,2),itiEventsSum);
h = errorbar(x,y,yerr,'ko-');
h.MarkerFaceColor = [1 1 1];
figXAxis([],'Session Start Time From Drug Delivery (min)',...
    [0 nrun+1],x,drugTimes)
figYAxis([],'N Events',[])
figAxForm
print(fullfile(fnout,'itiEventsOverSessions'),'-dpdf','-fillpage')

% figure
% fLim = [0 10000];
[~,sortedCellsInd] = sort(itiEventsSum{1});
exampleCells = sortedCellsInd(end-4:end);

nExampleCells = 5;

% exampleCells = randsample(1:nCells,nExampleCells);

setFigParams4Print('landscape')
figure
for icell = 1:nExampleCells
    xc = exampleCells(icell);
    subplot(nExampleCells,1,icell)
    for irun = 1:nrun
        tcRun = itiTC{irun}(:,xc);
        
        drugTime = itiTimeFromDrugMin{irun};
%         nTrials = size(tcRun,2);
%         for itrial = 1:nTrials
            hold on
            plot(drugTime,tcRun,'k-')
%         end
    end
    title(sprintf('Cell %s',num2str(xc)))
    xLabel = sprintf('Time From %s Delivery (min)',expt(slct_expt).drug);
    figXAxis([],xLabel,exptTimeLim)
    figYAxis([],'F',[])
    figAxForm([],0)
end

print(fullfile(fnout,'itiTCExampleCells'),'-dpdf','-fillpage')

%% plot all tc
[~,sortedCellsInd] = sort(eventsSum{1});
top5Cells = sortedCellsInd(end-4:end);

figure
suptitle({sprintf('%s-%s-%s',mouse,expDate,expt(iexp).img_loc{1});...
    'Entire Session'})
for i = 1:5
    subplot(5,1,i)
    c = top5Cells(i);
    y = tc{1}(:,c);
    offset = max(y)+min(y)+300;
    x = 1:size(y,1);
    h = plot(x,y,'k-');
    hold on
    y = tc{2}(:,c)+offset;
    h = plot(x,y,'-');
    h.Color = colors(2,:);
    y = tc{end}(:,c)+(offset*2);
    h = plot(x,y,'-');
    h.Color = colors(end,:);
    fLim = [0 max(y)+300];
    figXAxis([],'Session Frame Number',[1 x(end)])
    figYAxis([],'F',fLim)
    figAxForm([],0)
    title(sprintf('Cell %s',num2str(c)))
%     vline(nOff:(nOn+nOff):size(y,1),'k-')
    if i == 1
        legend(drugTimes([1 2 end]))
    end
end
print(fullfile(fnout,'TCExampleCells'),'-dpdf','-fillpage')

%% plot stimulus responses
[~,sortedCellsInd] = sort(eventsSum{1});
top6Cells = sortedCellsInd(end-5:end);
tt = double(-nBaselineFr+1:nOn)./params.frameRate.*1000;
ttLabel = -500:250:max(tt);

exStimSize = nSize;
figure;
suptitle({sprintf('%s-%s-%s',mouse,expDate,expt(iexp).img_loc{1});...
    sprintf('%s deg stim',num2str(sizes(exStimSize)))})
colormap(brewermap([],'Blues'))
for i = 1:6
    subplot(2,3,i)
    c = top6Cells(i);
    for irun = 1:nrun
        hold on
        y = sizeTC{irun}(exStimSize,:,c);
        h = plot(tt,y,'-');
        if irun == 1
            h.Color = 'k';
        else
            h.Color = colors(irun-1,:);
        end
    end
    figXAxis([],'Time From Vis Stim (ms)',[min(tt) max(tt)],ttLabel,ttLabel)
    figYAxis([],'dF/F',[])
    figAxForm
    title(sprintf('Cell %s',num2str(c)))
    if i == 1
        L = legend(drugTimes,'location','northwest');
        title(L,'Time From Drug')
    end
end
fileName = sprintf('%sdegTCExampleCells',num2str(sizes(exStimSize)));
print(fullfile(fnout,fileName),'-dpdf','-fillpage')

respwin = (nBaselineFr+1):(nOn+nBaselineFr);
sizeResp = cellfun(@(x) squeeze(mean(x(:,respwin,:),2)),sizeTC,'unif',0);
figure
suptitle(sprintf('%s-%s-%s',mouse,expDate,expt(iexp).img_loc{1}))
for i = 1:6
    subplot(2,3,i);
    c = top6Cells(i);
    for irun = 1:nrun
        hold on
        y = sizeResp{irun}(:,c);
        h = plot(sizes,y,'-');
        if irun == 1
            h.Color = 'k';
        else
            h.Color = colors(irun-1,:);
        end
    end
    figXAxis([],'Stim Size (deg)',[0 max(sizes)+10],sizes,sizes)
    figYAxis([],'dF/F',[])
    figAxForm
    title(sprintf('Cell %s',num2str(c)))
    if i == 1
        L = legend(drugTimes,'location','northwest');
        title(L,'Time From Drug')
    end
end
print(fullfile(fnout,'sizeTuningExampleCells'),'-dpdf','-fillpage')

figure;
exStimSize = 4;
suptitle({sprintf('%s-%s-%s',mouse,expDate,expt(iexp).img_loc{1});...
    sprintf('%s deg stim',num2str(sizes(exStimSize)))})
colormap(brewermap([],'Blues'))
for i = 1:6
    subplot(2,3,i)
    c = top6Cells(i);
    for irun = 1:nrun
        hold on
        y = sizeTC{irun}(exStimSize,:,c);
        h = plot(tt,y,'-');
        if irun == 1
            h.Color = 'k';
        else
            h.Color = colors(irun-1,:);
        end
    end
    figXAxis([],'Time From Vis Stim (ms)',[min(tt) max(tt)],ttLabel,ttLabel)
    figYAxis([],'dF/F',[])
    figAxForm
    title(sprintf('Cell %s',num2str(c)))
    if i == 1
        L = legend(drugTimes,'location','northwest');
        title(L,'Time From Drug')
    end
end
fileName = sprintf('%sdegTCExampleCells',num2str(sizes(exStimSize)));
print(fullfile(fnout,fileName),'-dpdf','-fillpage')


