clear all
close all
ds = 'szTuning_dreadds_PM';
rc = behavConstsAV;
eval(ds)
slct_expt = 1;
%%
iexp = slct_expt;

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

%% extract iti activity
itiEventsSum = cell(1,nrun);
itiTC = cell(1,nrun);
for irun = 1:nrun
    mwRun = mw{irun};
    tcRun = tc{irun};
    nOn = mwRun.nScansOn;
    nOff = mwRun.nScansOff;
    [nFr,nCells] = size(tcRun);
    nTrials = nFr./(nOn+nOff);
    tcTrials = reshape(tcRun,[nOn+nOff, nTrials, nCells]);
    itiFrameInd = (nOff-nItiFr+1):nOff;
    itiTC{irun} = reshape(tcTrials(itiFrameInd,:,:),[nTrials*nItiFr, nCells]);
    
    meanCell = mean(itiTC{irun},1);
    stdCell = std(itiTC{irun},1);
    twoStdThreshold = (2.*stdCell)+meanCell;
    itiEvents = itiTC{irun} > twoStdThreshold;
    itiEventsSum{irun} = sum(itiEvents,1);
end
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

%% stimulus-driven activity
sizeTC = cell(1,nrun);
for irun = 1:nrun
    mwRun = mw{irun};
    tcRun = tc{irun};
    nOn = mwRun.nScansOn;
    nOff = mwRun.nScansOff;
    [nFr,nCells] = size(tcRun);
    nTrials = nFr./(nOn+nOff);    
    tcTrials = reshape(tcRun,[nOn+nOff, nTrials, nCells]);
    f = mean(tcTrials((nOff-nBaselineFr+1):nOff,:,:),1);
    dff = (tcTrials-f)./f;
    tSize = round(celleqel2mat_padded(mwRun.tGratingDiameterDeg),2,'significant');
    sizes = unique(tSize);
    nSize = length(sizes);
    stimFrInd = (nOff-nBaselineFr+1):(nOff+nOn);
    sizeTCRun = nan(nSize,nOn+nBaselineFr,nCells);
    for i = 1:nSize
        ind = tSize == sizes(i);
        sizeTCRun(i,:,:) = squeeze(mean(dff(stimFrInd,ind,:),2));
    end
    sizeTC{irun} = sizeTCRun;
end
%% plot iti activity
drugTimes = cellfun(@num2str,expt(iexp).sizeTuningTimeFromDrugMin,'unif',0);


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

setFigParams4Print('landscape')
figure
fLim = [0 10000];
[~,sortedCellsInd] = sort(itiEventsSum{1});
top5Cells = sortedCellsInd(end-4:end);
for i = 1:5
    subplot(nrun,1,i)
    c = top5Cells(i);
    y = itiTC{1}(:,c);
    fLim = [0 max(y)+300];
    x = 1:size(y,1);
    h = plot(x,y,'k-');
    hold on
    y = itiTC{2}(:,c);
    h = plot(x,y,'-');
    h.Color = colors(2,:);
    y = itiTC{end}(:,c);
    h = plot(x,y,'-');
    h.Color = colors(end,:);
    figXAxis([],'ITI Frame Number',[1 x(end)])
    figYAxis([],'F',fLim)
    figAxForm([],0)
    title(sprintf('Cell %s',num2str(c)))
    if i == 1
        legend(drugTimes([1 2 end]))
    end
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

