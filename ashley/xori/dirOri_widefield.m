clear all
close all
ds = 'dirOri_emx_widefield';
rc = behavConstsAV;
eval(ds);
nexp = size(expt,2);
fnout = fullfile(rc.ashleyAnalysis,'Expt Summaries',ds);
if exist(fnout,'dir') == 0
    mkdir(fnout)
end
%% load data
data = struct;
for iexp = 1:nexp

    mouse = expt(iexp).mouse;
    subnum = mouse;
    expDate = expt(iexp).date;
    fn = fullfile(rc.ashleyAnalysis,mouse,'widefield imaging',expDate);
    load(fullfile(fn,'exptSum_struct'))
    if iexp == 1
        data = wfExpt;
    else
        data(iexp) = wfExpt;
    end
end

areas = [];
for iexp = 1:nexp
    areas = cat(1,data.roiLabels);
end
areas = unique(areas);
nAreas = length(areas);

stimInfo = data(1).stimInfo.stimulusTypes;
nStimType = length(stimInfo);
nBaselineFr = data(1).stimInfo.baseline_frames;
%% combine expts
areaSummary = struct;
for iexp = 1:nexp
    respwin = (nBaselineFr+1):(nBaselineFr+2);
    exptAreas = data(iexp).roiLabels;
    for iarea = 1:nAreas
        ind = strcmp(exptAreas,areas(iarea));
        if any(ind)
            tc = cellfun(@(x) squeeze(x(:,ind,:)),...
                data(iexp).tc_roiXstimType,'unif',0);
            if iexp == 1
                areaSummary(iarea).tc = cellfun(@(x) mean(x,2),tc,'unif',0);
                areaSummary(iarea).tcErr = cellfun(@(x) ste(x,2),tc,'unif',0);
                areaSummary(iarea).resp = cellfun(@(x) mean(mean(x(respwin,:),2)),...
                    tc,'unif',0);
                areaSummary(iarea).respErr = cellfun(@(x) ste(mean(x(respwin,:),1),2),...
                    tc,'unif',0);
            else
                areaSummary(iarea).tc = cat(1,areaSummary(iarea).tc,...
                    cellfun(@(x) mean(x,2),tc,'unif',0));
                areaSummary(iarea).tcErr = cat(1,areaSummary(iarea).tcErr,...
                    cellfun(@(x) ste(x,2),tc,'unif',0));
                areaSummary(iarea).resp = cat(1,areaSummary(iarea).resp,...
                    cellfun(@(x) mean(mean(x(respwin,:),2)),...
                    tc,'unif',0));
                areaSummary(iarea).respErr = cat(1,areaSummary(iarea).tc,...
                    cellfun(@(x) ste(mean(x(respwin,:),1),2),...
                    tc,'unif',0));
            end
        end
    end
end

%% plot responses
[nRow,nCol] = optimizeSubplotDim(nAreas);
figure
suptitle('Mean Response in Each Area Across Expt')
for iarea = 1:nAreas
    d = cell2mat(areaSummary(iarea).resp);
    y = mean(d,1);
    yerr = ste(d,1);
    subplot(nRow,nCol,iarea)
    h = errorbar(1:nStimType,y,yerr,'ko');
    title(areas{iarea})
    figXAxis([],'',[0 nStimType+1],1:nStimType,stimInfo)
    xtickangle(45)
    figYAxis([],'dF/F',[0 0.04])
    figAxForm([],0)
end

print(fullfile(fnout,'resp_allExpt'),'-dpdf','-fillpage')
%% plot time-courses
nTrialFr = length(areaSummary(1).tc{1,1});
timeFromStimMs = ((-nBaselineFr+1):(nTrialFr-nBaselineFr))./params.frameRate.*1000;
% respwinMs = ([respwin(1) respwin(2)]-nBaselineFr)./params.frameRate.*1000;

figure
suptitle('90 deg, static-blk, drift-red')
for iarea = 1:nAreas
    d = areaSummary(iarea).tc;
    subplot(nRow,nCol,iarea)
    
    ind = strcmp(stimInfo,'S-90');
    y = mean(cell2mat(d(:,ind)'),2);
    yerr = ste(cell2mat(d(:,ind)'),2);
    if all(~isnan(yerr))
        h = shadedErrorBar_chooseColor(timeFromStimMs,y,yerr,[0 0 0]);
        hold on
        ind = strcmp(stimInfo,'D-90');
        y = mean(cell2mat(d(:,ind)'),2);
        yerr = ste(cell2mat(d(:,ind)'),2);
        h = shadedErrorBar_chooseColor(timeFromStimMs,y,yerr,[1 0 0]);
    else
        h = plot(timeFromStimMs,y);
        h.Color = [0 0 0];
        hold on
        ind = strcmp(stimInfo,'D-90');
        y = mean(cell2mat(d(:,ind)'),2);
        yerr = ste(cell2mat(d(:,ind)'),2);
        h = plot(timeFromStimMs,y);
        h.Color = [1 0 0];
    end
    title(areas{iarea})
    figXAxis([],'Time (ms)',[])
    figYAxis([],'dF/F',[])
    figAxForm([],0)
end

figure
suptitle('0 deg, static-blk, drift-red')
for iarea = 1:nAreas
    d = areaSummary(iarea).tc;
    subplot(nRow,nCol,iarea)
    
    ind = strcmp(stimInfo,'S-0');
    y = mean(cell2mat(d(:,ind)'),2);
    yerr = ste(cell2mat(d(:,ind)'),2);
    if all(~isnan(yerr))
        h = shadedErrorBar_chooseColor(timeFromStimMs,y,yerr,[0 0 0]);
        hold on
        ind = strcmp(stimInfo,'D-0');
        y = mean(cell2mat(d(:,ind)'),2);
        yerr = ste(cell2mat(d(:,ind)'),2);
        h = shadedErrorBar_chooseColor(timeFromStimMs,y,yerr,[1 0 0]);
    else
        h = plot(timeFromStimMs,y);
        h.Color = [0 0 0];
        hold on
        ind = strcmp(stimInfo,'D-0');
        y = mean(cell2mat(d(:,ind)'),2);
        yerr = ste(cell2mat(d(:,ind)'),2);
        h = plot(timeFromStimMs,y);
        h.Color = [1 0 0];
    end
    title(areas{iarea})
    figXAxis([],'Time (ms)',[])
    figYAxis([],'dF/F',[])
    figAxForm([],0)
end
figure
suptitle('0/180 drifting, optic flow?,0-red,180-dark red')
for iarea = 1:nAreas
    d = areaSummary(iarea).tc;
    subplot(nRow,nCol,iarea)
    
    ind = strcmp(stimInfo,'D-0');
    y = mean(cell2mat(d(:,ind)'),2);
    yerr = ste(cell2mat(d(:,ind)'),2);
    if all(~isnan(yerr))
        h = shadedErrorBar_chooseColor(timeFromStimMs,y,yerr,[1 0 0]);
        hold on
        ind = strcmp(stimInfo,'D-180');
        y = mean(cell2mat(d(:,ind)'),2);
        yerr = ste(cell2mat(d(:,ind)'),2);
        h = shadedErrorBar_chooseColor(timeFromStimMs,y,yerr,[0.5 0 0]);
    else
        h = plot(timeFromStimMs,y);
        h.Color = [1 0 0];
        hold on
        ind = strcmp(stimInfo,'D-180');
        y = mean(cell2mat(d(:,ind)'),2);
        yerr = ste(cell2mat(d(:,ind)'),2);
        h = plot(timeFromStimMs,y);
        h.Color = [0.5 0 0];
    end
    title(areas{iarea})
    figXAxis([],'Time (ms)',[])
    figYAxis([],'dF/F',[])
    figAxForm([],0)
end