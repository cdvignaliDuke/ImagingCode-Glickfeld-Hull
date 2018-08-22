clear all
close all
ds = 'szTuning_dreadds_widefield';
rc = behavConstsAV;
eval(ds)
slct_expt = 6;
doFullFieldStim = true;
%%
nBaselineFr = round(params.nBaselineMs.*params.frameRate./1000);
nRespFr = round(params.nPostStimMs.*params.frameRate./1000);

%%
iexp = slct_expt;

%%
mouse = expt(iexp).mouse;
subnum = mouse;
expDate = expt(iexp).date;
runs = expt(iexp).visStimFolder;
nrun = length(runs);

fnout = fullfile(rc.isilon,'Analysis',mouse,'widefield imaging',expDate);
fnin = fullfile(rc.isilon,'data',mouse,'widefield imaging',expDate);
%% get outline of injection site
redImage = readtiff(fullfile(fnout, sprintf('AVG_%s_MMStack.ome.tif',expt(iexp).labelFolder)));
redLabel = imCellPolyEditInteractive(redImage);
redMask = bwlabel(redLabel);
figure;imagesc(redImage);colormap gray
% windowOutline = logical(imCellPolyEditInteractive(redImage));
close all

%% get outline of stimulation in injection site
folderName = expt(1).visStimFolder_preDrug;
imageName = sprintf('%s_MMStack.ome.tif',folderName);
data = readtiff(fullfile(fnin,folderName,imageName));
[ypix,xpix,nFr] = size(data);
% for i = 1:nFr
%     img = data(:,:,i);
%     img(~windowOutline) = 0;
%     data(:,:,i) = img;
% end

mw = loadMworksFile(subnum,expDate,expt(iexp).visStimTime_preDrug);
nOn = mw.nScansOn;
nOff = mw.nScansOff;
nTrialFr = nOn+nOff;
nTrials = nFr./nTrialFr;
data_trials = reshape(data,[ypix, xpix, nTrialFr, nTrials]);
f0 = uint8(mean(data_trials(:,:,(nOff-nBaselineFr+1):nOff,:),3));
df = (data_trials - f0);
dff = df./f0;
respwin = (nOff+2):(nOff+1+nRespFr);
meanRespImage = mean(mean(df(:,:,(nOff+2):(nOff+1+nRespFr),:),4),3);
figure;imagesc(meanRespImage);

injSiteTC_preDrug_red =  stackGetTimeCourses(data,redMask);
if doFullFieldStim
    retImage = readtiff(fullfile(fnout, sprintf('retROIguide.tif',expt(iexp).labelFolder)));
    fprintf('Select injection site vis stim ROI')
    visStimLabel =redLabel;
    fprintf('Select V1 ROI')
    V1Label = imCellPolyEditInteractive(retImage);
    fprintf('Select other ROI')
    otherHVALabel = imCellPolyEditInteractive(retImage);
    otherHVAName = input('Which control HVA did you select?');

    overlapMask = bwlabel((visStimLabel + redLabel) == 2);
    V1Mask = bwlabel(V1Label);
    otherHVAMask = bwlabel(otherHVALabel);
else
    fprintf('Select injection site vis stim ROI')
    visStimLabel = imCellPolyEditInteractive(meanRespImage);
    fprintf('Select V1 ROI')
    V1Label = imCellPolyEditInteractive(meanRespImage);
    fprintf('Select other ROI')
    otherHVALabel = imCellPolyEditInteractive(meanRespImage);
    otherHVAName = input('Which control HVA did you select?');

    overlapMask = bwlabel((visStimLabel + redLabel) == 2);
    V1Mask = bwlabel(V1Label);
    otherHVAMask = bwlabel(otherHVALabel);
end
injSiteTC_preDrug_visStim = stackGetTimeCourses(data,overlapMask);
V1_preDrug = stackGetTimeCourses(data,V1Mask);
other_preDrug = stackGetTimeCourses(data,otherHVAMask);
% 
% dff = df./f0;
% meanRespImage = mean(mean(dff(:,:,(nOff+2):(nOff+1+nRespFr),:),4),3);
wheelTicksTpS = cellfun(@(x) x./0.5,mw.wheelSpeedValues,'unif',0);
alignTrials = find(cellfun(@(x) length(x) ~= nTrialFr,wheelTicksTpS));
for i = 1:length(alignTrials)
    ind = alignTrials(i);
    if length(wheelTicksTpS{ind}) > nTrialFr
        chopTrialInd = (length(wheelTicksTpS{ind}) - nTrialFr + 1):length(wheelTicksTpS{ind});
        wheelTicksTpS{ind} = wheelTicksTpS{ind}(chopTrialInd);
    else
        addTrialN = (nTrialFr - length(wheelTicksTpS{ind}));
        wheelTicksTpS{ind} = cat(2,wheelTicksTpS{ind},nan(1,addTrialN));
    end
end
wheelTicksTpS = cell2mat(wheelTicksTpS');

clear data mw dff df f0
%%
respImg = nan(ypix,xpix,nrun+2);
respImg(:,:,1) = meanRespImage;


%% load all post-drug data

injSiteTC_red = cell(1,nrun);
injSiteTC_visStim = cell(1,nrun);
% mw_mat = cell(1,nrun);
V1TC = cell(1,nrun);
otherTC = cell(1,nrun);

for irun = 1:nrun
    folderName = sprintf('stim_drug_%s',expt(iexp).visStimFolder{irun});
    imageName = sprintf('%s_MMStack.ome.tif',folderName);
    data = readtiff(fullfile(fnin,folderName,imageName));
    tc = stackGetTimeCourses(data,redMask);
    injSiteTC_red{irun} = tc;
    
    injSiteTC_visStim{irun} = stackGetTimeCourses(data,overlapMask);
    V1TC{irun} = stackGetTimeCourses(data,V1Mask);
    otherTC{irun} = stackGetTimeCourses(data,otherHVAMask); 
    
    data_trials = reshape(data,[ypix, xpix, nTrialFr, nTrials]);
    f0 = uint8(mean(data_trials(:,:,(nOff-nBaselineFr+1):nOff,:),3));
    df = (data_trials - f0);
%     dff = df./f0;
    meanRespImage = mean(mean(df(:,:,(nOff+2):(nOff+1+nRespFr),:),4),3); 
        
    respImg(:,:,irun+1) = meanRespImage;     
    
end
clear data dff df f0

%% speed data
speedTpS = nan(nrun+1,nTrialFr);
speedErrTpS = nan(nrun+1,nTrialFr);
% speedTpS(1,:) = nanmean(wheelTicksTpS,1);
% speedErrTpS(1,:) = ste(wheelTicksTpS,1);

speedTpS_allTrials = cell(1,nrun);
for irun = 1:(nrun+1)      
    if irun == 1
        mw = loadMworksFile(subnum,expDate,expt(iexp).visStimTime_preDrug);
    else
        mw = loadMworksFile(subnum,expDate,expt(iexp).visStimTime{irun-1});
    end
    
    wheelTicksTpS = cellfun(@(x) x./0.5,mw.wheelSpeedValues,'unif',0);
    if length(wheelTicksTpS) > nTrials
        wheelTicksTpS = wheelTicksTpS(1:nTrials);
    end
    alignTrials = find(cellfun(@(x) length(x) ~= nTrialFr,wheelTicksTpS));
    for i = 1:length(alignTrials)
        ind = alignTrials(i);
        if length(wheelTicksTpS{ind}) > nTrialFr
            chopTrialInd = (length(wheelTicksTpS{ind}) - nTrialFr + 1):length(wheelTicksTpS{ind});
            wheelTicksTpS{ind} = wheelTicksTpS{ind}(chopTrialInd);
        else
            addTrialN = (nTrialFr - length(wheelTicksTpS{ind}));
            wheelTicksTpS{ind} = cat(2,wheelTicksTpS{ind},nan(1,addTrialN));
        end
    end
    speedTpS(irun,:) = nanmean(cell2mat(wheelTicksTpS'),1);
    speedErrTpS(irun,:) = ste(cell2mat(wheelTicksTpS'),1);
    speedTpS_allTrials{irun} = cell2mat(wheelTicksTpS');
end
%% stim-evoked mean dF
drugTimes = cellfun(@num2str,expt(iexp).visStimTimeFromDrugMin,'unif',0);
respImg(:,:,nrun+2) = redImage;
setFigParams4Print('landscape')
[nRows,nCols] = optimizeSubplotDim(nrun+2);
figure;colormap gray
for i = 1:(nrun+2)
    subplot(nRows,nCols,i)
    imagesc(respImg(:,:,i))
    if i == nrun+2
        title('mCherry')
    else
        title(drugTimes{i})
        clim([0 10])
    end
    axis square
end
print(fullfile(fnout,'respImagesAcrossDrugTimes'),'-dpdf','-fillpage')
writetiff(respImg,fullfile(fnout,'respImagesAcrossDrugTimes'));
%% stim-evoked activity
tt = double((-nOff+1):nOn)./params.frameRate;

basewin = (nOff-nBaselineFr+1):nOff;
respwin = (nOff+2):(nOff+1+nRespFr);

colors = brewermap(nrun+1,'Blues');
colors = colors(2:end,:);

injSiteTC_allTimes = cat(2,{injSiteTC_preDrug_visStim},injSiteTC_visStim );
V1TC_allTimes = cat(2, {V1_preDrug}, V1TC);
otherTC_allTimes = cat(2,{other_preDrug}, otherTC);

%% stim aligned running speed
setFigParams4Print('portrait')
figure
subplot 211
for i = 1:(nrun+1)
    hold on
    y = speedTpS(i,:);
    yerr = speedErrTpS(i,:);
    h = shadedErrorBar(tt,y,yerr,'k-');
    if i == 1
        leg = [];
        h.mainLine.Color = 'k';
        h.mainLine.LineWidth = 2;
        h.edge(1).LineStyle = 'none';
        h.edge(2).LineStyle = 'none';
    else
        h.mainLine.Color = colors(i-1,:);
        h.mainLine.LineWidth = 1;
        h.patch.FaceColor = h.mainLine.Color + ([1 1 1] - h.mainLine.Color)*0.5;
        h.edge(1).LineStyle = 'none';
        h.edge(2).LineStyle = 'none';
    end
    leg(i) = h.mainLine;
end
legend(leg,drugTimes,'location','northeastoutside')
figXAxis([],'Time (s)',[])
figYAxis([],'Speed (Ticks/S)',[])
figAxForm([],0)
title({[mouse '-' expDate '-' expt(iexp).drug];...
    'Stimulus Aligned Running'})

tc_trials = cellfun(@(x) reshape(x,[nTrialFr,nTrials]),injSiteTC_allTimes,'unif',0);
dff = cellfun(@(x) (x - mean(x(basewin,:),1))./x,tc_trials,'unif',0);
visStimResp = cellfun(@(x) mean(mean(x,2)),dff);

speedResp = mean(speedTpS(:,respwin),2);

subplot 212
for i = 1:(nrun+1)
    hold on
    h = plot(speedResp(i),visStimResp(i),'o');
    if i == 1
        h.Color = 'k';
        h.MarkerFaceColor = 'k';
        h.LineWidth = 2;
    else
        h.Color = colors(i-1,:);
        h.MarkerFaceColor = colors(i-1,:);
        h.LineWidth = 1;
    end
end
legend(drugTimes,'location','northeastoutside')
figXAxis([],'Speed (Ticks/S)',[])
figYAxis([],'Visual Resp (dF/F)',[])
figAxForm

print(fullfile(fnout,'running_drugResp_injSite'),'-dpdf','-fillpage')

visStimResp = cellfun(@(x) mean(x(respwin,:),1),dff,'unif',0);
speedResp = cellfun(@(x) mean(x(:,respwin),2),speedTpS_allTrials,'unif',0);
figure
for i = 1:(nrun+1)
    hold on
    h = plot(speedResp{i},visStimResp{i},'o');
    if i == 1
        h.Color = 'k';
        h.MarkerFaceColor = 'k';
        h.LineWidth = 2;
    else
        h.Color = colors(i-1,:);
        h.MarkerFaceColor = colors(i-1,:);
        h.LineWidth = 1;
    end
end
legend(drugTimes,'location','northeastoutside')
figXAxis([],'Speed (Ticks/S)',[])
figYAxis([],'Visual Resp (dF/F)',[])
figAxForm
%% plot responses
% in injection site 
tc_trials = cellfun(@(x) reshape(x,[nTrialFr,nTrials]),injSiteTC_allTimes,'unif',0);
dff = cellfun(@(x) (x - mean(x(basewin,:),1))./x,tc_trials,'unif',0);
visStimResp = cellfun(@(x) mean(x,2),dff,'unif',0);
visStimRespErr = cellfun(@(x) ste(x,2),dff,'unif',0);

setFigParams4Print('portrait')
figure
subplot 211
for i = 1:(nrun+1)
    hold on
    y = visStimResp{i};
    h = plot(tt,y,'-');
    if i == 1
        h.Color = 'k';
        h.LineWidth = 2;
    else
        h.Color = colors(i-1,:);
        h.LineWidth = 1;
    end
end
legend(drugTimes,'location','northeastoutside')
figXAxis([],'Time (s)',[])
figYAxis([],'dF/F',[])
figAxForm([],0)
title({[mouse '-' expDate '-' expt(iexp).inj_loc '-' expt(iexp).drug];...
    'Vis Stim Response At Injection Site'})

subplot 212
for i = 1:(nrun+1)
    hold on
    y = visStimResp{i};
    yerr = visStimRespErr{i};
    h = shadedErrorBar(tt,y,yerr,'k-');
    if i == 1
        h.mainLine.Color = 'k';
        h.mainLine.LineWidth = 2;
        h.edge(1).LineStyle = 'none';
        h.edge(2).LineStyle = 'none';
    else
        h.mainLine.Color = colors(i-1,:);
        h.mainLine.LineWidth = 1;
        h.patch.FaceColor = h.mainLine.Color + ([1 1 1] - h.mainLine.Color)*0.5;
        h.edge(1).LineStyle = 'none';
        h.edge(2).LineStyle = 'none';
    end
end
legend(drugTimes,'location','northeastoutside')
figXAxis([],'Time (s)',[])
figYAxis([],'dF/F',[])
figAxForm([],0)
title({[mouse '-' expDate '-' expt(iexp).inj_loc '-' expt(iexp).drug];...
    'Vis Stim Response At Injection Site'})
print(fullfile(fnout,'tc_drugResp_injSite'),'-dpdf','-fillpage')

% in V1
tc_trials = cellfun(@(x) reshape(x,[nTrialFr,nTrials]),V1TC_allTimes,'unif',0);
dff = cellfun(@(x) (x - mean(x(basewin,:),1))./x,tc_trials,'unif',0);
visStimResp = cellfun(@(x) mean(x,2),dff,'unif',0);
visStimRespErr = cellfun(@(x) ste(x,2),dff,'unif',0);

colors = brewermap(nrun+1,'Blues');
colors = colors(2:end,:);

setFigParams4Print('portrait')
figure
subplot 211
for i = 1:(nrun+1)
    hold on
    y = visStimResp{i};
    h = plot(tt,y,'-');
    if i == 1
        h.Color = 'k';
        h.LineWidth = 2;
    else
        h.Color = colors(i-1,:);
    h.LineWidth = 1;
    end
end
legend(drugTimes,'location','northeastoutside')
figXAxis([],'Time (s)',[])
figYAxis([],'dF/F',[])
figAxForm([],0)
title({[mouse '-' expDate '-V1-' expt(iexp).drug];...
    'Vis Stim Response'})

subplot 212
for i = 1:(nrun+1)
    hold on
    y = visStimResp{i};
    yerr = visStimRespErr{i};
    h = shadedErrorBar(tt,y,yerr,'k-');
    if i == 1
        h.mainLine.Color = 'k';
        h.mainLine.LineWidth = 2;
        h.edge(1).LineStyle = 'none';
        h.edge(2).LineStyle = 'none';
    else
        h.mainLine.Color = colors(i-1,:);
        h.mainLine.LineWidth = 1;
        h.patch.FaceColor = h.mainLine.Color + ([1 1 1] - h.mainLine.Color)*0.5;
        h.edge(1).LineStyle = 'none';
        h.edge(2).LineStyle = 'none';
    end
end
legend(drugTimes,'location','northeastoutside')
figXAxis([],'Time (s)',[])
figYAxis([],'dF/F',[])
figAxForm([],0)
title({[mouse '-' expDate '-V1-' expt(iexp).drug];...
    'Vis Stim Response'})
print(fullfile(fnout,'tc_drugResp_V1'),'-dpdf','-fillpage')

% in HVA
tc_trials = cellfun(@(x) reshape(x,[nTrialFr,nTrials]),otherTC_allTimes,'unif',0);
dff = cellfun(@(x) (x - mean(x(basewin,:),1))./x,tc_trials,'unif',0);
visStimResp = cellfun(@(x) mean(x,2),dff,'unif',0);
visStimRespErr = cellfun(@(x) ste(x,2),dff,'unif',0);

colors = brewermap(nrun+1,'Blues');
colors = colors(2:end,:);

setFigParams4Print('portrait')
figure
subplot 211
for i = 1:(nrun+1)
    hold on
    y = visStimResp{i};
    h = plot(tt,y,'-');
    if i == 1
        h.Color = 'k';
        h.LineWidth = 2;
    else
        h.Color = colors(i-1,:);
    h.LineWidth = 1;
    end
end
legend(drugTimes,'location','northeastoutside')
figXAxis([],'Time (s)',[])
figYAxis([],'dF/F',[])
figAxForm([],0)
title({[mouse '-' expDate '-' otherHVAName '-' expt(iexp).drug];...
    'Vis Stim Response'})

subplot 212
for i = 1:(nrun+1)
    hold on
    y = visStimResp{i};
    yerr = visStimRespErr{i};
    h = shadedErrorBar(tt,y,yerr,'k-');
    if i == 1
        h.mainLine.Color = 'k';
        h.mainLine.LineWidth = 2;
        h.edge(1).LineStyle = 'none';
        h.edge(2).LineStyle = 'none';
    else
        h.mainLine.Color = colors(i-1,:);
        h.mainLine.LineWidth = 1;
        h.patch.FaceColor = h.mainLine.Color + ([1 1 1] - h.mainLine.Color)*0.5;
        h.edge(1).LineStyle = 'none';
        h.edge(2).LineStyle = 'none';
    end
end
legend(drugTimes,'location','northeastoutside')
figXAxis([],'Time (s)',[])
figYAxis([],'dF/F',[])
figAxForm([],0)
title({[mouse '-' expDate '-' otherHVAName '-' expt(iexp).drug];...
    'Vis Stim Response'})
print(fullfile(fnout,'tc_drugResp_otherHVA'),'-dpdf','-fillpage')
