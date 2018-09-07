clear all
close all
ds = 'xori_emx_widefield';
rc = behavConstsAV;
eval(ds)
iexp = 3;
%%
nBaselineFr = round(...
    params.nBaselineMs.*params.frameRate./1000);
nRespFr = round(...
    params.nPostStimMs.*params.frameRate./1000);
nTrialFr = nBaselineFr+nRespFr;

%% experiment info

mouse = expt(iexp).mouse;
subnum = mouse;
expDate = expt(iexp).date;

fnout = fullfile(rc.ashleyAnalysis,mouse,'widefield imaging',expDate);
fnin = fullfile(rc.ashleyData,mouse,'widefield imaging',expDate);

%% load data, get important experiment info

images = readtiff(fullfile(fnin,expt(iexp).visStimFolder, ...
    sprintf('%s_MMStack.ome.tif',expt(iexp).visStimFolder)));
exptInfo = loadMworksFile(subnum,expDate,expt(iexp).visStimTime);

tStimOnFrames = cell2mat(exptInfo.cStimOneOn);
stimOnTimeMs = double(exptInfo.stimOneGratingOnTimeMs);
stimOnTimeFr = round(stimOnTimeMs.*(1/1000).*params.frameRate);

respwin = (nBaselineFr+1):(nBaselineFr+stimOnTimeFr);
basewin = 1:nBaselineFr;

[yPix,xPix,nFrames] = size(images);
nTrials = length(tStimOnFrames);

stimOri = cell2mat(exptInfo.tStimOneGratingDirectionDeg);
maskOri = cell2mat(exptInfo.tMaskOneGratingDirectionDeg);
stimCon = cell2mat(cellfun(@double,exptInfo.tStimOneGratingContrast,'unif',0));
maskCon = cell2mat(cellfun(@double,exptInfo.tMaskOneGratingContrast,'unif',0));

%% trial types - get indexes for each unique condition

orientations = unique(stimOri);
nOri = length(orientations);
contrasts = unique(stimCon+maskCon);
nCon = length(contrasts);
nX = nchoosek(nOri,nOri);

nTrialTypes = 5;
trialTypeName = {'40%-0'; '40%-90';'20%-0';'20%-90';'x-Ori'};
trialTypeInd{1} = stimOri == 0 & stimOri == maskOri & stimCon == 0.2 & maskCon == 0.2;
trialTypeInd{2} = stimOri == 90 & stimOri == maskOri & stimCon == 0.2 & maskCon == 0.2;
trialTypeInd{3}= stimCon == 0.2 & maskCon == 0 & stimOri == 0;
trialTypeInd{4} = stimCon == 0.2 & maskCon == 0 & stimOri == 90;
trialTypeInd{5} = stimOri ~= maskOri & stimCon == 0.2 & maskCon == 0.2;

%% get frames 2s before and 2s after trial starts

images_trialAligned = nan(yPix, xPix, nTrialFr, nTrials);
for itrial = 1:nTrials
    ind = (tStimOnFrames(itrial) - nBaselineFr +1):(tStimOnFrames(itrial) + nRespFr);
    images_trialAligned(:,:,:,itrial) = images(:,:,ind);
end

images_baseF = mean(images_trialAligned(:,:,basewin,:),3);
images_dFF = (images_trialAligned - images_baseF)./images_baseF;
clear images_baseF images_trialAligned

%%
images_trialType = cell(1,5);
for itype = 1:nTrialTypes
    ind = trialTypeInd{itype};
    images_trialType{itype} = images_dFF(:,:,:,ind);
end

imageTC_allType = cellfun(@(x) mean(x,4),images_trialType,'unif',0);
imageResp_allType = cellfun(@(x) mean(mean(x(:,:,respwin,:),4),3),...
    images_trialType,'unif',0);

for itype = 1:nTrialTypes
    writetiff(imageTC_allType{itype},fullfile(fnout,...
        ['tc_' trialTypeName{itype}]));
    writetiff(imageResp_allType{itype},fullfile(fnout,...
        ['resp_' trialTypeName{itype}]));
end

%% select ROIs for visual cortex areas and extract time-courses
ROIimage = imadjust(imageResp_allType{1});
drawnROIS = imCellPolyEditInteractive(ROIimage);
roiMask = bwlabel(drawnROIS);
figure;imagesc(roiMask);colorbar
roiLabels = {'V1';'AL';'RL';'PM';'AM'};
nRoi = length(roiLabels);

tc = stackGetTimeCourses(images,roiMask);
tc_trialAligned = nan(nTrialFr,nRoi,nTrials);
for itrial = 1:nTrials
    ind = (tStimOnFrames(itrial) - nBaselineFr +1):(tStimOnFrames(itrial) + nRespFr);
    tc_trialAligned(:,:,itrial) = tc(ind,:);
end

tc_F = mean(tc_trialAligned(basewin,:,:),1);
tc_dFF = (tc_trialAligned - tc_F)./tc_F;

%% get mean time-courses
tc_allTypes = cell(1,5);
for itype = 1:nTrialTypes
    tc_allTypes{itype} = tc_dFF(:,:,trialTypeInd{itype});
end
tc_mean_all = cellfun(@(x) mean(x,3),tc_allTypes,'unif',0);

resp_allTypes = cellfun(@(x) mean(mean(x(respwin,:,:),3),1),...
    tc_allTypes,'unif',0);
resp_err_allTypes = cellfun(@(x) mean(ste(x(respwin,:,:),3),1),...
    tc_allTypes,'unif',0);

%% plot time-courses for each stimulus condition
timeFromStimMs = ((-nBaselineFr+1):nRespFr)./params.frameRate.*1000;
respwinMs = ([respwin(1) respwin(2)]-nBaselineFr)./params.frameRate.*1000;
colors = [0 0 0;...
          0.5 0.5 0.5;...          
          0.5 1 0.5;...          
          0 1 0;...
          1 0 0;...
          1 0.5 0.5];

setFigParams4Print('portrait')
figure
suptitle('everything')
for iroi = 1:nRoi
    for itype = 1:nTrialTypes+1
        if itype == nTrialTypes+1
            subplot(3,2,iroi)
            h = plot(itype,...
            resp_allTypes{3}(iroi)+resp_allTypes{4}(iroi),...
            'ko');
            hold on            
        else
            subplot(3,2,iroi)
            h = errorbar(itype,...
            resp_allTypes{itype}(iroi),resp_err_allTypes{itype}(iroi),...
            'ko');
            hold on
        end
    end
    title(roiLabels{iroi})
    figXAxis([],'',[0 nTrialTypes+2],1:nTrialTypes+1,cat(1,trialTypeName,{'sum'}))
    xtickangle(30)
    figYAxis([],'dF/F',[0 0.03])
    figAxForm([],0)
end
print(fullfile(fnout,'xoriWFresults2'),'-dpdf','-fillpage')
      

for iroi = 1:nRoi
    legend_label = [];
    figure
    suptitle(roiLabels{iroi})
    subplot 211
    for itype = 1:nType
    h = shadedErrorBar_chooseColor(timeFromStimMs,...
        tc_100con_0(:,iroi),tc_100con_0_err(:,iroi),colors(1,:));
    legend_label(1) = h.mainLine;
    hold on
    vline(respwinMs,'k--')
    legend(legend_label,{'100% 0deg';'100% 90deg';'50% 0deg';'50% 90deg';'x-Ori'},'location','northeastoutside')
    figXAxis([],'Time (ms)',[])
    figYAxis([],'dF/F',[])
    figAxForm([],0)
    
    subplot 212
    h = errorbar(1,resp_100con(iroi),resp_100con_err(iroi),'o');
    h.Color = colors(1,:);
    hold on
    h = errorbar(2,resp_50con_0(iroi),resp_50con_0_err(iroi),'o');
    h.Color = colors(2,:);
    h = errorbar(2,resp_50con_90(iroi),resp_50con_90_err(iroi),'o');
    h.Color = colors(3,:);
    h = errorbar(3,resp_xori(iroi),resp_xori_err(iroi),'o');
    h.Color = colors(4,:);
    legend({'100%';'50% 0deg';'50% 90deg';'x-Ori'},'location','northeastoutside')
    figXAxis([],'',[0 4])
    figYAxis([],'dF/F',[])
    figAxForm
end
      
%%
      
      
tc_100con_0 = mean(tc_dFF(:,:,trInd_100con_0),3);
tc_100con_90 = mean(tc_dFF(:,:,trInd_100con_90),3);
tc_50con_0 = mean(tc_dFF(:,:,trInd_50con_0),3);
tc_50con_90 = mean(tc_dFF(:,:,trInd_50con_90),3);
tc_xori = mean(tc_dFF(:,:,trInd_xori),3);

tc_100con_0_err = ste(tc_dFF(:,:,trInd_100con_0),3);
tc_100con_90_err = ste(tc_dFF(:,:,trInd_100con_90),3);
tc_50con_0_err = ste(tc_dFF(:,:,trInd_50con_0),3);
tc_50con_90_err = ste(tc_dFF(:,:,trInd_50con_90),3);
tc_xori_err = ste(tc_dFF(:,:,trInd_xori),3);

resp_100con_0 = mean(mean(tc_dFF(respwin,:,trInd_100con_0),1),3);
resp_100con_90 = mean(mean(tc_dFF(respwin,:,trInd_100con_90),1),3);
resp_50con_0 = mean(mean(tc_dFF(respwin,:,trInd_50con_0),1),3);
resp_50con_90 = mean(mean(tc_dFF(respwin,:,trInd_50con_90),1),3);
resp_xori = mean(mean(tc_dFF(respwin,:,trInd_xori),1),3);

resp_100con_0_err = ste(mean(tc_dFF(respwin,:,trInd_100con_0),1),3);
resp_100con_90_err = ste(mean(tc_dFF(respwin,:,trInd_100con_90),1),3);
resp_50con_0_err = ste(mean(tc_dFF(respwin,:,trInd_50con_0),1),3);
resp_50con_90_err = ste(mean(tc_dFF(respwin,:,trInd_50con_90),1),3);
resp_xori_err = ste(mean(tc_dFF(respwin,:,trInd_xori),1),3);

for iroi = 1:nRoi
    legend_label = [];
    figure
    suptitle(roiLabels{iroi})
    subplot 211
    h = shadedErrorBar_chooseColor(timeFromStimMs,...
        tc_100con_0(:,iroi),tc_100con_0_err(:,iroi),colors(1,:));
    legend_label(1) = h.mainLine;
    hold on
    h = shadedErrorBar_chooseColor(timeFromStimMs,...
        tc_100con_90(:,iroi),tc_100con_90_err(:,iroi),colors(2,:));
    legend_label(2) = h.mainLine;
    h = shadedErrorBar_chooseColor(timeFromStimMs,...
        tc_50con_0(:,iroi),tc_50con_0_err(:,iroi),colors(3,:));
    legend_label(3) = h.mainLine;
    h = shadedErrorBar_chooseColor(timeFromStimMs,...
        tc_50con_90(:,iroi),tc_50con_90_err(:,iroi),colors(4,:));
    legend_label(4) = h.mainLine;
    h = shadedErrorBar_chooseColor(timeFromStimMs,...
        tc_xori(:,iroi),tc_xori_err(:,iroi),colors(5,:));
    legend_label(5) = h.mainLine;
    vline(respwinMs,'k--')
    legend(legend_label,{'100% 0deg';'100% 90deg';'50% 0deg';'50% 90deg';'x-Ori'},'location','northeastoutside')
    figXAxis([],'Time (ms)',[])
    figYAxis([],'dF/F',[])
    figAxForm([],0)
    
    subplot 212
    h = errorbar(1,resp_100con(iroi),resp_100con_err(iroi),'o');
    h.Color = colors(1,:);
    hold on
    h = errorbar(2,resp_50con_0(iroi),resp_50con_0_err(iroi),'o');
    h.Color = colors(2,:);
    h = errorbar(2,resp_50con_90(iroi),resp_50con_90_err(iroi),'o');
    h.Color = colors(3,:);
    h = errorbar(3,resp_xori(iroi),resp_xori_err(iroi),'o');
    h.Color = colors(4,:);
    legend({'100%';'50% 0deg';'50% 90deg';'x-Ori'},'location','northeastoutside')
    figXAxis([],'',[0 4])
    figYAxis([],'dF/F',[])
    figAxForm
end
colors = [0 0 0;...         
          0.5 1 0.5;...          
          0 1 0;...
          1 0.5 0.5;...
          1 0 0];
figure
subplot 121
iroi = 1;
h = errorbar(1,resp_100con_0(iroi),resp_100con_0_err(iroi),'o');
h.Color = colors(1,:);
hold on
h = errorbar(2,resp_50con_0(iroi),resp_50con_0_err(iroi),'o');
h.Color = colors(2,:);
h = errorbar(3,resp_50con_90(iroi),resp_50con_90_err(iroi),'o');
h.Color = colors(3,:);
h = errorbar(4,resp_50con_0(iroi)+resp_50con_90(iroi),...
    resp_50con_0_err(iroi)+resp_50con_90_err(iroi),'o');
h.Color = colors(4,:);
h = errorbar(5,resp_xori(iroi),resp_xori_err(iroi),'o');
h.Color = colors(5,:);

legend({'100% 0deg';'50% 0deg';'50% 90deg';'50%+50%';'x-Ori'},'location','northeastoutside')
figXAxis([],'',[0 6])
figYAxis([],'dF/F',[])
figAxForm
title(roiLabels{iroi})

subplot 122
iroi = 3;
h = errorbar(1,resp_100con_0(iroi),resp_100con_0_err(iroi),'o');
h.Color = colors(1,:);
hold on
h = errorbar(2,resp_50con_0(iroi),resp_50con_0_err(iroi),'o');
h.Color = colors(2,:);
h = errorbar(3,resp_50con_90(iroi),resp_50con_90_err(iroi),'o');
h.Color = colors(3,:);
h = errorbar(4,resp_50con_0(iroi)+resp_50con_90(iroi),...
    resp_50con_0_err(iroi)+resp_50con_90_err(iroi),'o');
h.Color = colors(4,:);
h = errorbar(5,resp_xori(iroi),resp_xori_err(iroi),'o');
h.Color = colors(5,:);

legend({'100% 0deg';'50% 0deg';'50% 90deg';'50%+50%';'x-Ori'},'location','northeastoutside')
figXAxis([],'',[0 6])
figYAxis([],'dF/F',[])
figAxForm
title(roiLabels{iroi})

