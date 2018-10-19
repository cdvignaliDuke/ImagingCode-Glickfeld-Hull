clear all
close all
ds = 'dirOri_emx_widefield';
rc = behavConstsAV;
eval(ds)
iexp = 2;
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
%% get stimulus-driven image for each trial type

images_trialAligned = nan(yPix, xPix, nTrialFr, nTrials);
for itrial = 1:nTrials
    ind = (tStimOnFrames(itrial) - nBaselineFr +1):(tStimOnFrames(itrial) + nRespFr);
    images_trialAligned(:,:,:,itrial) = images(:,:,ind);
end

images_baseF = mean(images_trialAligned(:,:,basewin,:),3);
images_dFF = (images_trialAligned - images_baseF)./images_baseF;
clear images_baseF images_trialAligned

trInd_100con = stimOri == maskOri & stimCon == 0.5 & maskCon == 0.5;
trInd_100con_0 = stimOri == 0 & stimOri == maskOri & stimCon == 0.5 & maskCon == 0.5;
trInd_100con_90 = stimOri == 90 & stimOri == maskOri & stimCon == 0.5 & maskCon == 0.5;
trInd_50con_0 = stimCon == 0.5 & maskCon == 0 & stimOri == 0;
trInd_50con_90 = stimCon == 0.5 & maskCon == 0 & stimOri == 90;
trInd_xori = stimOri ~= maskOri & stimCon == 0.5 & maskCon == 0.5;

images_100con = images_dFF(:,:,:,...
    trInd_100con);
images_50con_0 = images_dFF(:,:,:,...
    trInd_50con_0);
images_50con_90 = images_dFF(:,:,:,...
    trInd_50con_90);
images_xori = images_dFF(:,:,:,...
    trInd_xori);

imgTC_100con = mean(images_100con,4);
imgTC_50con_0 = mean(images_50con_0,4);
imgTC_50con_90 = mean(images_50con_90,4);
imgTC_xori = mean(images_xori,4);

img_100con = mean(mean(images_100con(:,:,respwin,:),4),3);
img_50con_0 = mean(mean(images_50con_0(:,:,respwin,:),4),3);
img_50con_90 = mean(mean(images_50con_90(:,:,respwin,:),4),3);
img_xori = mean(mean(images_xori(:,:,respwin,:),4),3);

if ~exist(fnout,'dir')
    mkdir(fnout)
end

writetiff(imgTC_100con,fullfile(fnout,'tc_100con'));
writetiff(imgTC_50con_0,fullfile(fnout,'tc_50con_0'));
writetiff(imgTC_50con_90,fullfile(fnout,'tc_50con_90'));
writetiff(imgTC_xori,fullfile(fnout,'tc_xori'));

writetiff(img_100con,fullfile(fnout,'resp_100con'));
writetiff(img_50con_0,fullfile(fnout,'resp_50con_0'));
writetiff(img_50con_90,fullfile(fnout,'resp_50con_90'));
writetiff(img_xori,fullfile(fnout,'resp_xori'));

figure; colormap gray
subplot 221
imagesc(img_100con);
axis square
title('100% Contrast')
subplot 222
imagesc(img_50con_0);
axis square
title('50% Contrast')
subplot 223
imagesc(img_xori);
axis square
subplot 222
imagesc(img_50con_90);
axis square
title('Cross Orientation 50/50% Contrast')

%% select ROIs for visual cortex areas and extract time-courses
ROIimage = imadjust(img_100con);
drawnROIS = imCellPolyEditInteractive(ROIimage);
roiMask = bwlabel(drawnROIS);
figure;imagesc(roiMask);colorbar
roiLabels = {'V1';'M';'PM';'AL';'AM'};
nRoi = length(roiLabels);

tc = stackGetTimeCourses(images,roiMask);
tc_trialAligned = nan(nTrialFr,nRoi,nTrials);
for itrial = 1:nTrials
    ind = (tStimOnFrames(itrial) - nBaselineFr +1):(tStimOnFrames(itrial) + nRespFr);
    tc_trialAligned(:,:,itrial) = tc(ind,:);
end

tc_F = mean(tc_trialAligned(basewin,:,:),1);
tc_dFF = (tc_trialAligned - tc_F)./tc_F;

%% plot time-courses for each stimulus condition
timeFromStimMs = ((-nBaselineFr+1):nRespFr)./params.frameRate.*1000;
respwinMs = ([respwin(1) respwin(2)]-nBaselineFr)./params.frameRate.*1000;
colors = [0 0 0;...
          0.5 0.5 0.5;...          
          0.5 1 0.5;...          
          0 1 0;...
          1 0 0];

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

