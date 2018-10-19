clear all
close all
ds = 'dirOri_emx_widefield';
rc = behavConstsAV;
eval(ds);
iexp = 1;
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

%% get frames 2s before and 2s after trial starts

images_trialAligned = nan(yPix, xPix, nTrialFr, nTrials);
for itrial = 1:nTrials
    ind = (tStimOnFrames(itrial) - nBaselineFr +1):(tStimOnFrames(itrial) + nRespFr);
    images_trialAligned(:,:,:,itrial) = images(:,:,ind);
end

images_baseF = mean(images_trialAligned(:,:,basewin,:),3);
images_dFF = (images_trialAligned - images_baseF)./images_baseF;
clear images_baseF images_trialAligned

%% get stimulus types and trial indexes for each type
movementInd = logical(cell2mat(exptInfo.tBlock2TrialNumber));
tDir = cell2mat(exptInfo.tStimOneGratingDirectionDeg);

nTrialTypes = 6;
trialTypeNames = {'S-0';'S-90';'D-0';'D-90';...
    'D-180';'D-270'};
trialTypeInd{1} = movementInd == 0 & tDir == 0;
trialTypeInd{2} = movementInd == 0 & tDir == 90;
trialTypeInd{3} = movementInd == 1 & tDir == 0;
trialTypeInd{4} = movementInd == 1 & tDir == 90;
trialTypeInd{5} = movementInd == 1 & tDir == 180;
trialTypeInd{6} = movementInd == 1 & tDir == 270;

%%
images_types = cell(1,nTrialTypes);
for itype = 1:nTrialTypes
    images_types{itype} = ...
        mean(mean(images_dFF(:,:,respwin,trialTypeInd{itype}),4),3);
end

for itype = 1:nTrialTypes
    writetiff(images_types{itype},fullfile(fnout,...
        ['resp_' trialTypeNames{itype}]));
end

%% which image to find ROIs?
image_mat_max = max(reshape(cell2mat(images_types),[yPix xPix nTrialTypes])...
    ,[],3);
image_mat_mean = mean(reshape(cell2mat(images_types),[yPix xPix nTrialTypes]),3);

ROIimage = imadjust(image_mat_mean);
figure;imagesc(ROIimage)
%% select ROIs for visual cortex areas and extract time-courses
ROIimage = imadjust(image_mat_mean);
drawnROIS = imCellPolyEditInteractive(ROIimage);
roiMask = bwlabel(drawnROIS);

setFigParams4Print('landscape')
figure;subplot 121
imagesc(roiMask);colorbar
axis square;title('ROI mask')
subplot 122
imagesc(imadjust(image_mat_mean))
colorbar
axis square;title('adjusted mean stim-driven response')
print(fullfile(fnout,'maskAndSelectionImage'),'-dpdf','-fillpage')
%% enter manually
roiLabels = {'V1';'LM';'AL';'PM';'RL';'AM';'Mys'};
%%
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
tc_allTypes = cell(1,nTrialTypes);
for itype = 1:nTrialTypes
    tc_allTypes{itype} = tc_dFF(:,:,trialTypeInd{itype});
end
tc_mean_all = cellfun(@(x) mean(x,3),tc_allTypes,'unif',0);
tc_err_all = cellfun(@(x) ste(x,3),tc_allTypes,'unif',0);

resp_allTypes = cellfun(@(x) mean(mean(x(respwin,:,:),3),1),...
    tc_allTypes,'unif',0);
resp_err_allTypes = cellfun(@(x) mean(ste(x(respwin,:,:),3),1),...
    tc_allTypes,'unif',0);

%% plot time-courses for each stimulus condition
timeFromStimMs = ((-nBaselineFr+1):nRespFr)./params.frameRate.*1000;
respwinMs = ([respwin(1) respwin(2)]-nBaselineFr)./params.frameRate.*1000;

setFigParams4Print('landscape')
figure
suptitle('90 deg, static-blk, drift-red')
for iroi = 1:nRoi
    subplot(3,3,iroi)
    h = shadedErrorBar_chooseColor(timeFromStimMs,...
        tc_mean_all{2}(:,iroi),tc_err_all{1}(:,iroi),[0 0 0]);
    hold on
    h = shadedErrorBar_chooseColor(timeFromStimMs,...
        tc_mean_all{4}(:,iroi),tc_err_all{3}(:,iroi),[1 0 0]);
    title(roiLabels{iroi})
end
print(fullfile(fnout,'tc_statDrift_90'),'-dpdf','-fillpage')
figure
suptitle('0 deg, static-blk, drift-red')
for iroi = 1:nRoi
    subplot(3,3,iroi)
    h = shadedErrorBar_chooseColor(timeFromStimMs,...
        tc_mean_all{1}(:,iroi),tc_err_all{1}(:,iroi),[0 0 0]);
    hold on
    h = shadedErrorBar_chooseColor(timeFromStimMs,...
        tc_mean_all{3}(:,iroi),tc_err_all{3}(:,iroi),[1 0 0]);
    title(roiLabels{iroi})
    figXAxis([],'Time (ms)',[])
    figYAxis([],'dF/F',[])
    figAxForm([],0)
end
print(fullfile(fnout,'tc_statDrift_0'),'-dpdf','-fillpage')

figure
suptitle('0/180 drifting, optic flow?,0-red,180-dark red')
for iroi = 1:nRoi
    subplot(3,3,iroi)
    h = shadedErrorBar_chooseColor(timeFromStimMs,...
        tc_mean_all{5}(:,iroi),tc_err_all{1}(:,iroi),[0.5 0 0]);
    hold on
    h = shadedErrorBar_chooseColor(timeFromStimMs,...
        tc_mean_all{3}(:,iroi),tc_err_all{3}(:,iroi),[1 0 0]);
    title(roiLabels{iroi})
    figXAxis([],'Time (ms)',[])
    figYAxis([],'dF/F',[])
    figAxForm([],0)
end
print(fullfile(fnout,'tc_drift_vert'),'-dpdf','-fillpage')

figure
suptitle('everything')
for iroi = 1:nRoi;
    for itype = 1:nTrialTypes
        iplot = iroi;
        subplot(3,3,iplot)
        h = errorbar(itype,...
        resp_allTypes{itype}(iroi),resp_err_allTypes{itype}(iroi),...
        'ko');
        hold on
    end
    title(roiLabels{iroi})
    figXAxis([],'',[0 nTrialTypes+1],1:nTrialTypes,trialTypeNames)
    xtickangle(45)
    figYAxis([],'dF/F',[0 0.04])
    figAxForm([],0)
end
print(fullfile(fnout,'resp_allCon'),'-dpdf','-fillpage')

%%
wfExpt = struct;
wfExpt.mouse = subnum;
wfExpt.date = expDate;
wfExpt.roiLabels = roiLabels;
wfExpt.stimInfo.stimulusTypes = trialTypeNames;
wfExpt.stimInfo.stimOn_frames = stimOnTimeFr;
wfExpt.stimInfo.stimOn_ms = stimOnTimeMs;
wfExpt.stimInfo.baseline_frames = nBaselineFr;
wfExpt.tc_roiXstimType = tc_allTypes;

save(fullfile(fnout,'exptSum_struct'),'wfExpt')
 