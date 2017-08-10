clear all
close all
ds = '_V1gad';
rc = behavConstsAV;
eval(['awFSAVdatasets' ds])
slct_expt = 2;
%%
iexp = slct_expt;
subnum = expt(iexp).SubNum;
mouse = expt(iexp).mouse;
expDate = expt(iexp).date;
retFolder = expt(iexp).rettuning{1};
retTime = expt(iexp).rettuning{2};
downSampleRate = 10;
%% load data & register
visStimData = loadMworksFile(subnum,expDate,retTime);
fName = [retFolder '_000_000'];
retData = loadsbx_choosepmt(1,mouse,expDate,retFolder,fName);

fn = fullfile(rc.ashleyAnalysis,mouse,'two-photon imaging', expDate);
load(fullfile(fn,'regOuts&Img.mat'))

retDataDownSampled = stackGroupProject(retData,downSampleRate);
[~,retDataRegistered] = stackRegister(retDataDownSampled,data_corr_img);

%% get trial types
tAzimuth = cell2mat_padded(visStimData.tGratingAzimuthDeg);
tElevation = cell2mat_padded(visStimData.tGratingElevationDeg);
ntrials = length(tAzimuth);

azimuths = unique(tAzimuth);
nazimuths = length(azimuths);
elevations = unique(tElevation);
nelevations = length(elevations);

npositions = nazimuths*nelevations;

positionInd = zeros(1,ntrials);
azimuthAtPosInd = nan(1,npositions);
elevationAtPosInd = nan(1,npositions);
for iaz = 1:nazimuths
    if iaz == 1
        ipos = 1;
    end
    for iel = 1:nelevations
        azInd = tAzimuth == azimuths(iaz);
        elInd = tElevation == elevations(iel);
        positionInd(azInd & elInd) = ipos;
        azimuthAtPosInd(ipos) = azimuths(iaz);
        elevationAtPosInd(ipos) = elevations(iel);
        ipos = ipos+1;
    end    
end
%% image params
frameRateHz = expt(iexp).frame_rate/downSampleRate;
[ypix,xpix,nframes] = size(retDataRegistered);
nFramesOn = visStimData.nScansOn/downSampleRate;
nFramesOff = visStimData.nScansOff/downSampleRate;

nBaselineFrames = ceil(nFramesOff/2);
nTrialFrames = nFramesOn;
respWindowFrames = nBaselineFrames+1:nBaselineFrames+1+ceil(0.5*frameRateHz);
trialStartInd = nFramesOff+1:nFramesOn+nFramesOff:nframes;
%% get image dF/F for each trial
retDFF = getTrialDffImages(retDataRegistered,trialStartInd,...
    nBaselineFrames,nTrialFrames);

resp2StimDFF = squeeze(mean(retDFF(:,:,respWindowFrames,:),3));
if size(resp2StimDFF,3) < ntrials
    ntrials = size(resp2StimDFF,3);
    positionInd = positionInd(1:ntrials);
end

%% for each stimulus type
positionRespDFF = nan(ypix,xpix,npositions);
for ipos = 1:npositions
    thisPositionInd = positionInd == ipos;
    positionRespDFF(:,:,ipos) = mean(resp2StimDFF(:,:,thisPositionInd),3);
end

azimuthRespDFF = nan(ypix,xpix,nazimuths);
for iaz = 1:nazimuths
    thisAzimuthInd = find(azimuthAtPosInd == azimuths(iaz));
    thisAzimuthTrialInd = ismember(positionInd,thisAzimuthInd);
    azimuthRespDFF(:,:,iaz) = mean(resp2StimDFF(:,:,thisAzimuthTrialInd),3);
end

elevationRespDFF = nan(ypix,xpix,nelevations);
for iel = 1:nelevations
    thisElevationInd = find(elevationAtPosInd == elevations(iel));
    thisElevationTrialInd = ismember(positionInd,thisElevationInd);
    elevationRespDFF(:,:,iel) = mean(resp2StimDFF(:,:,thisElevationTrialInd),3);
end
%% normalize each image in position stack
positionRespDFFNorm = normalizeImageStack(positionRespDFF);
azimuthRespDFFNorm = normalizeImageStack(azimuthRespDFF);
elevationRespDFFNorm = normalizeImageStack(elevationRespDFF);

%% plot images
if ~exist(fullfile(fn,retFolder),'dir')
    mkdir(fullfile(fn,retFolder))
end
positionTitles = arrayfun(@(x,y) sprintf('Az: %s, El: %s',num2str(x),...
       num2str(y)),azimuthAtPosInd,elevationAtPosInd,'unif',0);
azimuthTitles = arrayfun(@(x) sprintf('Az: %s',num2str(x)),azimuths,'unif',0);
elevationTitles = arrayfun(@(x) sprintf('El: %s',num2str(x)),elevations,'unif',0);

respPositionFig = plotImageStackAsSubplots(positionRespDFFNorm,positionTitles);
print(fullfile(fn,retFolder,'positionMeanDFF'),'-dpdf','-fillpage')
respAzimuthFig = plotImageStackAsSubplots(azimuthRespDFFNorm,azimuthTitles);
print(fullfile(fn,retFolder,'azimuthMeanDFF'),'-dpdf','-fillpage')
respElevationFig = plotImageStackAsSubplots(elevationRespDFFNorm,elevationTitles);
print(fullfile(fn,retFolder,'elevationMeanDFF'),'-dpdf','-fillpage')

plotRGBImage(azimuthRespDFFNorm);
title(sprintf('azimuth: %s (RGB)',num2str(azimuths)))
print(fullfile(fn,retFolder,'azimuthRGB'),'-dpdf','-fillpage')
plotRGBImage(elevationRespDFFNorm)
title(sprintf('elevation: %s (RGB)',num2str(elevations)))
print(fullfile(fn,retFolder,'elevationRGB'),'-dpdf','-fillpage')