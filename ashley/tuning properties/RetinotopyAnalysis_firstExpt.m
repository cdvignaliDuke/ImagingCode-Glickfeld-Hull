clear all
close all
ds = 'awFSAVdatasets_temp';
rc = behavConstsAV;
eval(ds)
iexp = 6;

downSampleRate = 10;

SubNum = expt(iexp).SubNum;
mouse = expt(iexp).mouse;
expDate = expt(iexp).date;
runFolder = expt(iexp).rettuning{1};
expTime = expt(iexp).rettuning{2};
fName = [runFolder '_000_000'];

data = loadsbx_choosepmt(1,mouse,expDate,runFolder,fName);
visStimData = loadMworksFile(SubNum,expDate,expTime);
data_sub = data-min(min(min(data,[],1),[],2),[],3);
data = data_sub;
clear data_sub

nfr = size(data,3);

fnout = fullfile(rc.ashleyAnalysis,mouse,'two-photon imaging',expDate,runFolder);

nMeanImages = 9;
[nRows,nCols] = optimizeSubplotDim(nMeanImages+1);

randStarterFrames = sort(randsample(nfr-100,nMeanImages));
setFigParams4Print('landscape')
figure
suptitle([SubNum '-' expDate])
for iimg = 1:nMeanImages
    subplot(nRows,nCols,iimg)
    imagesc(mean(data(:,:,randStarterFrames(iimg):randStarterFrames(iimg)+100),3))
    title(sprintf('%s:%s',num2str(randStarterFrames(iimg)),num2str(randStarterFrames(iimg)+100)))
end

print(['Z:\Analysis\FSAV Summaries\' ds '\rand image samples_' SubNum '-' expDate],'-dpdf','-fillpage')
savefig(['Z:\Analysis\FSAV Summaries\' ds '\rand image samples_' SubNum '-' expDate])

stFrame = expt(iexp).regImgStartFrame;
retImg = mean(data(:,:,stFrame:stFrame+99),3);
retDataDownSampled = stackGroupProject(data,downSampleRate);
[~,retDataRegistered] = stackRegister(retDataDownSampled,retImg);

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
if ~exist(fullfile(fnout,runFolder),'dir')
    mkdir(fullfile(fnout,runFolder))
end
positionTitles = arrayfun(@(x,y) sprintf('Az: %s, El: %s',num2str(x),...
       num2str(y)),azimuthAtPosInd,elevationAtPosInd,'unif',0);
azimuthTitles = arrayfun(@(x) sprintf('Az: %s',num2str(x)),azimuths,'unif',0);
elevationTitles = arrayfun(@(x) sprintf('El: %s',num2str(x)),elevations,'unif',0);

respPositionFig = plotImageStackAsSubplots(positionRespDFFNorm,positionTitles);
print(fullfile(fnout,runFolder,'positionMeanDFF'),'-dpdf','-fillpage')
respAzimuthFig = plotImageStackAsSubplots(azimuthRespDFFNorm,azimuthTitles);
print(fullfile(fnout,runFolder,'azimuthMeanDFF'),'-dpdf','-fillpage')
respElevationFig = plotImageStackAsSubplots(elevationRespDFFNorm,elevationTitles);
print(fullfile(fnout,runFolder,'elevationMeanDFF'),'-dpdf','-fillpage')

plotRGBImage(azimuthRespDFFNorm);
title(sprintf('azimuth: %s (RGB)',num2str(azimuths)))
print(fullfile(fnout,runFolder,'azimuthRGB'),'-dpdf','-fillpage')
plotRGBImage(elevationRespDFFNorm)
title(sprintf('elevation: %s (RGB)',num2str(elevations)))
print(fullfile(fnout,runFolder,'elevationRGB'),'-dpdf','-fillpage')