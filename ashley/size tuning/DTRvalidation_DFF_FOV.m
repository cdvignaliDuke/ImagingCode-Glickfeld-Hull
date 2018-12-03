clear all
close all
ds = 'DTRvalidation_EMX';
rc = behavConstsAV;
eval(ds)
slct_expt = 2;
%% load imaging data

iexp = slct_expt

mouse = expt(iexp).mouse;
expDate = expt(iexp).date;
dirFolder = expt(iexp).dirTuningFolder{1};
dirTime = expt(iexp).dirTuningTime{1};

if strcmp(rc.name, 'ashle')
    fn = fullfile(rc.ashleyAnalysis,mouse,'two-photon imaging',expDate);
    if ~exist(fullfile(fn,'data processing'),'dir')
        mkdir(fn,'data processing')
    end
else
    error('you do not belong here')
end

fnout = fullfile(fn,'data processing');

data = loadsbx_choosepmt(1,mouse,expDate,dirFolder,[dirFolder '_000_000'],[],...
    fullfile(rc.ashleyData,mouse,'two-photon imaging',expDate,dirFolder));
mw = loadMworksFile(mouse,expDate,dirTime);
%% down-sample and choose registration image
data_down = stackGroupProject(data,params.downSampleRate);
[ypix,xpix,nfr] = size(data_down);

nMeanImages = 9;
[nRows,nCols] = optimizeSubplotDim(nMeanImages+1);

randStarterFrames = sort(randsample(nfr-10,nMeanImages));
setFigParams4Print('landscape')
figure
suptitle([mouse '-' expDate])
for iimg = 1:nMeanImages
    subplot(nRows,nCols,iimg)
    imagesc(mean(data_down(:,:,randStarterFrames(iimg):randStarterFrames(iimg)+9),3))
    title(sprintf('%s:%s',num2str(randStarterFrames(iimg)),num2str(randStarterFrames(iimg)+100)))
end

print(fullfile(fnout,...
    ['rand image samples_' mouse '-' expDate]),'-dpdf','-fillpage')
savefig(fullfile(fnout,...
    ['rand image samples_' mouse '-' expDate]))

%% register
regImg = mean(data_down(:,:,expt(iexp).regImgStartFrame:(expt(iexp).regImgStartFrame+9)),3);
figure;imagesc(regImg);colormap gray
[out_down,data_down_reg] = stackRegister(data_down,regImg);
save(fullfile(fnout,'regOuts&Img.mat'),'out_down', 'regImg')
figure;imagesc(mean(data_down_reg,3));colormap gray
writetiff(data_down_reg,fullfile(fnout,'downSampleData'));
%% align images to stim on
nBaselineFr = ceil((params.nBaselineMs./1000).*params.frameRate./params.downSampleRate);
nStimFr = ceil((params.nPostStimMs./1000).*params.frameRate./params.downSampleRate);
nTrialFr = nBaselineFr+nStimFr;
basewin = 1:nBaselineFr;
respwin = (nBaselineFr+params.nFramesVisDelay_VSR+1):(nBaselineFr+params.nFramesVisDelay_VSR+2);

on = mw.nScansOn./params.downSampleRate;
off = mw.nScansOff./params.downSampleRate;

tDirection = cell2mat(mw.tGratingDirectionDeg);
directions = unique(tDirection);
nDirection = length(directions);
nTrials = length(tDirection);
tStimOn = off+1:(on+off):nfr;
tStimOn = tStimOn(1:nTrials);

F = mean(data_down_reg,3);
dF= data_down_reg-F;
dFF = dF./F;

maxDFF = max(dFF,[],3);
figure;imagesc(maxDFF);colormap gray

dFF_trialTC = nan(ypix,xpix,nTrialFr,nTrials);
for i = 1:nTrials
    ind = (tStimOn(i)-nBaselineFr):(tStimOn(i)+nStimFr-1);
    dFF_trialTC(:,:,:,i) = data_down_reg(:,:,ind);
end

img_trialResp = squeeze(mean(dFF_trialTC(:,:,respwin,:),3));
maxTrialResp = max(img_trialResp,[],3);
writetiff(maxTrialResp,fullfile(fnout,'maxRespAcrossTrials'));
