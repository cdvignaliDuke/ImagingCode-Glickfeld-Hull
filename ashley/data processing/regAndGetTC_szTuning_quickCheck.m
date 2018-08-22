clear all
close all
ds = 'szTuning_dreadds_PM';
rc = behavConstsAV;
eval(ds)
doRedChannel = 1;
doOtherRedChannel = 0;
iexp = 1;
%%
down = 10;

nBaselineFr = round((params.nBaselineMs/1000)*params.frameRate);
nRespwinFr_sz = round(0.5*params.frameRate);

respwin_sz = (1+nBaselineFr+params.nFramesVisDelay_VSR):...
    (nBaselineFr+params.nFramesVisDelay_VSR+nRespwinFr_sz);

basewin_sz = (nBaselineFr - nRespwinFr_sz +1):nBaselineFr;
%%

mouse = expt(iexp).mouse;
subnum = mouse;
expDate = expt(iexp).date;

fn = fullfile(rc.ashleyAnalysis,mouse,'two-photon imaging',expDate);
fnout = fullfile(fn,'data processing');

%% load and register one dataset

imgFolder = expt(iexp).sizeTuningFolder{1};

    fName = [imgFolder '_000_000'];
    data_temp = loadsbx_choosepmt(1,mouse,expDate,imgFolder,fName);

        regImgFrames = ...
            expt(iexp).regImgStartFrame:(expt(iexp).regImgStartFrame+99);
        regImg = mean(data_temp(:,:,regImgFrames),3);

    [outs, data_reg] = stackRegister(data_temp,regImg);

    if ~exist(fullfile(fn,imgFolder{irun}),'dir')
        mkdir(fullfile(fn,imgFolder{irun}))
    end
    save(fullfile(fn,imgFolder,'regOuts&Img'),'regImg','outs')



%%
mw = loadMworksFile(subnum,expDate,expt(iexp).sizeTuningTime{1});
tSize = round(celleqel2mat_padded(mw.tGratingDiameterDeg),2,'significant');

nOn = mw.nScansOn;
nOff = mw.nScansOff;
nFr = size(data_reg,3);

trialStart = nOn:(nOn+nOff):nFr;
trialType = findgroups(tSize);
maxDFF = getStimMaxDFF(data_reg,...
    trialStart,trialType,basewin_sz,respwin_sz);

[ypix,xpix,~] = size(maxDFF);
%crop it
tun_img = max(maxDFF,[],3);
% tuning image
figure;colormap gray; imagesc(tun_img)

%**enter vals here***
xcrop = [1:2 794:xpix];
ycrop = [1:10 510:ypix];

tun_crop = tun_img;
tun_crop(:,xcrop) = 0;
tun_crop(ycrop,:) = 0;

imagesc(tun_crop)

% crop orig max images
maxDFF(:,xcrop,:) = 0;
maxDFF(ycrop,:,:) = 0;


%%
green_mask = maskFromMultiMaxDFFStack(maxDFF);
close all
imagesc(green_mask)
title(sprintf('%s cells selected',num2str(length(unique(green_mask(:)))-1)))
%%
tc = stackGetTimeCourses(data_reg,green_mask);
dataTC = getWeightedNeuropilTimeCourse(...
    data_reg,tc,green_mask,4,6);

%%
nt = length(tSize);
sizes = unique(tSize);
nSize = length(sizes);

nCells = size(dataTC,2);
nStimFr = nRespwinFr_sz+nBaselineFr;
F = nan(nBaselineFr+nStimFr,nt,nCells);
for i = 1:nt
   ind = trialStart(i);
    F(:,i,:) = dataTC((ind-nBaselineFr+1):(ind+nStimFr),:);
end
F0 = mean(F(basewin_sz,:,:),1);
dFF = (F-F0)./F0;

sizeTC = nan(nBaselineFr+nStimFr,nSize,nCells);
for i = 1:nSize
    ind = tSize == sizes(i);
    sizeTC(:,i,:) = mean(dFF(:,ind,:),2);
end

%% do global max dF/F if no signal
F0 = uint16(mean(data_reg,3));
dF = data_reg-F0;
clear F0
maxDF = max(dF,[],3);
imagesc(maxDF)
%% check retinotopy if no signal
imgFolder = expt(iexp).retinotopyFolder{1};

fName = [imgFolder '_000_000'];
data_temp = loadsbx_choosepmt(1,mouse,expDate,imgFolder,fName);

[outs, data_reg] = stackRegister(data_temp,regImg);
clear data_temp
mw = loadMworksFile(subnum,expDate,expt(iexp).retinotopyTime{1});

[ypix,xpix,nFr] = size(data_reg);
nOn = mw.nScansOn;
nOff = mw.nScansOff;
basewin = round(nOff/2):nOff;
respwin = (nOff+2):(nOff+10);
tAz = cell2mat_padded(mw.tGratingAzimuthDeg);
tEl = cell2mat_padded(mw.tGratingElevationDeg);
azimuths = unique(tAz);
elevations = unique(tEl);
nAz = length(azimuths);
nEl = length(elevations);
nStim = nAz.*nEl;
nt = length(tAz);

data_tr = double(reshape(data_reg,ypix,xpix,nOn+nOff,nt));
F0 = mean(data_tr(:,:,basewin,:),3);
dFF = (data_tr - F0)./F0;

clear data_reg data_tr F0

data_stims = nan(ypix,xpix,nStim);
stims = cell(1,nStim);
for iEl = 1:nEl
    if iEl == 1
        stimInd = 1;
    end
    indE = tEl == elevations(iEl);
    for iAz = 1:nAz
        indA = tAz == azimuths(iAz);
        ind = indE & indA;
        
        data_stims(:,:,stimInd) = mean(mean(dFF(:,:,respwin,ind),4),3);
        stims{stimInd} = [elevations(iEl), azimuths(iAz)];
        stimInd = stimInd+1;
    end
end

[nRow,nCol] = optimizeSubplotDim(nStim);

figure
colormap gray
for i = 1:nStim
    subplot(nRow,nCol,i)
    imagesc(data_stims(:,:,i))
    title(stims{i})
end

ret_mask = maskFromMultiMaxDFFStack(data_stims);

tc = stackGetTimeCourses(data_reg,green_mask);
dataTC = getWeightedNeuropilTimeCourse(...
    data_reg,tc,green_mask,4,6);