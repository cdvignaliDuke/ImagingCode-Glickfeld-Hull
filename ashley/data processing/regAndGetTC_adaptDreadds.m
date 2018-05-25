clear all
close all
ds = 'adaptDREADDS_V1_SOM_control';
rc = behavConstsAV;
eval(ds)
doRedChannel = 1;
doOtherRedChannel = 0;
iexp = 5;
%%
down = 10;

nBaselineFr = (params.nBaselineMs/1000)*params.frameRate;
nRespwinFr_adapt = 0.1*params.frameRate;
nRespwinFr_dir = 0.5*params.frameRate;

respwin_adapt = (1+nBaselineFr+params.nFramesVisDelay_FS):...
    (nBaselineFr+params.nFramesVisDelay_FS+nRespwinFr_adapt);
respwin_dir = (1+nBaselineFr+params.nFramesVisDelay_VSR):...
    (nBaselineFr+params.nFramesVisDelay_FS+nRespwinFr_dir);

basewin_adapt = (nBaselineFr - nRespwinFr_adapt+1):nBaselineFr;
basewin_dir = (nBaselineFr - nRespwinFr_dir +1):nBaselineFr;
%%

mouse = expt(iexp).mouse;
subnum = mouse;
expDate = expt(iexp).date;

fn = fullfile(rc.ashleyAnalysis,mouse,'two-photon imaging',expDate);
fnout = fullfile(fn,'data processing');

%% load and register all  data

imgFolders = cat(1,expt(iexp).adapt(1),expt(iexp).retinotopy(1),...
    expt(iexp).sizeTuning(1),expt(iexp).dirTuning(1));
nrun = length(imgFolders);

regImg = [];
regImgFrames = [];
registeredData = cell(1,nrun);
for irun = 1:nrun
    fName = [imgFolders{irun} '_000_000'];
    data_temp = loadsbx_choosepmt(1,mouse,expDate,imgFolders{irun},fName);
    if irun == 1
        regImgFrames = ...
            expt(iexp).regImgStartFrame:(expt(iexp).regImgStartFrame+99);
        regImg = mean(data_temp(:,:,regImgFrames),3);
    end
    [outs, data_reg] = stackRegister(data_temp,regImg);
    registeredData{irun} = data_reg;
    if ~exist(fullfile(fn,imgFolders{irun}),'dir')
        mkdir(fullfile(fn,imgFolders{irun}))
    end
    save(fullfile(fn,imgFolders{irun},'regOuts&Img'),'regImg','outs')
    clear data_temp data_reg 
end

%%
if doRedChannel
    
    if doOtherRedChannel
        otherExpDate = expt(iexp).matchExptDate;
        otherExpInd = strcmp({expt.date},otherExpDate);
        
        redData = loadsbx_choosepmt(2,mouse,otherExpDate,...
            expt(otherExpInd).labelFolder,...
            [expt(otherExpInd).labelFolder '_000_000']);
        greenData = loadsbx_choosepmt(1,mouse,otherExpDate,...
            expt(otherExpInd).labelFolder,...
            [expt(otherExpInd).labelFolder '_000_000']);
        redData = stackGroupProject(redData,down);
        greenData = stackGroupProject(greenData,down);
        greenImg = mean(greenData(:,:,1:10),3);
        
        
        greenImgNorm = greenImg./max(greenImg(:));
        regImgNorm = regImg./max(regImg(:));
        [input_points, base_points] = ...
            cpselect(greenImgNorm,regImgNorm,'Wait', true); 
        sz_target  = size(regImgNorm);
        mytform    = maketform('affine',input_points(1:3,:), ...
            base_points(1:3,:));
        redDataShift = imtransform(redData,mytform,...
            'XData',[1 sz_target(2)],'YData',[1 sz_target(1)]);
        greenDataShift = imtransform(greenData,mytform,...
            'XData',[1 sz_target(2)],'YData',[1 sz_target(1)]);
        
        outs = stackRegister(greenDataShift,regImg);
        [~,redData_reg] = stackRegister_MA(redDataShift,[],[],outs);
        redImage = mean(redData_reg,3);
        figure
        imagesc(redImage)
        colormap gray
        writetiff(redImage,fullfile(fnout,'redImage'));
        save(fullfile(fnout,'redImage'),'redImage','outs','mytform');
        
        bwout = imCellEditInteractive(redImage);
        red_mask = bwlabel(bwout);
    else    
        redData = loadsbx_choosepmt(2,mouse,expDate,expt(iexp).labelFolder,...
            [expt(iexp).labelFolder '_000_000']);
        greenData = loadsbx_choosepmt(1,mouse,expDate,expt(iexp).labelFolder,...
            [expt(iexp).labelFolder '_000_000']);
        redData = stackGroupProject(redData,down);
        greenData = stackGroupProject(greenData,down);
        outs = stackRegister(greenData,regImg);
        [~,redData_reg] = stackRegister_MA(redData,[],[],outs);
        redImage = mean(redData_reg,3);
        figure
        imagesc(redImage)
        colormap gray
        writetiff(redImage,fullfile(fnout,'redImage'));
        save(fullfile(fnout,'redImage'),'redImage','outs');

        bwout = imCellEditInteractive(redImage);
        red_mask = bwlabel(bwout);
        save(fullfile(fnout,'redmask'),'red_mask')
    end
end

%% 
if doRedChannel
    redTC = cell(1,nrun);
    for i = 1:nrun
        tc = stackGetTimeCourses(registeredData{i},red_mask);
        redTC{i} = getWeightedNeuropilTimeCourse(registeredData{i},tc,red_mask,4,6);
    end
    save(fullfile(fnout,'redTC'),'redTC')
end

%%
mwAdapt = loadMworksFile(subnum,expDate,expt(iexp).adapt{2});
trialStart = cell2mat(mwAdapt.cStimOn);
trialType = ones(1,length(trialStart));
adaptMaxDFF = getStimMaxDFF(registeredData{1},...
    trialStart,trialType,basewin_adapt,respwin_adapt);

mwDir = loadMworksFile(subnum,expDate,expt(iexp).dirTuning{2});
nOn = mwDir.nScansOn;
nOff = mwDir.nScansOff;
trialStart = nOff+1:(nOn+nOff):size(registeredData{4},3);
tDirection = cell2mat(mwDir.tGratingDirectionDeg);
trialType = findgroups(tDirection);
dirMaxDFF = getStimMaxDFF(registeredData{4},...
    trialStart,trialType,basewin_dir,respwin_dir);

[ypix,xpix,~] = size(adaptMaxDFF);
%crop it
tun_img = max(dirMaxDFF,[],3);
bx_img = max(adaptMaxDFF,[],3);
% tuning image
figure;colormap gray; imagesc(tun_img)

%**enter vals here***
xcrop = [1:2 794:xpix];
ycrop = [1:10 258:ypix];

tun_crop = tun_img;
tun_crop(:,xcrop) = 0;
tun_crop(ycrop,:) = 0;

imagesc(tun_crop)

% check if behavior image still needs cropping
figure;colormap gray; imagesc(bx_img)
bx_crop = bx_img;
bx_crop(:,xcrop) = 0;
bx_crop(ycrop,:) = 0;

imagesc(bx_crop)

% crop orig max images
dirMaxDFF(:,xcrop,:) = 0;
dirMaxDFF(ycrop,:,:) = 0;

adaptMaxDFF(:,xcrop,:) = 0;
adaptMaxDFF(ycrop,:,:) = 0;

if doRedChannel
    adaptMaxDFF(red_mask > 0) = 0;
    nStim = size(dirMaxDFF,3);
    for i = 1:nStim
        img = dirMaxDFF(:,:,i);
        img(red_mask > 0) = 0;
        dirMaxDFF(:,:,i) = img;
    end
end

save(fullfile(fnout,'max_images_crop.mat'),'dirMaxDFF', 'adaptMaxDFF','xcrop','ycrop');
close all

%%
dFF_stack = cat(3,dirMaxDFF,adaptMaxDFF);
green_mask = maskFromMultiMaxDFFStack(dFF_stack);
save(fullfile(fnout,'greenmask'),'green_mask')

close all
imagesc(green_mask)
title(sprintf('%s cells selected',num2str(length(unique(green_mask(:)))-1)))
%%
dataTC = cell(1,nrun);
for i = 1:nrun
    tc = stackGetTimeCourses(registeredData{i},green_mask);
    dataTC{i} = getWeightedNeuropilTimeCourse(...
        registeredData{i},tc,green_mask,4,6);
end
save(fullfile(fnout,'dataTC'),'dataTC')