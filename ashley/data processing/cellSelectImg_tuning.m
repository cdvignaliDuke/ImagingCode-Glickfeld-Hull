clear all
close all
ds = 'movDotsSpeedTun_V13trialtypes';
%%
rc = behavConstsAV;
eval(ds)
%%
iexp = 1;
irun = 1;
    %%
%% expt specs
SubNum = expt(iexp).SubNum;
mouse = expt(iexp).mouse;
expDate = expt(iexp).date;
imgFolder = expt(iexp).regImg;

fnout = fullfile(rc.ashleyAnalysis,mouse,'two-photon imaging',expDate,imgFolder);

%% load max images
load(fullfile(fnout,'tun_max_images.mat'))

tun_img = max(dFF_stimmax,[],3);

%% image specs
xpix = size(tun_img,2);
ypix = size(tun_img,1);

%% ******choose crop parameters*******

% tuning image
figure;colormap gray; imagesc(tun_img)

%**enter vals here***
xcrop = [1:2 794:xpix];
ycrop = [1:8 527:ypix];

tun_crop = tun_img;
tun_crop(:,xcrop) = 0;
tun_crop(ycrop,:) = 0;

imagesc(tun_crop)

%% crop orig max images
stim_crop = dFF_stimmax;
stim_crop(:,xcrop,:) = 0;
stim_crop(ycrop,:,:) = 0;

%% save cropped images
save(fullfile(fnout,'max_images_crop.mat'),'stim_crop', 'xcrop','ycrop');

