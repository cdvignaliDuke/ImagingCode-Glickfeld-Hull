clear all
close all
rc = behavConstsAV;
awFSAVdatasets_V1

%%
iexp = 20;
    %%
%% expt specs
SubNum = expt(iexp).SubNum;
mouse = expt(iexp).mouse;
expDate = expt(iexp).date;
dirFolder = expt(iexp).dirtuning;
% % % down = 10;

fnbx = fullfile(rc.ashleyAnalysis,mouse,'two-photon imaging',expDate);
fntun = fullfile(rc.ashleyAnalysis,mouse,'two-photon imaging',expDate,dirFolder);

%% load max images
load(fullfile(fnbx,'bx_max_images.mat'))
load(fullfile(fntun,'tun_max_images.mat'))

tun_img = max(dFF_dirmax,[],3);
bx_img = max(cat(3, start_max, long_max, tar_max),[],3);

%% image specs
xpix = size(tun_img,2);
ypix = size(tun_img,1);

%% ******choose crop parameters*******

% tuning image
figure;colormap gray; imagesc(tun_img)

%**enter vals here***
xcrop = [1:2 794:xpix];
ycrop = [1:2 262:ypix];

tun_crop = tun_img;
tun_crop(:,xcrop) = 0;
tun_crop(ycrop,:) = 0;

imagesc(tun_crop)

% check if behavior image still needs cropping
% tuning image
figure;colormap gray; imagesc(bx_img)

bx_crop = bx_img;
bx_crop(:,xcrop) = 0;
bx_crop(ycrop,:) = 0;

imagesc(bx_crop)

% % % if needs more crop, change vals here
% % xcrop = [1:2 794:xpix];
% % ycrop = [1:14 254:ypix];
% % 
% % bx_crop = bx_img;
% % bx_crop(:,xcrop) = 0;
% % bx_crop(ycrop,:) = 0;
% % 
% % imagesc(bx_crop)

%% crop orig max images
dir_crop = dFF_dirmax;
dir_crop(:,xcrop,:) = 0;
dir_crop(ycrop,:,:) = 0;

bx_crop = cat(3, start_max, long_max, tar_max);
bx_crop(:,xcrop,:) = 0;
bx_crop(ycrop,:,:) = 0;

%% save cropped images
save(fullfile(fnbx,'max_images_crop.mat'),'dir_crop', 'bx_crop','xcrop','ycrop');

