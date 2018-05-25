clear all
close all
rc = behavConstsAV;
FSAV_V1_GAD
iexp = 8;

% expt specs
SubNum = expt(iexp).SubNum;
mouse = expt(iexp).mouse;
expDate = expt(iexp).date;
dirFolder = expt(iexp).dirtuning;
% % % down = 10;

fnbx = fullfile(rc.ashleyAnalysis,mouse,'two-photon imaging',expDate,'data processing');

fntun = fullfile(rc.ashleyAnalysis,mouse,'two-photon imaging',expDate,'data processing');

% load max images
load(fullfile(fnbx,'bx_max_images.mat'))
load(fullfile(fntun,'tun_max_images.mat'))

tun_img = max(dFF_dirmax,[],3);
try
    bx_img = max(cat(3, start_max, long_max, tar_max),[],3);
catch
    bx_img = max(dFF_bxMax,[],3);
end

% image specs
xpix = size(tun_img,2);
ypix = size(tun_img,1);

% tuning image
figure;colormap gray; imagesc(tun_img)
%% ******choose crop parameters*******


%**enter vals here***
xcrop = [1:5 794:xpix];
ycrop = [1:12 260:ypix];
%%
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
dir_crop = dFF_dirmax;
dir_crop(:,xcrop,:) = 0;
dir_crop(ycrop,:,:) = 0;

try
    bx_crop = cat(3, start_max, long_max, tar_max);
catch
    bx_crop = dFF_bxMax;
end
bx_crop(:,xcrop,:) = 0;
bx_crop(ycrop,:,:) = 0;

% save cropped images
save(fullfile(fnbx,'max_images_crop.mat'),'dir_crop', 'bx_crop','xcrop','ycrop');
