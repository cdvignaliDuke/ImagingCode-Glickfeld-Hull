clear all
close all
ds = 'awData_audMod_V13trialtypes';
%%
rc = behavConstsAV;
eval(ds)
%%
iexp = 2;
    %%
%% expt specs
SubNum = expt(iexp).SubNum;
mouse = expt(iexp).mouse;
expDate = expt(iexp).date;

fnout = fullfile(rc.ashleyAnalysis,mouse,'two-photon imaging',expDate);

%% load max dF/F for each run in expt

for irun = 1:expt(iexp).nrun
    runFolder = expt(iexp).runs(irun,:);
   temp_max = readtiff(fullfile(fnout, runFolder,'max images.tif'));
   if irun == 1
       maxDFFstack = zeros(size(temp_max,1), size(temp_max,2), expt(iexp).nrun);
   end
       maxDFFstack(:,:,irun) = temp_max;   
end 
max_img = max(maxDFFstack,[],3);

figure;
for i = 1:size(maxDFFstack,3)
    imagesc(maxDFFstack(:,:,i))
    drawnow
end
%% image specs
xpix = size(max_img,2);
ypix = size(max_img,1);

%% ******choose crop parameters*******

% tuning image
figure;colormap gray; imagesc(max_img)

%**enter vals here***
xcrop = [1:2 793:xpix];
ycrop = [1:2 649:ypix];

max_crop = max_img;
max_crop(:,xcrop) = 0;
max_crop(ycrop,:) = 0;

imagesc(max_crop)

%% crop orig max images
stim_crop = maxDFFstack;
stim_crop(:,xcrop,:) = 0;
stim_crop(ycrop,:,:) = 0;

%% save cropped images
save(fullfile(fnout,'max_images_crop.mat'),'stim_crop', 'xcrop','ycrop');

