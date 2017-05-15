clear all
close all
rc = behavConstsAV;
awFSAVdatasets_V1
for iexp = 16:size(expt,2)

SubNum = expt(iexp).SubNum;
mouse = expt(iexp).mouse;
expDate = expt(iexp).date;

fnout = fullfile(rc.ashleyAnalysis,mouse,expt(iexp).folder,expDate);

load(fullfile(fnout,'final_mask.mat'));

%% load all max dF/F images, crop if necessary

load(fullfile(fnout,'max_images_crop.mat'))    
 
dFF_stack = cat(3,dir_crop,bx_crop);

%% subtract off cell mask 
xpix = size(dFF_stack,2);
ypix = size(dFF_stack,1);

dFF_2D = reshape(dFF_stack,xpix*ypix,[]);
dFF_2D(mask_cell(:) > 0,:) = 0;

dFF_stack_den = reshape(dFF_2D,ypix,xpix,[]);

%% ****select cells
nstim = size(dir_crop,3);

mask_cell = maskFromMultiMaxDFFStack(dFF_stack_den);

figure; setFigParams4Print('portrait')
imagesc(mask_cell);
title({[num2str(length(unique(mask_cell(:)))-1) ' dendrites with behavior'];[mouse '-' expDate]})
print(fullfile(fnout,'dendrite_mask'),'-dpdf')

save(fullfile(fnout,'dendrite_mask.mat'),'mask_cell');

end