clear all
close all
rc = behavConstsAV;
awFSAVdatasets_audControl
for iexp = 4:size(expt,2)

SubNum = expt(iexp).SubNum;
mouse = expt(iexp).mouse;
expDate = expt(iexp).date;

fnout = fullfile(rc.ashleyAnalysis,mouse,expt(iexp).folder,expDate);

%% load all max dF/F images, crop if necessary

load(fullfile(fnout,'max_images_crop.mat'))    
 
dFF_stack = cat(3,dir_crop,bx_crop);
%% ****select cells
nstim = size(dir_crop,3);

mask_cell = maskFromMultiMaxDFFStack(dFF_stack);

figure; setFigParams4Print('portrait')
imagesc(mask_cell);
title({[num2str(length(unique(mask_cell(:)))-1) ' cells with behavior'];[mouse '-' expDate]})
print(fullfile(fnout,'final_mask'),'-dpdf')

save(fullfile(fnout,'final_mask.mat'),'mask_cell');
close all
end