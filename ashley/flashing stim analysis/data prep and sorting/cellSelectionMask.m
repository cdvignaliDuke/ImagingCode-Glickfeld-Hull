clear all
close all
rc = behavConstsAV;
ds = 'FSAV_V1_GAD';
eval(ds)
doFOVsegment = 1;
for iexp = 1:size(expt,2)
SubNum = expt(iexp).SubNum;
mouse = expt(iexp).mouse;
expDate = expt(iexp).date;

fnout = fullfile(rc.ashleyAnalysis,mouse,expt(iexp).folder,expDate);

%% load all max dF/F images, crop if necessary

load(fullfile(fnout,'data processing','max_images_crop.mat'))    
 
dFF_stack = cat(3,dir_crop,bx_crop);
%% ****select cells
nstim = size(dir_crop,3);

if doFOVsegment
    if expt(iexp).areaBorders    
        borders = readtiff(fullfile(fnout,'FOVborders.tif'));
        dFFstackWithBorders = zeros(size(dFF_stack));
        for iimage = 1:size(dFF_stack,3)
            thisImage = dFF_stack(:,:,iimage);
            imageWithBorders = mean(cat(3,thisImage,borders),3);
            dFFstackWithBorders(:,:,iimage) = imageWithBorders;
        end
        mask_cell = maskFromMultiMaxDFFStack(dFFstackWithBorders);
        figure; setFigParams4Print('portrait')
        imagesc(mask_cell);
        title({[num2str(length(unique(mask_cell(:)))-1) ' cells with behavior'];[mouse '-' expDate]})
        print(fullfile(fnout,['final_mask' ds]),'-dpdf')
        save(fullfile(fnout,['final_mask' ds '.mat']),'mask_cell');
        close all
    else
        mask_cell = maskFromMultiMaxDFFStack(dFF_stack);
        figure; setFigParams4Print('portrait')
        imagesc(mask_cell);
        title({[num2str(length(unique(mask_cell(:)))-1) ' cells with behavior'];[mouse '-' expDate]})
        print(fullfile(fnout,'final_mask'),'-dpdf')

        save(fullfile(fnout,'final_mask.mat'),'mask_cell');
        close all
    end
else
    mask_cell = maskFromMultiMaxDFFStack(dFF_stack);
    figure; setFigParams4Print('portrait')
    imagesc(mask_cell);
    title({[num2str(length(unique(mask_cell(:)))-1) ' cells with behavior'];[mouse '-' expDate]})
    print(fullfile(fnout,'final_mask'),'-dpdf')

    save(fullfile(fnout,'final_mask.mat'),'mask_cell');
    close all
end

end