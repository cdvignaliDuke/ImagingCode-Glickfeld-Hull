clear all
close all
ds = 'awData_audMod_V13trialtypes';
slct_exp = 2;
%%
rc = behavConstsAV;
eval(ds)
for iexp = slct_exp
    SubNum = expt(iexp).SubNum;
    mouse = expt(iexp).mouse;
    expDate = expt(iexp).date;
    fnout = fullfile(rc.ashleyAnalysis,mouse,expt(iexp).folder,expDate);

    %% load all max dF/F images, crop if necessary

    load(fullfile(fnout,'max_images_crop.mat'))    

    %% ****select cells

    mask_cell = maskFromMultiMaxDFFStack(stim_crop);

    close all
    figure; setFigParams4Print('portrait')
    imagesc(mask_cell);
    title({[num2str(length(unique(mask_cell(:)))-1) ' cells'];[mouse '-' expDate]})
    print(fullfile(fnout,'final_mask'),'-dpdf')

    save(fullfile(fnout,'final_mask.mat'),'mask_cell');
end