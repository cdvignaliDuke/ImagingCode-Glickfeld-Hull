clear all
close all
CRP_expt_list_Crus
for id = 3
lg_out = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lindsey\Analysis\2P\Jake';
nexp = size(expt(id).date,1);
fprintf(['Day ' num2str(id) '\n'])
    for iexp = 1:nexp
        mouse = expt(id).mouse(iexp,:);
        date = expt(id).date(iexp,:);
        run = expt(id).run(iexp,:);
        fprintf([date ' ' mouse '\n'])
        img_fn = [date '_' mouse];

        load(fullfile(lg_out,img_fn, [img_fn '_reg.mat']))
        sz = size(img_ref);
        figure; 
        imagesc(img_ref)
        [x y] = ginput;
        
        
        xvals = 1:sz(2);
        yvals = 1:sz(1);
        
        y(1) = 1;
        y(end) = sz(1);
        x_int = interp1(y,x,yvals);
        
        split_img = zeros(size(img_ref));
        for i = yvals
            split_img(i,round(x_int(i)):sz(2)) = 1;
        end
        figure; imagesc(split_img)
        
        load(fullfile(lg_out,img_fn, [img_fn '_ROI_TCs.mat']))
        nmask = size(mask3D,3);
        [maskCat maskCat_map] = splitMasks(split_img, mask3D, nmask);
        figure; imagesc(maskCat_map)

        save(fullfile(lg_out,img_fn, [img_fn '_splitImage.mat']), 'split_img', 'x', 'y','maskCat','maskCat_map');
    end
end