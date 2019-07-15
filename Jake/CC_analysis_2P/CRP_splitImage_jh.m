clear all
close all
CRP_expt_list_Crus_jh
for id = 2
jake_out = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\jake\Analysis\Cue_reward_pairing_analysis\CC_analysis_2P';
nexp = size(expt(id).date,1);
fprintf(['Day ' num2str(id) '\n'])
disp('Draw line from top right to bottom left by clicking points in the FoV');
disp('Points will not appear in FoV');
    for iexp = 6%:nexp
        mouse = expt(id).mouse(iexp,:);
        date = expt(id).date(iexp,:);
        run = expt(id).run(iexp,:);
        fprintf([date ' ' mouse '\n'])
        img_fn = [date '_' mouse];

        load(fullfile(jake_out,img_fn, [img_fn '_reg.mat']))
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
        
        load(fullfile(jake_out,img_fn, [img_fn '_ROI_TCs.mat']))
        nmask = size(mask3D,3);
        [maskCat maskCat_map] = splitMasks(split_img, mask3D, nmask);
        figure; imagesc(maskCat_map)

        save(fullfile(jake_out,img_fn, [img_fn '_splitImage.mat']), 'split_img', 'x', 'y','maskCat','maskCat_map');
    end
end