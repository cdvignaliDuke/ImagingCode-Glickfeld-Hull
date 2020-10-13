ds = 'CrossOriRandDir_ExptList';
eval(ds)
nexp = length(expt);
LG_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\lindsey';

for iexp = 1:nexp
    if size(expt(iexp).img_loc,2)>1
        if ~exist(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_splitImage.mat']))
            mouse = expt(iexp).mouse;
            date = expt(iexp).date;
            fprintf([mouse ' ' date '\n'])
            ImgFolder = expt(iexp).coFolder;
            nrun = length(ImgFolder);
            run_str = catRunName(cell2mat(ImgFolder), nrun);
            load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_mask_cell.mat']))

            img = max(data_dfof,[],3);
            sz = size(img);

            figure; imagesc(img)
            [x y] = ginput;

            xvals = 1:sz(2);
            yvals = 1:sz(1);

            y(1) = 1;
            y(end) = sz(1);
            x_int = interp1(y,x,yvals);

            split_img = zeros(size(img));
            for i = yvals
                split_img(i,round(x_int(i)):sz(2)) = 1;
            end
            figure; imagesc(split_img)

            nmask = max(mask_cell(:));
            mask3D = zeros(sz(1),sz(2),nmask);
            for i = 1:nmask
                ind = find(mask_cell == i);
                mask_temp = zeros(sz(1),sz(2));
                mask_temp(ind) = 1;
                mask3D(:,:,i) = mask_temp;
            end
            [maskCat maskCat_map] = splitMasks(split_img, mask3D, nmask);
            figure; imagesc(maskCat_map)

            save(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_splitImage.mat']), 'split_img', 'x', 'y','maskCat','maskCat_map');
        end
    end
end