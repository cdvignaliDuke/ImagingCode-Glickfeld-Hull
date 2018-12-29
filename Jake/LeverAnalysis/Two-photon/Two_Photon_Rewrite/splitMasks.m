function [maskCat maskCat_map] = splitMasks(split_img, mask, nmask)

maskCat = zeros(nmask,1);
maskCat_map = zeros(size(mask(:,:,1)));
for ic=1:nmask
    [i,j] = find(mask(:,:,ic)==1);                             
    tc = zeros(size(i,1),1);
    for a = 1:size(i,1)
        tc(a,:) = split_img(i(a),j(a));
    end
    if mean(tc,1) == 0
        maskCat(ic) = 1;
        maskCat_map(find(mask(:,:,ic)==1)) = 1;
    elseif mean(tc,1) == 1
        maskCat(ic) = 2;
        maskCat_map(find(mask(:,:,ic)==1)) = 2;
    else
        maskCat(ic) = 0;
        maskCat_map(find(mask(:,:,ic)==1)) = 3;
    end
end
end