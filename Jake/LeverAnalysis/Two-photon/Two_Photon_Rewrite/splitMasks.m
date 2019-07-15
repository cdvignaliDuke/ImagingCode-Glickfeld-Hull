function [maskCat maskCat_map] = splitMasks(split_img, mask, nmask)

maskCat = zeros(nmask,1);
maskCat_map = zeros(size(mask(:,:,1)));
for ic=1:nmask
    [i,j] = find(mask(:,:,ic)==1); %find coords of each pixel in mask                         
    tc = zeros(size(i,1),1);  %tc lenght = # of pixels in mask
    for a = 1:size(i,1)   %for each pixel...
        tc(a,:) = split_img(i(a),j(a));  %determine which side of hte line each pixel in the mask belongs to. 0 for top left, 1 for bottom right
    end
    %assign each mask to one of three categories. Allow for 5% of error
    if mean(tc,1) < 0.05  
        maskCat(ic) = 1;
        maskCat_map(find(mask(:,:,ic)==1)) = 1;
    elseif mean(tc,1) > 0.95
        maskCat(ic) = 2;
        maskCat_map(find(mask(:,:,ic)==1)) = 2;
    else
        maskCat(ic) = 0;
        maskCat_map(find(mask(:,:,ic)==1)) = 3;
    end
end
end