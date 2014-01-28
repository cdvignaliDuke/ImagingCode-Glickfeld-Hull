function [dist,ctrs] = mask2dist(mask,fov_siz)

props = regionprops(mask,'Centroid');

nCells= max(mask(:));

for i1 = 1:nCells
    for i2 = 1:nCells
        dist_pix(i1,i2) = sqrt(sum((props(i1).Centroid - props(i2).Centroid).^2));
        dist(i1,i2) = dist_pix(i1,i2)/256*fov_siz;
    end;
end;

ctrs = reshape([props(:).Centroid],2,nCells)'/256*fov_siz;

return;
