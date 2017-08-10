[x,y,z] = size(img);
sm_logical = zeros(x,y);
sm = zeros(x,y,z);
Ithreshold = 95;
for i = 1: size(img,3)
    sm_logical (img(:,:,i) > 0.6*mean([max(prctile(img(:,:,i),Ithreshold,1)) max(prctile(img(:,:,i),Ithreshold,2))])) = 1;
    sm(:,:,i) = sm_logical;
end
sm_logical = zeros(x,y);
sm_logical(sum(sm,3) > round(0.9*z)) = 1;
bwmask = bwlabel(sm_logical);
bwRP =  regionprops(bwmask, 'Area'); minRoiPixArea = 200;
regionIdx = find([bwRP.Area] > minRoiPixArea);

bkMask = 0*sm_logical;
for i = 1:length(regionIdx)
    bkMask(bwmask == regionIdx(i)) = 1;
end
bkMask = stackFilter(bkMask);
bkMask(bkMask > 0) = 1;
figure;imagesc(bkMask);

% apply mask
imgMasked = bsxfun(@times, img, cast(bkMask, 'like', img));
figure;imagesc(imgMasked(:,:,1));