function mask_final = getCells(img1, Ithreshold)
I = img1(:,:,1:4:end);
% figure;imshow(mat2gray(I));
[x,y,nt] = size(I);

% Ithreshold = 96;
sm_logical = zeros(x,y);
% cellMask = zeros(size(I));
cellMask = [];
se = strel('disk',2);
minRoiPixArea = 50; maxRoiPixArea = 500; eRatio = 1.5;
for i = 1:nt
    i
    sm_logical (I(:,:,i) > mean([max(prctile(I(:,:,i),Ithreshold,1)) max(prctile(I(:,:,i),Ithreshold,2))])) = 1;
    ele = bwlabel(sm_logical);
    bw = regionprops(ele, 'MajorAxisLength', 'MinorAxisLength', 'Area', 'Orientation');
     
    regionIdx = find([bw.Area] >= minRoiPixArea & [bw.Area] <= maxRoiPixArea & [bw.MajorAxisLength]./[bw.MinorAxisLength] > eRatio...
        & [bw.Orientation] < 90 & [bw.Orientation] >= 0);
    if ~isempty(regionIdx)
        bkMask = 0*sm_logical;
        for j = 1:length(regionIdx)
            bkMask(ele == regionIdx(j)) = 1;
        end
        masktmp = imerode( stackFilter( imfill(bkMask, 'holes'), 1), se);
        masktmp (masktmp < 0.1) = 0;
        masktmp = bwlabel(masktmp);
        for j = 1 : max(max(masktmp))
            maskSep = zeros(x,y);
            maskSep(masktmp == j) = 1;
            cellMask = cat(3, cellMask, maskSep);
        end
    end
    if mod(i,500) == 0
        cellNewMask = bwlabel(WF_processMask(cellMask));
        cellMask = [];
        for k = 1:max(max(cellNewMask))
            maskSep = zeros(x,y);
            maskSep(cellNewMask == k) = 1;
            cellMask = cat(3, cellMask, maskSep);
        end
    end
end

mask_final = bwlabel(WF_processMask(cellMask));
bw = regionprops(mask_final, 'MajorAxisLength', 'MinorAxisLength', 'Orientation');

regionIdx = find([bw.MajorAxisLength]./[bw.MinorAxisLength] > eRatio...
    & [bw.Orientation] < 90 & [bw.Orientation] >= 0);
if ~isempty(regionIdx)
    bkMask = 0*sm_logical;
    for j = 1:length(regionIdx)
        bkMask(ele == regionIdx(j)) = 1;
    end
end
mask_final = bwlabel(bkMask);
end