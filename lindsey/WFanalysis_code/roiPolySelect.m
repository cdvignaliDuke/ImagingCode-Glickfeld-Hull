bdata_avg = nanmean(mean(bdata_trials_gsub(:,:,preframes+5:preframes+15,:),3),4);
figure; H = imagesc(bdata_avg); colormap(gray);
clim([-.005 .03])
%adjust number of ROIs to choose
for i = 1:nROI
    roi(i) = impoly;
end

mask = zeros(sz(1),sz(2),nROI);
for i = 1:nROI
    mask(:,:,i) = createMask(roi(i));
end

roi_cluster = sum(mask,3);
mask_cell = bwlabel(roi_cluster);
figure; imagesc(mask_cell);