function ROIsize = getROIsize(mask_cell)

n_pix = hist(mask_cell(:),unique(mask_cell(:)));
ROIsize = n_pix(:,2:end);

end