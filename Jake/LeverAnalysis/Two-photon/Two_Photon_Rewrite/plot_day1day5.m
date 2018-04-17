load([out_dir2, 'ROI_TCs.mat']);
day5_mask = reshape(mask_final, npw, nph);
day5_warp_mask = imwarp(day5_mask, D);
mask_temp = reshape(day5_warp_mask, npw, nph);
[ ~, mask3D, ~] = finalMask(img_pca, mask_temp, threshold, out_dir_comb);
mask_final = processMask(mask3D);
day5_mask = reshape(mask_final, npw,nph);
load([out_dir, 'ROI_TCs.mat']);
day1_mask = reshape(mask_final, npw,nph);

aa = mask1;
aa(aa>0) = 1;

bb = mask2;
bb(bb>0) = 0.8;
figure;imshow(mat2gray(aa));truesize
green = cat(3, zeros(size(aa)), ones(size(aa)), zeros(size(aa)));
hold on
h = imshow(green);
hold off
set(h, 'AlphaData', bb)

