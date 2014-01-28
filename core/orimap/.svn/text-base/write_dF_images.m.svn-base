function [max_dir_dF, max_dir_ratio_change, max_ori_dF, max_ori_ratio_change] = write_dF_images (dir_dF_sm, ori_dF_sm, dir_ratio, ori_ratio);

% write dF images & ratio change images as tif files
% auto scaled between zero and max value

max_dir_dF=max(max(max(dir_dF_sm)));
max_dir_ratio_change=max(max(max(dir_ratio)));
max_ori_dF=max(max(max(ori_dF_sm)));
max_ori_ratio_change=max(max(max(ori_ratio)));

dim=size(dir_dF_sm);
ndir=dim(3);
nori=ndir/2;

for i=1:ndir
    imwrite(dir_dF_sm(:,:,i)./max_dir_dF, ['dir_dF_', num2str(floor(i)), '.tif']);
    imwrite(dir_ratio(:,:,i)./max_dir_ratio_change, ['dir_ratio_', num2str(floor(i)), '.tif']);
end

for i=1:nori
    imwrite(ori_dF_sm(:,:,i)./max_ori_dF, ['ori_dF_', num2str(floor(i)), '.tif']);
    imwrite(ori_ratio(:,:,i)./max_ori_ratio_change, ['ori_ratio_', num2str(floor(i)), '.tif']);
end
