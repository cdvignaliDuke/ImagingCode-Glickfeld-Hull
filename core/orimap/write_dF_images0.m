function [max_dir_dF, max_dir_ratio_change] = write_dF_images0 (dir_dF_sm, dir_ratio, stim_names);

% write dF images & ratio change images as tif files
% auto scaled between zero and max value
dim=size(dir_dF_sm);

ndir=dim(3);

if nargin <3
    for i=1:ndir
        stim_names{i}=num2str(i);
    end
end

max_dir_dF=max(max(max(dir_dF_sm)));
max_dir_ratio_change=max(max(max(dir_ratio)));


for i=1:ndir
    imwrite(dir_dF_sm(:,:,i)./max_dir_dF, ['dF_', stim_names{i}, '.tif']);
    imwrite(dir_ratio(:,:,i)./max_dir_ratio_change, ['ratio_', stim_names{i}, '.tif']);
end
