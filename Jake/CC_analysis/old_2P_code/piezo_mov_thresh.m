%function for determining threshold for movement onset from 2P piezo data
% to be used in piezo_analysis

function [mov_means, mov_stds] = piezo_mov_thresh(piezo_data, window_size);

%find the moving mean and std
mov_means = movmean(piezo_data, window_size);
mov_stds = movstd(piezo_data, window_size);

%window size is truncated at the ends fix this by just using the mean/std
%over the first window_size datapoints for the first window_size/2 datapoints
half_win = ceil(window_size/2);
mov_means(1:half_win) = mean(piezo_data(1:window_size));
mov_means(end-half_win:end) = mean(piezo_data(end-window_size:end));
mov_stds(1:half_win) = std(piezo_data(1:window_size));
mov_stds(end-half_win:end) = std(piezo_data(end-window_size:end));

return



