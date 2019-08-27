%function for determining movement onsets
% to be used in piezo_analysis

function [all_mov_onset_inds, mov_inds, bout_buffer] = piezo_mov_onset(piezo_data, mov_means, mov_stds);

%define threshold for deciding if a movement bout is its own bout or part of a bout which came before it
bout_buffer = ceil(600/(1000/90)); %ASSUMES a 90hz piezo collection rate

%find all indeces where movement trace is above threshold
mov_thresh_pos = mov_means+mov_stds;
mov_thresh_neg = mov_means-mov_stds;
mov_inds = find(piezo_data > mov_thresh_pos | piezo_data<mov_thresh_neg);

%indicate the onset index of each bout by finding suprathreshold values
%with only subthreshold values in hte 600ms before that.
mov_inds_diff = diff(mov_inds); 
all_mov_onset_inds = [0, find(mov_inds_diff > bout_buffer)]+1; %the first suprathreshold movement will always be the onset. Add one to find the correct index due to diff. 

return





