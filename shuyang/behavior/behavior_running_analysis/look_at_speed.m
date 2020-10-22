% look at the running trials, what does the speeds look like, what does the
% speed of 1.5 seconds before and after running look like?

clear;
days = '1064-200321_1';
behav_dest = ['Z:\Analysis\WF_MovingDots_Analysis\behavioral_analysis\' days '\'];
behav_output = load([behav_dest days '_behavAnalysis.mat']);
frames_run_cell = behav_output.frames_run_cell;
speed = double(behav_output.speed);
% get the last and first frame number of each trial
trial_start_end = zeros(length(frames_run_cell),2);
for i = 1:length(frames_run_cell)
    trial_start_end(i,1) = frames_run_cell{i}(1);
    trial_start_end(i,2) = frames_run_cell{i}(end);
end

% put speeds of these running trials and speeds before and after running in
% a matrix to have an idea of how does it look like
trial_len = zeros(1,length(frames_run_cell));
for i = 1:length(frames_run_cell)
    trial_len(i) = length(frames_run_cell{i});
end
max_len = max(trial_len);

speeds_run_mat = NaN(length(frames_run_cell),max_len);
for i = 2:size(speeds_run_mat,1)
    speeds_run_mat(i,1:trial_len(i)) = speed(frames_run_cell{i});
end
speeds_befo_run_mat = zeros(length(frames_run_cell),30);
for i = 2:size(speeds_befo_run_mat,1)
    speeds_befo_run_mat(i,:) = speed(frames_run_cell{i}(1)-30:frames_run_cell{i}(1)-1);
end
speeds_aft_run_mat = zeros(length(frames_run_cell),45);
for i = 2:size(speeds_aft_run_mat,1)
    speeds_aft_run_mat(i,:) = speed(frames_run_cell{i}(end):frames_run_cell{i}(end)+44);
end
speed_run_mat = [speeds_befo_run_mat,speeds_run_mat,speeds_aft_run_mat];


% # of frames between running trials, when shorter than 1s (30 frames)
%from all running trials in session 190429-1021, 190430-1023, 190507-1024,
%190603-1025
frm_bt = [7,9,7,7,10,7,10,7,8,8,7,7,7,7,10,7,10,10,8,9,7,7,10,9,11,8,10,8,12,...
    20,9,7,7,9,8,7,10,8,11,10,7,7,23,7,8,7,7,9,12,9,13,9,11,7,7,13,12,17,11,...
    8,19,11,8,14,7,7,7,7,8,10,12,9,7,7,10,7,8,21,7,17,7,9,7,7,18];
figure; hist(frm_bt);title('number of frames between consecutive running trials');
xlabel('frames'); ylabel('#of trials');




