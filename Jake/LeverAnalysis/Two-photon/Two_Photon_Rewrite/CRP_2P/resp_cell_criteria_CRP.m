% New script for determining which neurons to include or exclude
% to be run within main_CRP2

% 1) load data
% 2) extract mean TCs
% 3) exclude trials which excede x stds from the mean
% 4) look for responsive neurons
%    -compare cue vs baseline
%    -compare rew vs pre-rew windows   -two stage t-test  
% 5) exclude ramping neurons
% 6) Make combinations of responsive subtypes for plotting later

%% 1) load data
load([dest '_cue_movies.mat']);  %NR_movie dim1=trials  dim2=cell# dim3=frame#
load([dest 'parse_behavior.mat']);
load([dest '_cue_movies_lick.mat']);
nCells = size(NR_movie,2);
% out_dir_TS = ['Z:\Analysis\Cue_reward_pairing_analysis\2P\PCA_ICA_recalibration\resp cell troubleshoot\', [session_date, '_', mouse_ID]];
% if ~exist(out_dir_TS)
%     mkdir(out_dir_TS)
% end

%% 2. calculate average timecourses for NR/release events
trial_cond = {'all'};
%find avg and sem across trials then across trials and cells
[avg_NR, sem_NR, avg_NR_all, sem_NR_all] = get_movie_mean_sem(NR_movie, trial_cond, 'all');
[avg_OR, sem_OR, avg_OR_all, sem_OR_all] = get_movie_mean_sem(OR_movie, trial_cond, 'all');
[avg_UR, sem_UR, avg_UR_all, sem_UR_all] = get_movie_mean_sem(UR_movie, trial_cond, 'all');

%average and sem across Trials for no lick trials
[avg_NR_nolick, sem_NR_nolick, avg_NR_all_nolick, sem_NR_all_nolick] = get_movie_mean_sem(NR_movie_nolick, trial_cond, 'nolick');
[avg_OR_nolick, sem_OR_nolick, avg_OR_all_nolick, sem_OR_all_nolick] = get_movie_mean_sem(OR_movie_nolick, trial_cond, 'nolick');
[avg_UR_nolick, sem_UR_nolick, avg_UR_all_nolick, sem_UR_all_nolick] = get_movie_mean_sem(UR_movie_nolick, trial_cond', 'nolick');

if strcmp(trial_cond, 'half')
    save([dest '_cell_TCs_half.mat'], 'avg_NR', 'sem_NR', 'avg_OR', 'sem_OR', 'avg_UR', 'sem_UR', ...
        'avg_NR_nolick', 'sem_NR_nolick', 'avg_OR_nolick', 'sem_OR_nolick', 'avg_UR_nolick', 'sem_UR_nolick');
elseif strcmp(trial_cond, 'all')
    save([dest '_cell_TCs.mat'], 'avg_NR', 'sem_NR', 'avg_OR', 'sem_OR', 'avg_UR', 'sem_UR', ...
        'avg_NR_nolick', 'sem_NR_nolick', 'avg_OR_nolick', 'sem_OR_nolick', 'avg_UR_nolick', 'sem_UR_nolick');
end

%% 3) exclude trials which excede x stds from the mean
pre_rew_frames = pre_cue_frames + round((input.RewardDelayDurationMs + mean(cell2mat(input.reactTimesMs)))/ifi);
artifact_rng = [1:pre_rew_frames+1500/double(ifi)];
if session_subset ==4
    artifact_rng = [1:pre_rew_frames+1000/double(ifi)];
end
outlier_dur = 1;
pre_cue_outlier_dur = 1;
NR_movie_avg = squeeze(mean(NR_movie,1));
NR_movie_std = squeeze(std(NR_movie,1));
NR_movie_thresh = NR_movie_avg + 4*NR_movie_std;
% NR_movie_thresh2 = NR_movie_avg + 5*NR_movie_std;
for cell_num = 1:nCells
    %NR trials
    outlier_frames = squeeze(NR_movie(:,cell_num,:)) > repmat(NR_movie_thresh(cell_num,:), size(NR_movie,1) ,1) ;
    outliers_per_trial = sum(outlier_frames(:,[artifact_rng]), 2);
    NR_ex_trials{cell_num} = find(outliers_per_trial > outlier_dur); 
    outliers_pre_cue = sum(outlier_frames(:,[1:pre_cue_frames]), 2);
    NR_ex_trials2{cell_num} = find(outliers_pre_cue > pre_cue_outlier_dur); 
    %OR trials
    if ~isempty(OR_movie)
        outlier_frames = squeeze(OR_movie(:,cell_num,:)) > repmat(NR_movie_thresh(cell_num,:), size(OR_movie,1) ,1) ;
        outliers_per_trial = sum(outlier_frames(:,[artifact_rng]), 2);
        OR_ex_trials{cell_num} = find(outliers_per_trial > outlier_dur); 
        outliers_pre_cue = sum(outlier_frames(:,[1:pre_cue_frames]), 2);
        OR_ex_trials2{cell_num} = find(outliers_pre_cue > pre_cue_outlier_dur);  
    else 
        OR_ex_trials = [];  OR_ex_trials2 = [];
    end
    %UR trials
    if ~isempty(UR_movie)
        outlier_frames = squeeze(UR_movie(:,cell_num,:)) > repmat(NR_movie_thresh(cell_num,:), size(UR_movie,1) ,1) ;
        outliers_per_trial = sum(outlier_frames(:,[artifact_rng]), 2);
        UR_ex_trials{cell_num} = find(outliers_per_trial > outlier_dur); 
        outliers_pre_cue = sum(outlier_frames(:,[1:pre_cue_frames]), 2);
        UR_ex_trials2{cell_num} = find(outliers_pre_cue > pre_cue_outlier_dur); 
    else 
        UR_ex_trials = [];  UR_ex_trials2 = [];
    end
end
ex_ops.artifact_rng = artifact_rng; ex_ops.outlier_val = outlier_dur;  ex_ops.pre_cue_outlier_dur = pre_cue_outlier_dur; ex_ops.NR_movie_thresh=NR_movie_thresh;
save([dest, '\', '_exclude_trials.mat'], 'NR_ex_trials', 'NR_ex_trials2', 'OR_ex_trials', 'OR_ex_trials2', 'UR_ex_trials', 'UR_ex_trials2', 'ex_ops');   %'exclude_trials', 'exclude_trials2', 'ex_ops');

%% 4) look for responsive neurons
ops.cue_frame = pre_cue_frames+1;
if session_subset == 4
    ops.rew_frame = pre_cue_frames+1 + round(1080/double(ifi));
else
    ops.rew_frame = pre_cue_frames+1 + round(580/double(ifi));
end

%remove trials from specific neurons
NR_3D_ex_nolick = NR_movie_nolick;
OR_3D_ex_nolick = OR_movie_nolick;
UR_3D_ex_nolick = UR_movie_nolick;
for cell_num = 1:nCells
    this_cell_ex_NR = unique([NR_ex_trials{cell_num}', NR_ex_trials2{cell_num}']);
    NR_3D_ex_nolick([this_cell_ex_NR],cell_num,:) = NaN;
    if ~isempty(OR_movie_nolick)
        this_cell_ex_OR = unique([OR_ex_trials{cell_num}', OR_ex_trials2{cell_num}']);
        OR_3D_ex_nolick([this_cell_ex_OR],cell_num,:) = NaN;
    end
    if ~isempty(UR_movie_nolick)
        this_cell_ex_UR = unique([UR_ex_trials{cell_num}', UR_ex_trials2{cell_num}']);
        UR_3D_ex_nolick([this_cell_ex_UR],cell_num,:) = NaN;
    end
end

%POSITIVE CUE responsive neurons 
ops.base_cue_buffer = 500;
ops.resp_cue_buffer = 400;
ops.effect_sign = 'pos';
ops.event_type = 'cue';
ops.ifi = ifi;
[NR_Cue_h_pos, NR_Cue_p_pos, NR_Cue_resp_cells_pos, NR_Cue_resp_avg_pos, NR_Cue_resp_sem_pos, NR_Cue_base_pos, NR_Cue_resp_pos] = findRespCell_CRP(NR_3D_ex_nolick, ops); %NR_movie_nolick
[OR_Cue_h_pos, OR_Cue_p_pos, OR_Cue_resp_cells_pos, OR_Cue_resp_avg_pos, OR_Cue_resp_sem_pos, OR_Cue_base_pos, OR_Cue_resp_pos] = findRespCell_CRP(OR_3D_ex_nolick, ops); % OR_movie_nolick
[UR_Cue_h_pos, UR_Cue_p_pos, UR_Cue_resp_cells_pos, UR_Cue_resp_avg_pos, UR_Cue_resp_sem_pos, UR_Cue_base_pos, UR_Cue_resp_pos] = findRespCell_CRP(UR_3D_ex_nolick, ops); % UR_movie_nolick

%NEGATIVE CUE responsive neurons
ops.effect_sign =  'neg';
[NR_Cue_h_neg, NR_Cue_p_neg, NR_Cue_resp_cells_neg, NR_Cue_resp_avg_neg, NR_Cue_resp_sem_neg, NR_Cue_base_neg, NR_Cue_resp_neg] = findRespCell_CRP(NR_3D_ex_nolick, ops);
[OR_Cue_h_neg, OR_Cue_p_neg, OR_Cue_resp_cells_neg, OR_Cue_resp_avg_neg, OR_Cue_resp_sem_neg, OR_Cue_base_neg, OR_Cue_resp_neg] = findRespCell_CRP(OR_3D_ex_nolick, ops);
[UR_Cue_h_neg, UR_Cue_p_neg, UR_Cue_resp_cells_neg, UR_Cue_resp_avg_neg, UR_Cue_resp_sem_neg, UR_Cue_base_neg, UR_Cue_resp_neg] = findRespCell_CRP(UR_3D_ex_nolick, ops);

%POSITIVE REWARD responsive neurons
ops.event_type = 'reward';
ops.effect_sign = 'pos';
ops.base_cue_buffer = 400;
ops.resp_cue_buffer = 1000;
%pre_rew_frames = pre_cue_frames + round((trial_outcome.normalReward(1) - trial_outcome.normalRewardCue(1))/ifi);
[NR_Rew_h_pos, NR_Rew_p_pos, NR_Rew_resp_cells_pos, NR_Rew_resp_avg_pos, NR_Rew_resp_sem_pos, NR_Rew_base_pos, NR_Rew_resp_pos] = findRespCell_CRP(NR_3D_ex_nolick, ops);
[OR_Rew_h_pos, OR_Rew_p_pos, OR_Rew_resp_cells_pos, OR_Rew_resp_avg_pos, OR_Rew_resp_sem_pos, OR_Rew_base_pos, OR_Rew_resp_pos] = findRespCell_CRP(OR_3D_ex_nolick, ops);
[UR_Rew_h_pos, UR_Rew_p_pos, UR_Rew_resp_cells_pos, UR_Rew_resp_avg_pos, UR_Rew_resp_sem_pos, UR_Rew_base_pos, UR_Rew_resp_pos] = findRespCell_CRP(UR_3D_ex_nolick, ops);

%NEGATIVE REWARD responsive neurons
ops.effect_sign =  'neg';
[NR_Rew_h_neg, NR_Rew_p_neg, NR_Rew_resp_cells_neg, NR_Rew_resp_avg_neg, NR_Rew_resp_sem_neg, NR_Rew_base_neg, NR_Rew_resp_neg] = findRespCell_CRP(NR_3D_ex_nolick, ops);
[OR_Rew_h_neg, OR_Rew_p_neg, OR_Rew_resp_cells_neg, OR_Rew_resp_avg_neg, OR_Rew_resp_sem_neg, OR_Rew_base_neg, OR_Rew_resp_neg] = findRespCell_CRP(OR_3D_ex_nolick, ops);
[UR_Rew_h_neg, UR_Rew_p_neg, UR_Rew_resp_cells_neg, UR_Rew_resp_avg_neg, UR_Rew_resp_sem_neg, UR_Rew_base_neg, UR_Rew_resp_neg] = findRespCell_CRP(UR_3D_ex_nolick, ops);

save([dest, '\', '_cell_resp.mat'], 'ops', ...
    'NR_Cue_base_pos', 'NR_Cue_resp_pos', 'OR_Cue_base_pos', 'OR_Cue_resp_pos', 'UR_Cue_base_pos', 'UR_Cue_resp_pos', ...
    'NR_Rew_base_pos', 'NR_Rew_resp_pos', 'OR_Rew_base_pos', 'OR_Rew_resp_pos', 'UR_Rew_base_pos', 'UR_Rew_resp_pos', ...
	'NR_Cue_base_neg', 'NR_Cue_resp_neg', 'OR_Cue_base_neg', 'OR_Cue_resp_neg', 'UR_Cue_base_neg', 'UR_Cue_resp_neg', ...
    'NR_Rew_base_neg', 'NR_Rew_resp_neg', 'OR_Rew_base_neg', 'OR_Rew_resp_neg', 'UR_Rew_base_neg', 'UR_Rew_resp_neg');

%% 5) exclude cells based on pre-cue activity
pre_cue_ex_cell = zeros(1,nCells);
pre_cue_win = NR_movie_avg(:,[1:pre_cue_frames]);
for cell_num = 1:nCells
    temp_diff = diff(pre_cue_win(cell_num,:)) > 0;
    for frame_num = 1:length(temp_diff)-10
        if temp_diff(frame_num) == 1 & sum(temp_diff(frame_num:frame_num+9))==10 &  mean(pre_cue_win(cell_num,[end-5:end]),2) - mean(pre_cue_win(cell_num,[1:5]),2) > 0.05
            pre_cue_ex_cell(cell_num) = 1;
            break
        end
    end
end
save([dest, '\', '_pre_cue_ex_cell.mat'], 'pre_cue_ex_cell');   

%% 6) Make combinations of responsive subtypes for plotting later
if ~isempty(OR_movie)
    allresp_cells = NR_Cue_h_pos | OR_Cue_h_pos | NR_Rew_h_pos | OR_Rew_h_pos | NR_Cue_h_neg | OR_Cue_h_neg | NR_Rew_h_neg | OR_Rew_h_neg;
    allresp_cells_pos = NR_Cue_h_pos | OR_Cue_h_pos | NR_Rew_h_pos | OR_Rew_h_pos;
    allresp_cells_neg = NR_Cue_h_neg | OR_Cue_h_neg | NR_Rew_h_neg | OR_Rew_h_neg;
    cue_cells = NR_Cue_h_pos | OR_Cue_h_pos | NR_Cue_h_neg | OR_Cue_h_neg;
    cue_cells_pos = NR_Cue_h_pos | OR_Cue_h_pos;
    cue_cells_neg = NR_Cue_h_neg | OR_Cue_h_neg;
    rew_cells = NR_Rew_h_pos | OR_Rew_h_pos | NR_Rew_h_neg | OR_Rew_h_neg;
    rew_cells_pos = NR_Rew_h_pos | OR_Rew_h_pos;
    rew_cells_neg = NR_Rew_h_neg | OR_Rew_h_neg;
else
    allresp_cells = NR_Cue_h_pos | NR_Rew_h_pos | NR_Cue_h_neg | NR_Rew_h_neg;
    allresp_cells_pos = NR_Cue_h_pos | NR_Rew_h_pos;
    allresp_cells_neg = NR_Cue_h_neg | NR_Rew_h_neg;
    cue_cells = NR_Cue_h_pos | NR_Cue_h_neg;
    cue_cells_pos = NR_Cue_h_pos;
    cue_cells_neg = NR_Cue_h_neg;
    rew_cells = NR_Rew_h_pos | NR_Rew_h_neg;
    rew_cells_pos = NR_Rew_h_pos;
    rew_cells_neg = NR_Rew_h_neg;
end

if ~isempty(UR_movie)
    allresp_cells = allresp_cells | UR_Cue_h_pos | UR_Rew_h_pos | UR_Cue_h_neg | UR_Rew_h_neg;
    allresp_cells_pos = allresp_cells_pos | UR_Cue_h_pos | UR_Rew_h_pos;
    allresp_cells_neg = allresp_cells_neg | UR_Cue_h_neg | UR_Rew_h_neg;
    cue_cells = cue_cells | UR_Cue_h_pos | UR_Cue_h_neg;
    cue_cells_pos = cue_cells_pos | UR_Cue_h_pos;
    cue_cells_neg = cue_cells_neg | UR_Cue_h_neg;
    rew_cells = rew_cells | UR_Rew_h_pos | UR_Rew_h_neg;
    rew_cells_pos = rew_cells_pos | UR_Rew_h_pos;
    rew_cells_neg = rew_cells_neg | UR_Rew_h_neg;
end

cue_cells_perc = sum(cue_cells) / nCells;
cue_cells_perc_pos = sum(cue_cells_pos) / nCells;
cue_cells_perc_neg = sum(cue_cells_neg) / nCells;
rew_cells_perc = sum(rew_cells) / nCells;
rew_cells_perc_pos = sum(rew_cells_pos) / nCells;
rew_cells_perc_neg = sum(rew_cells_neg) / nCells;

save([dest '\_cell_categories.mat'], 'NR_Cue_resp_cells_pos', 'NR_Cue_resp_cells_neg', 'OR_Cue_resp_cells_pos', 'OR_Cue_resp_cells_neg', ...
    'UR_Cue_resp_cells_pos', 'UR_Cue_resp_cells_neg', 'NR_Rew_resp_cells_pos', 'NR_Rew_resp_cells_neg', 'OR_Rew_resp_cells_pos', 'OR_Rew_resp_cells_neg', 'NR_Rew_resp_cells_pos', 'NR_Rew_resp_cells_neg',...
    'allresp_cells', 'allresp_cells_pos', 'allresp_cells_neg', 'nCells', 'cue_cells', 'cue_cells_pos', 'cue_cells_neg',  'rew_cells', 'rew_cells_pos', 'rew_cells_neg', ...
    'cue_cells_perc', 'cue_cells_perc_pos', 'cue_cells_perc_neg', 'rew_cells_perc', 'rew_cells_perc_pos', 'rew_cells_perc_neg');
save([dest '\_pvals.mat'], 'NR_Cue_p_pos', 'NR_Cue_p_neg', 'NR_Cue_h_pos', 'NR_Cue_h_neg', 'NR_Rew_p_pos', 'NR_Rew_p_neg', 'NR_Rew_h_pos', 'NR_Rew_h_neg', ...
    'OR_Cue_p_pos', 'OR_Cue_p_neg', 'OR_Cue_h_pos', 'OR_Cue_h_neg', 'OR_Rew_p_pos', 'OR_Rew_p_pos', 'OR_Rew_h_pos', 'OR_Rew_h_neg', ...
    'UR_Cue_p_pos', 'UR_Cue_p_neg', 'UR_Cue_h_pos', 'UR_Cue_h_neg', 'UR_Rew_p_pos', 'UR_Rew_p_neg', 'UR_Rew_h_pos', 'UR_Rew_h_neg');



