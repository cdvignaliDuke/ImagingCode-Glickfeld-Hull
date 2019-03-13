% triggers time-courses off of event times
% 1. finds frame and lever times based on events and outcomes
% 2. obtain df/f timecourse
% 3. create event triggered movies
% 4. correct for visual artifact in two datasets
% 5. save data

%% 1. find frame and lever times
if exist('frame_times') == 1
    ifi = (frame_times(end)-frame_times(1))/length(frame_times);
else
    ifi = mode(diff(cell2mat(cellfun(@int64,input.counterTimesUs,'UniformOutput',0))))/1000;
end

load([dest, 'img_reg.mat'], 'laser_on_ind_conserv');
Sampeling_rate = 1000/ifi;

% ---- parse behavior
holdT_min  = 500000;  %us
[lever, frame_info, trial_outcome, lick_data] = cleanBehav_CRP(input, ifi, holdT_min, laser_on_ind_conserv);
frame_info.ifi = ifi;

data_dest = [dest 'parse_behavior.mat'];
save(data_dest, 'lever', 'frame_info', 'trial_outcome', 'Sampeling_rate', 'holdT_min', 'ifi', 'lick_data');

%% 2. Obtain a df/f TC from baseline times

%adjust tc_avg to account for trimmed frames
load([dest, 'img_reg.mat'], 'img_nframes');
if isfield(frame_info, 'laser_on_ind_conserv') & ~strcmp(session, '171227_img067')
    laser_on_ind_conserv = frame_info.laser_on_ind_conserv;
end

%just in case some sessions do not have the laser power toggling.
if exist('laser_on_ind_conserv', 'var') == 1
    data_tc = zeros(frame_info.counter_by_time(end), size(tc_avg,2));
    if laser_on_ind_conserv(end) > img_nframes | laser_on_ind_conserv(end) > size(data_tc,1)
        laser_on_ind_conserv = laser_on_ind_conserv(find(laser_on_ind_conserv<=img_nframes));
    end
    assert(length(laser_on_ind_conserv) == size(tc_avg,1))
    data_tc([laser_on_ind_conserv],:) = tc_avg;
    assert(size(data_tc,1) == frame_info.counter_by_time(end));
    data_tc = data_tc';
else
    data_tc = tc_avg';
end

%quick check to make sure the matrix is oriented correctly
if size(data_tc,1) > 10000
    data_tc = data_tc';
end

%avoid choosing empty trial
empty_ind = find(~cellfun(@isempty, input.counterValues));
if input.counterValues{empty_ind(1)}(1) == 0
    startT = round(input.counterTimesUs{empty_ind(2)}(1)./1000);
    start_i  =  empty_ind(2) + 1;
else
    startT = round(input.counterTimesUs{empty_ind(1)}(1)./1000);
    start_i  =  empty_ind(1) + 1;
end

trial_len = cell2mat(cellfun(@length, input.counterValues, 'UniformOutput', 0));
valid_trial = find(trial_len > 2);
stop_i = valid_trial(end) - 1;

% extract F for entire trial and baseF from baseline_timesMs (500-300ms before press
tc_dfoverf = zeros(size(data_tc));
first_baseline = find(~isnan(lever.baseline_timesMs(1,:)),1, 'first');    %find the first trial / baseline_timesMs window that is not NaN
F_range = [];
for iT = start_i:stop_i;
    % labjack counter could be larger than total frames collected
    if ~isnan(lever.baseline_timesMs(1,iT)) && frame_info.counter(lever.baseline_timesMs(2,iT)) > size(data_tc,2)
        break
    else
        if ~isnan(lever.baseline_timesMs(1,iT));
            F_range = frame_info.counter(lever.baseline_timesMs(1,iT)):frame_info.counter(lever.baseline_timesMs(2,iT));
        elseif isempty(F_range)
            F_range = frame_info.counter(lever.baseline_timesMs(1,first_baseline)):frame_info.counter(lever.baseline_timesMs(2,first_baseline));
        end
        F_avg= mean(data_tc(:,F_range),2);
        if frame_info.counter(cell2mat(input.tThisTrialStartTimeMs(iT+1))-startT) > size(data_tc,2)
            t_range = frame_info.counter(cell2mat(input.tThisTrialStartTimeMs(iT))-startT):size(data_tc,2);
        else
            t_range = frame_info.counter(cell2mat(input.tThisTrialStartTimeMs(iT))-startT):frame_info.counter(cell2mat(input.tThisTrialStartTimeMs(iT+1))-startT);
        end
        t_df = bsxfun(@minus, double(data_tc(:,t_range)), F_avg);
        t_dfoverf = bsxfun(@rdivide, t_df, F_avg);
        tc_dfoverf(:,t_range) = t_dfoverf;
        zeros_in_baselines(iT) = length(find( data_tc(1,F_range)==0 ));
    end
end
data_dest2 = [dest '_dFOverF_TC.mat'];
save(data_dest2, 'tc_dfoverf', 'data_tc')

%% 3. create event triggered movies

%define some variables
pre_cue_ms = 1000;
post_cue_ms = 2500;
pre_cue_frames = round(pre_cue_ms./double(ifi));
post_cue_frames = round(post_cue_ms./double(ifi));
cue_2_rew_delay_ms = mode(trial_outcome.normalReward - trial_outcome.normalRewardCue);
do_lickAna = 1; do_alignLick = 0;

%use event times to extract df/f and licking traces/info ------------
use_ev_NR = trial_outcome.normalRewardCue;
use_ev_NR(isnan(use_ev_NR)) = [];
[NR_movie, ~, ~, ~] = trigger_movie_by_event_2P(tc_dfoverf, frame_info, ...  %NR_movie dim1=trials  dim2=cells  dim3=frames
    use_ev_NR, pre_cue_frames, post_cue_frames, lick_data, [], do_alignLick)
[trigger_licks, trigger_licks_10ms, lick] = trigger_movie_by_event_2P(frame_info, lick_data, use_ev_NR, ...
     pre_cue_frames, post_cue_frames, do_lickAna, cue_2_rew_delay_ms)
% [NR_movie, ~, ~, ~, lick_trace_NR, lick_trace_NR_10ms, NR_lick_info] = trigger_movie_by_event_2P(tc_dfoverf, frame_info, ...
%     use_ev_NR, pre_cue_frames, post_cue_frames, lick_data, [], do_lickAna, do_alignLick, cue_2_rew_delay_ms);  
% get rid of trials with licking bout
NR_movie_nolick = nan(size(NR_movie));
NR_movie_nolick(logical(NR_lick_info.no_lick_cue_to_500) ,:,:) = NR_movie(logical(NR_lick_info.no_lick_cue_to_500), :, :);

%use event times to extract df/f and licking traces/info ------------
use_ev_OR = trial_outcome.omitRewardCue;
use_ev_OR(isnan(use_ev_OR)) = [];
[OR_movie, ~, ~, ~, lick_trace_OR, lick_trace_OR_10ms, OR_lick_info] = trigger_movie_by_event_2P(tc_dfoverf, frame_info, ...
    use_ev_OR, pre_cue_frames, post_cue_frames, lick_data, [], do_lickAna, do_alignLick, cue_2_rew_delay_ms);
% get rid of trials with licking bout
OR_movie_nolick = nan(size(OR_movie));
OR_movie_nolick(logical(OR_lick_info.no_lick_cue_to_500),:,:) = OR_movie(logical(OR_lick_info.no_lick_cue_to_500), :, :);

%use event times to extract df/f and licking traces/info ------------
use_ev_UR = trial_outcome.unexpRewardCue;
use_ev_UR(isnan(use_ev_UR)) = [];
[UR_movie, ~, ~, ~, lick_trace_UR, lick_trace_UR_10ms, UR_lick_info] = trigger_movie_by_event_2P(tc_dfoverf, frame_info, ...
    use_ev_UR, pre_cue_frames, post_cue_frames, lick_data, [], do_lickAna, do_alignLick, cue_2_rew_delay_ms);
% get rid of trials with licking bout
UR_movie_nolick = nan(size(UR_movie));
UR_movie_nolick(logical(UR_lick_info.no_lick_cue_to_500),:,:) = UR_movie(logical(UR_lick_info.no_lick_cue_to_500), :, :);

%% 4. correct for visual artifact in two datasets
if strcmp(session, '180502_img084') | strcmp(session, '180403_img077') | strcmp(session, '180404_img077')
    if strcmp(session, '180502_img084')
        artifact_ind = [2:5]+pre_cue_frames+1;
    elseif strcmp(session, '180403_img077')
        artifact_ind = [2:4]+pre_cue_frames+1;
    elseif strcmp(session, '180404_img077')
        artifact_ind = [3:6]+pre_cue_frames+1;
    end
    for n_cell = 1:size(NR_movie,2)
        for n_trial = 1:size(NR_movie,1) %dim1=trials   dim2=cells   dim3=frames
            NR_movie(n_trial, n_cell, artifact_ind) = linspace(NR_movie(n_trial, n_cell, artifact_ind(1)-1), NR_movie(n_trial, n_cell, artifact_ind(end)+1), length(artifact_ind));
        end
        for n_trial = 1:size(OR_movie,1)
            OR_movie(n_trial, n_cell, artifact_ind) = linspace(OR_movie(n_trial, n_cell, artifact_ind(1)-1), OR_movie(n_trial, n_cell, artifact_ind(end)+1), length(artifact_ind));
        end
        for n_trial = 1:size(NR_movie_nolick,1)
            NR_movie_nolick(n_trial, n_cell, artifact_ind) = linspace(NR_movie_nolick(n_trial, n_cell, artifact_ind(1)-1), NR_movie_nolick(n_trial, n_cell, artifact_ind(end)+1), length(artifact_ind));
        end
        for n_trial = 1:size(OR_movie_nolick,1)
            OR_movie_nolick(n_trial, n_cell, artifact_ind) = linspace(OR_movie_nolick(n_trial, n_cell, artifact_ind(1)-1), OR_movie_nolick(n_trial, n_cell, artifact_ind(end)+1), length(artifact_ind));
        end
    end
end

%% 5. save data
%save([dest  'parse_behavior.mat'], 'lever', 'frame_info', 'trial_outcome', 'lick_data', 'Sampeling_rate', 'holdT_min', 'ifi')
save([dest '_cue_movies.mat'], 'NR_movie','OR_movie', 'UR_movie', 'NR_movie_nolick','OR_movie_nolick', 'UR_movie_nolick', 'pre_cue_frames', 'post_cue_frames', 'ifi', 'cue_2_rew_delay_ms');
save([dest '_cue_movies_lick.mat'], 'lick_trace_NR','lick_trace_OR', 'lick_trace_UR', 'NR_lick_info', 'OR_lick_info', 'UR_lick_info');




