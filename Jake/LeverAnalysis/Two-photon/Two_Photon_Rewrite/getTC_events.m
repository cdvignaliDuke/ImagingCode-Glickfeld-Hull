% triggers time-courses off of event times
% 1. finds frame and lever times based on events and outcomes
% 2. obtain df/f timecourse
% 3. create event triggered movies

%load frame and lever info

%% 1. find frame and lever times
if exist('frame_times') == 1
    ifi = (frame_times(end)-frame_times(1))/length(frame_times);
else
    ifi = mode(diff(cell2mat(cellfun(@int64,input.counterTimesUs,'UniformOutput',0))))/1000;
end

Sampeling_rate = 1000/ifi;

% ---- parse behavior
holdT_min  = 500000;  %us
[lever, frame_info, trial_outcome, lick_data] = cleanBehav(input, ifi, holdT_min);

data_dest = [dest 'parse_behavior.mat'];
save(data_dest, 'lever', 'frame_info', 'trial_outcome', 'Sampeling_rate', 'holdT_min', 'ifi', 'lick_data')

%% 2. Obtain a df/f TC from baseline times
data_tc = tc_avg';   %tc_avg is derived from the .sbx file and is the raw F not df/f

% avoid choosing empty trial
empty_ind = find(~cellfun(@isempty, input.counterValues));
if input.counterValues{empty_ind(1)}(1) == 0
    startT = round(input.counterTimesUs{empty_ind(2)}(1)./1000);
    start_i  =  empty_ind(2) + 1;
else
    startT = round(input.counterTimesUs{empty_ind(1)}(1)./1000);
    start_i  =  empty_ind(1) + 1;
end

%defines the trial length, which trials actually have frames, and which trial will be the last one included in analysis
trial_len = cell2mat(cellfun(@length, input.counterValues, 'UniformOutput', 0));
valid_trial = find(trial_len > 2);
stop_i = valid_trial(end) - 1;

% extrack F for entire trial and baseF from baseline_timesMs (500-300ms before  press
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
        if frame_info.counter(cell2mat(input.tThisTrialStartTimeMs(iT+1))-startT) > size(data_tc,2)  %if the first frame in this trial exceeds the number of frames collected by .sbx...
            t_range = frame_info.counter(cell2mat(input.tThisTrialStartTimeMs(iT))-startT):size(data_tc,2);
        else
            t_range = frame_info.counter(cell2mat(input.tThisTrialStartTimeMs(iT))-startT):frame_info.counter(cell2mat(input.tThisTrialStartTimeMs(iT+1))-startT);
        end
        t_df = bsxfun(@minus, double(data_tc(:,t_range)), F_avg);
        t_dfoverf = bsxfun(@rdivide, t_df, F_avg);
        tc_dfoverf(:,t_range) = t_dfoverf;
    end
end
data_dest2 = [dest '_dFOverF_TC.mat'];
save(data_dest2, 'tc_dfoverf')

%% 3. create event triggered movies   LEVER
if input.doLeverSolenoidAllTrials == 1
    func = @mean;
    pre_release_ms = 500;
    post_release_ms = 500;
    pre_release_frames = round(pre_release_ms./double(ifi));
    post_release_frames = round(post_release_ms./double(ifi));
    
    %successes
    use_ev_success = trial_outcome.success_time;
    if strcmp(input.trialOutcomeCell{1}, 'success')  %removing first and last trials from consideration if they were successful
        use_ev_success(1) = [];
    elseif strcmp(input.trialOutcomeCell{end}, 'success')
        use_ev_success(end) = [];
    end
    success_movie = trigger_movie_by_event(tc_dfoverf, frame_info, ...
        use_ev_success, pre_release_frames, post_release_frames);
    
    %failures
    use_ev_fail = trial_outcome.early_time;
    if strcmp(input.trialOutcomeCell{1}, 'failure')   %removing first and last trials from consideration if they were failures
        use_ev_fail(1) = [];
    elseif strcmp(input.trialOutcomeCell{end}, 'failure')
        use_ev_fail(end) = [];
    end
    fail_movie = trigger_movie_by_event(tc_dfoverf, frame_info, ...
        use_ev_fail, pre_release_frames, post_release_frames);
    
    %tooFast correct
    use_ev_tooFastCorrect = trial_outcome.tooFastCorrects;
    if strcmp(input.trialOutcomeCell{1}, 'success')  %removing first and last trials from consideration if they were successful (tooFastCorrects are coded as successes here)
        use_ev_tooFastCorrect(1) = [];
    elseif strcmp(input.trialOutcomeCell{end}, 'success')
        use_ev_tooFastCorrect(end) = [];
    end
    tooFastCorrect_movie = trigger_movie_by_event(tc_dfoverf, frame_info, ...
        use_ev_tooFastCorrect, pre_release_frames, post_release_frames);
    
    %fidgets
    use_ev_fidget = trial_outcome.fidget;
    if strcmp(input.trialOutcomeCell{1}, 'failure')   %removing first and last trials from consideration if they were failures
        use_ev_fidget(1) = [];
    elseif strcmp(input.trialOutcomeCell{end}, 'failure')
        use_ev_fidget(end) = [];
    end
    fidget_movie = trigger_movie_by_event(tc_dfoverf, frame_info, ...
        use_ev_fidget, pre_release_frames, post_release_frames);
    
    
    % ---- Trigger movie off all lever presses at trial start
    pre_press_ms = 500;
    post_press_ms = 500;
    pre_press_frames = round(pre_press_ms./double(ifi));
    post_press_frames = round(post_press_ms./double(ifi));
    
    pressTime = NaN(1,length(lever.baseline_timesMs));
    releaseTime = NaN(1,length(lever.baseline_timesMs));
    % for iT = 2:length(lever.baseline_timesMs)-1
    %     if ~isempty(input.leverTimesUs{iT})
    %         leverTimes = round((cell2mat(input.leverTimesUs(iT))-input.counterTimesUs{1}(1))./1000);
    %     end
    %     if ~isnan(trial_outcome.ind_press_prerelease(iT)) & size(input.leverTimesUs{iT},2)>1
    %         pressTime(1,iT) = leverTimes(trial_outcome.ind_press_prerelease(iT));
    %         releaseTime(1,iT) = leverTimes(trial_outcome.ind_press_prerelease(iT)+1);
    %     end
    % end
    
    pressTime(2:end-1) = lever.press(2:end-1);
    releaseTime(2:end-1) = lever.release(2:end-1);
    
    use_ev_press = pressTime;
    use_ev_press(isnan(use_ev_press)) = [];
    
    press_movie = trigger_movie_by_event(tc_dfoverf, frame_info, ...
        use_ev_press, pre_press_frames, post_press_frames);
    
    %break up presses by hold time
    holdTime = releaseTime-pressTime;
    ind_200 = find(holdTime<200);
    ind_500 = find(holdTime<500);
    ind_both = find(ismember(ind_500,ind_200));
    ind_500(ind_both) = [];
    ind_long = find(holdTime>=500);
    ind_both = find(ismember(ind_long,ind_500));
    ind_long(ind_both) = [];
    
    use_200_press = pressTime(ind_200);
    use_500_press = pressTime(ind_500);
    use_long_press = pressTime(ind_long);
    
    use_long_release = releaseTime(ind_long);
    
    release_long_movie = trigger_movie_by_event(tc_dfoverf, frame_info, ...
        use_long_release, pre_release_frames, post_release_frames);
    
    press_200_movie = trigger_movie_by_event(tc_dfoverf, frame_info, ...
        use_200_press, pre_press_frames, post_press_frames);
    press_500_movie = trigger_movie_by_event(tc_dfoverf, frame_info, ...
        use_500_press, pre_press_frames, post_press_frames);
    press_long_movie = trigger_movie_by_event(tc_dfoverf, frame_info, ...
        use_long_press, pre_press_frames, post_press_frames);
    
    save([dest '_release_movies.mat'],'fail_movie','success_movie', 'tooFastCorrect_movie', 'release_long_movie', 'use_ev_fidget','pre_release_frames','post_release_frames','ifi');
    
    save([dest '_press_movies.mat'],'press_200_movie','press_500_movie','press_long_movie','press_movie','pre_press_frames', 'post_press_frames', 'ifi');
    
else   %CRP
    pre_cue_ms = 2000;
    post_cue_ms = 2500;
    pre_cue_frames = round(pre_cue_ms./double(ifi));
    post_cue_frames = round(post_cue_ms./double(ifi));
    
    use_ev_NR = trial_outcome.normalReward;
    use_ev_NR(isnan(use_ev_NR)) = [];
    
    [NR_movie, ~, lick_trace_NR, lick_trace_NR_10ms] = trigger_movie_by_event(tc_dfoverf, frame_info, ...
        use_ev_NR, pre_cue_frames, post_cue_frames, lick_data);
    
    use_ev_OR = trial_outcome.omitReward;
    use_ev_OR(isnan(use_ev_OR)) = [];
    
    [OR_movie, ~, lick_trace_OR, lick_trace_OR_10ms] = trigger_movie_by_event(tc_dfoverf, frame_info, ...
        use_ev_OR, pre_cue_frames, post_cue_frames, lick_data);
    
    use_ev_UR = trial_outcome.unexpReward;
    use_ev_UR(isnan(use_ev_UR)) = [];
    
    [UR_movie, ~, lick_trace_UR, lick_trace_UR_10ms] = trigger_movie_by_event(tc_dfoverf, frame_info, ...
        use_ev_UR, pre_cue_frames, post_cue_frames, lick_data);
    
    save([dest '_cue_movies.mat'], 'NR_movie','OR_movie', 'UR_movie', 'pre_cue_frames', 'post_cue_frames', 'ifi');
    save([dest '_cue_movies_lick.mat'], 'lick_trace_NR','lick_trace_OR', 'lick_trace_UR');
    
end