% triggers time-courses off of event times
% 1. finds frame and lever times based on events and outcomes
% 2. obtain df/f timecourse
% 3. create event triggered movies



%load frame and lever info
% if exist([dest 'parse_behavior.mat']) == 2
%     load([dest 'parse_behavior.mat'])
%     [lever, frame_info, trial_outcome, lick_data] = cleanBehav(input, ifi, holdT_min);
%     load([dest '_dFOverF_TC.mat']);
% else
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
    frame_info.ifi = ifi;
    
    data_dest = [dest 'parse_behavior.mat'];
%     save(data_dest, 'lever', 'frame_info', 'trial_outcome', 'Sampeling_rate', 'holdT_min', 'ifi')
    
    %% 2. Obtain a df/f TC from baseline times
    data_tc = tc_avg;
    data_tc = data_tc';
    
    % avoid choosing empty trial
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
    
    
    % extrack F for entire trial and baseF from baseline_timesMs (500-300ms before
    % press
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
        end
    end
    data_dest2 = [dest '_dFOverF_TC.mat'];
    save(data_dest2, 'tc_dfoverf')
% end
%% 3. create event triggered movies
if input.doLeverSolenoidAllTrials == 1 || input.doLever == 1
    func = @mean;
    pre_release_ms = 500;
    post_release_ms = 1000;
    pre_release_frames = round(pre_release_ms./double(ifi));
    post_release_frames = round(post_release_ms./double(ifi));
    
    pre_release_ms2 = 500;
    post_release_ms2 = 800;
    pre_release_frames2 = round(pre_release_ms2./double(ifi));
    post_release_frames2 = round(post_release_ms2./double(ifi));
    
    %successes
    use_ev_success = trial_outcome.success_time;
    %     if strcmp(input.trialOutcomeCell{1}, 'success')  %removing first and last trials from consideration if they were successful
    %         use_ev_success(1) = [];
    %     elseif strcmp(input.trialOutcomeCell{end}, 'success')
    %         use_ev_success(end) = [];
    %     end
    do_lickAna = 0; do_alignLick = 0;
    [success_movie, ~, succ_hold_dur, ~, success_lick] = trigger_movie_by_event_2P(tc_dfoverf, frame_info, ...
        use_ev_success, pre_release_frames, post_release_frames, lick_data, trial_outcome.succ_hold_dur, do_lickAna, do_alignLick);
    
    [success_movie_long, ~, ~, ~, success_lick_long] = trigger_movie_by_event_2P(tc_dfoverf, frame_info, ...
        use_ev_success, pre_release_frames2, post_release_frames2, lick_data, [], do_lickAna, do_alignLick);
    
    use_ev_early_success = trial_outcome.early_correct_time;
    early_success_movie =                              trigger_movie_by_event_2P(tc_dfoverf, frame_info, ...
        use_ev_early_success, pre_release_frames, post_release_frames, lick_data, [], do_lickAna, do_alignLick);
    
    use_ev_late_success = trial_outcome.late_correct_time;
    late_success_movie =                               trigger_movie_by_event_2P(tc_dfoverf, frame_info, ...
        use_ev_late_success, pre_release_frames, post_release_frames, lick_data, [], do_lickAna, do_alignLick);
    
    %failures
    use_ev_fail = trial_outcome.early_time;
    
    [fail_movie, remove_idx, fail_hold_dur, ~, fail_lick] = trigger_movie_by_event_2P(tc_dfoverf, frame_info, ...
        use_ev_fail, pre_release_frames, post_release_frames, lick_data, trial_outcome.fail_hold_dur, do_lickAna, do_alignLick);
    
    use_ev_early_fail = trial_outcome.early_early_time;
    [early_fail_movie, ~, ~,~,  early_fail_lick] =           trigger_movie_by_event_2P(tc_dfoverf, frame_info, ...
        use_ev_early_fail, pre_release_frames, post_release_frames, lick_data, [], do_lickAna, do_alignLick);
    
    use_ev_late_fail = trial_outcome.late_early_time;
    [late_fail_movie, ~,~,~, late_fail_lick] =               trigger_movie_by_event_2P(tc_dfoverf, frame_info, ...
        use_ev_late_fail, pre_release_frames, post_release_frames, lick_data, [], do_lickAna, do_alignLick);
    
%     lick_data.lickRateF(remove_idx) = [];
    
    [fail_movie_long, ~, ~, ~, fail_lick_long] = trigger_movie_by_event_2P(tc_dfoverf, frame_info, ...
        use_ev_fail, pre_release_frames2, post_release_frames2, lick_data, [], do_lickAna, do_alignLick);
    
    %tooFast correct
    use_ev_tooFastCorrect = trial_outcome.tooFastCorrects;
   
    [tooFastCorrect_movie,~, tf_hold_dur] = trigger_movie_by_event_2P(tc_dfoverf, frame_info, ...
        use_ev_tooFastCorrect, pre_release_frames, post_release_frames, [], trial_outcome.tf_hold_dur, do_lickAna, do_alignLick);
    
    tooFastCorrect_movie_long = trigger_movie_by_event_2P(tc_dfoverf, frame_info, ...
        use_ev_tooFastCorrect, pre_release_frames2, post_release_frames2, [], [], do_lickAna, do_alignLick);
    
    use_ev_early_tooFast = trial_outcome.early_tooFast_time;
    early_tooFast_movie =               trigger_movie_by_event_2P(tc_dfoverf, frame_info, ...
        use_ev_early_tooFast, pre_release_frames, post_release_frames, lick_data, [], do_lickAna, do_alignLick);
    
    use_ev_late_tooFast = trial_outcome.late_tooFast_time;
    late_tooFast_movie =               trigger_movie_by_event_2P(tc_dfoverf, frame_info, ...
        use_ev_late_tooFast, pre_release_frames, post_release_frames, lick_data, [], do_lickAna, do_alignLick);
    
    %fidgets
    use_ev_fidget = trial_outcome.fidget;
    %     if strcmp(input.trialOutcomeCell{1}, 'failure')   %removing first and last trials from consideration if they were failures
    %         use_ev_fidget(1) = [];
    %     elseif strcmp(input.trialOutcomeCell{end}, 'failure')
    %         use_ev_fidget(end) = [];
    %     end
    fidget_movie = trigger_movie_by_event_2P(tc_dfoverf, frame_info, ...
        use_ev_fidget, pre_release_frames, post_release_frames, [], [], do_lickAna, do_alignLick);
    
    if isfield(input, 'rewardOmissionPercent') && input.rewardOmissionPercent ~= 0
        %Reward Omission
        use_ev_omitReward = trial_outcome.omitReward;
        
        [omitReward_movie, ~, ~, ~, omitR_lick] = trigger_movie_by_event_2P(tc_dfoverf, frame_info, ...
            use_ev_omitReward, pre_release_frames, post_release_frames, lick_data, [], do_lickAna, do_alignLick);
        
        [omitReward_movie_long, ~, ~, ~, omitR_lick_long] = trigger_movie_by_event_2P(tc_dfoverf, frame_info, ...
        use_ev_omitReward, pre_release_frames2, post_release_frames2, lick_data, [], do_lickAna, do_alignLick);
    else
        omitReward_movie = []; omitReward_movie_long = []; 
        omitR_lick = []; omitR_lick_long =[];
    end
    
    if isfield(input, 'itiRewardUnexpectPercent') && input.itiRewardUnexpectPercent ~= 0
        %Unexpected Iti Reward
        use_ev_itiReward = trial_outcome.ItiunexpReward;
        
        [itiReward_movie, ~, ~, ~, itiR_lick] = trigger_movie_by_event_2P(tc_dfoverf, frame_info, ...
            use_ev_itiReward, pre_release_frames, post_release_frames, lick_data, [], do_lickAna, 1);
        
        [itiReward_movie_long, ~, ~, ~, itiR_lick_long] = trigger_movie_by_event_2P(tc_dfoverf, frame_info, ...
        use_ev_itiReward, pre_release_frames2, post_release_frames2, lick_data, [], do_lickAna, 1);
    else
        itiReward_movie = []; itiR_lick= [];
        itiReward_movie_long = []; itiR_lick_long= [];
    end
    
    % ---- Trigger movie off all lever presses at trial start
    pre_press_ms = 500;
    post_press_ms = 1000;
    pre_press_frames = round(pre_press_ms./double(ifi));
    post_press_frames = round(post_press_ms./double(ifi));
    
    pre_press_ms2 = 500;
    post_press_ms2 = 1000;
    pre_press_frames2 = round(pre_press_ms2./double(ifi));
    post_press_frames2 = round(post_press_ms2./double(ifi));
    
    %     pressTime = NaN(1,length(lever.baseline_timesMs));
    %     releaseTime = NaN(1,length(lever.baseline_timesMs));
    % for iT = 2:length(lever.baseline_timesMs)-1
    %     if ~isempty(input.leverTimesUs{iT})
    %         leverTimes = round((cell2mat(input.leverTimesUs(iT))-input.counterTimesUs{1}(1))./1000);
    %     end
    %     if ~isnan(trial_outcome.ind_press_prerelease(iT)) & size(input.leverTimesUs{iT},2)>1
    %         pressTime(1,iT) = leverTimes(trial_outcome.ind_press_prerelease(iT));
    %         releaseTime(1,iT) = leverTimes(trial_outcome.ind_press_prerelease(iT)+1);
    %     end
    % end
    
    pressTime = lever.press(2:end-1);
    releaseTime = lever.release(2:end-1);
    
    use_ev_press = pressTime;
    use_ev_press(isnan(use_ev_press)) = [];
    
    press_movie = trigger_movie_by_event_2P(tc_dfoverf, frame_info, ...
        use_ev_press, pre_press_frames, post_press_frames, [], [], do_lickAna, do_alignLick);
    
    use_ev_press_success = trial_outcome.success_ptime;
    press_success_movie = trigger_movie_by_event_2P(tc_dfoverf, frame_info, ...
        use_ev_press_success, pre_press_frames, post_press_frames, [], [], do_lickAna, do_alignLick);
    
     use_ev_press_fail = trial_outcome.early_ptime;
    press_fail_movie = trigger_movie_by_event_2P(tc_dfoverf, frame_info, ...
        use_ev_press_fail, pre_press_frames, post_press_frames, [], [], do_lickAna, do_alignLick);
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
    
    release_long_movie = trigger_movie_by_event_2P(tc_dfoverf, frame_info, ...
        use_long_release, pre_release_frames, post_release_frames, [], [], do_lickAna, do_alignLick);
    
    press_200_movie = trigger_movie_by_event_2P(tc_dfoverf, frame_info, ...
        use_200_press, pre_press_frames, post_press_frames, [], [], do_lickAna, do_alignLick);
    press_500_movie = trigger_movie_by_event_2P(tc_dfoverf, frame_info, ...
        use_500_press, pre_press_frames, post_press_frames, [], [], do_lickAna, do_alignLick);
    press_long_movie = trigger_movie_by_event_2P(tc_dfoverf, frame_info, ...
        use_long_press, pre_press_frames, post_press_frames, [], [], do_lickAna, do_alignLick);
    press_long_movie_long = trigger_movie_by_event_2P(tc_dfoverf, frame_info, ...
        use_long_press, pre_press_frames2, post_press_frames2, [], [], do_lickAna, do_alignLick);
    
    
    % Trigger movie for release 
    pre_cuerelease_frames = round(400/double(ifi));
    post_cuerelease_frames = round(400/double(ifi));
    
    use_ev_cuerelease = lever.cue_release;
    [cuerelease_movie, cue_remove_idx] = trigger_movie_by_event_2P(tc_dfoverf, frame_info, ...
        use_ev_cuerelease, pre_cuerelease_frames, post_cuerelease_frames, [], [], do_lickAna, do_alignLick);
    
    % Trigger movie for licking
    if ~isempty(lick_data)
        pre_lick_frames = round(2500/double(ifi));%lick_data.buffer_lick;
        post_lick_frames = round(2500/double(ifi));
        use_ev_lick = lick_data.bout_onsetTime;
        
        [lick_movie, ~,~,~, lickb_lick] = trigger_movie_by_event_2P(tc_dfoverf, frame_info, ...
            use_ev_lick, pre_lick_frames, post_lick_frames, lick_data, [], do_lickAna, do_alignLick);
        
        use_ev_single_lick = lick_data.single_lickTime;
        
        [single_lick_movie, ~,~,~, licks_lick] = trigger_movie_by_event_2P(tc_dfoverf, frame_info, ...
            use_ev_single_lick, pre_lick_frames, post_lick_frames, lick_data, [], do_lickAna, do_alignLick);
    else
        lick_movie = []; single_lick_movie = []; lickb_lick = []; licks_lick = [];
    end
    
    trial_outcome.succ_hold_dur = succ_hold_dur;
    trial_outcome.fail_hold_dur = fail_hold_dur;
    trial_outcome.tf_hold_dur = tf_hold_dur;
    
    save([dest  'parse_behavior.mat'], 'lever', 'frame_info', 'trial_outcome', 'lick_data', 'Sampeling_rate', 'holdT_min', 'ifi')
    
    save([dest '_release_movies.mat'],'fail_movie','success_movie', 'tooFastCorrect_movie', 'cuerelease_movie','fail_movie_long','success_movie_long', ...
        'tooFastCorrect_movie_long', 'release_long_movie', 'omitReward_movie', 'omitReward_movie_long', 'itiReward_movie', 'itiReward_movie_long', 'use_ev_fidget','pre_release_frames','post_release_frames','ifi',...
        'pre_release_frames2', 'post_release_frames2', 'pre_cuerelease_frames', 'post_cuerelease_frames', 'cue_remove_idx', 'early_fail_movie', 'late_fail_movie', ...
        'early_success_movie', 'late_success_movie', 'early_tooFast_movie', 'late_tooFast_movie');
    
    save([dest '_press_movies.mat'],'press_200_movie','press_500_movie','press_long_movie','press_long_movie_long','press_movie',...
        'pre_press_frames', 'post_press_frames', 'ifi', 'press_success_movie', 'press_fail_movie');
    
    save([dest '_lick_movies.mat'], 'lick_movie', 'single_lick_movie');
    
    save([dest '_lick_stats.mat'], 'success_lick', 'itiR_lick', 'omitR_lick', 'fail_lick', 'success_lick_long', 'fail_lick_long', 'omitR_lick_long', 'itiR_lick_long',...
        'lickb_lick', 'licks_lick', 'early_fail_lick', 'late_fail_lick');
else
    pre_cue_ms = 1000;
    post_cue_ms = 2500;
    pre_cue_frames = round(pre_cue_ms./double(ifi));
    post_cue_frames = round(post_cue_ms./double(ifi));
    
    use_ev_NR = trial_outcome.normalRewardCue;
    use_ev_NR(isnan(use_ev_NR)) = [];
    
    do_lickAna = 1; do_alignLick = 0;
    [NR_movie, ~, ~, ~, lick_trace_NR, lick_trace_NR_10ms, NR_lick_info] = trigger_movie_by_event_2P(tc_dfoverf, frame_info, ...
        use_ev_NR, pre_cue_frames, post_cue_frames, lick_data, [], do_lickAna, do_alignLick);
    
    % get rid of trials with licking bout
    NR_movie_nolick = nan(size(NR_movie));
    NR_movie_nolick(~NR_lick_info.lickTrial,:,:) = NR_movie(~NR_lick_info.lickTrial, :, :);
    
    use_ev_OR = trial_outcome.omitRewardCue;
    use_ev_OR(isnan(use_ev_OR)) = [];
    
    [OR_movie, ~, ~, ~, lick_trace_OR, lick_trace_OR_10ms, OR_lick_info] = trigger_movie_by_event_2P(tc_dfoverf, frame_info, ...
        use_ev_OR, pre_cue_frames, post_cue_frames, lick_data, [], do_lickAna, do_alignLick);
    
    % get rid of trials with licking bout
    OR_movie_nolick = nan(size(OR_movie));
    OR_movie_nolick(~OR_lick_info.lickTrial,:,:) = OR_movie(~OR_lick_info.lickTrial, :, :);
    
    use_ev_UR = trial_outcome.unexpRewardCue;
    use_ev_UR(isnan(use_ev_UR)) = [];
    
    [UR_movie, ~, ~, ~, lick_trace_UR, lick_trace_UR_10ms, UR_lick_info] = trigger_movie_by_event_2P(tc_dfoverf, frame_info, ...
        use_ev_UR, pre_cue_frames, post_cue_frames, lick_data, [], do_lickAna, do_alignLick);
    
    % get rid of trials with licking bout
    UR_movie_nolick = nan(size(UR_movie));
    UR_movie_nolick(~UR_lick_info.lickTrial,:,:) = UR_movie(~UR_lick_info.lickTrial, :, :);
    
    save([dest  'parse_behavior.mat'], 'lever', 'frame_info', 'trial_outcome', 'lick_data', 'Sampeling_rate', 'holdT_min', 'ifi')
    save([dest '_cue_movies.mat'], 'NR_movie','OR_movie', 'UR_movie', 'NR_movie_nolick','OR_movie_nolick', 'UR_movie_nolick', 'pre_cue_frames', 'post_cue_frames', 'ifi');
    save([dest '_cue_movies_lick.mat'], 'lick_trace_NR','lick_trace_OR', 'lick_trace_UR', 'NR_lick_info', 'OR_lick_info', 'UR_lick_info');
end