% Quantify df/f onset latency and licking onset latency for the WF data
% across days. 
%use the df/f onset latency code from the lever project. 
%adapt the licking onset latency code from the 2P CRP project, 

clear; 
WF_CRP_list_of_days;
F_TC_dir    = 'Z:\Analysis\WF Lever Analysis\LeverSummaryFolder\';
lick_TC_dir = 'Z:\Analysis\WF Lever Analysis\BxAndAnalysisOutputs\BxOutputs\'; 
out_dir = 'Z:\Analysis\Cue_reward_pairing_analysis\WF\onset_latencies\';
days = days_UR;
ROI_cell = days_UR_ROI;
cue_time_interp = 501;
reward_time_interp = 1101;
no_lick_window = [9:11];

rew_onset_times_all =  [];
rew_onset_times_all_sem = [];
rew_om_onset_times_all =  [];
rew_om_onset_times_all_sem = [];
unexp_onset_times_all =  [];
unexp_onset_times_all_sem = [];
rew_om_peak_to_base_diff_mean = [];
rew_om_peak_to_base_diff_sem = [];
rew_peak_to_base_diff_mean = [];
rew_peak_to_base_diff_sem = [];
unexp_peak_to_base_diff_mean = []; 
unexp_peak_to_base_diff_sem = [];

for ii = 1:length(days)
    if isempty(days{ii});
        continue
    end
    
    %load TCs and and licking -------------------------------------
    load([lick_TC_dir, days{ii}, '_bx_outputs'], 'licking_data');
    lick_trace_rew = licking_data.lick_trace_rew;
    load([F_TC_dir, days{ii}, '_reward_trials']);
    if exist([F_TC_dir, days{ii}, '_rew_om_trials.mat'])
        if strcmp(days{1}, days_1{1}) | strcmp(days{1}, days_post{1}) | strcmp(days{1}, days_1000{1});
            load([F_TC_dir, days{ii}, '_rew_om_trials']);
            lick_trace_rew_om = licking_data.lick_trace_rew_om;
        end
    end
    if exist([F_TC_dir, days{ii}, '_unexp_trials.mat']) && strcmp(days{1}, days_UR{1});
        load([F_TC_dir, days{ii}, '_unexp_trials']);
        lick_trace_unexp_rew = licking_data.lick_trace_unexp_rew;
    end
    
    %only use 1st half of trials from day 1 -------------------------
    if strcmp(days{1}, days_1{1})
        rewarded_roi = rewarded_roi([1:round(size(rewarded_roi,1)/2)], : ,:);
        rew_om_roi = rew_om_roi([1:round(size(rew_om_roi,1)/2)],: ,:);
        lick_trace_rew = lick_trace_rew([1:round(size(rewarded_roi,1)/2)] ,:);
        lick_trace_rew_om = lick_trace_rew_om([1:round(size(rew_om_roi,1)/2)] ,:);
    end
    
    %identify trials where lick bout starts before reward delivery  ------------
    %remove trials with licking between cue and reward
    rew_no_lick_inx = remove_trials_with_licks(lick_trace_rew, no_lick_window, days{ii});
    lick_trace_rew = lick_trace_rew([rew_no_lick_inx], :);
    rewarded_roi = rewarded_roi([rew_no_lick_inx], :, :);
    disp([days{ii}, ' rewarded n=', num2str(length(rew_no_lick_inx))]);
    if exist('rew_om_roi')
        rew_om_no_lick_inx = remove_trials_with_licks(lick_trace_rew_om, no_lick_window, days{ii});
        lick_trace_rew_om = lick_trace_rew_om([rew_om_no_lick_inx], :);
        rew_om_roi = rew_om_roi([rew_om_no_lick_inx], :, :);
        disp([days{ii}, ' omission n=', num2str(length(rew_om_no_lick_inx))]);
    elseif exist('unexp_rew_roi');
        unexp_no_lick_inx = remove_trials_with_licks(lick_trace_unexp_rew, no_lick_window, days{ii});
        lick_trace_unexp_rew = lick_trace_unexp_rew([unexp_no_lick_inx], :);
        unexp_rew_roi = unexp_rew_roi([unexp_no_lick_inx], :, :);
        disp([days{ii}, ' unexpected n=', num2str(length(unexp_no_lick_inx))]);
    end
    
    %interpolate data ---------------------------------
    XVector = [1:.01:size(rewarded_roi,3)];
    rew_roi_interp = nan(length(XVector), size(rewarded_roi,1), size(rewarded_roi,2)); %dims: 1=T2 2=trial# 3=ROI#
    for ROI_num = 1:size(rewarded_roi,2); %for each roi...
        rew_temp = squeeze(rewarded_roi(:,ROI_num,:))';  %dims 1=Time 2=trial number
        rew_temp = interp1(rew_temp, XVector);
        rew_roi_interp(:,:,ROI_num) = rew_temp;  %dims 1=Time 2=trialNumber  3=roi
    end
    rew_roi_interp = squeeze(mean(rew_roi_interp(:,:, ROI_cell{ii}),3));
    
    if exist('rew_om_roi')
        rew_om_roi_interp = nan(length(XVector), size(rew_om_roi,1), size(rew_om_roi,2)); %dims: 1=T2 2=trial# 3=ROI#
        for ROI_num = 1:size(rew_om_roi,2); %for each roi...
            rew_om_temp = squeeze(rew_om_roi(:,ROI_num,:))';  %dims 1=Time 2=trial number
            rew_om_temp = interp1(rew_om_temp, XVector);
            rew_om_roi_interp(:,:,ROI_num) = rew_om_temp;  %dims 1=Time 2=trialNumber  3=roi
        end
        rew_om_roi_interp = squeeze(mean(rew_om_roi_interp(:,:, ROI_cell{ii}),3));
    elseif exist('unexp_rew_roi')
        unexp_rew_roi_interp = nan(length(XVector), size(unexp_rew_roi,1), size(unexp_rew_roi,2)); %dims: 1=T2 2=trial# 3=ROI#
        for ROI_num = 1:size(unexp_rew_roi,2); %for each roi...
            rew_UR_temp = squeeze(unexp_rew_roi(:,ROI_num,:))';  %dims 1=Time 2=trial number
            rew_UR_temp = interp1(rew_UR_temp, XVector);
            unexp_rew_roi_interp(:,:,ROI_num) = rew_UR_temp;  %dims 1=Time 2=trialNumber  3=roi
        end
        unexp_rew_roi_interp = squeeze(mean(unexp_rew_roi_interp(:,:, ROI_cell{ii}),3));
    end
    
    if strcmp(days{1}, days_1{1});
        save(['Z:\Analysis\Cue_reward_pairing_analysis\WF\onset latencies and peak mags\mean interpolated TCs\day_1_vars_' ,days{ii}], 'rew_roi_interp', 'rew_om_roi_interp');
    elseif strcmp(days{1}, days_post{1});
        save(['Z:\Analysis\Cue_reward_pairing_analysis\WF\onset latencies and peak mags\mean interpolated TCs\day_N_vars_' ,days{ii}], 'rew_roi_interp', 'rew_om_roi_interp');
    elseif strcmp(days{1}, days_UR{1});
        save(['Z:\Analysis\Cue_reward_pairing_analysis\WF\onset latencies and peak mags\mean interpolated TCs\day_N+1_vars_' ,days{ii}], 'rew_roi_interp', 'unexp_rew_roi_interp');
    end
    
    %baseline TCs------------------------------------------
%     for trial_num = 1:size(rew_roi_interp,2);  %dim1=time(10ms)  dim2=trial#
%         shift = mean(rew_roi_interp([200:401], trial_num));
%         rew_roi_interp(:, trial_num) = rew_roi_interp(:,trial_num)-shift;
%     end
%     if exist('rew_om_roi_interp')
%         for trial_num = 1:size(rew_om_roi_interp,2);  %dim1=time(10ms)  dim2=trial#
%             shift = mean(rew_om_roi_interp([200:401], trial_num));
%             rew_om_roi_interp(:, trial_num) = rew_om_roi_interp(:,trial_num)-shift;
%         end
%     elseif exist('unexp_rew_roi_interp')
%         for trial_num = 1:size(unexp_rew_roi_interp,2);  %dim1=time(10ms)  dim2=trial#
%             shift = mean(unexp_rew_roi_interp([800:1001], trial_num));
%             unexp_rew_roi_interp(:, trial_num) = unexp_rew_roi_interp(:,trial_num)-shift;
%         end
%     end
    
    %find peak magnitude and latencey on a trial by trial basis-----------------------
    rew_tbyt_lat_peak = []; 
    rew_tbyt_peak_mag = []; 
    peak_val_window = [700:1701];
    for trial_num = 1:size(rew_roi_interp,2);
        temp_val = find(   rew_roi_interp(peak_val_window,trial_num)==max(rew_roi_interp(peak_val_window, trial_num))) + min(peak_val_window)-1;  %finds the peak time of each trials TC
        rew_tbyt_lat_peak = [rew_tbyt_lat_peak, temp_val]; %give peak time for each trial. Will occur in increments of 10 because the peak will always be on one of the actual frames. Not interpolated data points
        rew_tbyt_peak_mag = [rew_tbyt_peak_mag, rew_roi_interp(temp_val, trial_num)]; %collects all the trial by trial peak magnitudes
    end
    
    if exist('rew_om_roi_interp')
        rew_om_tbyt_lat_peak = [];
        rew_om_tbyt_peak_mag = [];
        for trial_num = 1:size(rew_om_roi_interp,2);
            temp_val = find(rew_om_roi_interp(peak_val_window,trial_num)==max(rew_om_roi_interp(peak_val_window, trial_num))) + min(peak_val_window)-1;  %finds the peak time of each trials TC
            rew_om_tbyt_lat_peak = [rew_om_tbyt_lat_peak, temp_val]; %give peak time for each trial. Will occur in increments of 10 because the peak will always be on one of the actual frames. Not interpolated data points
            rew_om_tbyt_peak_mag = [rew_om_tbyt_peak_mag, rew_om_roi_interp(temp_val, trial_num)]; %collects all the trial by trial peak magnitudes
        end
    elseif exist('unexp_rew_roi_interp')
        unexp_tbyt_lat_peak = [];
        unexp_tbyt_peak_mag = [];
        for trial_num = 1:size(unexp_rew_roi_interp,2);
            temp_val = find(unexp_rew_roi_interp(peak_val_window,trial_num)==max(unexp_rew_roi_interp(peak_val_window, trial_num))) + min(peak_val_window)-1;  %finds the peak time of each trials TC
            unexp_tbyt_lat_peak = [unexp_tbyt_lat_peak, temp_val]; %give peak time for each trial. Will occur in increments of 10 because the peak will always be on one of the actual frames. Not interpolated data points
            unexp_tbyt_peak_mag = [unexp_tbyt_peak_mag, unexp_rew_roi_interp(temp_val, trial_num)]; %collects all the trial by trial peak magnitudes
        end
    end
    
    %find baseline values and indeces for each trial's F response------------------------
    rew_TC_baselines = [];     
    rew_base_ind_vec = [];  % this is a vector which will collect the indexed value of the min value calculated for each trial
    baseline_window_dur = 500;   %baseline value will be the min value within 500ms before the peak 
    for trial_num = 1:size(rew_roi_interp,2)
        %rew_TC_baselines = [rew_TC_baselines, mean(rew_roi_interp(cue_time_interp-baseline_window_dur:cue_time_interp-1, trial_num))]; 
        if rew_tbyt_lat_peak(trial_num)-baseline_window_dur >= cue_time_interp
            rew_TC_baselines = [rew_TC_baselines, min(rew_roi_interp(    [rew_tbyt_lat_peak(trial_num)-baseline_window_dur:rew_tbyt_lat_peak(trial_num)-1],  trial_num))];
        else
            rew_TC_baselines = [rew_TC_baselines, min(rew_roi_interp(    [cue_time_interp:rew_tbyt_lat_peak(trial_num)-1],  trial_num))];
        end
        rew_base_ind = find(rew_roi_interp(cue_time_interp-baseline_window_dur:cue_time_interp-1, trial_num) == rew_TC_baselines(end), 1, 'last') +(cue_time_interp-baseline_window_dur);
        rew_base_ind_vec = [rew_base_ind_vec, rew_base_ind];
    end
    
    if exist('rew_om_roi_interp')
        rew_om_TC_baselines = [];   
        rew_om_base_ind_vec = []; 
        for trial_num = 1:size(rew_om_roi_interp,2)
            if rew_om_tbyt_lat_peak(trial_num)-baseline_window_dur >= cue_time_interp
                rew_om_TC_baselines = [rew_om_TC_baselines, min(rew_om_roi_interp(    [rew_om_tbyt_lat_peak(trial_num)-baseline_window_dur:rew_om_tbyt_lat_peak(trial_num)-1],  trial_num))];
            else
                rew_om_TC_baselines = [rew_om_TC_baselines, min(rew_om_roi_interp(    [cue_time_interp:rew_om_tbyt_lat_peak(trial_num)-1],  trial_num))];
            end
            rew_om_base_ind = find(rew_om_roi_interp([cue_time_interp-baseline_window_dur:cue_time_interp-1], trial_num) == rew_om_TC_baselines(end), 1, 'last') +(cue_time_interp-baseline_window_dur);
            rew_om_base_ind_vec = [rew_om_base_ind_vec, rew_om_base_ind];
        end
    elseif exist('unexp_rew_roi_interp')
        unexp_TC_baselines = []; 
        unexp_base_ind_vec = []; 
        for trial_num = 1:size(unexp_rew_roi_interp,2)
            if unexp_tbyt_lat_peak(trial_num)-baseline_window_dur >= cue_time_interp
                unexp_TC_baselines = [unexp_TC_baselines, min(unexp_rew_roi_interp(    [unexp_tbyt_lat_peak(trial_num)-baseline_window_dur:unexp_tbyt_lat_peak(trial_num)-1],  trial_num))];
            else
                unexp_TC_baselines = [unexp_TC_baselines, min(unexp_rew_roi_interp(    [cue_time_interp:unexp_tbyt_lat_peak(trial_num)-1],  trial_num))];
            end
            unexp_base_ind =  find(unexp_rew_roi_interp(cue_time_interp-baseline_window_dur:cue_time_interp-1, trial_num) == unexp_TC_baselines(end), 1, 'last') +(cue_time_interp-baseline_window_dur);
            unexp_base_ind_vec = [unexp_base_ind_vec, unexp_base_ind];
        end
    end
    
    % find the 10% and 90% values and times. baseline=0% peak=100% need to take into account that some values could cross zero.
    % Code needs to find the difference between the peak and baseline values. Then it will look the left of the peak time and log the 1st time the
    % value hits 90% or less of the peak and the 1st time at which it hits 10% or less of the peak.
    
    %rise time data for rew ect trials -----------------------------------
    rew_peak_to_base_diff = rew_tbyt_peak_mag-rew_TC_baselines; %find the difference between the peak magnitude and the magnitude of the baseline F
    good_diff_ind = find(rew_peak_to_base_diff>0.02); 
    rew_ninety_val = (0.9*rew_peak_to_base_diff)+rew_TC_baselines;
    rew_ten_val    = (0.1*rew_peak_to_base_diff)+rew_TC_baselines;
    rew_ninety_time = [];
    rew_ten_time = [];
    for trial_num = good_diff_ind;
        rew_ninety_time = [rew_ninety_time, find(rew_roi_interp([1:rew_tbyt_lat_peak(trial_num)],trial_num)<=rew_ninety_val(trial_num),1,'last')];  %find the last value from 1 to the time of the peak where the TC is less than the 90% value.
        rew_ten_time    = [rew_ten_time, find(rew_roi_interp([1:rew_tbyt_lat_peak(trial_num)],trial_num)<=rew_ten_val(trial_num),1,'last')];  %find the last value from 1 to the time of the peak where the TC is less than the 10% value.
        assert(rew_ninety_time(end)<=rew_tbyt_lat_peak(trial_num)); % | rew_ninety_time(trial_num)>=rew_base_ind_vec(trial_num));
        assert(rew_ten_time(end)<=rew_tbyt_lat_peak(trial_num)); % | rew_ten_time(trial_num)>=rew_base_ind_vec(trial_num));
    end
    rew_rise_times = rew_ninety_time-rew_ten_time;
    rew_onset_times = rew_ten_time;
    [rew_peak_to_base_diff_mean(end+1), rew_peak_to_base_diff_sem(end+1)] = get_mean_and_sem(rew_peak_to_base_diff);
    
    if exist('rew_om_roi_interp')
        rew_om_peak_to_base_diff = rew_om_tbyt_peak_mag-rew_om_TC_baselines; %find the difference between the peak magnitude and the magnitude of the baseline F
        good_diff_ind = find(rew_om_peak_to_base_diff>0.02); 
        rew_om_ninety_val = (0.9*rew_om_peak_to_base_diff)+rew_om_TC_baselines;
        rew_om_ten_val    = (0.1*rew_om_peak_to_base_diff)+rew_om_TC_baselines;
        rew_om_ninety_time = [];
        rew_om_ten_time = [];
        for trial_num = good_diff_ind;
            rew_om_ninety_time = [rew_om_ninety_time, find(rew_om_roi_interp([1:rew_om_tbyt_lat_peak(trial_num)],trial_num)<=rew_om_ninety_val(trial_num),1,'last')];  %find the last value from 1 to the time of the peak where the TC is less than the 90% value.
            rew_om_ten_time    = [rew_om_ten_time, find(rew_om_roi_interp([1:rew_om_tbyt_lat_peak(trial_num)],trial_num)<=rew_om_ten_val(trial_num),1,'last')];  %find the last value from 1 to the time of the peak where the TC is less than the 10% value.
            assert(rew_om_ninety_time(end)<=rew_om_tbyt_lat_peak(trial_num)); % | rew_om_ninety_time(trial_num)>=rew_om_base_ind_vec(trial_num));
            assert(rew_om_ten_time(end)<=rew_om_tbyt_lat_peak(trial_num)); % | rew_om_ten_time(trial_num)>=rew_om_base_ind_vec(trial_num));
        end
        rew_om_rise_times = rew_om_ninety_time-rew_om_ten_time;
        rew_om_onset_times = rew_om_ten_time;
        [rew_om_peak_to_base_diff_mean(end+1), rew_om_peak_to_base_diff_sem(end+1)] = get_mean_and_sem(rew_om_peak_to_base_diff);
        
    elseif exist('unexp_rew_roi_interp')
        unexp_peak_to_base_diff = unexp_tbyt_peak_mag-unexp_TC_baselines; %find the difference between the peak magnitude and the magnitude of the baseline F
        good_diff_ind = find(unexp_peak_to_base_diff>0.02); 
        unexp_ninety_val = (0.9*unexp_peak_to_base_diff)+unexp_TC_baselines;
        unexp_ten_val    = (0.1*unexp_peak_to_base_diff)+unexp_TC_baselines;
        unexp_ninety_time = [];
        unexp_ten_time = [];
        for trial_num = good_diff_ind;
            unexp_ninety_time = [unexp_ninety_time, find(unexp_rew_roi_interp([1:unexp_tbyt_lat_peak(trial_num)],trial_num)<=unexp_ninety_val(trial_num),1,'last')];  %find the last value from 1 to the time of the peak where the TC is less than the 90% value.
            unexp_ten_time    = [unexp_ten_time, find(unexp_rew_roi_interp([1:unexp_tbyt_lat_peak(trial_num)],trial_num)<=unexp_ten_val(trial_num),1,'last')];  %find the last value from 1 to the time of the peak where the TC is less than the 10% value.
            assert(unexp_ninety_time(end)<=unexp_tbyt_lat_peak(trial_num)); % | unexp_ninety_time(trial_num)>=unexp_base_ind_vec(trial_num));
            assert(unexp_ten_time(end)<=unexp_tbyt_lat_peak(trial_num)); % | unexp_ten_time(trial_num)>=unexp_base_ind_vec(trial_num));
        end
        unexp_rise_times = unexp_ninety_time-unexp_ten_time;
        unexp_onset_times = unexp_ten_time;
        [unexp_peak_to_base_diff_mean(end+1), unexp_peak_to_base_diff_sem(end+1)] = get_mean_and_sem(unexp_peak_to_base_diff);
    end
    
    %store mean and sem across animals ------------------------------------
    if length(rew_onset_times)<= 2
        rew_onset_times_all =  [rew_onset_times_all, NaN];
        rew_onset_times_all_sem = [rew_onset_times_all_sem, NaN];
    else
        rew_onset_times= rew_onset_times-cue_time_interp;
        [rew_onset_times_mean, rew_onset_times_sem] = get_mean_and_sem(rew_onset_times);
        rew_onset_times_all =  [rew_onset_times_all, rew_onset_times_mean];
        rew_onset_times_all_sem = [rew_onset_times_all_sem, rew_onset_times_sem];
    end
    if exist('rew_om_roi_interp') && length(rew_om_onset_times) <= 2
        rew_om_onset_times_all =  [rew_om_onset_times_all, NaN];
        rew_om_onset_times_all_sem = [rew_om_onset_times_all_sem, NaN];
    elseif exist('rew_om_roi_interp') && length(rew_om_onset_times) > 2
        rew_om_onset_times = rew_om_onset_times-cue_time_interp;
        [rew_om_onset_times_mean, rew_om_onset_times_sem] = get_mean_and_sem(rew_om_onset_times);
        rew_om_onset_times_all =  [rew_om_onset_times_all, rew_om_onset_times_mean];
        rew_om_onset_times_all_sem = [rew_om_onset_times_all_sem, rew_om_onset_times_sem];
    end
    if exist('unexp_rew_roi_interp') && length(unexp_onset_times) <= 2
        unexp_onset_times_all =  [unexp_onset_times_all, NaN];
        unexp_onset_times_all_sem = [unexp_onset_times_all_sem, NaN];
    elseif exist('unexp_rew_roi_interp') &&length(unexp_onset_times) > 2
        unexp_onset_times= unexp_onset_times-cue_time_interp;
        [unexp_onset_times_mean, unexp_onset_times_sem] = get_mean_and_sem(unexp_onset_times);
        unexp_onset_times_all =  [unexp_onset_times_all, unexp_onset_times_mean];
        unexp_onset_times_all_sem = [unexp_onset_times_all_sem, unexp_onset_times_sem];
    end
end

if strcmp(days{1}, days_1{1});
    disp(['DAY 1: rewarded trials df/f onset time =', num2str(nanmean(rew_onset_times_all)), '  rewarded trials df/f onset sem =' num2str(nanstd(rew_onset_times_all)/sqrt(sum(~isnan(rew_onset_times_all))))]);
    disp(['DAY 1: omission trials df/f onset time =', num2str(nanmean(rew_om_onset_times_all)), '  omission trials df/f onset sem =' num2str(nanstd(rew_om_onset_times_all)/sqrt(length(rew_om_onset_times_all)))]);
elseif strcmp(days{1}, days_post{1})
    disp(['DAY N: rewarded trials df/f onset time =', num2str(nanmean(rew_onset_times_all)), '  rewarded trials df/f onset sem =' num2str(nanstd(rew_onset_times_all)/sqrt(sum(~isnan(rew_onset_times_all))))]);
    disp(['DAY N: omission trials df/f onset time =', num2str(nanmean(rew_om_onset_times_all)), '  omission trials df/f onset sem =' num2str(nanstd(rew_om_onset_times_all)/sqrt(sum(~isnan(rew_om_onset_times_all))))]);
elseif strcmp(days{1}, days_UR{1})
    disp(['DAY N+1: rewarded trials df/f onset time =', num2str(nanmean(rew_onset_times_all)), '  rewarded trials df/f onset sem =' num2str(nanstd(rew_onset_times_all)/sqrt(sum(~isnan(rew_onset_times_all))))]);
    disp(['DAY N+1: unexpected trials df/f onset time =', num2str(nanmean(unexp_onset_times_all)), '  unexpected trials df/f onset sem =' num2str(nanstd(unexp_onset_times_all)/sqrt(sum(~isnan(unexp_onset_times_all))))]);
end

if strcmp(days{1}, days_post{1})
    figure; hold on;
    for ii = 1:length(rew_onset_times_all)
        if ~isnan(rew_onset_times_all(ii)) & ~isnan(rew_om_onset_times_all(ii))
            %errorbar(rew_onset_times_all(ii), rew_om_onset_times_all(ii), rew_om_onset_times_all_sem(ii), 'k');
            scatter(rew_onset_times_all(ii), rew_om_onset_times_all(ii), 'k');
            plot([rew_onset_times_all(ii), rew_onset_times_all(ii)], [rew_om_onset_times_all(ii)-rew_om_onset_times_all_sem(ii),rew_om_onset_times_all(ii)+rew_om_onset_times_all_sem(ii)], 'k')
            plot([ rew_onset_times_all(ii)-rew_onset_times_all_sem(ii), rew_onset_times_all(ii)+rew_onset_times_all_sem(ii) ], [rew_om_onset_times_all(ii),rew_om_onset_times_all(ii)], 'k')
        end
    end
    yy = [1:100:1400]; plot(yy,yy, 'k');
    xlim([0 1200]); ylim([0 1200]);
    title('day N onset latency WF'); xlabel('rewarded trials (ms relative to cue onset)'); ylabel('omission trials (ms relative to cue onset)');
elseif strcmp(days{1}, days_UR{1});
    figure; hold on;
    for ii = 1:length(rew_onset_times_all)
        if ~isnan(rew_onset_times_all(ii)) & ~isnan(unexp_onset_times_all(ii))
            scatter(rew_onset_times_all(ii), unexp_onset_times_all(ii), 'k');
            plot([rew_onset_times_all(ii), rew_onset_times_all(ii)], [unexp_onset_times_all(ii)-unexp_onset_times_all_sem(ii), unexp_onset_times_all(ii)+unexp_onset_times_all_sem(ii)], 'k');
            plot([ rew_onset_times_all(ii)-rew_onset_times_all_sem(ii), rew_onset_times_all(ii)+rew_onset_times_all_sem(ii) ], [unexp_onset_times_all(ii),unexp_onset_times_all(ii)], 'k')
        end
    end
    yy = [1:100:1400]; plot(yy,yy, 'k');
    xlim([0 1200]); ylim([0 1200]);
    title('day N+1 onset latency WF'); xlabel('rewarded trials (ms relative to cue onset)'); ylabel('unexpected trials (ms relative to cue onset)');
end

if strcmp(days{1}, days_post{1})
    save('Z:\Analysis\Cue_reward_pairing_analysis\WF\onset latencies and peak mags\dayN_variables', 'rew_om_peak_to_base_diff_mean', 'rew_om_peak_to_base_diff_sem', 'rew_peak_to_base_diff_mean', 'rew_peak_to_base_diff_sem', 'rew_onset_times_all', 'rew_onset_times_all_sem', 'rew_om_onset_times_all', 'rew_om_onset_times_all_sem', 'days')
elseif strcmp(days{1}, days_1{1});
    save('Z:\Analysis\Cue_reward_pairing_analysis\WF\onset latencies and peak mags\day1_variables', 'rew_om_peak_to_base_diff_mean', 'rew_om_peak_to_base_diff_sem', 'rew_peak_to_base_diff_mean', 'rew_peak_to_base_diff_sem', 'rew_onset_times_all', 'rew_onset_times_all_sem', 'rew_om_onset_times_all', 'rew_om_onset_times_all_sem', 'days')
elseif strcmp(days{1}, days_UR{1});
    %save('Z:\Analysis\Cue_reward_pairing_analysis\WF\onset latencies and peak mags\dayUR_variables', 'rew_peak_to_base_diff', 'rew_onset_times_all', 'rew_om_onset_times_all', 'days')
end

%% 
clear; 
dayN_vars = load('Z:\Analysis\Cue_reward_pairing_analysis\WF\onset latencies and peak mags\dayN_variables');
day1_vars = load('Z:\Analysis\Cue_reward_pairing_analysis\WF\onset latencies and peak mags\day1_variables');

%onset latencies for OMISSION  day1 vs dayN 
% figure; subplot(1,2,1);  hold on;
% for ii = 1:length(day1_vars.rew_onset_times_all)
%     if ~isnan(day1_vars.rew_om_onset_times_all(ii)) & ~isnan(dayN_vars.rew_om_onset_times_all(ii))
%         errorbar(dayN_vars.rew_om_onset_times_all(ii), day1_vars.rew_om_onset_times_all(ii), day1_vars.rew_om_onset_times_all_sem(ii));
%         plot( [dayN_vars.rew_om_onset_times_all(ii)-dayN_vars.rew_om_onset_times_all_sem(ii), dayN_vars.rew_om_onset_times_all(ii)+dayN_vars.rew_om_onset_times_all_sem(ii)] ...
%             ,[day1_vars.rew_om_onset_times_all(ii), day1_vars.rew_om_onset_times_all(ii)]);
%     end
% end
% yy = [1:100:1400]; plot(yy,yy);
% plot(yy,yy); xlim([0 1200]); ylim([0 1200]);
% title('omission: df/f onset latencies (ms relative to cue onset) across days'); xlabel('day N'); ylabel('day 1');

%onset latencies for REWARD  day1 vs dayN 
%subplot(1,2,2); 
figure; hold on;
for ii = 1:length(day1_vars.rew_onset_times_all)
    if ~isnan(day1_vars.rew_onset_times_all(ii)) & ~isnan(dayN_vars.rew_onset_times_all(ii))
        %errorbar(dayN_vars.rew_onset_times_all(ii), day1_vars.rew_onset_times_all(ii), day1_vars.rew_onset_times_all_sem(ii));
        scatter(dayN_vars.rew_onset_times_all(ii), day1_vars.rew_onset_times_all(ii), 'k');
        plot([dayN_vars.rew_onset_times_all(ii), dayN_vars.rew_onset_times_all(ii)], ...
            [day1_vars.rew_onset_times_all(ii)-day1_vars.rew_onset_times_all_sem(ii), day1_vars.rew_onset_times_all(ii)+day1_vars.rew_onset_times_all_sem(ii)], 'k');
        plot( [dayN_vars.rew_onset_times_all(ii)-dayN_vars.rew_onset_times_all_sem(ii), dayN_vars.rew_onset_times_all(ii)+dayN_vars.rew_onset_times_all_sem(ii)] ...
            ,[day1_vars.rew_onset_times_all(ii), day1_vars.rew_onset_times_all(ii)], 'k');
    end
end
yy = [1:100:1400]; plot(yy,yy, 'k');
xlim([0 1200]); ylim([0 1200]);
title('Reward: df/f onset latencies (ms relative to cue onset) across days'); xlabel('day N'); ylabel('day 1');

%% plotpeak magnitudes for omission  day1 vs dayN 
% figure; subplot(1,2,1);  hold on;
% for ii = 1:length(day1_vars.rew_onset_times_all)
%     if ~isnan(day1_vars.rew_om_onset_times_all(ii)) & ~isnan(dayN_vars.rew_om_onset_times_all(ii))
%         errorbar(dayN_vars.rew_om_peak_to_base_diff_mean(ii), day1_vars.rew_om_peak_to_base_diff_mean(ii), day1_vars.rew_om_peak_to_base_diff_sem(ii));
%         plot( [dayN_vars.rew_om_peak_to_base_diff_mean(ii)-dayN_vars.rew_om_peak_to_base_diff_sem(ii), dayN_vars.rew_om_peak_to_base_diff_mean(ii)+dayN_vars.rew_om_peak_to_base_diff_sem(ii)] ...
%             ,[day1_vars.rew_om_peak_to_base_diff_mean(ii), day1_vars.rew_om_peak_to_base_diff_mean(ii)]);
%     end
% end
% yy = [-1:.01:1]; %plot(yy,yy);
% plot(yy,yy); 
% xlim([-0.02 .22]); ylim([-0.02 .22]);
% title('omission: df/f peak df/f across days'); xlabel('day N'); ylabel('day 1');
% 
% % plotpeak magnitudes for reward  day1 vs dayN 
% subplot(1,2,2);  hold on;
% for ii = 1:length(day1_vars.rew_onset_times_all)
%     if ~isnan(day1_vars.rew_onset_times_all(ii)) & ~isnan(dayN_vars.rew_onset_times_all(ii))
%         errorbar(dayN_vars.rew_peak_to_base_diff_mean(ii), day1_vars.rew_peak_to_base_diff_mean(ii), day1_vars.rew_peak_to_base_diff_sem(ii));
%         plot( [dayN_vars.rew_peak_to_base_diff_mean(ii)-dayN_vars.rew_peak_to_base_diff_sem(ii), dayN_vars.rew_peak_to_base_diff_mean(ii)+dayN_vars.rew_peak_to_base_diff_sem(ii)] ...
%             ,[day1_vars.rew_peak_to_base_diff_mean(ii), day1_vars.rew_peak_to_base_diff_mean(ii)]);
%     end
% end
% yy = [-1:.01:1]; %plot(yy,yy);
% plot(yy,yy); 
% xlim([-0.02 .22]); ylim([-0.02 .22]);
% title('rewarded: df/f peak df/f across days'); xlabel('day N'); ylabel('day 1');








