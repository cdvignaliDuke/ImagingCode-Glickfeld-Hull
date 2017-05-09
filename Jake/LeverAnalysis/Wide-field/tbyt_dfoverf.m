function tc_dfoverf = tbyt_dfoverf(b_data,bx_out_dir);
%Uses the trial by trial baseline times to create a df/f timecourse

%load variables from bx_outputs
load(bx_out_dir, 'lever', 'frame_info', 'data_tc');

baseline_timesMs = lever.baseline_timesMs;
first_baseline = find(~isnan(baseline_timesMs(1,:)),1, 'first');    %find the first trial / baseline_timesMs window that is not NaN
StartT = frame_info.imaging_start_MW_T; %time of imaging onset in MWorks time. 

tc_dfoverf = nan(size(data_tc));
F_range = [];
for iT=frame_info.f_frame_trial_num+1: frame_info.l_frame_trial_num-1;    %only looks at the first and last fully imaged trials
    %finding F_range
    %F_range is the # of each frame which we will use to generate f. It uses an iti based f
    if ~isnan(baseline_timesMs(1,iT));   %if there is a valid baseline interval then make that the new F
        F_range = frame_info.counter(baseline_timesMs(1,iT)):frame_info.counter(baseline_timesMs(2,iT));
    elseif isempty(F_range)   %if these are trials before there was ever a valid F_range then use the first valid F_range as this trials F_range also.
        F_range = frame_info.counter(baseline_timesMs(1,first_baseline)):frame_info.counter(baseline_timesMs(2,first_baseline));
    end    %if there was no valid F_range but there was previously a valid F_range then F_range will remain unchanged and the most recent one will be use.
    
    %Find t_range (the time interval encompasing the whole trial. Used to find the frames to which we will apply this df/f).
    if iT == frame_info.f_frame_trial_num+1; %if this is the first fully imaged trial then t_range includes all frames up to this point
        t_range = 1:(frame_info.counter(round(cell2mat(b_data.tThisTrialStartTimeMs(iT+1)))-StartT)-1);
    else
        t_range = frame_info.counter(round(cell2mat(b_data.tThisTrialStartTimeMs(iT))-StartT)):(frame_info.counter(round(cell2mat(b_data.tThisTrialStartTimeMs(iT+1))-StartT))-1);
    end
    if iT == frame_info.l_frame_trial_num-1;
        t_range = (t_range(1)+4):length(data_tc);
    elseif iT == frame_info.f_frame_trial_num+1;
        t_range = 1:(t_range(end)+4);
    else
        t_range = t_range + 4;    %added this shift because we have a 1s anaylsis window post release but trial ends 600ms after release.
    end
    
    %assign the frame numbers which correspond to this trial
    for i = 1:size(data_tc,1); %find the avf_f for F_range and apply it to all the frames in that trial
        F_avg= mean(data_tc(i,F_range));
        t_df = bsxfun(@minus, data_tc(i, t_range), F_avg);   %use bsxfun because sometimes F_avg is a vector
        t_dfoverf = bsxfun(@rdivide, t_df, F_avg);
        tc_dfoverf(i,t_range) = t_dfoverf;
    end
end

