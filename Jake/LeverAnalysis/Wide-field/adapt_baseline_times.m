function baseline_timesMs = adapt_baseline_times(baseline_timesMs, imaging_start_MW_T, l_frame_MWorks_time, f_frame_trial_num, l_frame_trial_num);
%proceesses baseline_timesMs to remove unimaged trials and subtract off the MW time of imaging start. 
%used in parse_behavior_for_HAD as a part of the WF lever analysis

%remove events from baseline_times which do not have associated frames. JH
for t = 1:length(baseline_timesMs)
    if baseline_timesMs(1,t)<imaging_start_MW_T;
        baseline_timesMs(:,t)=NaN;
    elseif baseline_timesMs(2,t)>l_frame_MWorks_time;
        baseline_timesMs(:,t)=NaN;
    end
end

%convert 1st and last trials' baseline_times to NaNs  
baseline_timesMs(:,f_frame_trial_num)=NaN;    %the first and last trials with camera pulses will be incompletely imaged. Therefore they are unusable baseline_times
baseline_timesMs(:,l_frame_trial_num)=NaN;

%subtract the MWtime of camera start in order to align baseline_timesMS to frame.counter
baseline_timesMs = baseline_timesMs - imaging_start_MW_T;
end