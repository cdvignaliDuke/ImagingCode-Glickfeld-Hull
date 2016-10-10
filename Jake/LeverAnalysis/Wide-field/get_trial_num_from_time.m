function event_times_trial_num = get_trial_num_from_time(trial_start_times, event_times);
%small function to use event_times to find the trial number those events
%correspond to. 
    event_times_trial_num = [];
    for iii = 1:length(event_times)
        trial_num = find(trial_start_times<event_times(iii),1,'last');  %finds the final trial start time which happened before the time of the succesful lever release. This IDs the trial number of that lever release.
        event_times_trial_num  = [event_times_trial_num, trial_num];
    end
    assert(size(event_times_trial_num,2)==size(event_times,2));
end
