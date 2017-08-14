function rew_no_lick_inx = remove_trials_with_licks(lick_trace_rew, no_lick_window, date_mouse);

rew_no_lick_inx = [];
for trial_num = 1:size(lick_trace_rew,1);
    if sum(lick_trace_rew(trial_num, no_lick_window)) == 0
        rew_no_lick_inx = [rew_no_lick_inx, trial_num];
    end
end
    


