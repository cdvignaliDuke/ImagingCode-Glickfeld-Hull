%script for plotting lick rate vs hold time
%pulls data from bxOutputs variables licking_data and trial_outcome

clear; 
WF_plotting_lists_of_days;
bxOutputs_dir = 'Z:\Analysis\WF Lever Analysis\BxAndAnalysisOutputs\BxOutputs\';
output_dir = 'Z:\Analysis\WF Lever Analysis\licking_investigation\lick_rate_vs_hold_duration\';
lick_rate_window = [-2:2];
lever_release_frame = 6;

for session_num = 1:length(days)
    %load data
    load([bxOutputs_dir, days{session_num}, '_bx_outputs'], 'licking_data', 'trial_outcome');
    succ_hold_dur = trial_outcome.succ_hold_dur;
    fail_hold_dur = trial_outcome.fail_hold_dur;
    lick_trace_early = licking_data.lick_trace_fail;
    lick_trace_corr = licking_data.lick_trace_succ;
    
    %determine lick rate just after lever release 
    lick_rate_corr = sum(lick_trace_corr(:,[lick_rate_window+lever_release_frame]),2)' / length(lick_rate_window) * 10;  %converts lick rate to Hz 
    lick_rate_early = sum(lick_trace_early(:,[lick_rate_window+lever_release_frame]),2)' / length(lick_rate_window) * 10;
    
    %plot lick rate against hold duration
    figure;
    subplot(2,1,1);
    scatter(succ_hold_dur, lick_rate_corr, 'g');
    title(['correct', ]);
    xlabel('hold duration (ms)');
    ylabel('lick rate (Hz)');
    xlim([0 6000]);
    ylim([0 20]);
    
    subplot(2,1,2);
    scatter(fail_hold_dur, lick_rate_early, 'r');
    title('early');
    xlabel('hold duration (ms)');
    ylabel('lick rate (Hz)');
    xlim([0 6000]);
    ylim([0 20]);
    
    %title and save
    suptitle([[days{session_num}, 'lick rate vs hold duration. '], ['Rate window = frames ', num2str(lick_rate_window), ' relative to release']]);
    savefig([output_dir, days{session_num}, '_-2_2']);
     
end


