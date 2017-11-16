%look for sig diff between peaks of early trials with and without licks
clear;

WF_plotting_lists_of_days;
baseline_window = [1:3];
peak_window = [5:7];
fail_TC_dir = ['Z:\Analysis\WF Lever Analysis\licking_investigation\early_trials_lick_v_no_lick\']; 

for session_num = 1:length(days)
    %load TCs
    if exist([fail_TC_dir, days{session_num}, '.mat'])
        load([fail_TC_dir, days{session_num}, '.mat']);
    else
        continue;
    end
    %only look at LS
    no_lick_TCs = early_TCs.no_lick_TCs(:,[LS_ROIs{session_num}],:);
    lick_TCs = early_TCs.lick_TCs(:,[LS_ROIs{session_num}],:);
    
    %only include sessions with 4 or more trials in each condition
    if size(no_lick_TCs,1) < 4 
        continue
    elseif size(lick_TCs,1) < 4 
        continue
    end
    
    %baseline all trials to the first 3 frames
    shift = mean(no_lick_TCs(:,:,baseline_window),3);
    shift = repmat(shift, 1,1,size(no_lick_TCs,3));
    no_lick_TCs = no_lick_TCs - shift;
    
    shift = mean(lick_TCs(:,:,baseline_window),3);
    shift = repmat(shift, 1,1,size(lick_TCs,3));
    lick_TCs = lick_TCs - shift;
    
    %find the peak value for each trial/ROI 
    peak_vals_no_lick = max(no_lick_TCs(:,:,peak_window),[],3);
    peak_vals_lick = max(lick_TCs(:,:,peak_window), [], 3);
    
    %testing the null hypothesis that trials with licks have a larger peak
    %0 failed to reject null   1 rejection of null. licking is not greater
    ttest_results_mat = NaN(2,size(no_lick_TCs,2));
    for ROI_num = 1:size(no_lick_TCs,2)
        [ttest_results_mat(1,ROI_num), ttest_results_mat(2,ROI_num)] = ttest2(peak_vals_lick(:,ROI_num), peak_vals_no_lick(:,ROI_num), 'tail', 'right');
    end
    ttest_results{session_num} = ttest_results_mat;
    
end

%save([fail_TC_dir, 'ttest_results'], 'ttest_results', 'days', 'LS_ROIs', 'V_ROIs', 'C1_ROIs');





