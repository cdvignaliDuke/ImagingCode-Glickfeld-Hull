% ADAPTED from WF_lever_plotting_TCs and fail_lick_vs_no_lick to plot the
% TCs and licking of only LS ROIs and only trials which have a lick in the
% analysis window or do not have a lick in the window. 

%PLOTTING GRAPHS 
%plots correlation coefficient between the ROIs
%plots TCs of the df/f for success, earlies, tooFast, fidget and press.  
clear
BEHAVE_DIR = 'Z:\Data\WidefieldImaging\GCaMP\behavior\';
ANALYSIS_DIR ='Z:\Analysis\WF Lever Analysis\';
CLUSTER_DIR  ='Z:\Analysis\WF Lever Analysis\BxAndAnalysisOutputs\'; 
TC_dir = 'Z:\Analysis\WF Lever Analysis\LeverSummaryFolder\';
WF_plotting_lists_of_days;

%define variables to help with plotting
func = @mean; %func = @std; @median
pre_frames = 5;
post_frames = 10;
sampling_rate = 10;
tot_frame = pre_frames + post_frames+1;
colors = [1,0,0; 0,1,0; 0,0,1; 0.5,0.5,0.5; 1,0,1; 1,1,0; 0,1,1]; %sets up the color scheme for plotting multiple ROIs with errorbar

for session_num=1:length(days)
    %load files
    bx_out_dir  = [CLUSTER_DIR 'BxOutputs\' days{session_num} '_bx_outputs'];
    load(bx_out_dir);
    load([TC_dir, days{session_num}, '_success']);
    load([TC_dir, days{session_num}, '_fail']);
    lick_trace_succ = licking_data.lick_trace_succ;
    lick_trace_fail = licking_data.lick_trace_fail;
    ts = (-pre_frames:post_frames)*1000/double(round(sampling_rate));
    ts = repmat(ts,[length(LS_ROIs{session_num}) 1]);
    
    %add in a condition so it only plots the datasets with >=4 trials in correct and earlies
    load(['Z:\Analysis\WF Lever Analysis\licking_investigation\early_trials_lick_v_no_lick\no_licks_-2_to_2\', days{session_num}, '_trial_indeces'], 'lick_fail_trials', 'no_lick_fail_trials');
    early_trials_no_lick_window = load(['Z:\Analysis\WF Lever Analysis\licking_investigation\early_trials_lick_v_no_lick\no_licks_-2_to_2\', days{session_num}, '_trial_indeces'], 'no_lick_window');
    load(['Z:\Analysis\WF Lever Analysis\licking_investigation\correct_trials_lick_v_no_lick\no_licks_-2_to_2\', days{session_num}, '_trial_indeces'], 'lick_corr_trials', 'no_lick_corr_trials');
    corr_trials_no_lick_window =  load(['Z:\Analysis\WF Lever Analysis\licking_investigation\correct_trials_lick_v_no_lick\no_licks_-2_to_2\', days{session_num}, '_trial_indeces'], 'no_lick_window');
    assert(isequal(early_trials_no_lick_window.no_lick_window, corr_trials_no_lick_window.no_lick_window)); %assert that both data sets come from the same analysis window
    
    %determine if this session has at least 4 corr/early trials with licks        (or alternatively without licks)
    %     if length(no_lick_corr_trials) < 4 | length(no_lick_fail_trials) < 4 %modify days, LS_ROI, colors, summary_succ/fail
    %         continue
    %     else  %remove the no lick trials and non LS ROIs
    %         success_roi = success_roi([no_lick_corr_trials], LS_ROIs{session_num}, :);
    %         lick_trace_succ = lick_trace_succ([no_lick_corr_trials], :);
    %         fail_roi = fail_roi([no_lick_fail_trials], LS_ROIs{session_num}, :);
    %         lick_trace_fail = lick_trace_fail([no_lick_fail_trials], :);
    %     end
    if length(lick_corr_trials) < 4 | length(lick_fail_trials) < 4 %modify days, LS_ROI, colors, summary_succ/fail
        continue
    else  %remove the no lick trials and non LS ROIs
        success_roi = success_roi([lick_corr_trials], LS_ROIs{session_num}, :);
        lick_trace_succ = lick_trace_succ([lick_corr_trials], :);
        fail_roi = fail_roi([lick_fail_trials], LS_ROIs{session_num}, :);
        lick_trace_fail = lick_trace_fail([lick_fail_trials], :);
    end
    
    % PLOT SUCCESSFUL TRIALS----------------------------------------------
    avg_success_roi = squeeze(func(success_roi,1));
    if length(LS_ROIs{session_num}) == 1
        avg_success_roi = avg_success_roi';
    end
    std_success = squeeze(std(squeeze(success_roi),1));
    sm_success = std_success./sqrt(size(success_roi,1));
    for ROI_num = 1:length(LS_ROIs{session_num})  %baseline each curve so it passes through zero
        shift = (-1)*mean(avg_success_roi(ROI_num,[1:2]),2);
        avg_success_roi(ROI_num,:) = avg_success_roi(ROI_num,:)+shift;
    end
    figure; suptitle([days{session_num}, ' LS ROIs only. Only trials with a lick -2:2']);
    subplot(1,2,1); lickBars = bar(ts(1,:), mean(lick_trace_succ)/10); hold on
    alpha(.25);
    for i = 1:size(ts,1);
        subplot(1,2,1); errorbar(ts(i,:), avg_success_roi(i,:), sm_success(i,:), 'Color', colors(i,:)); hold on;
    end
    xlabel('Time from release (ms)');
    ylabel('dF/F');
    title(['success n=', num2str(size(success_roi,1))]);
    axis tight;
    hold off
    
    %PLOT FAILED TRIALS---------------------------------------------------
    avg_fail_roi = squeeze(func(fail_roi,1));
    if length(LS_ROIs{session_num}) == 1
        avg_fail_roi = avg_fail_roi';
    end
    std_fail = squeeze(std(squeeze(fail_roi),1));
    sm_fail = std_fail./sqrt(size(fail_roi,1));
    for ROI_num = 1:length(LS_ROIs{session_num})  %baseline each curve so it passes through zero
        shift = mean((-1)*avg_fail_roi(ROI_num,[1:2]),2);
        avg_fail_roi(ROI_num,:) = avg_fail_roi(ROI_num,:)+shift;
    end
    subplot(1,2,2); lickBars = bar(ts(1,:), mean(lick_trace_fail)/10); hold on
    alpha(.25);
    for i = 1:size(ts,1);
        subplot(1,2,2); errorbar(ts(i,:), avg_fail_roi(i,:), sm_fail(i,:), 'Color', colors(i,:)); hold on;
    end
    xlabel('Time from release (ms)');
    ylabel('dF/F');
    title(['fail n=', num2str(size(fail_roi,1))]);
    axis tight;
    hold off
    
    %set y scale to be equal for all 3 subplots----
    YL = [];
    for i =1:2
        subplot(1,2,i); YL(i,:) = ylim;
    end
    YL(YL ==1) = 0; %if there are zero trials in a catagory then it will default to values of -1 and 1
    YL(YL ==-1)= 0;
    
    %STANDARDIZE YLIMS, SAVE VARIABLES, REPORT Ns---------------------------------
    subplot(1,2,1); ylim([min(YL(:,1)) max(YL(:,2))]);
    subplot(1,2,2); ylim([min(YL(:,1)) max(YL(:,2))]);

    %save the figure
    savefig(['Z:\Analysis\WF Lever Analysis\licking_investigation\TCs_with_licks\', days{session_num}]);
end




