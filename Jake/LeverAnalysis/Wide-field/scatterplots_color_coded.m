%scatterplots color coded by lick triggered average magnitude relative to
%diff in corr-early df/f response 
%or color coded according to the ratio of early trials w/ vs w/o licks
%For LM
%or color coded accoring to the diff of peak df/f for early trials with licks - early trials without licks
%or diff in corr - earlies with licks 

% of failed trials with a lick vs lick triggered average magnitude 
clear;

%define some variables and pathnames
WF_plotting_lists_of_days;
early_TC_dir = 'Z:\Analysis\WF Lever Analysis\licking_investigation\early_trials_lick_v_no_lick\no_licks_-2_to_2\';
corr_TC_dir = 'Z:\Analysis\WF Lever Analysis\licking_investigation\correct_trials_lick_v_no_lick\no_licks_-2_to_2\';
corr_early_TCs_dir = ['Z:\Analysis\WF Lever Analysis\licking_investigation\corr_early_diff_by_ROI\'];
output_dir = 'Z:\Analysis\WF Lever Analysis\licking_investigation\percent licks on earlies scatterplots\';
LTA_dir = 'Z:\Analysis\WF Lever Analysis\licking_investigation\lick_triggered_averages\lick_trig_mags_bout';
colors = {'k','k', 'b','b', 'g','g', 'c','c', 'm','m', 'r','r', 'y','y', 'k', 'b', 'r'}; 
plot_symbols = {'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 's', 's', 's'};
baseline_window = [1:20];
peak_window = [51:110];
min_ratio = 0.6;

%generate some empty cell arrays to be filled in
percent_early_with_lick{1} = [];
early_with_lick_mags{1} = [];
early_no_lick_mags{1} = [];

%load LTA data
load(LTA_dir);

%load tc data and determine differences
for session_num = [[5:length(days)]]
    %load early trials w/ vs w/o licks
    load([early_TC_dir, days{session_num}, '_trial_indeces']);
    load([corr_TC_dir, days{session_num}, '_trial_indeces']);
    load([corr_early_TCs_dir, days{session_num}, '_all']);
    
    %baseline individual TCs
    shift = mean(success_roi_interp(baseline_window,:,:),1);
    shift = repmat(shift, size(success_roi_interp,1),1,1);
    success_roi_interp = success_roi_interp-shift;
    shift = mean(fail_roi_interp(baseline_window,:,:),1);
    shift = repmat(shift, size(fail_roi_interp,1),1,1);
    fail_roi_interp = fail_roi_interp-shift;
    
    %find peak of each trial for all ROIs 51:110
    corr_peak_mags = squeeze(max(success_roi_interp(peak_window,:,:),[],1));
    early_peak_mags = squeeze(max(fail_roi_interp(peak_window,:,:),[],1));
    
    %get peak df/f for early 
    corr_peak_mags_mean{session_num} = mean(corr_peak_mags([lick_corr_trials],:),1);
    early_peak_mags_mean{session_num} = mean(early_peak_mags([lick_fail_trials],:),1);
    corr_peak_mags_sem{session_num} = std(corr_peak_mags([lick_corr_trials],:),[],1)/ sqrt(length(lick_corr_trials));
    early_peak_mags_sem{session_num} = std(early_peak_mags([lick_fail_trials],:),[],1)/ sqrt(length(lick_fail_trials));
    
    %Find the ratio of the peak magnitude for earlies no-licks/licks
    if length(no_lick_fail_trials) >=4 &  length(lick_fail_trials) >=4
        early_by_lick_ratio{session_num} = mean(early_peak_mags([no_lick_fail_trials], :),1) ./ mean(early_peak_mags([lick_fail_trials], :),1);
    else 
        early_by_lick_ratio{session_num} = [];
    end
end

%plot the peak df/f for corr vs early. But Select ROIs based on the ratio
%of the peak df/f for early trials no-lick/lick
figure; subplot(1,2,2);
for session_num = [[5:length(days)]]
    if ~isempty(early_by_lick_ratio{session_num})
        for ROI_num = 1:length(corr_peak_mags_mean{session_num})
            if early_by_lick_ratio{session_num}(ROI_num) >= min_ratio
                if ismember(ROI_num, LS_ROIs{session_num})
                    plot(corr_peak_mags_mean{session_num}(ROI_num), early_peak_mags_mean{session_num}(ROI_num), strcat(plot_symbols{session_num}, colors{session_num})); hold on; 
                end
            end
        end
    end
end
diag_line = [0:0.1:1];
plot(diag_line, diag_line, '-k');
ylim([min([corr_peak_mags_mean{:}, early_peak_mags_mean{:}])-0.02 max([corr_peak_mags_mean{:}, early_peak_mags_mean{:}])+0.02]);
xlim([min([corr_peak_mags_mean{:}, early_peak_mags_mean{:}])-0.02 max([corr_peak_mags_mean{:}, early_peak_mags_mean{:}])+0.02]);
xlabel('correct trials with licks');
ylabel('early trials with licks');
title({'peak df/f for correct vs early trials'; ['only for LS ROIs where']; ['peak df/f for early trials without licks / with licks >= ', num2str(min_ratio)]});









