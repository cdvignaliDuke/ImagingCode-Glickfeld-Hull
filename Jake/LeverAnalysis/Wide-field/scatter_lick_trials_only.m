

clear;
%set directory and days to be analyzed
F_TC_dir    = 'Z:\Analysis\WF Lever Analysis\LeverSummaryFolder\';
lick_TC_dir = 'Z:\Analysis\WF Lever Analysis\BxAndAnalysisOutputs\BxOutputs\'; 
corr_w_lick_dir = ['Z:\Analysis\WF Lever Analysis\licking_investigation\correct_trials_lick_v_no_lick\no_licks_-2_to_2\'];
early_w_lick_dir = ['Z:\Analysis\WF Lever Analysis\licking_investigation\early_trials_lick_v_no_lick\no_licks_-2_to_2\'];
fig_dir = ['Z:\Analysis\WF Lever Analysis\licking_investigation\ROI_criteria_folder\'];

%set variables
WF_plotting_lists_of_days;
lever_release = 6; 
lick_window = [lever_release-2:lever_release+2];
peak_window = [lever_release-2:lever_release+2];
corr_tbyt_peak_mag_mean{1} = [];
fail_tbyt_peak_mag_mean{1} = [];
colors_final_plot = {'k','k', 'b','b', 'g','g', 'c','c', 'm','m', 'r','r', 'y','y', 'k', 'b', 'r'}; 
plot_symbols = {'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 's', 's', 's'};

for session_num = [1,2, 5:length(days)] %img30 has some confound in the licking data
    %load files
    load([F_TC_dir, days{session_num}, '_success']);
    load([F_TC_dir, days{session_num}, '_fail']);
    load([lick_TC_dir, days{session_num}, '_bx_outputs'], 'licking_data');
    load([corr_w_lick_dir, days{session_num}, '_trial_indeces']);
    assert(isequal(no_lick_window,[4:8]));
    load([early_w_lick_dir, days{session_num}, '_trial_indeces']);
    assert(isequal(no_lick_window,[4:8]));
    lick_trace_succ = licking_data.lick_trace_succ;
    lick_trace_fail = licking_data.lick_trace_fail;
    
    %get rid of non-LS ROIs 
    success_roi = success_roi(:,[LS_ROIs{session_num}],:);
    fail_roi = fail_roi(:,[LS_ROIs{session_num}],:);
    
    %remove trials without licks
    if length(lick_corr_trials) >5 & length(lick_fail_trials) > 5
        success_roi = success_roi([lick_corr_trials],:,:);
        fail_roi = fail_roi([lick_fail_trials],:,:);
        lick_trace_succ = lick_trace_succ([lick_corr_trials],:);
        lick_trace_fail = lick_trace_fail([lick_fail_trials],:);
    else
        continue
    end
    
    %find peak mags for remaining trials
    corr_tbyt_peak_mag = max(success_roi(:,:,peak_window), [], 3);
    fail_tbyt_peak_mag = max(fail_roi(:,:,peak_window), [], 3);
    
    %get mean peak mag for each ROI and store these values in a cell
    corr_tbyt_peak_mag_mean{session_num} = squeeze(mean(corr_tbyt_peak_mag,1));
    fail_tbyt_peak_mag_mean{session_num} = squeeze(mean(fail_tbyt_peak_mag,1));
    
end

%plot each ROI corr vs early
figure; 
for session_num = [1,2, 5:length(days)] 
    if ~isempty(corr_tbyt_peak_mag_mean{session_num})
        plot([corr_tbyt_peak_mag_mean{session_num}],[fail_tbyt_peak_mag_mean{session_num}], strcat( plot_symbols{session_num}, colors_final_plot{session_num})); hold on;
    end
end
x = [0:.1:1]; y =x;
plot(x,y);
xlabel('correct trials');
ylabel('early trials');
title('trial by trial peak df/f for only trials with licks');
%savefig([fig_dir, 'main_scatter_only_trials_with_licks_ROI']);

%select the ROIs with a minimum ratio of peak df/f corr:early
peak_mag_ratios{1} = [];
min_ratio = 1.15;
valid_ROIs{1} = [];
for session_num = [1,2, 5:length(days)] 
    if ~isempty(corr_tbyt_peak_mag_mean{session_num})
        peak_mag_ratios{session_num} = [corr_tbyt_peak_mag_mean{session_num}]./[fail_tbyt_peak_mag_mean{session_num}];
        valid_ROIs{session_num} = find(peak_mag_ratios{session_num}>min_ratio); %(6 animals 6 sessions 11 ROIs > 1.15)   (7 animals 7 sessions 14ROIs >1.10)
    end
end
save([fig_dir, 'corr_early_ratio_trials_with_licks_ROI'], 'peak_mag_ratios', 'min_ratio', 'valid_ROIs');

%determine number of valid animals
min_ratio_mat = [];
for session_num = 1:length(valid_ROIs)
    if ~isempty(valid_ROIs{session_num})
        min_ratio_mat = [min_ratio_mat, 1];
    else
        min_ratio_mat = [min_ratio_mat, 0];
    end
end
days{find(min_ratio_mat)}



