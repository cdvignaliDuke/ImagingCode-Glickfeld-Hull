%WF lever task: measure lick rate post release vs magnitude of the df/f
%response. Testing if lick rate correlates with df/f across trials within
%animals to determine if CFs signal a licking motor signal which should
%occur with a ~1:1 ratio with licks.

clear; 
TC_dir = ['Z:\Analysis\WF Lever Analysis\LeverSummaryFolder\'];
bx_outputs_dir = ['Z:\Analysis\WF Lever Analysis\BxAndAnalysisOutputs\BxOutputs\'];
%struct_dir = 'Z:\Analysis\WF Lever Analysis\StructuresForPlotting\'; 
figure_output = ['Z:\Analysis\WF Lever Analysis\licking_investigation\within_session_lick_rate_vs_f\'];
WF_plotting_lists_of_days;
lick_rate_window = [0:2];
lever_release_frame = 6;
peak_window = [4:9];
early_Scorr{1}=[];
corr_Scorr{1} = [];
mean_mags{1} = [];
plot_colors = {'r', 'b', 'g', 'k', 'c', 'm'};
colors_final_plot = {'k','k', 'b','b', 'g','g', 'c','c', 'm','m', 'r','r', 'y','y', 'k', 'b', 'r'}; 
plot_symbols = {'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 's', 's', 's'};

for session_num = [1,2, 5:length(days)] %img30 has some confound in the licking data
    %load data
    bx_out_dir  = [bx_outputs_dir days{session_num} '_bx_outputs'];
    load(bx_out_dir);
    lick_trace_early = licking_data.lick_trace_fail;
    lick_trace_corr = licking_data.lick_trace_succ;
    load([TC_dir, days{session_num}, '_fail']);   %already baselined
    load([TC_dir, days{session_num}, '_success']);
    
    %remove non-LS ROIs
    fail_roi = fail_roi(:,[LS_ROIs{session_num}],:);
    success_roi = success_roi(:,[LS_ROIs{session_num}],:);
    
    %determine peak magnitude for each ROI  borrowed code from struct maker
    assert(length(size(success_roi))==3);
    corr_tbyt_peak_mag = max(success_roi(:,:,peak_window), [], 3);
    fail_tbyt_peak_mag = max(fail_roi(:,:,peak_window), [], 3);
    
    %ensure the peak magnitudes and lick traces have the same number of trials
    assert(size(corr_tbyt_peak_mag,1) == size(lick_trace_corr, 1));
    assert(size(fail_tbyt_peak_mag,1) == size(lick_trace_early, 1));
    
    %determine lick rate just after lever release 
    lick_rate_corr  = sum(lick_trace_corr(:,[lick_rate_window+lever_release_frame]),2)'  / length(lick_rate_window) * 10;  %converts lick rate to Hz 
    lick_rate_early = sum(lick_trace_early(:,[lick_rate_window+lever_release_frame]),2)' / length(lick_rate_window) * 10;
    
    %remove trials with no licks 
    corr_no_lick_inx = find(lick_rate_corr==0);
    early_no_lick_inx = find(lick_rate_early==0);
    lick_rate_corr([corr_no_lick_inx]) = [];
    lick_rate_early([early_no_lick_inx]) = [];
    corr_tbyt_peak_mag([corr_no_lick_inx],:) = [];
    fail_tbyt_peak_mag([early_no_lick_inx],:) = [];
    
    %OPTIONAL - match lick rates-------
    figure_output = ['Z:\Analysis\WF Lever Analysis\licking_investigation\within_session_lick_rate_vs_f\match_lick_rates\']; 
    plot_mean_mags = 1;
    %calculate lick rate for corr vs early and verify corrLR > early LR
    if length(lick_rate_corr)<6 | length(lick_rate_early)<6 
        continue
    end
    mean_LR_corr = mean(lick_rate_corr);
    mean_LR_early = mean(lick_rate_early);
    if mean_LR_corr>mean_LR_early;
        ind_order = 'descend';
    elseif mean_LR_corr < mean_LR_early;
        ind_order = 'ascend';
    elseif mean_LR_corr == mean_LR_early;
        mean_mags{session_num} = [mean(corr_tbyt_peak_mag); mean(fail_tbyt_peak_mag)];
    end
    %Remove corr trial with highest lick rate
    [sorted_corr_LR, sorted_corr_LR_inx] = sort(lick_rate_corr, ind_order);
    sorted_corr_peaks = corr_tbyt_peak_mag([sorted_corr_LR_inx],:);
    for trial_num = 2:length(lick_rate_corr)
        temp_mean = mean(sorted_corr_LR(trial_num:end));
        prev_temp_mean = mean(sorted_corr_LR(trial_num-1:end));
        if temp_mean == mean_LR_early
            break;
        elseif (temp_mean<mean_LR_early & mean_LR_early<prev_temp_mean) | (temp_mean>mean_LR_early & mean_LR_early>prev_temp_mean)
            if abs(temp_mean - mean_LR_early) >= abs(prev_temp_mean - mean_LR_early)
                trial_num = trial_num-1; %using trial_num as there reference for which trials to remove
            end
            break;
        end
        prev_temp_mean = temp_mean; %set prev temp mean before entering the next iteration of the forloop
    end
    lick_rate_corr = sorted_corr_LR(trial_num:end);
    corr_tbyt_peak_mag = sorted_corr_peaks([trial_num:end],:);
    mean_mags{session_num} = [mean(corr_tbyt_peak_mag); mean(fail_tbyt_peak_mag)];
    
    %----------
    
    %determine spearman's correlation
    if length(lick_rate_corr)<6 | length(lick_rate_early)<6 
        continue
    end
    corr_Scorr_temp = [];
    early_Scorr_temp = [];
    for ROI_num = 1:size(corr_tbyt_peak_mag,2)
        [corr_rho, corr_p] = corr(corr_tbyt_peak_mag(:,ROI_num), lick_rate_corr', 'type', 'Spearman');
        [early_rho, early_p] = corr(fail_tbyt_peak_mag(:,ROI_num), lick_rate_early', 'type', 'Spearman');
        corr_Scorr_temp(:,ROI_num) = [corr_rho; corr_p];
        early_Scorr_temp(:,ROI_num) = [early_rho, early_p];
    end
    corr_Scorr{session_num} = corr_Scorr_temp;
    early_Scorr{session_num} = early_Scorr_temp;
    
    %plot lick rate vs peak magnitude for each trial
    figure;    plot_colors;
    x_max = max(max(corr_tbyt_peak_mag));
    for ROI_num = 1:size(corr_tbyt_peak_mag,2)
        subplot(size(corr_tbyt_peak_mag,2), 1, ROI_num);
        scatter(corr_tbyt_peak_mag(:,ROI_num)', lick_rate_corr, plot_colors{ROI_num}); hold on;
        ylim([0 20]);
        xlabel('peak df/f');
        xlim([0 x_max]);
    end
    ylabel('lick rate (Hz)');
    suptitle([days{session_num}, 'correct trials']);
    savefig([figure_output, days{session_num}, '_corr_', num2str(length(lick_rate_window))]);
    
    figure;
    x_max = max(max(fail_tbyt_peak_mag));
    for ROI_num = 1:size(corr_tbyt_peak_mag,2)
        subplot(size(corr_tbyt_peak_mag,2), 1, ROI_num);
        scatter(fail_tbyt_peak_mag(:,ROI_num), lick_rate_early); hold on;
        xlim([0 x_max]);
        ylim([0 20]);
        ylabel('lick rate (Hz)');
    end
    xlabel('peak df/f'); 
    suptitle([days{session_num}, 'early trials']);
    savefig([figure_output, days{session_num}, '_fail_', num2str(length(lick_rate_window))]);
end

save([figure_output, 'Spearman_correlation_', num2str(length(lick_rate_window))], 'corr_Scorr', 'early_Scorr');

mean_correlation_corr = [];
mean_correlation_fail = [];
all_correlation_corr = [];
all_correlation_fail = [];
for session_num = 1:length(corr_Scorr) %img30 has some confound in the licking data
    if isempty(corr_Scorr{session_num})
        continue
    end
    
    mean_correlation_corr = [mean_correlation_corr, mean(corr_Scorr{session_num}(1,:))];
    mean_correlation_fail = [mean_correlation_fail, mean(early_Scorr{session_num}(1,:))];
    all_correlation_corr = [all_correlation_corr, corr_Scorr{session_num}(1,:)];
    all_correlation_fail = [all_correlation_fail, early_Scorr{session_num}(1,:)];
end

figure;
scatter([ones(1,length(mean_correlation_corr))], mean_correlation_corr); hold on;
scatter(1, mean(mean_correlation_corr),'k', 'filled');
errorbar(1, mean(mean_correlation_corr), std(mean_correlation_corr)/sqrt(length(mean_correlation_corr)), 'k');
scatter([ones(1,length(find(~isnan(mean_correlation_fail))))]+1, mean_correlation_fail(find(~isnan(mean_correlation_fail))), 'r');
scatter(2, nanmean(mean_correlation_fail),'k', 'filled')
errorbar(2, mean(mean_correlation_fail), std(mean_correlation_fail)/sqrt(length(mean_correlation_fail)), 'k');
title(['within session correlations between lick rate and df/f. Window=' num2str(length(lick_rate_window))]);
ylabel('Rho');
xlabel('blue: correct trials      red: early trials')
xlim([-1 4]);
hline(0, 'k');
savefig([figure_output, 'ROI_correlations_', num2str(length(lick_rate_window))]);

if exist('plot_mean_mags')
    mean_mags_mat_corr = [];
    mean_mags_mat_early = [];
    figure;
    for session_num = 1:length(mean_mags)
        if ~isempty(mean_mags{session_num})
            %for ROI_num = 1:size(mean_mags{session_num},2)
            mean_mags_mat_corr = [mean_mags_mat_corr mean(mean_mags{session_num}(1,:),2)];
            mean_mags_mat_early = [mean_mags_mat_early, mean(mean_mags{session_num}(2,:),2)];
            plot([1,2], [mean(mean_mags{session_num}(1,:),2)/mean(mean_mags{session_num}(1,:),2), mean(mean_mags{session_num}(2,:),2)/mean(mean_mags{session_num}(1,:),2)], strcat('-', plot_symbols{session_num}, colors_final_plot{session_num}));   hold on;
            
            %end
        end
    end
    errorbar(1,mean(mean_mags_mat_corr), std(mean_mags_mat_corr)/sqrt(length(mean_mags_mat_corr)), 'k');
    errorbar(2,mean(mean_mags_mat_early), std(mean_mags_mat_early)/sqrt(length(mean_mags_mat_early)), 'k');
    scatter([1,2], [mean(mean_mags_mat_corr), mean(mean_mags_mat_early)], 'filled', 'k');
    [hh,pp] = ttest(mean_mags_mat_corr, mean_mags_mat_early)
    xlabel('corr trials     early trials');
    ylabel('peak magnitude df/f');
    title('comparison of df/f peak magnitudes for correct vs early trials with matched lick rates. By ROI.');
    savefig([figure_output, 'compare_peak_mag_session_base_LR_matching_', num2str(length(lick_rate_window))]);
end



