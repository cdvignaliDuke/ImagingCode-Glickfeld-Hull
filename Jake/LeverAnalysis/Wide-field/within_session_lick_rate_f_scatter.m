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
lick_rate_window = [0:4];
lever_release_frame = 6;
peak_window = [5:10];
baseline_window = [1:3];
early_Scorr{1}=[];
corr_Scorr{1} = [];
mean_mags{1} = [];
lick_rate_corr_mat = [];
lick_rate_early_mat = [];
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
    
    %baseline each trial's TC
    for ROI_num = 1:size(success_roi,2);
        for trial_num=1:size(success_roi,1)
            shift = squeeze(mean(success_roi(trial_num, ROI_num, baseline_window),3));
            success_roi(trial_num, ROI_num, :) = success_roi(trial_num, ROI_num, :)-shift;
        end
        for trial_num=1:size(fail_roi,1)
            shift = squeeze(mean(fail_roi(trial_num, ROI_num, baseline_window),3));
            fail_roi(trial_num, ROI_num, :) = fail_roi(trial_num, ROI_num, :)-shift;
        end
    end
    
    %determine peak magnitude for each ROI  borrowed code from struct maker
    assert(length(size(success_roi))==3);
    peak_mag_window = 0; %set to zero to get a single value. Otherwise it well get the mean of peak_ind-peak_mag_window : peak_ind+peak_mag_window
    [~, corr_tbyt_peak_mag] = get_peak_lat_and_mag(success_roi, peak_window, lever_release_frame, peak_mag_window);
    [~, fail_tbyt_peak_mag] = get_peak_lat_and_mag(fail_roi, peak_window, lever_release_frame, peak_mag_window);
    
    %ensure the peak magnitudes and lick traces have the same number of trials
    corr_tbyt_peak_mag = corr_tbyt_peak_mag';
    fail_tbyt_peak_mag = fail_tbyt_peak_mag';
    assert(size(corr_tbyt_peak_mag,1) == size(lick_trace_corr, 1));
    assert(size(fail_tbyt_peak_mag,1) == size(lick_trace_early, 1));
    
    %determine lick rate just after lever release 
    lick_rate_corr  = sum(lick_trace_corr(:,[lick_rate_window+lever_release_frame]),2)'  / length(lick_rate_window) * 10;  %converts lick rate to Hz 
    lick_rate_early = sum(lick_trace_early(:,[lick_rate_window+lever_release_frame]),2)' / length(lick_rate_window) * 10;
    
    %remove trials with no licks 
    corr_no_lick_inx = find(lick_rate_corr==0);
    early_no_lick_inx = find(lick_rate_early==0);
    lick_rate_corr([corr_no_lick_inx]) = [];
    corr_tbyt_peak_mag([corr_no_lick_inx],:) = [];
    lick_rate_early([early_no_lick_inx]) = [];
    fail_tbyt_peak_mag([early_no_lick_inx],:) = [];
    
    %plot individual TCs
%     figure; 
%     subplot(1,2,1); 
%     for ROI_num = 1:size(success_roi,2)
%         for trial_num = 1:size(success_roi,1)
%             if ismember(trial_num, corr_no_lick_inx)
%                 continue
%             end
%             plot(squeeze(success_roi(trial_num, ROI_num, :)), plot_colors{ROI_num}); hold on;
%         end
%     end
%     title('corrects');
%     subplot(1,2,2); 
%     for ROI_num = 1:size(fail_roi,2)
%         for trial_num = 1:size(fail_roi,1)
%             if ismember(trial_num, early_no_lick_inx)
%                 continue
%             end
%             plot(squeeze(fail_roi(trial_num, ROI_num, :)), plot_colors{ROI_num}); hold on;
%         end
%     end
%     title('earlies');
%      suptitle(days{session_num});
    
    %OPTIONAL - match lick rates-------
%     figure_output = ['Z:\Analysis\WF Lever Analysis\licking_investigation\within_session_lick_rate_vs_f\match_lick_rates\']; 
%     plot_mean_mags = 1;
%     %calculate lick rate for corr vs early and verify corrLR > early LR
%     if length(lick_rate_corr)<6 | length(lick_rate_early)<6 
%         continue
%     end
%     mean_LR_corr = mean(lick_rate_corr);
%     mean_LR_early = mean(lick_rate_early);
%     if mean_LR_corr>mean_LR_early;
%         ind_order = 'descend';
%     elseif mean_LR_corr < mean_LR_early;
%         ind_order = 'ascend';
%     end
%     if mean_LR_corr == mean_LR_early;
%         mean_mags{session_num} = [mean(corr_tbyt_peak_mag); mean(fail_tbyt_peak_mag)];
%     else
%         %Remove corr trial with highest lick rate
%         [sorted_corr_LR, sorted_corr_LR_inx] = sort(lick_rate_corr, ind_order);
%         sorted_corr_peaks = corr_tbyt_peak_mag([sorted_corr_LR_inx],:);
%         for trial_num = 2:length(lick_rate_corr)
%             temp_mean = mean(sorted_corr_LR(trial_num:end));
%             prev_temp_mean = mean(sorted_corr_LR(trial_num-1:end));
%             if temp_mean == mean_LR_early
%                 break;
%             elseif (temp_mean<mean_LR_early & prev_temp_mean>mean_LR_early) | (temp_mean>mean_LR_early & prev_temp_mean<mean_LR_early)
%                 if abs(temp_mean - mean_LR_early) >= abs(prev_temp_mean - mean_LR_early) %Go with which ever trial num has the smallest abs diff between LRs for corr and early
%                     trial_num = trial_num-1; %using trial_num as there reference for which trials to remove
%                 end
%                 break;
%             end
%             prev_temp_mean = temp_mean; %set prev temp mean before entering the next iteration of the forloop
%             if length(trial_num:length(sorted_corr_LR)) < 5 
%                 continue
%             end
%         end
%         if length(trial_num:length(sorted_corr_LR)) < 5
%             continue
%         end
%         lick_rate_corr = sorted_corr_LR(trial_num:end);
%         corr_tbyt_peak_mag = sorted_corr_peaks([trial_num:end],:);
%         mean_mags{session_num} = [mean(corr_tbyt_peak_mag); mean(fail_tbyt_peak_mag)];  %==============================problem is that the algorithm can never get lick_rate_corr as low as lick_rate_fail. Keeps removing trials
%     days{session_num}
%     size(corr_tbyt_peak_mag)
%     size(fail_tbyt_peak_mag)
%     end
    %----------
    
    %determine spearman's correlation
    if length(lick_rate_corr)<6 | length(lick_rate_early)<6 
        continue
    end
    lick_rate_corr_mat = [lick_rate_corr_mat, mean(lick_rate_corr)];
    lick_rate_early_mat = [lick_rate_early_mat, mean(lick_rate_early)];
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
    
    %plot lick rate vs peak magnitude for each trial - correct trials
%     figure;    plot_colors;
%     x_max = max(max(corr_tbyt_peak_mag));
%     for ROI_num = 1:size(corr_tbyt_peak_mag,2)
%         subplot(size(corr_tbyt_peak_mag,2), 1, ROI_num);
%         scatter(corr_tbyt_peak_mag(:,ROI_num)', lick_rate_corr, plot_colors{ROI_num}); hold on;
%         ylim([0 20]);
%         xlabel('peak df/f');
%         xlim([0 x_max]);
%     end
%     ylabel('lick rate (Hz)');
%     suptitle([days{session_num}, 'correct trials']);
    %savefig([figure_output, days{session_num}, '_corr_', num2str(length(lick_rate_window))]);
    
    %plot lick rate vs peak magnitude for each trial - early trials
%     figure;
%     x_max = max(max(fail_tbyt_peak_mag));
%     for ROI_num = 1:size(corr_tbyt_peak_mag,2)
%         subplot(size(corr_tbyt_peak_mag,2), 1, ROI_num);
%         scatter(fail_tbyt_peak_mag(:,ROI_num), lick_rate_early); hold on;
%         xlim([0 x_max]);
%         ylim([0 20]);
%         ylabel('lick rate (Hz)');
%     end
%     xlabel('peak df/f'); 
%     suptitle([days{session_num}, 'early trials']);
    %savefig([figure_output, days{session_num}, '_fail_', num2str(length(lick_rate_window))]);
end
%save([figure_output, 'Spearman_correlation_', num2str(length(lick_rate_window))], 'corr_Scorr', 'early_Scorr');

%store values for spearmans correlations in matrices to plot
mean_correlation_corr = [];
mean_correlation_fail = [];
all_correlation_corr = [];
all_correlation_fail = [];
correlation_P_corr = [];
correlation_P_fail = [];
for session_num = 1:length(corr_Scorr) %img30 has some confound in the licking data
    if isempty(corr_Scorr{session_num})
        continue
    end
    mean_correlation_corr = [mean_correlation_corr, mean(corr_Scorr{session_num}(1,:))];
    mean_correlation_fail = [mean_correlation_fail, mean(early_Scorr{session_num}(1,:))];
    correlation_P_corr = [correlation_P_corr, corr_Scorr{session_num}(2,:)];
    correlation_P_fail = [correlation_P_fail, early_Scorr{session_num}(2,:)];
    all_correlation_corr = [all_correlation_corr, corr_Scorr{session_num}(1,:)];
    all_correlation_fail = [all_correlation_fail, early_Scorr{session_num}(1,:)];
end

%plot all the spearman's correlations
figure;
scatter([ones(1,length(mean_correlation_corr))], mean_correlation_corr); hold on;
scatter(1, mean(mean_correlation_corr),'k', 'filled');
errorbar(1, mean(mean_correlation_corr), std(mean_correlation_corr)/sqrt(length(mean_correlation_corr)), 'k');
scatter([ones(1,length(find(~isnan(mean_correlation_fail))))]+1, mean_correlation_fail(find(~isnan(mean_correlation_fail))), 'r');
scatter(2, nanmean(mean_correlation_fail),'k', 'filled')
errorbar(2, mean(mean_correlation_fail), std(mean_correlation_fail)/sqrt(length(mean_correlation_fail)), 'k');
title(['within session correlations between lick rate and df/f. Window=' num2str(length(lick_rate_window))]);
ylabel('Rho');
xlabel('blue: correct trials      red: early trials');
xlim([-1 4]);
hline(0, 'k');
[t_corr,p_corr, ~, stats_corr] = ttest(mean_correlation_corr);
[t_fail,p_fail, ~, stats_fail] = ttest(mean_correlation_fail);
%savefig([figure_output, 'ROI_correlations_', num2str(length(lick_rate_window))]);

%plot spearman's correlations in the form of a scatterplot with confidence intervals
figure; hold on;
corr_v_early_model = fitlm(mean_correlation_corr, mean_correlation_fail(find(~isnan(mean_correlation_fail))));
plot(corr_v_early_model);
ylim([-0.65 0.65]); xlim([-0.2 0.2])
hline(0, '--k'); vline(0, '--k');
title('spearmans correlation between lick rate and peak df/f on a trial by trial basis');
ylabel('early trials');
xlabel('correct trials');



if exist('plot_mean_mags')
    mean_mags{11} = (mean_mags{11} + mean_mags{12})/2;  %averaging together both sessions for img38 so it does not skew the average. 
    mean_mags{12} = [];
    mean_mags_mat_corr = [];
    mean_mags_mat_early = [];
    figure;
    for session_num = 1:length(mean_mags)
        if ~isempty(mean_mags{session_num})
            %for ROI_num = 1:size(mean_mags{session_num},2)
                mean_mags_mat_corr = [mean_mags_mat_corr mean(mean_mags{session_num}(1,:),2)];
                mean_mags_mat_early = [mean_mags_mat_early, mean(mean_mags{session_num}(2,:),2)];
                plot([1,2], [mean(mean_mags{session_num}(1,:),2)/mean(mean_mags{session_num}(1,:),2), mean(mean_mags{session_num}(2,:),2)/mean(mean_mags{session_num}(1,:),2)], strcat('-', plot_symbols{session_num}, colors_final_plot{session_num}));   hold on;
                %plot([1,2], [mean_mags{session_num}(1,ROI_num)/mean_mags{session_num}(1,ROI_num), mean_mags{session_num}(2,ROI_num)/mean_mags{session_num}(1,ROI_num)], strcat('-', plot_symbols{session_num}, colors_final_plot{session_num}));   hold on;
            %end
        end
    end
    %errorbar(1,mean(mean_mags_mat_corr), std(mean_mags_mat_corr)/sqrt(length(mean_mags_mat_corr)), 'k');
    %errorbar(2,mean(mean_mags_mat_early), std(mean_mags_mat_early)/sqrt(length(mean_mags_mat_early)), 'k');
    %scatter([1,2], [mean(mean_mags_mat_corr), mean(mean_mags_mat_early)], 'filled', 'k');
    plot(2, mean(mean_mags_mat_early./mean_mags_mat_corr), 'ok', 'MarkerFaceColor', 'k');
    errorbar(2, mean(mean_mags_mat_early./mean_mags_mat_corr), std(mean_mags_mat_early./mean_mags_mat_corr)/sqrt(length(mean_mags_mat_corr)), 'k');
    [hh,pp] = ttest(mean_mags_mat_corr, mean_mags_mat_early);
    xlabel('corr trials     early trials');
    ylabel('peak magnitude df/f');
    title('comparison of df/f peak magnitudes for correct vs early trials with matched lick rates. By ROI.');
    ylim([0 1.2]); xlim([0.5 2.5]);
    %savefig([figure_output, 'compare_peak_mag_session_base_LR_matching_', num2str(length(lick_rate_window))]);
end



