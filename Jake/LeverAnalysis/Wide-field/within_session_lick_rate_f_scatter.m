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
early_Scorr = {}; early_Scorr{1}=[];
corr_Scorr = {}; corr_Scorr{1} = [];
plot_colors = {'r', 'b', 'g', 'k', 'c', 'm'};

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
    lick_rate_corr = sum(lick_trace_corr(:,[lick_rate_window+lever_release_frame]),2)' / length(lick_rate_window) * 10;  %converts lick rate to Hz 
    lick_rate_early = sum(lick_trace_early(:,[lick_rate_window+lever_release_frame]),2)' / length(lick_rate_window) * 10;
    
    %determine spearman's correlation
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
for session_num = 1:length(days)
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
scatter(1, mean(mean_correlation_corr),'k', 'filled')
scatter([ones(1,length(find(~isnan(mean_correlation_fail))))]+1, mean_correlation_fail(find(~isnan(mean_correlation_fail))));
scatter(2, nanmean(mean_correlation_fail),'k', 'filled')
title(['within session correlations between lick rate and df/f. Window=' num2str(length(lick_rate_window))]);
ylabel('Rho');
xlabel('blue: correct trials      red: early trials')
xlim([-1 4]);
hline(0, 'k');
savefig([figure_output, 'ROI_correlations_', num2str(length(lick_rate_window))]);

