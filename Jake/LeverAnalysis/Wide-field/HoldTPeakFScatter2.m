%make scatter plots of the hold time duration vs the peak df/f for corr or
%early trials (trial by trial)
%binning in the x axis for each session. Then plotting the mean across
%sessions and fitting a line to that plot. 

clear;
WF_plotting_lists_of_days;
struct_dir = 'Z:\Analysis\WF Lever Analysis\StructuresForPlotting\';
TC_dir = 'Z:\Analysis\WF Lever Analysis\LeverSummaryFolder\';
bx_data_dir = 'Z:\Data\WidefieldImaging\GCaMP\behavior\';
output_dir = 'Z:\Analysis\WF Lever Analysis\hold time vs peak f\';
bx_outputs_dir = ['Z:\Analysis\WF Lever Analysis\BxAndAnalysisOutputs\BxOutputs\']; 
colors = {'k','k', 'b','b', 'g','g', 'c','c', 'm','m', 'r','r', 'y','y', 'k', 'b', 'r'}; 
plot_symbols = {'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 's', 's', 's'};
colors_roi = {'b', 'r', 'g', 'k', 'c', 'm'};
peak_window = [5:9];
%peak_window = [1:2]; %modified version used to find the pre-release df/f
peak_window = [6:10]; %consistent with the scatterplots
shift_window = [1:3];
bin_duration = 250;
lick_rate_window = [6:11];
peak_diff{1} = [];
peak_corr_long_cell{1} = [];
peak_corr_short_cell{1} = [];
peak_early_long_cell{1} = [];
peak_early_short_cell{1} = [];
corr_hold_dur_cell{1} = []; 
early_hold_dur_cell{1} = [];
corr_values_correct{1} = [];
corr_values_early{1} = [];
sample_sizes{1} = [];
across_animals_fail_peak_binned = [];
across_animals_corr_peak_binned = [];
across_animals_corr_lick_binned = [];
across_animals_fail_lick_binned = [];
peak_corr_long_mat= [];
peak_early_long_mat= [];
peak_corr_short_mat= [];
peak_early_short_mat=[];
across_animals_corr_bin_N = [];
across_animals_early_bin_N = [];
corr_trial_num{1} = [];
early_trial_num{1} = [];

for session_num = 1:length(days)
    %load data
    bx_data = get_bx_data(bx_data_dir, days{session_num}); %use the bx data file to determine max Random hold duration 1500 vs 5000
    load([TC_dir, days{session_num}, '_success']);
    load([TC_dir, days{session_num}, '_fail']);
    load([struct_dir, days{session_num}]);
    corr_hold_dur = curr_struct.corr_hold_dur;
    early_hold_dur = curr_struct.early_hold_dur;
    bx_out_dir  = [bx_outputs_dir days{session_num} '_bx_outputs'];
    load(bx_out_dir, 'licking_data');
    corr_licks = licking_data.lick_trace_succ;
    fail_licks = licking_data.lick_trace_fail;
    lick_rate_corr = mean(corr_licks(:,[lick_rate_window]),2);
    lick_rate_fail = mean(fail_licks(:,[lick_rate_window]),2);
    
    %determine max rand hold
    if bx_data.randReqHoldMaxMs > 4000
        long_hold = 1;
    elseif bx_data.randReqHoldMaxMs < 2000
        long_hold = 0;
    else
        disp('error in max rand hold duration value'); 
        pause
    end
      
    %remove non-LS ROIs
    success_roi = success_roi(:,[LS_ROIs{session_num}],:);
    fail_roi = fail_roi(:,[LS_ROIs{session_num}],:);
    
    %apply shift to TCs and get peak mags 
    shift = mean(success_roi(:,:,shift_window),3);
    shift = repmat(shift, [1,1,size(success_roi,3)]);
    success_roi = success_roi - shift;
    shift = mean(fail_roi(:,:,shift_window),3);
    shift = repmat(shift, [1,1,size(fail_roi,3)]);
    fail_roi = fail_roi - shift;
    
    %find the peak values and get the difference
    peak_corr_roi = max(success_roi(:,:,peak_window), [], 3); %use max when looking for the peak
    peak_early_roi = max(fail_roi(:,:,peak_window), [], 3);
%     peak_corr_roi = mean(success_roi(:,:,peak_window), 3);  %used the mean when looking at the pre-release
%     peak_early_roi = mean(fail_roi(:,:,peak_window), 3);
    peak_diff{session_num} = mean(peak_corr_roi)-mean(peak_early_roi);
    
%     %plot each trials hold time vs peak df/f 
%     figure;  
%     subplot(1,2,1);
%     for ROI_num = 1:size(success_roi,2)
%         scatter(corr_hold_dur, peak_corr(:,ROI_num), colors_roi{ROI_num}); hold on;
%         model_1 = fitlm(corr_hold_dur, peak_corr(:,ROI_num));
%         plot(model_1)
%     end
%     ylabel('correct trials peak df/f'); xlabel('hold duration (ms)');
%     subplot(1,2,2);
%     for ROI_num = 1:size(fail_roi,2)
%         scatter(early_hold_dur, peak_early(:,ROI_num), colors_roi{ROI_num}); hold on;
%     end
%     ylabel('early trials peak df/f'); xlabel('hold duration (ms)');
%     savefig([output_dir, days{session_num}]);
    
    %collapse across ROIs
    peak_corr = mean(peak_corr_roi,2);
    peak_early = mean(peak_early_roi,2);
    
    %normalize df/f values to the mean
    peak_corr_mean = mean(peak_corr);
    peak_early_mean = mean(peak_early);
    peak_corr_norm = peak_corr/peak_corr_mean;
    peak_early_norm = peak_early/peak_corr_mean;
    %normalize df/f for each roi
    peak_corr_roi_mean = mean(peak_corr_roi);
    peak_corr_roi_mean = repmat(peak_corr_roi_mean, size(peak_corr_roi,1), 1);
    peak_corr_roi_norm = peak_corr_roi./peak_corr_roi_mean;
    %peak_early_roi_mean = mean(peak_early_roi);
    %peak_early_roi_mean = repmat(peak_early_roi_mean, size(peak_early_roi,1), 1);
    peak_corr_roi_mean = repmat(mean(peak_corr_roi), size(peak_early_roi,1), 1); %normalize earlies to the corrects mean
    peak_early_roi_norm = peak_early_roi./peak_corr_roi_mean;
    
    %bin the hold durations and lick rates
    corr_peak_binned = [];
    fail_peak_binned = [];
    corr_lick_binned = [];
    fail_lick_binned = [];
    curr_edge_min = 0;
    curr_edge_max = bin_duration;
    num_bins = 5500/bin_duration;
    ROI_dim_size = size(success_roi,2);
    for this_bin = 1:num_bins
        holds_this_bin = find(corr_hold_dur>curr_edge_min & corr_hold_dur<=curr_edge_max);
        if isempty(holds_this_bin)
            corr_peak_binned(this_bin, [1:ROI_dim_size]) = NaN;
            corr_lick_binned(this_bin) = NaN;
            session_holds_corr(this_bin) = 0;
        else
            peaks_this_bin = peak_corr_roi([holds_this_bin], :);
            %peaks_this_bin = peak_corr_roi_norm([holds_this_bin], :);
            lick_rates_this_bin = lick_rate_corr([holds_this_bin]);
            corr_peak_binned(this_bin, [1:ROI_dim_size]) = mean(peaks_this_bin,1);
            corr_lick_binned(this_bin) = mean(lick_rates_this_bin);  %dims = 1 by N with N being the number of bins
            session_holds_corr(this_bin) = length(holds_this_bin);
        end
        holds_this_bin = find(early_hold_dur>curr_edge_min & early_hold_dur<=curr_edge_max);
        if isempty(holds_this_bin)
            fail_peak_binned(this_bin, [1:ROI_dim_size]) = NaN;
            fail_lick_binned(this_bin) = NaN;
            session_holds_early(this_bin) = 0;
        else
            peaks_this_bin = peak_early_roi([holds_this_bin], :);
            %peaks_this_bin = peak_early_roi_norm([holds_this_bin], :);
            lick_rates_this_bin = lick_rate_fail([holds_this_bin]);
            fail_peak_binned(this_bin, [1:ROI_dim_size]) = mean(peaks_this_bin,1);
            fail_lick_binned(this_bin) = mean(lick_rates_this_bin);
            session_holds_early(this_bin) = length(holds_this_bin);
        end
        curr_edge_min = curr_edge_max;
        curr_edge_max = curr_edge_max + bin_duration;
    end
    across_animals_corr_bin_N(session_num,:)= session_holds_corr;
    across_animals_early_bin_N(session_num,:)= session_holds_early;
    across_animals_corr_peak_binned(:,end+1) = mean(corr_peak_binned,2);
    across_animals_fail_peak_binned(:,end+1) = mean(fail_peak_binned,2);
    corr_peak_binned_cell{session_num} = corr_peak_binned;
    fail_peak_binned_cell{session_num} = fail_peak_binned;
    if session_num == 2 | session_num == 3
    else
        across_animals_corr_lick_binned(:,end+1) = corr_lick_binned'*10; %converting to Hz since lick rate was collected as the MEAN over a 500ms window
        across_animals_fail_lick_binned(:,end+1) = fail_lick_binned'*10;
    end
    
    %log sample sizes
    sample_sizes{session_num, :} = [length(corr_hold_dur), length(early_hold_dur)];
    corr_trial_num{session_num} = length(peak_corr);
    early_trial_num{session_num} = length(peak_early);
    
    %store normalized df/f for corr/early in long_dur or short_dur
    if long_hold == 1
        peak_corr_long_cell{session_num} = peak_corr_norm;
        peak_early_long_cell{session_num} = peak_early_norm;
        peak_corr_long_mat(:,end+1) = mean(corr_peak_binned,2);
        peak_early_long_mat(:,end+1) = mean(fail_peak_binned,2);
    else
        peak_corr_short_cell{session_num} = peak_corr_norm;
        peak_early_short_cell{session_num} = peak_early_norm;
        peak_corr_short_mat(:,end+1) = mean(corr_peak_binned,2);
        peak_early_short_mat(:,end+1) = mean(fail_peak_binned,2);
    end
    
    %store hold dur corr/early no need for seperate long or short. Will
    %just use indeces to extract the right ones
    corr_hold_dur_cell{session_num} = corr_hold_dur; 
    early_hold_dur_cell{session_num} = early_hold_dur;
end
x_axis = [1:22].*bin_duration;
figure; %subplot(2,1,1);
errorbar(x_axis,nanmean(across_animals_corr_peak_binned,2), nanstd(across_animals_corr_peak_binned, [],2)./sqrt(sum(~isnan(across_animals_corr_peak_binned),2)), '-ok'); hold on;
%errorbar(x_axis,nanmean(across_animals_fail_peak_binned,2), nanstd(across_animals_fail_peak_binned, [],2)./sqrt(sum(~isnan(across_animals_fail_peak_binned),2)), '-or');
errorbar(x_axis(1:end-5),nanmean(across_animals_fail_peak_binned([1:end-5],:),2), nanstd(across_animals_fail_peak_binned([1:end-5],:), [],2)./sqrt(sum(~isnan(across_animals_fail_peak_binned([1:end-5],:)),2)), '-or'); %excluded two data points to match the 2P data
plot(x_axis,nanmean(across_animals_corr_peak_binned,2),'ok');
%plot(x_axis,nanmean(across_animals_fail_peak_binned,2),'or');
plot(x_axis(1:end-5),nanmean(across_animals_fail_peak_binned([1:end-5],:),2),'or');
xlabel('hold duration (250ms bins)'); ylabel('normalized df/f'); 
%ylim([0 3.6]); 
xlim([0 5700]);
%title('hold duration vs pre-release df/f: prerelease = frames -5 to -4 relative to lever release');
title('hold duration vs peak df/f: all sessions n=17');
% corr_fit = fitlm(x_axis,nanmean(across_animals_corr_peak_binned,2));
% early_fit = fitlm(x_axis,nanmean(across_animals_fail_peak_binned,2));
[corr_fit, corr_gof] = fit(x_axis(3:end)',nanmean(across_animals_corr_peak_binned([3:end],:),2), 'poly1');
[early_fit, early_gof] = fit(x_axis(1:end-5)',nanmean(across_animals_fail_peak_binned([1:end-5],:),2), 'poly1');
%subplot(2,1,2);
%plot(corr_fit, x_axis(3:end)',nanmean(across_animals_corr_peak_binned([3:end],:),2)); hold on;
plot(early_fit); hold on;
plot(corr_fit);
[Rho_corr, p_corr] = corr(x_axis(3:end)',nanmean(across_animals_corr_peak_binned([3:end],:),2), 'type', 'Spearman');
[Rho_early, p_early] = corr(x_axis(1:end-5)',nanmean(across_animals_fail_peak_binned([1:end-5],:),2), 'type', 'Spearman');

modifier_mat = nanmean(across_animals_fail_peak_binned([1:end-5],:),2)- nanmean(across_animals_fail_peak_binned(1,:),2)
[Rho_early, p_early] = corr(x_axis(1:end-5)',nanmean(across_animals_fail_peak_binned([1:end-5],:),2)-modifier_mat, 'type', 'Spearman');

corr_mdl = fitlm(x_axis(3:end)',nanmean(across_animals_corr_peak_binned([3:end],:),2));

%plot only the long duration sessions
figure; 
errorbar(x_axis, nanmean(peak_corr_long_mat,2), nanstd(peak_corr_long_mat, [],2)./sqrt(sum(~isnan(peak_corr_long_mat),2)), '-ok'); hold on;
errorbar(x_axis,nanmean(peak_early_long_mat,2), nanstd(peak_early_long_mat, [],2)./sqrt(sum(~isnan(peak_early_long_mat),2)), '-or');
%  corr_fit_lin = fitlm(x_axis,nanmean(peak_corr_long_mat,2));
%  early_fit_lin = fitlm(x_axis,nanmean(peak_early_long_mat,2));
%[corr_fit_lin, corr_fit_lin_gof, corr_fit_lin_output] = fit(x_axis(3:end)', nanmean(peak_corr_long_mat([3:end],:),2), 'poly2');
early_fit_lin = fitlm(x_axis,nanmean(peak_early_long_mat,2));
plot(corr_fit_lin, x_axis(3:end)', nanmean(peak_corr_long_mat([3:end],:),2))
plot(early_fit_lin);
[corr_fit, corr_fit_gof, corr_fit_output] = fit(x_axis(3:end)', nanmean(peak_corr_long_mat([3:end],:),2), 'exp2');
plot(corr_fit, x_axis(3:end)', nanmean(peak_corr_long_mat([3:end],:),2));
[early_fit, early_gof] = fit(x_axis(1:end-2)', nanmean(peak_early_long_mat([1:end-2],:),2), 'poly1');
%plot(early_fit, x_axis(1:end-2)', nanmean(peak_early_long_mat([1:end-2],:),2))
title('long duration sessions');
xlabel('hold duration (250ms bins)');
ylabel('df/f');

%plot only the short duration sessions
figure; 
errorbar(x_axis,nanmean(peak_corr_short_mat,2), nanstd(peak_corr_short_mat, [],2)./sqrt(sum(~isnan(peak_corr_short_mat),2)), 'b'); hold on;
errorbar(x_axis,nanmean(peak_early_short_mat,2), nanstd(peak_early_short_mat, [],2)./sqrt(sum(~isnan(peak_early_short_mat),2)), 'r');
corr_fit = fitlm(x_axis,nanmean(peak_corr_long_mat,2));
early_fit = fitlm(x_axis,nanmean(peak_early_short_mat,2));
title('short duration sessions');

% %plot binned lick rate
% figure; 
% plot(x_axis, nanmean(across_animals_corr_lick_binned,2), 'ok');
% plot(x_axis, nanmean(across_animals_fail_lick_binned,2), 'or');
% errorbar(x_axis, nanmean(across_animals_corr_lick_binned,2), nanstd(across_animals_corr_lick_binned,[],2)./sqrt(sum(~isnan(across_animals_corr_lick_binned),2)), 'k'); hold on;
% errorbar(x_axis, nanmean(across_animals_fail_lick_binned,2), nanstd(across_animals_fail_lick_binned, [],2)./sqrt(sum(~isnan(across_animals_fail_lick_binned),2)), 'r');
% title(['lick rate vs hold duration: n=', num2str(size(across_animals_fail_lick_binned,2))]);
% xlabel('hold durations (250ms bins)');
% ylabel('lick rate (Hz)');







