% probability of a significant df/f response on lever events
clear;

%assign some variables
lever_frame = 6;
peak_window = [lever_frame-1:lever_frame+3];
shift_window = [1:3];
pre_frames = 5;
post_frames = 10;
percent_significant_ROI_corr{1} = [];
percent_significant_ROI_fail{1} = [];
percent_significant_corr{1} = [];
percent_significant_fail{1} = [];
trial_num_corr{1} = [];
trial_num_fail{1} = [];

%set directories and load experiment labels
bx_source      = ['Z:\Data\WidefieldImaging\GCaMP\behavior\'];
bx_outputs_dir = ['Z:\Analysis\WF Lever Analysis\BxAndAnalysisOutputs\BxOutputs\'];  
baseline_output = ['Z:\Analysis\WF Lever Analysis\licking_investigation\baselines_ITI_no_licks'];
TC_dir = 'Z:\Analysis\WF Lever Analysis\LeverSummaryFolder\';
old_cd = cd; %save old cd so I can restore it later
WF_plotting_lists_of_days;

%datasets with good licking traces
days = {'151021_img29', '151022_img29', [], [], '151211_img32', '151212_img32', '160129_img35', '160131_img35', '160129_img36', '160131_img36', '160314_img38', '160315_img38', '160319_img41', '160320_img41', '160606_img46', '160722_img53', '160904_img55'}; %'150718_img27', '150719_img27', '150716_img28', '150717_img28', 
colors_roi = {'b', 'r', 'g', 'k', 'c', 'm'};
colors = {'k','k', 'b','b', 'g','g', 'c','c', 'm','m', 'r','r', 'y','y', 'k', 'b', 'r'}; 
plot_symbols = {'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 's', 's', 's'};

for session_num = 1:length(days);
    if isempty(days{session_num})
        lick_trig_mags{session_num} = [];
        continue
    end
    b_data = get_bx_data(bx_source, days{session_num});  %find the correct behavior file and loads it.
    bx_out_dir  = [bx_outputs_dir days{session_num} '_bx_outputs'];
    load(bx_out_dir);
    load([baseline_output, '\', days{session_num}]);
    load([TC_dir, days{session_num}, '_success']);
    load([TC_dir, days{session_num}, '_fail']);
        
    %assign variables
    imaging_start_MW_T = frame_info.imaging_start_MW_T;
    licksByFrame = licking_data.licksByFrame; 
    leverPress = cell2mat(b_data.holdStartsMs)-imaging_start_MW_T;
    holdDur = cell2mat(b_data.holdTimesMs);
    leverRelease = leverPress + holdDur;
    reqHoldTime = cell2mat(b_data.tTotalReqHoldTimeMs);
    changeOrientation = trial_outcome.change_orientation;
    
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
    
    %get the mean and sem for each ROI across trials
    success_roi_roi_avg = squeeze(mean(success_roi,1));
    success_roi_roi_sem = squeeze(std(success_roi, [], 1)./sqrt(size(success_roi,1)));
    fail_roi_roi_avg = squeeze(mean(success_roi,1));
    fail_roi_roi_sem = squeeze(std(success_roi, [], 1)./sqrt(size(success_roi,1)));
    if size(success_roi,2) == 1
        success_roi_roi_avg = success_roi_roi_avg'
        success_roi_roi_sem = success_roi_roi_sem';
        fail_roi_roi_avg = fail_roi_roi_avg';
        fail_roi_roi_sem = fail_roi_roi_sem';
    end
    
    %find the peak values
    peak_corr = max(success_roi(:,:,peak_window), [], 3);
    peak_fail = max(fail_roi(:,:,peak_window), [], 3);
    
    %determine what percent of trials exceed 2 stdevs first by ROI
    this_session_sig_corr = [];
    this_session_sig_fail = [];
    for ROI_num = 1:size(peak_corr,2);
        this_session_sig_corr = [this_session_sig_corr, length(find(peak_corr(:,ROI_num)>2*baseline_frames_std(ROI_num)))/size(peak_corr,1)];
        this_session_sig_fail = [this_session_sig_fail, length(find(peak_fail(:,ROI_num)>2*baseline_frames_std(ROI_num)))/size(peak_fail,1)];
    end
    percent_significant_ROI_corr{session_num} = this_session_sig_corr;
    percent_significant_ROI_fail{session_num} = this_session_sig_fail;
        
    %then averaged across ROIs
    success_roi_avg = squeeze(mean(success_roi,2));
    peak_mag_avg_corr = max(success_roi_avg(:,peak_window),[],2);
    percent_significant_corr{session_num} = length(find(peak_mag_avg_corr>2*mean(baseline_frames_std)))/length(peak_mag_avg_corr);
    trial_num_corr{session_num} = length(peak_mag_avg_corr);
    fail_roi_avg = squeeze(mean(fail_roi,2));
    peak_mag_avg_fail = max(fail_roi_avg(:,peak_window),[],2);
    percent_significant_fail{session_num} = length(find(peak_mag_avg_fail>2*mean(baseline_frames_std)))/length(peak_mag_avg_fail);
    trial_num_fail{session_num} = length(peak_mag_avg_fail);
    
    %generate index of which trials have a significant response: for plotting
    sig_corr_resp_mat = zeros(size(peak_corr));
    for ROI_num = 1:size(peak_corr,2)
        sig_corr_resp_mat([find(peak_corr(:,ROI_num)>2*baseline_frames_std(ROI_num))], ROI_num) = 1; %find all trials with significant response and log them in this matrix on an ROI basis
    end
    sig_fail_resp_mat = zeros(size(peak_fail));
    for ROI_num = 1:size(peak_fail,2)
        sig_fail_resp_mat([find(peak_fail(:,ROI_num)>2*baseline_frames_std(ROI_num))], ROI_num) = 1; %find all trials with significant response and log them in this matrix on an ROI basis
    end
        
    %make a color coded plot: trial by trial TCs of sig vs non-sig trials
    figure;
    for ROI_num = 1:size(peak_corr,2)
        subplot(1,size(peak_corr,2),ROI_num); hold on;
        this_roi_sig_frac = sum(sig_corr_resp_mat(:,ROI_num))/size(sig_corr_resp_mat,1);
        title(['ROI #', num2str(ROI_num), ' % of sig. trials=',num2str(this_roi_sig_frac)]);
        ylabel('df/f'); xlabel('frames relative to release');
        for trial_num = 1:size(peak_corr,1)
            if sig_corr_resp_mat(trial_num, ROI_num) == 1
                plot([-pre_frames:post_frames], squeeze(success_roi(trial_num, ROI_num, :)), 'b');
            elseif sig_corr_resp_mat(trial_num, ROI_num) == 0
                plot([-pre_frames:post_frames], squeeze(success_roi(trial_num, ROI_num, :)), 'k');
            end
        end
        errorbar([-pre_frames:post_frames], success_roi_roi_avg(ROI_num,:), success_roi_roi_sem(ROI_num,:), 'm', 'LineWidth', 1.5);
        xlim([-5 10]);
    end
    suptitle([days(session_num), ' correct trials: blue=significant black=not']);
    
    figure;
    for ROI_num = 1:size(peak_fail,2)
        subplot(1,size(peak_fail,2),ROI_num); hold on;
        this_roi_sig_frac = sum(sig_fail_resp_mat(:,ROI_num))/size(sig_fail_resp_mat,1);
        title(['ROI #', num2str(ROI_num), ' % of sig. trials=',num2str(this_roi_sig_frac)]);
        ylabel('df/f'); xlabel('frames relative to release');
        for trial_num = 1:size(peak_fail,1)
            if sig_fail_resp_mat(trial_num, ROI_num) == 1
                plot([-pre_frames:post_frames], squeeze(fail_roi(trial_num, ROI_num, :)), 'b');
            elseif sig_fail_resp_mat(trial_num, ROI_num) == 0
                plot([-pre_frames:post_frames], squeeze(fail_roi(trial_num, ROI_num, :)), 'k');
            end
        end
        errorbar([-pre_frames:post_frames], fail_roi_roi_avg(ROI_num,:), fail_roi_roi_sem(ROI_num,:), 'm', 'LineWidth', 1.5);
        xlim([-5 10]);
    end
    suptitle([days(session_num), ' early trials: blue=significant black=not']);
end

%plot the fraction of significant trials for each ROI 
percent_sig_mat_corr = [];
percent_sig_mat_fail = [];
figure;
for session_num = 1:length(days)
    if ~isempty(percent_significant_ROI_corr{session_num})  | trial_num_corr{session_num} > 10
        plot(zeros(1, length(percent_significant_ROI_corr{session_num})), percent_significant_ROI_corr{session_num}, strcat(plot_symbols{session_num}, colors{session_num})); hold on;
        percent_sig_mat_corr = [percent_sig_mat_corr, percent_significant_ROI_corr{session_num}];
    end
    if ~isempty(percent_significant_ROI_fail{session_num})  | trial_num_fail{session_num} > 10
        plot(ones(1, length(percent_significant_ROI_fail{session_num})), percent_significant_ROI_fail{session_num}, strcat(plot_symbols{session_num}, colors{session_num})); hold on;
        percent_sig_mat_fail = [percent_sig_mat_fail, percent_significant_ROI_fail{session_num}];
    end
end
errorbar(0, mean(percent_sig_mat_corr), std(percent_sig_mat_corr)/sqrt(length(percent_sig_mat_corr)), 'o', 'MarkerFaceColor', 'k');
errorbar(1, mean(percent_sig_mat_fail), std(percent_sig_mat_fail)/sqrt(length(percent_sig_mat_fail)), 'o', 'MarkerFaceColor', 'k');
xlim([-0.5 1.5]);
ylim([0 1]); ylabel('percent of trials with a significant positive df/f response');
title('trial peak df/f > 2 standard deviations above baseline.');
xlabel('corr     early');
ylabel('# significant trials / # fraction of trials');







