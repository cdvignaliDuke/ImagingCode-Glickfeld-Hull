%WF_CRP_lever_plotting_TCs

%PLOTTING GRAPHS 
%plots correlation coefficient between the ROIs
%plots TCs of the df/f for success, earlies, tooFast, fidget and press.  
clear
BEHAVE_DIR = 'Y:\home\jake\Data\WidefieldImaging\GCaMP\behavior\';
ANALYSIS_DIR = 'Y:\home\jake\Analysis\WF Lever Analysis\';
CLUSTER_DIR  ='Y:\home\jake\Analysis\WF Lever Analysis\BxAndAnalysisOutputs\'; 
WF_CRP_list_of_days; 
days = days_1;%([1,4,5,6,7]);  %'151212_img32',  
save_data = 0;

for kk= 1%[1,3:7]%:length(days)
    %set directories and load bxOutputs and cluter data. 
    ROI_name  =  days{kk};
    b_data = get_bx_data(BEHAVE_DIR, days{kk});  %find the correct behavior file and loads it.
    load([ANALYSIS_DIR 'BxAndAnalysisOutputs\BxOutputs\', days{kk}, '_bx_outputs']);
    load([CLUSTER_DIR, days{kk}, '\', days{kk}, '_cluster']); 
    
    %define variables to help with plotting
    func = @mean; %func = @std; @median;
    pre_frames = 5;
    post_frames = 20;
    sampling_rate = round(1000/mode(diff(frame_info.times)));
    ts = (-pre_frames:post_frames)*1000/double(round(sampling_rate));
    ts = repmat(ts,[cluster.num_cluster 1]);
    tot_frame = pre_frames + post_frames+1;
    colors = [1,0,0; 0,1,0; 0,0,1; 0.5,0.5,0.5; 1,0,1; 1,1,0; 0,1,1]; %sets up the color scheme for plotting multiple ROIs with errorbar  
    
    %bin licking times according to frame number.
    licking_data = bin_licking_by_frame(lickTimes, frame_info);
    
    %PLOT the ROIs frome the cluter data----------------------------------
    figure;
    plot_ROIs_on_heatmap(avg_img, sz, b_data, days{kk}, cluster); 
    
    %----PLOT CUE REWARD PAIRING SESSIONS------------------------------------------
    %plot rewarded trials
    use_ev_rewarded = trial_outcome.rewarded_trial_time;
    use_ev_rewarded = use_ev_rewarded(1:round(length(use_ev_rewarded)/2));
    [rewarded_roi, use_times_rew, lick_trace_rew, lick_trace_rew_10ms] = trigger_movie_by_event_licks(tc_dfoverf, frame_info, licking_data, use_ev_rewarded, pre_frames, post_frames);
    avg_rewarded_roi = squeeze(func(rewarded_roi,1));
    if cluster.num_cluster == 1
        avg_rewarded_roi = avg_rewarded_roi';
    end
    std_rewarded = squeeze(std(squeeze(rewarded_roi),1));
    sm_rewarded = std_rewarded./sqrt(size(rewarded_roi,1));
    rewarded_roi_mean = squeeze(mean(rewarded_roi,2));
    sm_rewarded_avg = std(rewarded_roi_mean,1)/size(rewarded_roi_mean,1);
    for i = 1:cluster.num_cluster  %baseline each curve so it passes through zero
        shift = (-1)*avg_rewarded_roi(i,1);
        avg_rewarded_roi(i,:) = avg_rewarded_roi(i,:)+shift;
    end
    subplot(1,3,2); lickBars = bar(ts(1,:), mean(lick_trace_rew)/10); hold on;
    errorbar(ts(1,:), mean(lick_trace_rew)/10, [std(lick_trace_rew/10)/sqrt(size(lick_trace_rew,1))], 'LineStyle', 'none');
    %alpha(.25);
    for i = 1:size(ts,1);
        subplot(1,3,2); errorbar(ts(i,:), avg_rewarded_roi(i,:), sm_rewarded(i,:), 'Color', colors(i,:)); hold on;
    end
    %subplot(1,3,2); errorbar(ts(i,:), mean(rewarded_roi_mean), sm_rewarded_avg, 'k'); hold on;
    
    if b_data.rewardDelayPercent >0
        assert(b_data.rewardDelayPercent == 100);
        rewardDelay = b_data.RewardDelayDurationMs + round(mean(cell2mat(b_data.reactTimesMs)));
        xlabel(['Time from cue onset (ms): reward at ' num2str(rewardDelay), '(ms)']);
        vert_lines(rewardDelay);
    elseif b_data.rewardDelayPercent ==0;
        xlabel('Time from reward and cue onset (ms)');
    end
    ylabel('dF/F');
    title(['Reward trials n=', num2str(size(rewarded_roi,1)), ' Cue duration = ', num2str(round(mean(cell2mat(b_data.reactTimesMs))))]);
    axis tight;
    hold off
    
    %plot reward omission trials
    use_ev_rew_omission = trial_outcome.rew_omission_trial_time;
    if length(use_ev_rew_omission) >2
        use_ev_rew_omission = use_ev_rew_omission(1:round(length(use_ev_rew_omission)/2));
    
        [rew_om_roi, use_times_rew_om, lick_trace_rew_om, lick_trace_rew_om_10ms] = trigger_movie_by_event_licks(tc_dfoverf, frame_info, licking_data, use_ev_rew_omission, pre_frames, post_frames);
        avg_rew_om_roi = squeeze(func(rew_om_roi,1));
        if cluster.num_cluster == 1
            avg_rew_om_roi = avg_rew_om_roi';
        end
        std_rew_om = squeeze(std(squeeze(rew_om_roi),1));
        sm_rewarded = std_rew_om./sqrt(size(rew_om_roi,1));
        rew_om_roi_mean = squeeze(mean(rew_om_roi,2));
        sm_om_avg = std(rew_om_roi_mean,1)/size(rew_om_roi_mean,1);
        for i = 1:cluster.num_cluster  %baseline each curve so it passes through zero
            shift = (-1)*avg_rew_om_roi(i,1);
            avg_rew_om_roi(i,:) = avg_rew_om_roi(i,:)+shift;
        end
        subplot(1,3,3); lickBars = bar(ts(1,:), mean(lick_trace_rew_om)/10); hold on;
        errorbar(ts(1,:), mean(lick_trace_rew_om)/10, [std(lick_trace_rew_om/10)/sqrt(size(lick_trace_rew_om,1))], 'LineStyle', 'none');
        %alpha(.25);
        for i = 1:size(ts,1);
            subplot(1,3,3); errorbar(ts(i,:), avg_rew_om_roi(i,:), sm_rewarded(i,:), 'Color', colors(i,:)); hold on;
        end
       % subplot(1,3,3); errorbar(ts(i,:), mean(rew_om_roi_mean), sm_om_avg, 'k'); hold on;

        if b_data.rewardDelayPercent >0 %set axes for reward delay condition
            assert(b_data.rewardDelayPercent == 100);
            rewardDelay = b_data.RewardDelayDurationMs + round(mean(cell2mat(b_data.reactTimesMs)));
            xlabel(['Time from cue onset (ms)']);
            vert_lines(rewardDelay);
        elseif b_data.rewardDelayPercent ==0;
            xlabel('Time from cue onset (ms)');
            rewardDelay = b_data.RewardDelayDurationMs;
        end
        ylabel('dF/F');
        title(['Reward omission trials n=', num2str(size(rew_om_roi,1))]);
        legend(['avg licks/frame'; cellstr(num2str([1:cluster.num_cluster]'))]);
        axis tight;
    end
    hold off
    
    %STANDARDIZE YLIMS, SAVE VARIABLES, REPORT Ns---------------------------------
    %set y scale to be equal for all 3 subplots----
    YL = [];
    for i =2:3
        subplot(1,3,i); YL(i-1,:) = ylim;
    end
    YL(YL ==1) = 0; %if there are zero trials in a catagory then it will default to values of -1 and 1
    YL(YL ==-1)= 0;
    
    subplot(1,3,2); ylim([min(YL(:,1)) max(YL(:,2))]);
    subplot(1,3,3); ylim([min(YL(:,1)) max(YL(:,2))]);
    disp(['day/animal: ' num2str(ROI_name)])
    disp(['reward delay of ', num2str(rewardDelay), 'ms on ', num2str(b_data.rewardDelayPercent), '% of trials'])
    disp(['# of rewarded trials = ' num2str(size(rewarded_roi,1))])
    if exist('rew_om_roi'); disp(['# of reward omission trials = ' num2str(size(rew_om_roi,1))]); end
    
    %Save figure before opening next one
    if save_data == 1
        destyFig = strcat(ANALYSIS_DIR, 'LeverFigureFolder\', days{kk}, '_fig');
        savefig([destyFig]);
    end
    
    %PLOT UNEXPECTED REWARD TRIALS -------------------------------------------------
    if b_data.rewardUnexpectPercent >0
        figure;
        %aligned to time of reward
        use_ev_unexp_rew = trial_outcome.unexpected_rew_time;
        disp(['# of unexpected reward trials ', num2str(length(use_ev_unexp_rew))]);
        [unexp_rew_roi, use_times_unexp_rew, lick_trace_unexp_rew, lick_trace_unexp_rew_10ms] = trigger_movie_by_event_licks(tc_dfoverf, frame_info, licking_data, use_ev_unexp_rew, pre_frames, post_frames);
        avg_unexp_rew_roi = squeeze(func(unexp_rew_roi,1));
        if cluster.num_cluster == 1
            avg_unexp_rew_roi = avg_unexp_rew_roi';
        end
        std_unexp_rew = squeeze(std(squeeze(unexp_rew_roi),1));
        sm_unexp_rew = std_unexp_rew./sqrt(size(unexp_rew_roi,1));
        for i = 1:cluster.num_cluster  %baseline each curve so it passes through zero
            shift = (-1)*avg_unexp_rew_roi(i,1);
            avg_unexp_rew_roi(i,:) = avg_unexp_rew_roi(i,:)+shift;
        end
        subplot(1,2,1); lickBars = bar(ts(1,:), mean(lick_trace_unexp_rew)/10); hold on
        alpha(.25);
        for i = 1:size(ts,1);
            errorbar(ts(i,:), avg_unexp_rew_roi(i,:), sm_unexp_rew(i,:), 'Color', colors(i,:)); hold on;
        end
        if b_data.rewardDelayPercent >0 %set axes for reward delay condition
            assert(b_data.rewardDelayPercent == 100);
            rewardDelay = b_data.RewardDelayDurationMs;
            xlabel(['reward delivery at ', num2str(rewardDelay), '(ms)']);
            vert_lines(rewardDelay);
        elseif b_data.rewardDelayPercent ==0;
            xlabel('Time from reward delivery (ms)');
        end
        ylabel('dF/F');
        suptitle('Unexpected reward trials');
        title(['reward aligned trials n=', num2str(size(unexp_rew_roi,1))]);
        legend(['avg licks/frame'; cellstr(num2str([1:cluster.num_cluster]'))]);
        axis tight;
        ylim([min(YL(:,1)) max(YL(:,2))]);
        hold off
        
        %align unexpected reward trials to first lick after reward delivery
        unexp_rew_lick_align = [];
        for i = 1:length(use_ev_unexp_rew)
            unexp_rew_lick_align = [unexp_rew_lick_align, find(lickTimes>use_ev_unexp_rew(i), 1, 'first')]; %gives the indeces within lickTimes for the 1st lick after the reward delivery
        end
        unexp_rew_lick_align = lickTimes(unexp_rew_lick_align);
        if ~isempty(unexp_rew_lick_align)
            cue_lick_diff = unexp_rew_lick_align - use_ev_unexp_rew;
        end
        [unexp_rew_roi_lick_align, use_times_unexp_rew_lick_align, lick_trace_unexp_rew_lick_align, lick_trace_unexp_rew_lick_align_10ms] = trigger_movie_by_event_licks(tc_dfoverf, frame_info, licking_data, unexp_rew_lick_align, pre_frames, post_frames);
        avg_unexp_rew_roi_lick_align = squeeze(func(unexp_rew_roi_lick_align,1));
        if cluster.num_cluster == 1
            avg_unexp_rew_roi_lick_align = avg_unexp_rew_roi_lick_align';
        end
        std_unexp_rew_lick_align = squeeze(std(squeeze(unexp_rew_roi_lick_align),1));
        sm_unexp_rew_lick_align = std_unexp_rew_lick_align./sqrt(size(unexp_rew_roi_lick_align,1));
        for i = 1:cluster.num_cluster  %baseline each curve so it passes through zero
            shift = (-1)*avg_unexp_rew_roi_lick_align(i,1);
            avg_unexp_rew_roi_lick_align(i,:) = avg_unexp_rew_roi_lick_align(i,:)+shift;
        end
        subplot(1,2,2); lickBars = bar(ts(1,:), mean(lick_trace_unexp_rew_lick_align)/10); hold on
        alpha(.25);
        for i = 1:size(ts,1);
            errorbar(ts(i,:), avg_unexp_rew_roi_lick_align(i,:), sm_unexp_rew_lick_align(i,:), 'Color', colors(i,:)); hold on;
        end
        xlabel('Time from first lick after reward delivery (ms)');
        ylabel('dF/F');
        suptitle('Unexpected reward trials');
        title(['licking aligned trials. avg time to first lick = ', num2str(mean(cue_lick_diff)), ' +/- ', num2str(std(cue_lick_diff))]);
        axis tight
        ylim([min(YL(:,1)) max(YL(:,2))]);
        hold off
        if save_data == 1
            destyFig = strcat(ANALYSIS_DIR, 'LeverFigureFolder\', days{kk}, '_unexp_fig');
            savefig([destyFig]);
        end
    end
    
    %SAVE matfiles  ----------------------------------------------
    destyNR = strcat(ANALYSIS_DIR, 'LeverSummaryFolder\', days{kk}, '_reward_trials');
    rewarded_roi = squeeze(rewarded_roi);
    %save([destyNR], 'rewarded_roi');
    licking_data.lick_trace_rew = lick_trace_rew;
    licking_data.lick_trace_rew_10ms = lick_trace_rew_10ms;
    if exist('rew_om_roi') && size(rew_om_roi,1) >2
        destyOR = strcat(ANALYSIS_DIR, 'LeverSummaryFolder\', days{kk}, '_rew_om_trials');
        rew_om_roi = squeeze(rew_om_roi);
        %save([destyOR], 'rew_om_roi');
        licking_data.lick_trace_rew_om = lick_trace_rew_om;
        licking_data.lick_trace_rew_om_10ms = lick_trace_rew_om_10ms;
    end
    if b_data.rewardUnexpectPercent >0
        destyUR = strcat(ANALYSIS_DIR, 'LeverSummaryFolder\', days{kk}, '_unexp_trials');
        unexp_rew_roi = squeeze(unexp_rew_roi);
        %save([destyUR], 'unexp_rew_roi');
        licking_data.lick_trace_unexp_rew = lick_trace_unexp_rew;
        licking_data.lick_trace_unexp_rew_10ms=lick_trace_unexp_rew_10ms;
    end
    
    %SAVE lick traces
    if save_data == 1
        save([ANALYSIS_DIR 'BxAndAnalysisOutputs\BxOutputs\', days{kk}, '_bx_outputs'], 'licking_data', '-append');
    end
end


