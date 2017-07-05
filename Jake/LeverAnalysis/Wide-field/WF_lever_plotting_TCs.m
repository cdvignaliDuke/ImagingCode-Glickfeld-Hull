%PLOTTING GRAPHS 
%plots correlation coefficient between the ROIs
%plots TCs of the df/f for success, earlies, tooFast, fidget and press.  
clear
BEHAVE_DIR = 'Z:\Data\WidefieldImaging\GCaMP\behavior\';
ANALYSIS_DIR ='Z:\Analysis\WF Lever Analysis\';
CLUSTER_DIR  ='Z:\Analysis\WF Lever Analysis\BxAndAnalysisOutputs\'; 
days = {'170628_img98', '170704_img98'}; % '170321_img86'
 
for kk=1:length(days)
    %set directories and load bxOutputs and cluter data. 
    ROI_name  =  days{kk};
    b_data = get_bx_data(BEHAVE_DIR, days{kk});  %find the correct behavior file and loads it.
    load([ANALYSIS_DIR 'BxAndAnalysisOutputs\BxOutputs\', days{kk}, '_bx_outputs']);
    load([CLUSTER_DIR, days{kk}, '\', days{kk}, '_cluster']); 
    
    %check to see if FakeMouse exists
    if isfield(b_data, 'doFakeMouseSuccessOnly');
        fake_mouse = b_data.doFakeMouseSuccessOnly;
    else
        fake_mouse = 0;
    end
    
    %define variables to help with plotting
    func = @median; %func = @mean; %func = @std;
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
    
    % PLOT SUCCESSFUL TRIALS----------------------------------------------
    if fake_mouse == 0;
        use_ev_success = trial_outcome.success_time;
        [success_roi, use_times_succ, lick_trace_succ, lick_trace_succ_10ms] = trigger_movie_by_event_licks(tc_dfoverf, frame_info, licking_data, use_ev_success, pre_frames, post_frames);
        avg_success_roi = squeeze(func(success_roi,1));
        if cluster.num_cluster == 1
            avg_success_roi = avg_success_roi';
        end
        std_success = squeeze(std(squeeze(success_roi),1));
        sm_success = std_success./sqrt(size(success_roi,1));
        for i = 1:cluster.num_cluster  %baseline each curve so it passes through zero
            shift = (-1)*avg_success_roi(i,1);
            avg_success_roi(i,:) = avg_success_roi(i,:)+shift;
        end
        subplot(2,3,2); lickBars = bar(ts(1,:), mean(lick_trace_succ)/10); hold on
        alpha(.25);
        for i = 1:size(ts,1);
            subplot(2,3,2); errorbar(ts(i,:), avg_success_roi(i,:), sm_success(i,:), 'Color', colors(i,:)); hold on;
        end
        xlabel('Time from release (ms)');
        ylabel('dF/F');
        title(['success n=', num2str(size(success_roi,1))]);
        axis tight;
        hold off
    end
    
    %PLOT FAILED TRIALS---------------------------------------------------
    if fake_mouse == 0;
        use_ev_fail = trial_outcome.early_time;
        [fail_roi, use_times_fail, lick_trace_fail, lick_trace_fail_10ms] = trigger_movie_by_event_licks(tc_dfoverf, frame_info, licking_data, use_ev_fail, pre_frames, post_frames);
        avg_fail_roi = squeeze(func(fail_roi,1));
        if cluster.num_cluster == 1
            avg_fail_roi = avg_fail_roi';
        end
        if size(fail_roi,1) == 1
            sm_fail = zeros(size(fail_roi,2), size(fail_roi,3));
        else
            std_fail = squeeze(std(squeeze(fail_roi),[],1));
            sm_fail = std_fail./sqrt(size(fail_roi,1));
        end
        for i = 1:cluster.num_cluster  %baseline each curve so it passes through zero
            shift = (-1)*avg_fail_roi(i,1);
            avg_fail_roi(i,:) = avg_fail_roi(i,:)+shift;
        end
        subplot(2,3,3); bar(ts(1,:), mean(lick_trace_fail)/10); hold on
        alpha(.25);
        for i = 1:size(ts,1);
            hold on; subplot(2,3,3); errorbar(ts(i,:), avg_fail_roi(i,:), sm_fail(i,:), 'Color', colors(i,:));
        end
        xlabel('Time from release (ms)');
        ylabel('dF/F');
        title(['fail n=' num2str(size(fail_roi,1))]);
        axis tight;
        legend(['avg licks/frame'; cellstr(num2str([1:cluster.num_cluster]'))]);
        hold off
    end
    
    %PLOT FIDGETS-----------------------------------------------------
    if fake_mouse == 0;
        use_ev_fidget = trial_outcome.fidget;
        [fidget_roi, use_times_fidget, lick_trace_fidget, lick_trace_fidget_10ms] = trigger_movie_by_event_licks(tc_dfoverf, frame_info, licking_data, use_ev_fidget, pre_frames, post_frames);
        avg_fidget_roi = squeeze(func(fidget_roi,1));
        if cluster.num_cluster == 1
            avg_fidget_roi = avg_fidget_roi';
        end
        std_fidget = squeeze(std(squeeze(fidget_roi),1));
        sm_fidget = std_fidget./sqrt(size(fidget_roi,1));
        for i = 1:cluster.num_cluster  %baseline each curve so it passes through zero
            shift = (-1)*avg_fidget_roi(i,1);
            avg_fidget_roi(i,:) = avg_fidget_roi(i,:)+shift;
        end
        subplot(2,3,4); bar(ts(1,:), mean(lick_trace_fidget)/10); hold on
        alpha(.25);
        for i = 1:size(ts,1);
            hold on; subplot(2,3,4); errorbar(ts(i,:), avg_fidget_roi(i,:), sm_fidget(i,:), 'Color', colors(i,:));
        end
        xlabel('Time from release (ms)');
        ylabel('dF/F');
        title(['fidget n=' num2str(size(fidget_roi,1))]);
        axis tight;
        hold off
    end
     
    %PLOT SUCCESS-FAIL------------------------------------------------
    if fake_mouse == 0;
        sub_sm = sqrt(sm_fail.^2+sm_success.^2);
        for i = 1:size(ts,1);
            hold on; subplot(2,3,5); errorbar(ts(i,:), avg_success_roi(i,:) - avg_fail_roi(i,:), sm_fail(i,:), 'Color', colors(i,:));
        end
        xlabel('Time from release (ms)');
        ylabel('dF/F');
        title('success - fail');
        axis tight;
        hold off
    end
     
    %PLOT TOOFAST SUCCESSES---------------------------------------
    if fake_mouse == 0;
        use_ev_tooFast = trial_outcome.tooFastCorrects;
        [tooFast_roi, use_times_tooFast, lick_trace_tooFast, lick_trace_tooFast_10ms] = trigger_movie_by_event_licks(tc_dfoverf, frame_info, licking_data, use_ev_tooFast, pre_frames, post_frames);
        avg_tooFast_roi = squeeze(func(tooFast_roi,1));
        if cluster.num_cluster == 1
            avg_tooFast_roi = avg_tooFast_roi';
        end
        std_tooFast = squeeze(std(squeeze(tooFast_roi),1));
        sm_tooFast = std_tooFast./sqrt(size(squeeze(tooFast_roi),1));
        if size(tooFast_roi,1) == 1
            for i = 1:cluster.num_cluster  %baseline each curve so it passes through zero
                shift = (-1)*avg_tooFast_roi(i,1);
                avg_tooFast_roi(i,:) = avg_tooFast_roi(i,:)+shift;
            end
            for i = 1:size(avg_tooFast_roi,1);
                hold on; subplot(2,3,6); plot(ts(1,:), avg_tooFast_roi(i,:), 'Color', colors(i,:));
                axis tight;
            end
        elseif size(tooFast_roi,1) == 0
            subplot(2,3,6); plot(ts(1,:), zeros(length(ts))); ylim([-0.01 0.01]);
        else
            subplot(2,3,6); bar(ts(1,:), mean(lick_trace_tooFast)/10); hold on
            alpha(.25);
            for i = 1:cluster.num_cluster  %baseline each curve so it passes through zero
                shift = (-1)*avg_tooFast_roi(i,1);
                avg_tooFast_roi(i,:) = avg_tooFast_roi(i,:)+shift;
            end
            for i = 1:size(avg_tooFast_roi,1);
                hold on; subplot(2,3,6); errorbar(ts(1,:), avg_tooFast_roi(i,:), sm_tooFast(i,:), 'Color', colors(i,:));
                axis tight;
            end
        end
        xlabel('Time from release (ms)');
        ylabel('dF/F');
        title(['tooFast success n=' num2str(size(tooFast_roi,1))]);
        axis tight;
        
        %set y scale to be equal for all 3 subplots----
        YL = [];
        for i =2:6
            subplot(2,3,i); YL(i-1,:) = ylim;
        end
        YL(YL ==1) = 0; %if there are zero trials in a catagory then it will default to values of -1 and 1
        YL(YL ==-1)= 0;
    end
    
    %STANDARDIZE YLIMS, SAVE VARIABLES, REPORT Ns---------------------------------
    if fake_mouse == 0;
        subplot(2,3,2); ylim([min(YL(:,1)) max(YL(:,2))]);
        subplot(2,3,3); ylim([min(YL(:,1)) max(YL(:,2))]);
        subplot(2,3,4); ylim([min(YL(:,1)) max(YL(:,2))]);
        subplot(2,3,5); ylim([min(YL(:,1)) max(YL(:,2))]);
        subplot(2,3,6); ylim([min(YL(:,1)) max(YL(:,2))]);
        disp(['day/animal: ' num2str(ROI_name)])
        disp(['# of successful trials = ' num2str(size(success_roi,1))])
        disp(['# of failed trials = ' num2str(size(fail_roi,1))])
        disp(['# of fidget trials = ' num2str(size(fidget_roi,1))])
        disp(['# of tooFast_successes = ' num2str(size(tooFast_roi,1))])
        if exist('lapse_roi')
            if length(size(lapse_roi))==3
                disp(['# of lapse trials = ' num2str(size(lapse_roi,1))])
            else
                disp(['# of lapse trials = 1']);
            end
        end
        
        %Save figure before opening next one
        destyFig = strcat(ANALYSIS_DIR, 'LeverFigureFolder\', days{kk}, '_fig');
        savefig([destyFig]);
    end

    
    %PLOT SUCCESSES TRIGGERED OFF CUE CHANGE. 
    if fake_mouse == 0;
        cueTimes = trial_outcome.change_orientation;
        corr_inx = find(trial_outcome.corr_inx);
        cueTimes = cueTimes .* (trial_outcome.corr_inx~=0);
        use_ev_cueChange = round(cueTimes(frame_info.f_frame_trial_num+1:frame_info.l_frame_trial_num-1));
        use_ev_cueChange(isnan(use_ev_cueChange) ) = [];
        use_ev_cueChange(use_ev_cueChange<0) = [];
        use_ev_cueChange(use_ev_cueChange==0) = [];
        [cue_roi, use_times_cue, lick_trace_cue, lick_trace_cue_10ms] = trigger_movie_by_event_licks(tc_dfoverf, frame_info, licking_data, use_ev_cueChange, pre_frames, post_frames);
        if size(cue_roi,1) == size(success_roi,1)+1;
            cue_roi = cue_roi([1:end-1],:,:);
            use_times_cue = use_times_cue([1:end-1]);
            lick_trace_cue = lick_trace_cue([1:end-1],:);
            lick_trace_cue_10ms = lick_trace_cue_10ms([1:end-1],:);
            assert(use_times_cue(1)<use_times_succ(1));
            if length(trial_outcome.success_time) == size(success_roi,1)+1;
                trial_outcome.success_time(end) =[];
                trial_outcome.succ_hold_dur(end) = [];
                save([ANALYSIS_DIR 'BxAndAnalysisOutputs\BxOutputs\', days{kk}, '_bx_outputs'], 'trial_outcome', '-append')
            end
        end
        avg_cue_roi = squeeze(func(cue_roi,1));
        if cluster.num_cluster == 1
            avg_cue_roi = avg_cue_roi';
        end
        std_cue = squeeze(std(squeeze(cue_roi),1));
        sm_cue = std_cue./sqrt(size(cue_roi,1));
        for i = 1:cluster.num_cluster  %baseline each curve so it passes through zero
            shift = (-1)*avg_cue_roi(i,3);
            avg_cue_roi(i,:) = avg_cue_roi(i,:)+shift;
        end
        figure; bar(ts(1,:), mean(lick_trace_cue)/10); hold on
        alpha(.25);
        for i = 1:size(ts,1);
            errorbar(ts(i,:), avg_cue_roi(i,:), sm_cue(i,:), 'Color', colors(i,:)); hold on;
        end
        xlabel('Time from cue (ms)');
        ylabel('dF/F');
        title([days{kk}, ' cue change n=', num2str(size(cue_roi,1))]);
        axis tight;
        ylim([min(YL(:,1)) max(YL(:,2))]);
        legend(['avg licks/frame'; cellstr(num2str([1:cluster.num_cluster]'))])
        destyFig2 = strcat(ANALYSIS_DIR, 'LeverFigureFolder\', days{kk}, '_cue');
        savefig([destyFig2]);
        hold off
    end
    
    %main figure saved before plotting cue triggered fig
    if fake_mouse == 0;
        success_roi = squeeze(success_roi);
        fail_roi = squeeze(fail_roi);
        fidget_roi = squeeze(fidget_roi);
        tooFast_roi = squeeze(tooFast_roi);
        cue_roi = squeeze(cue_roi);
        destySucc = strcat(ANALYSIS_DIR, 'LeverSummaryFolder\', days{kk}, '_success');
        destyFail = strcat(ANALYSIS_DIR, 'LeverSummaryFolder\', days{kk}, '_fail');
        destyFidget = strcat(ANALYSIS_DIR, 'LeverSummaryFolder\', days{kk}, '_fidget');
        destyTooFast = strcat(ANALYSIS_DIR, 'LeverSummaryFolder\', days{kk}, '_tooFast');
        destyCue = strcat(ANALYSIS_DIR, 'LeverSummaryFolder\', days{kk}, '_cue');
        save([destySucc], 'success_roi');
        save([destyFail], 'fail_roi');
        save([destyFidget], 'fidget_roi');
        save([destyTooFast], 'tooFast_roi');
        save([destyCue], 'cue_roi');
    end
    
    %save lapsed trials 
    if fake_mouse == 0;
        if isempty(trial_outcome.late_time)==0
            use_ev_lapse = trial_outcome.late_time;
            [lapse_roi, use_times_lapse, lick_trace_lapse, lick_trace_lapse_10ms] = trigger_movie_by_event_licks(tc_dfoverf, frame_info, licking_data, use_ev_lapse, pre_frames, post_frames);
            lapse_roi = squeeze(lapse_roi);
            destyLapse = strcat(ANALYSIS_DIR, 'LeverSummaryFolder\', days{kk}, '_lapse');
            save([destyLapse], 'lapse_roi');
        end
    end
    
    %save lick traces
    if fake_mouse == 0;
        licking_data.lick_trace_succ = lick_trace_succ;
        licking_data.lick_trace_succ_10ms = lick_trace_succ_10ms;
        licking_data.lick_trace_fail = lick_trace_fail;
        licking_data.lick_trace_fail_10ms = lick_trace_fail_10ms;
        licking_data.lick_trace_fidget = lick_trace_fidget;
        licking_data.lick_trace_fidget_10ms = lick_trace_fidget_10ms;
        licking_data.lick_trace_tooFast = lick_trace_tooFast;
        licking_data.lick_trace_tooFast_10ms = lick_trace_tooFast_10ms;
        licking_data.lick_trace_cue = lick_trace_cue;
        licking_data.lick_trace_cue_10ms = lick_trace_cue_10ms;
        save([ANALYSIS_DIR 'BxAndAnalysisOutputs\BxOutputs\', days{kk}, '_bx_outputs'], 'licking_data', '-append');
    end
    
    %----PLOT CUE REWARD PAIRING SESSIONS------------------------------------------
    if fake_mouse == 1;
        %plot rewarded trials
        use_ev_rewarded = trial_outcome.rewarded_trial_time;
        [rewarded_roi, use_times_rew, lick_trace_rew, lick_trace_rew_10ms] = trigger_movie_by_event_licks(tc_dfoverf, frame_info, licking_data, use_ev_rewarded, pre_frames, post_frames);
        avg_rewarded_roi = squeeze(func(rewarded_roi,1));
        if cluster.num_cluster == 1
            avg_rewarded_roi = avg_rewarded_roi';
        end
        std_rewarded = squeeze(std(squeeze(rewarded_roi),1));
        sm_rewarded = std_rewarded./sqrt(size(rewarded_roi,1));
        for i = 1:cluster.num_cluster  %baseline each curve so it passes through zero
            shift = (-1)*avg_rewarded_roi(i,1);
            avg_rewarded_roi(i,:) = avg_rewarded_roi(i,:)+shift;
        end
        subplot(1,3,2); lickBars = bar(ts(1,:), mean(lick_trace_rew)/10); hold on
        alpha(.25);
        for i = 1:size(ts,1);
            subplot(1,3,2); errorbar(ts(i,:), avg_rewarded_roi(i,:), sm_rewarded(i,:), 'Color', colors(i,:)); hold on;
        end
        if b_data.rewardDelayPercent >0
            assert(b_data.rewardDelayPercent == 100);
            rewardDelay = b_data.RewardDelayDurationMs; 
            xlabel(['Time from cue termination (ms): reward at ' num2str(rewardDelay), '(ms)']);
            vert_lines(rewardDelay);
        elseif b_data.rewardDelayPercent ==0; 
            xlabel('Time from reward and cue termination (ms)');
        end
        ylabel('dF/F');
        title(['Reward trials n=', num2str(size(rewarded_roi,1)), ' Cue duration = ', num2str(round(mean(cell2mat(b_data.reactTimesMs))))]);
        axis tight;
        hold off
        
        %plot reward omission trials
        use_ev_rew_omission = trial_outcome.rew_omission_trial_time;
        [rew_om_roi, use_times_rew_om, lick_trace_rew_om, lick_trace_rew_om_10ms] = trigger_movie_by_event_licks(tc_dfoverf, frame_info, licking_data, use_ev_rew_omission, pre_frames, post_frames);
        avg_rew_om_roi = squeeze(func(rew_om_roi,1));
        if cluster.num_cluster == 1
            avg_rew_om_roi = avg_rew_om_roi';
        end
        std_rew_om = squeeze(std(squeeze(rew_om_roi),1));
        sm_rewarded = std_rew_om./sqrt(size(rew_om_roi,1));
        for i = 1:cluster.num_cluster  %baseline each curve so it passes through zero
            shift = (-1)*avg_rew_om_roi(i,1);
            avg_rew_om_roi(i,:) = avg_rew_om_roi(i,:)+shift;
        end
        subplot(1,3,3); lickBars = bar(ts(1,:), mean(lick_trace_rew_om)/10); hold on
        alpha(.25);
        for i = 1:size(ts,1);
            subplot(1,3,3); errorbar(ts(i,:), avg_rew_om_roi(i,:), sm_rewarded(i,:), 'Color', colors(i,:)); hold on;
        end
        if b_data.rewardDelayPercent >0 %set axes for reward delay condition
            assert(b_data.rewardDelayPercent == 100);
            rewardDelay = b_data.RewardDelayDurationMs; 
            xlabel(['Time from cue termination (ms)']);
            vert_lines(rewardDelay);
        elseif b_data.rewardDelayPercent ==0; 
            xlabel('Time from cue termination (ms)');
            rewardDelay = b_data.RewardDelayDurationMs; 
        end
        ylabel('dF/F');
        title(['Reward omission trials n=', num2str(size(rew_om_roi,1))]);
        legend(['avg licks/frame'; cellstr(num2str([1:cluster.num_cluster]'))]);
        axis tight;
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
        disp(['# of reward omission trials = ' num2str(size(rew_om_roi,1))])
       
        %Save figure before opening next one
        destyFig = strcat(ANALYSIS_DIR, 'LeverFigureFolder\', days{kk}, '_fig');
        savefig([destyFig]);
        
        %PLOT UNEXPECTED REWARD TRIALS -------------------------------------------------
        if b_data.rewardUnexpectPercent >0
            figure;
            %aligned to time of reward
            use_ev_unexp_rew = trial_outcome.unexpected_rew_time;
            disp(['# of unexpected reward trials ', num2str(use_ev_unexp_rew)]);
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
            destyFig = strcat(ANALYSIS_DIR, 'LeverFigureFolder\', days{kk}, '_unexp_fig');
            savefig([destyFig]);
        end
        
        %SAVE matfiles  ----------------------------------------------
        rewarded_roi = squeeze(rewarded_roi);
        rew_om_roi = squeeze(rew_om_roi);
        
        destySucc = strcat(ANALYSIS_DIR, 'LeverSummaryFolder\', days{kk}, '_reward_trials');
        destyFail = strcat(ANALYSIS_DIR, 'LeverSummaryFolder\', days{kk}, '_rew_om_trials');
        save([destySucc], 'rewarded_roi');
        save([destyFail], 'rew_om_roi');
        if b_data.rewardUnexpectPercent >0
            unexp_rew_roi = squeeze(unexp_rew_roi);
            save([destyFail], 'unexp_rew_roi');
            destyUnexp = strcat(ANALYSIS_DIR, 'LeverSummaryFolder\', days{kk}, '_unexp_rew_trials');
        end
        
        %SAVE lick traces
        licking_data.lick_trace_rew = lick_trace_rew;
        licking_data.lick_trace_rew_10ms = lick_trace_rew_10ms;
        licking_data.lick_trace_rew_om = lick_trace_rew_om;
        licking_data.lick_trace_rew_om_10ms = lick_trace_rew_om_10ms;
        save([ANALYSIS_DIR 'BxAndAnalysisOutputs\BxOutputs\', days{kk}, '_bx_outputs'], 'licking_data', '-append');
    end
end


%  Moveed to the end of the script for storage. Should go insdie the forloop for ues. 
    %PLOT ALIGNED TO LEVER PRESS  
%     if lever.press(1)<lever.release(1)
%         holdTime = lever.release-lever.press(1:length(lever.release));  %sometimes there is one more press than there are releases.
%     else
%         firstRelease = find(lever.release>lever.press(1),1,'first');
%         holdTime = lever.release(firstRelease:end)-lever.press(1:length(lever.release)-1);
%     end
%     nonfidgets = find(holdTime>500);
%     use_ev_press = round(lever.press(nonfidgets));
%     [press_roi, use_times_press, lick_trace_press] = trigger_movie_by_event_licks(tc_dfoverf, frame_info, licksByFrame, use_ev_press, pre_frames, post_frames+35);
%     ts_press = (-500:100:4500);
%     
%     figure;
%     avg_press_roi = squeeze(func(press_roi,1));
%     std_success = squeeze(std(press_roi,1));
%     sm_press = std_success./sqrt(size(press_roi,1));
%     for i = 1:cluster.num_cluster  %baseline each curve so it passes through zero
%         shift = (-1)*avg_press_roi(i,3);
%         avg_press_roi(i,:) = avg_press_roi(i,:)+shift;
%     end
%     bar(ts_press, mean(lick_trace_press)/10); hold on
%     for i = 1:size(ts,1);
%         errorbar(ts_press, avg_press_roi(i,:), sm_press(i,:), 'Color', colors(i,:)); hold on;
%     end
%     xlabel('Time from press (ms)');
%     ylabel('dF/F');
%     title([days{kk} ' presses n=', num2str(length(use_ev_press))]);
%     axis tight;
     
%     press_roi = squeeze(press_roi);
%     destyPress = strcat(ANALYSIS_DIR, 'LeverSummaryNoFolder\', days{kk}, '_press');
%     save([destyPress], 'press_roi');
%     destyFig = strcat(ANALYSIS_DIR, 'LeverFigureFolder\', days{kk}, '_press');
%     savefig([destyFig]);


%====================================================================
%also moved to the end of the script for storage
    
%      %plot the corr coef of licking and df/f  as well as between ROIs 
%     allLeverTimes = round(sort([lever.press, lever.release]));
%     leverFrameNums = frame_info.counter(allLeverTimes); 
%     leverFramesExc = [];
%     for i = 1:length(leverFrameNums)
%         leverFramesExc = [leverFramesExc, (leverFrameNums(i)-3):(leverFrameNums(i)+5)];   %isolate windows 3frames before to 5 frames after any lever event in order to exclude those frames and lick values
%     end
%     leverFramesExc = unique(leverFramesExc); %get rid of repeat frames due to lever events in quick succession 
%     licksByFrame2 = licksByFrame;
%     tc_dfoverf2 = tc_dfoverf;
%     if length(tc_dfoverf2) > length(licksByFrame2); %ensure that dFoF TC and licking TC are the same length. 
%         tc_dfoverf2(:,(length(licksByFrame2)+1):end)=[];
%     elseif length(tc_dfoverf2) < length(licksByFrame2);
%         licksByFrame2((length(tc_dfoverf2)+1):end)=[];
%     end
%     leverFramesExc(find(leverFramesExc > length(tc_dfoverf2))) = [];
%     leverFramesExc(find(leverFramesExc < 1)) = [];
%     licksByFrame2(leverFramesExc) = [];
%     tc_dfoverf2(:, leverFramesExc) = [];
%     
%     cnames = cell([1,size(tc_dfoverf,1)+1]);
%     rnames = cell([1,size(tc_dfoverf,1)+1]);
%     coefMat = corrcoef([tc_dfoverf2; licksByFrame2]');
%     coefMat2 = coefMat(end,1:size(tc_dfoverf,1));
%     for i = 1:size(tc_dfoverf,1);
%         cnames{1,i} = ['ROI' mat2str(i)];
%         rnames{1,i} = ['ROI' mat2str(i)];
%     end
%     cnames{1,i+1} = ['licking'];
%     rnames{1,i+1} = ['licking'];
% %     f = figure;
% %     title(days{kk});
% %     t = uitable(f, 'Data', coefMat,...
% %         'ColumnName', cnames,...
% %         'RowName', rnames);
% %     destyFig = strcat(ANALYSIS_DIR, 'LeverFigureFolder\', days{kk}, '_table');
% %     %savefig([destyFig]);
% %     
% %     destySucc = strcat(ANALYSIS_DIR, 'CorrCoefSummary\', days{kk}, '_corrCoefLick');
%     %save([destySucc], 'coefMat2');