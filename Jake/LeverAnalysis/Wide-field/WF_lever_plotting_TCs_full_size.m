%PLOTTING GRAPHS 
%ALTERED VERSION OF WF_PLOTTER MAINLY CHANGING PATHNAMES TO TROUBLE SHOOT
%FULL SIZED MOVIE ANALYSIS
%plots correlation coefficient between the ROIs
%plots TCs of the df/f for success, earlies, tooFast, fidget and press.  
clear

BEHAVE_DIR = 'Z:\Data\WidefieldImaging\GCaMP\behavior\';
ANALYSIS_DIR ='Z:\Misc\full_size_analysis\';

days = {'151021_img29', '151022_img29', '151009_img30', '151011_img30', '151211_img32', '151212_img32', '160129_img35', '160131_img35', '160129_img36', '160131_img36', '160314_img38', '160315_img38', '160319_img41', '160320_img41', '160606_img46', '160722_img53', }; %'150718_img27', '150719_img27', '150716_img28', '150717_img28', \
days = {'160129_img35'};
for kk=1:length(days)
    ROI_name  =  days{kk};
    bfile = dir([BEHAVE_DIR 'data-*i9' days{kk}(end-1:end) '-' days{kk}(1:6) '*' ]);
    behave_dest = [BEHAVE_DIR bfile.name];
    b_data = load(behave_dest);
    load([ANALYSIS_DIR, days{kk}, '_bx_outputs']);
    % ---- do simple movie analysis.
    func = @median;
    %func = @mean;
    %func = @std;
    pre_frames = 5;
    post_frames = 10;
    sampling_rate = 1000/mode(diff(frame_times));
    
    %  --- first plot  clustering
    colors = [1,0,0; 0,1,0; 0,0,1; 0.5,0.5,0.5; 1,0,1; 1,1,0; 0,1,1]; %sets up the color scheme for plotting multiple ROIs with errorbar     
    f2 = figure;
    subplot(2,3,1); imagesc(avg_img);
    if b_data.input.doLever == 0;
        title(['CONTROL: no lever ' days{kk}]);
    else
        title(days{kk});
    end
    shading flat; hold on;
    for i=1:cluster.num_cluster
        line(  cluster.roi_position{i}(:,1),   cluster.roi_position{i}(:,2) ,'color', 'w', 'linewidth', 2)
        text(mean(cluster.roi_position{i}(:,1)),mean(cluster.roi_position{i}(:,2)), ...
            [num2str(i)], 'color', 'k', 'FontSize', 30);
    end  
    
    ts = (-pre_frames:post_frames)*1000/double(round(sampling_rate));
    tot_frame = pre_frames + post_frames+1;
    use_ev_success = trial_outcome.success_time;
    
    %bin licking times according to frame number.
    frameTimes = frame_info.times; %frame times is monotonic by the time it gets here.    
    licksByFrame = [];
    for i = 1:length(frameTimes);
        if i == 1; %if this is the first frameTime 
            licksThisFrame = sum(lickTimes>=1 & lickTimes<=frameTimes(i));
            licksByFrame = [licksByFrame licksThisFrame];
        else
            licksThisFrame = sum(lickTimes>frameTimes(i-1) & lickTimes<=frameTimes(i));
            licksByFrame = [licksByFrame licksThisFrame];
        end
    end
    licking_data = [];
    licking_data.lickTimes = lickTimes;
    licking_data.licksByFrame = licksByFrame;
    
    % PLOT SUCCESSFUL TRIALS----------------------------------------------
    time_before = 500; % in ms, time before event w/o release
    time_after = 1000; % in ms, time after event w/o press
    %success_roi = trigger_movie_by_event(tc_dfoverf, frame_info, ...
    % use_ev_success, pre_frames, post_frames);
    [success_roi, use_times_succ, lick_trace_succ, lick_trace_succ_10ms] = trigger_movie_by_event_licks(tc_dfoverf, frame_info, licking_data, use_ev_success, pre_frames, post_frames);
    ts = repmat(ts,[cluster.num_cluster 1]);
    avg_success_roi = squeeze(func(success_roi,1));
    if cluster.num_cluster == 1
        avg_success_roi = avg_success_roi';
    end
    std_success = squeeze(std(squeeze(success_roi),1));
    sm_success = std_success./sqrt(size(success_roi,1));
    for i = 1:cluster.num_cluster  %baseline each curve so it passes through zero
        shift = (-1)*avg_success_roi(i,3);
        avg_success_roi(i,:) = avg_success_roi(i,:)+shift;
    end
    subplot(2,3,2); bar(ts(1,:), mean(lick_trace_succ)/10); hold on
    for i = 1:size(ts,1);
        subplot(2,3,2); errorbar(ts(i,:), avg_success_roi(i,:), sm_success(i,:), 'Color', colors(i,:)); hold on;
    end
    xlabel('Time from release (ms)');
    ylabel('dF/F');
    title(['success n=', num2str(size(success_roi,1))]);
    axis tight;
    hold off
    
    %PLOT FAILED TRIALS---------------------------------------------------
    use_ev_fail = trial_outcome.early_time;
    [fail_roi, use_times_fail, lick_trace_fail, lick_trace_fail_10ms] = trigger_movie_by_event_licks(tc_dfoverf, frame_info, licking_data, use_ev_fail, pre_frames, post_frames);
    %fail_roi = trigger_movie_by_event(tc_dfoverf, frame_info, ...
    %use_ev_fail, pre_frames, post_frames);
    avg_fail_roi = squeeze(func(fail_roi,1));
    if cluster.num_cluster == 1
        avg_fail_roi = avg_fail_roi';
    end
    std_fail = squeeze(std(squeeze(fail_roi),1));
    sm_fail = std_fail./sqrt(size(fail_roi,1));
    for i = 1:cluster.num_cluster  %baseline each curve so it passes through zero
        shift = (-1)*avg_fail_roi(i,3);
        avg_fail_roi(i,:) = avg_fail_roi(i,:)+shift;
    end
    subplot(2,3,3); bar(ts(1,:), mean(lick_trace_fail)/10); hold on
    for i = 1:size(ts,1);
        hold on; subplot(2,3,3); errorbar(ts(i,:), avg_fail_roi(i,:), sm_fail(i,:), 'Color', colors(i,:));
    end
    xlabel('Time from release (ms)');
    ylabel('dF/F');
    title(['fail n=' num2str(size(fail_roi,1))]);
    axis tight;
    legend(['avg licks/frame'; cellstr(num2str([1:cluster.num_cluster]'))]);
    hold off
    
    %PLOT FIDGETS-----------------------------------------------------
    use_ev_fidget = trial_outcome.fidget;
    [fidget_roi, use_times_fidget, lick_trace_fidget, lick_trace_fidget_10ms] = trigger_movie_by_event_licks(tc_dfoverf, frame_info, licking_data, use_ev_fidget, pre_frames, post_frames);
    %fail_roi = trigger_movie_by_event(tc_dfoverf, frame_info, ...
    %use_ev_fail, pre_frames, post_frames);
    avg_fidget_roi = squeeze(func(fidget_roi,1));
    if cluster.num_cluster == 1
        avg_fidget_roi = avg_fidget_roi';
    end
    std_fidget = squeeze(std(squeeze(fidget_roi),1));
    sm_fidget = std_fidget./sqrt(size(fidget_roi,1));
    for i = 1:cluster.num_cluster  %baseline each curve so it passes through zero
        shift = (-1)*avg_fidget_roi(i,3);
        avg_fidget_roi(i,:) = avg_fidget_roi(i,:)+shift;
    end
    subplot(2,3,4); bar(ts(1,:), mean(lick_trace_fidget)/10); hold on
    for i = 1:size(ts,1);
        hold on; subplot(2,3,4); errorbar(ts(i,:), avg_fidget_roi(i,:), sm_fidget(i,:), 'Color', colors(i,:));
    end
    xlabel('Time from release (ms)');
    ylabel('dF/F');
    title(['fidget n=' num2str(size(fidget_roi,1))]);
    axis tight;
    hold off
    
    %PLOT SUCCESS-FAIL------------------------------------------------
    sub_sm = sqrt(sm_fail.^2+sm_success.^2);
    for i = 1:size(ts,1);
        hold on; subplot(2,3,5); errorbar(ts(i,:), avg_success_roi(i,:) - avg_fail_roi(i,:), sm_fail(i,:), 'Color', colors(i,:));
    end
    xlabel('Time from release (ms)');
    ylabel('dF/F');
    title('success - fail');
    axis tight;
    hold off
    
    %PLOT TOOFAST SUCCESSES---------------------------------------
    use_ev_tooFast = trial_outcome.tooFastCorrects;
    [tooFast_roi, use_times_tooFast, lick_trace_tooFast, lick_trace_tooFast_10ms] = trigger_movie_by_event_licks(tc_dfoverf, frame_info, licking_data, use_ev_tooFast, pre_frames, post_frames);
    %tooFast_roi = trigger_movie_by_event(tc_dfoverf, frame_info, ...
    %use_ev_tooFast, pre_frames, post_frames);
    avg_tooFast_roi = squeeze(func(tooFast_roi,1));
    if cluster.num_cluster == 1
        avg_tooFast_roi = avg_tooFast_roi';
    end
    std_tooFast = squeeze(std(squeeze(tooFast_roi),1));
    sm_tooFast = std_tooFast./sqrt(size(squeeze(tooFast_roi),1));
    if size(tooFast_roi,1) == 1
        for i = 1:cluster.num_cluster  %baseline each curve so it passes through zero
            avg_tooFast_roi = avg_tooFast_roi;
            shift = (-1)*avg_tooFast_roi(i,3);
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
        for i = 1:cluster.num_cluster  %baseline each curve so it passes through zero
            shift = (-1)*avg_tooFast_roi(i,3);
            avg_tooFast_roi(i,:) = avg_tooFast_roi(i,:)+shift;
        end
        for i = 1:size(avg_tooFast_roi,1);
            hold on; subplot(2,3,6); errorbar(ts(1,:), avg_tooFast_roi(i,:), sm_tooFast(i,:), 'Color', colors(i,:));
            axis tight;
        end
    end
    %         for i = 1:size(ts,1);
    %             subplot(2,3,5); errorbar(ts(i,:), avg_tooFast_roi(i,:), sm_tooFast(i,:), 'Color', colors(i,:)); hold on;
    %         end
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
    
    
    %STANDARDIZE YLIMS, SAVE VARIABLES, REPORT Ns---------------------------------
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
    destyFig = strcat(ANALYSIS_DIR, days{kk}, '_fig');
    savefig([destyFig]);
    
    %PLOT SUCCESSES TRIGGERED OFF CUE CHANGE. 
    cueTimes = trial_outcome.change_orientation;
    corr_inx = find(trial_outcome.corr_inx);
    cueTimes = cueTimes .* (trial_outcome.corr_inx~=0);
    use_ev_cueChange = round(cueTimes(frame_info.f_frame_trial_num+1:frame_info.l_frame_trial_num-1));
    use_ev_cueChange(isnan(use_ev_cueChange) ) = [];
    use_ev_cueChange(use_ev_cueChange<0) = [];
    [cue_roi, use_times_cue, lick_trace_cue, lick_trace_cue_10ms] = trigger_movie_by_event_licks(tc_dfoverf, frame_info, licking_data, use_ev_cueChange, pre_frames, post_frames);
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
    for i = 1:size(ts,1);
        errorbar(ts(i,:), avg_cue_roi(i,:), sm_cue(i,:), 'Color', colors(i,:)); hold on;
    end
    xlabel('Time from cue (ms)');
    ylabel('dF/F');
    title([days{kk}, ' cue change n=', num2str(size(cue_roi,1))]);
    axis tight;
    ylim([min(YL(:,1)) max(YL(:,2))]);
    destyFig2 = strcat(ANALYSIS_DIR, days{kk}, '_cue');
    savefig([destyFig2]);
    hold off
    
    %main figure saved before plotting cue triggered fig
    success_roi = squeeze(success_roi);
    fail_roi = squeeze(fail_roi);
    fidget_roi = squeeze(fidget_roi);
    tooFast_roi = squeeze(tooFast_roi);
    cue_roi = squeeze(cue_roi);
    destySucc = strcat(ANALYSIS_DIR, days{kk}, '_success');
    destyFail = strcat(ANALYSIS_DIR, days{kk}, '_fail');
    destyFidget = strcat(ANALYSIS_DIR, days{kk}, '_fidget');
    destyTooFast = strcat(ANALYSIS_DIR, days{kk}, '_tooFast');
    destyCue = strcat(ANALYSIS_DIR, days{kk}, '_cue');
%     save([destySucc], 'success_roi');
%     save([destyFail], 'fail_roi');
%     save([destyFidget], 'fidget_roi');
%     save([destyTooFast], 'tooFast_roi');
%     save([destyCue], 'cue_roi');


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
    
     %plot the corr coef of licking and df/f  as well as between ROIs 
    allLeverTimes = round(sort([lever.press, lever.release]));
    leverFrameNums = frame_info.counter(allLeverTimes); 
    leverFramesExc = [];
    for i = 1:length(leverFrameNums)
        leverFramesExc = [leverFramesExc, (leverFrameNums(i)-3):(leverFrameNums(i)+5)];   %isolate windows 3frames before to 5 frames after any lever event in order to exclude those frames and lick values
    end
    leverFramesExc = unique(leverFramesExc); %get rid of repeat frames due to lever events in quick succession 
    licksByFrame2 = licksByFrame;
    tc_dfoverf2 = tc_dfoverf;
    if length(tc_dfoverf2) > length(licksByFrame2); %ensure that dFoF TC and licking TC are the same length. 
        tc_dfoverf2(:,(length(licksByFrame2)+1):end)=[];
    elseif length(tc_dfoverf2) < length(licksByFrame2);
        licksByFrame2((length(tc_dfoverf2)+1):end)=[];
    end
    leverFramesExc(find(leverFramesExc > length(tc_dfoverf2))) = [];
    leverFramesExc(find(leverFramesExc < 1)) = [];
    licksByFrame2(leverFramesExc) = [];
    tc_dfoverf2(:, leverFramesExc) = [];
    
    cnames = cell([1,size(tc_dfoverf,1)+1]);
    rnames = cell([1,size(tc_dfoverf,1)+1]);
    coefMat = corrcoef([tc_dfoverf2; licksByFrame2]');
    coefMat2 = coefMat(end,1:size(tc_dfoverf,1));
    for i = 1:size(tc_dfoverf,1);
        cnames{1,i} = ['ROI' mat2str(i)];
        rnames{1,i} = ['ROI' mat2str(i)];
    end
    cnames{1,i+1} = ['licking'];
    rnames{1,i+1} = ['licking'];
%     f = figure;
%     title(days{kk});
%     t = uitable(f, 'Data', coefMat,...
%         'ColumnName', cnames,...
%         'RowName', rnames);
%     destyFig = strcat(ANALYSIS_DIR, 'LeverFigureFolder\', days{kk}, '_table');
%     %savefig([destyFig]);
%     
%     destySucc = strcat(ANALYSIS_DIR, 'CorrCoefSummary\', days{kk}, '_corrCoefLick');
    %save([destySucc], 'coefMat2');
    
    %save lapsed trials 
    if isempty(trial_outcome.late_time)==0
        use_ev_lapse = trial_outcome.late_time;
        [lapse_roi, use_times_lapse, lick_trace_lapse, lick_trace_lapse_10ms] = trigger_movie_by_event_licks(tc_dfoverf, frame_info, licking_data, use_ev_lapse, pre_frames, post_frames);
        lapse_roi = squeeze(lapse_roi);
        destyLapse = strcat(ANALYSIS_DIR, days{kk}, '_lapse');
        %save([destyLapse], 'lapse_roi');
    end
    
    %save lick traces
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
    save([ANALYSIS_DIR days{kk}, '_bx_outputs'], 'licking_data', '-append');
end