%%plot all calclium traces on one graph
%PLOTTING GRAPHS 
clear

DATA_DIR =  'C:\Users\jake\TempData\';
FRAME_TIME_DIR = 'C:\Users\jake\TempData\';
BEHAVE_DIR = 'Z:\Data\WidefieldImaging\GCaMP\behavior\';
ANALYSIS_DIR ='Z:\Analysis\LeverAnalysis\';

days = {'150717_img28', '151021_img29', '160131_img36', '160131_img35', '160315_img38', '160320_img41'}; %rand=1000 
%days = {'160129_img36', '151009_img30', '151011_img30', '160314_img38'};  %rand=4500

%days = {'150716_img28', '150717_img28', '151021_img29', '151022_img29', '151009_img30', '151011_img30', '151211_img32', '151212_img32', '160129_img35', '160131_img35', '160129_img36','160131_img36', '160314_img38', '160315_img38', '160319_img41', '160320_img41'}; %'150718_img27', '150719_img27',
colors = {'r', 'r', 'b', 'b', 'b', 'b', 'm', 'm', 'g', 'g', 'k', 'k', 'c', 'c', 'y', 'y', 'r', 'r', 'b', 'b'};
ROIcell = {[1], [1:5], [2], [2], [1:3], [1,3], [1:2], [1:2], [1:2], [1:2], [1:2], [1:2], [3:6], [2,3,5], [4], [1:2]};
f2 = figure;
for kk=1:length(days)
    ROI_name  =  days{kk};
    currROI = days_roi_matcher.(strcat('i', days{kk}));
    %currROI = cell2mat(ROIcell(kk));
    bfile = dir([BEHAVE_DIR 'data-*i9' days{kk}(end-1:end) '-' days{kk}(1:6) '*' ]);
    behave_dest = [BEHAVE_DIR bfile.name];
    b_data = load(behave_dest);
    load([ANALYSIS_DIR 'BxAndAnalysisOutputs\BxOutputs\', days{kk}, '_bx_outputs']);
    % ---- do simple movie analysis.
    func = @median;
    %func = @mean;
    %func = @std;
    pre_frames = 5;
    post_frames = 10;
    
    ts = (-pre_frames:post_frames)*1000/round(sampling_rate);
    tot_frame = pre_frames + post_frames+1;
    use_ev_success = trial_outcome.success_time;
    
    %bin licking times according to frame number
    frameTimes = frame_info.times; %start by adjusting frameTimes so it accounts for missing frames
    frameMisses = find(diff(frameTimes)>150);
    frameShift=0;
    for i = frameMisses;
        missedFrame = ((frameTimes(i+frameShift+1)-frameTimes(i+frameShift))/2) + frameTimes(i+frameShift); %approximate the time of the missed frame 
        frameTimes = [frameTimes(1:i+frameShift), missedFrame, frameTimes(i+frameShift+1:end)];
        frameShift=frameShift+1; %necessary in order to account for the fact that frames are being inserted into the TC as the forloop progresses
    end
    licksByFrame = [];
    for i = 1:length(frameTimes);
        if i == length(frameTimes);
            licksThisFrame = sum(lickTimes>=frameTimes(i) & lickTimes<frameTimes(i)+diff(fliplr(frameTimes(i-1:i))));
            licksByFrame = [licksByFrame licksThisFrame];
            break
        end
        licksThisFrame = sum(lickTimes>=frameTimes(i) & lickTimes<frameTimes(i+1));
        licksByFrame = [licksByFrame licksThisFrame];
    end
    
    % PLOT SUCCESSFUL TRIALS----------------------------------------------
    time_before = 500; % in ms, time before event w/o release
    time_after =1000; % in ms, time after event w/o press
    [success_roi, use_times_succ, lick_trace_succ] = trigger_movie_by_event_licks(tc_dfoverf, frame_info, licksByFrame, use_ev_success, pre_frames, post_frames);
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
    for i = currROI;
        subplot(2,1,1); errorbar(ts(i,:), avg_success_roi(i,:), sm_success(i,:), 'Color', colors{kk}); hold on;
    end
   
    
    %PLOT FAILED TRIALS---------------------------------------------------
    use_ev_fail = trial_outcome.early_time;
    [fail_roi, use_times_fail, lick_trace_fail] = trigger_movie_by_event_licks(tc_dfoverf, frame_info, licksByFrame, use_ev_fail, pre_frames, post_frames);
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
    for i = currROI;
        subplot(2,1,2); errorbar(ts(i,:), avg_fail_roi(i,:), sm_fail(i,:), 'Color', colors{kk}); hold on; 
    end
   
    
%     
%     %PLOT FIDGETS-----------------------------------------------------
%     use_ev_fidget = trial_outcome.fidget;
%     [fidget_roi, use_times_fidget, lick_trace_fidget] = trigger_movie_by_event_licks(tc_dfoverf, frame_info, licksByFrame, use_ev_fidget, pre_frames, post_frames);
%     avg_fidget_roi = squeeze(func(fidget_roi,1));
%     if cluster.num_cluster == 1
%         avg_fidget_roi = avg_fidget_roi';
%     end
%     std_fidget = squeeze(std(squeeze(fidget_roi),1));
%     sm_fidget = std_fidget./sqrt(size(fidget_roi,1));
%     for i = 1:cluster.num_cluster  %baseline each curve so it passes through zero
%         shift = (-1)*avg_fidget_roi(i,3);
%         avg_fidget_roi(i,:) = avg_fidget_roi(i,:)+shift;
%     end
%     subplot(2,3,4); bar(ts(1,:), mean(lick_trace_fidget)/10); hold on
%     for i = 1:size(ts,1);
%         hold on; subplot(2,3,4); errorbar(ts(i,:), avg_fidget_roi(i,:), sm_fidget(i,:), 'Color', colors(i,:));
%     end
%     xlabel('Time from release (ms)');
%     ylabel('dF/F');
%     title(['fidget n=' num2str(size(fidget_roi,1))]);
%     axis tight;
%     hold off
%     
%     %PLOT SUCCESS-FAIL------------------------------------------------
%     sub_sm = sqrt(sm_fail.^2+sm_success.^2);
%     for i = 1:size(ts,1);
%         hold on; subplot(2,3,5); errorbar(ts(i,:), avg_success_roi(i,:) - avg_fail_roi(i,:), sm_fail(i,:), 'Color', colors(i,:));
%     end
%     xlabel('Time from release (ms)');
%     ylabel('dF/F');
%     title('success - fail');
%     axis tight;
%     hold off
%     
%     %PLOT TOOFAST SUCCESSES---------------------------------------
%     use_ev_tooFast = trial_outcome.tooFastCorrects;
%     [tooFast_roi, use_times_tooFast, lick_trace_tooFast] = trigger_movie_by_event_licks(tc_dfoverf, frame_info, licksByFrame, use_ev_tooFast, pre_frames, post_frames);
%     avg_tooFast_roi = squeeze(func(tooFast_roi,1));
%     if cluster.num_cluster == 1
%         avg_tooFast_roi = avg_tooFast_roi';
%     end
%     std_tooFast = squeeze(std(squeeze(tooFast_roi),1));
%     sm_tooFast = std_tooFast./sqrt(size(squeeze(tooFast_roi),1));
%     if size(tooFast_roi,1) == 1
%         for i = 1:cluster.num_cluster  %baseline each curve so it passes through zero
%             avg_tooFast_roi = avg_tooFast_roi;
%             shift = (-1)*avg_tooFast_roi(i,3);
%             avg_tooFast_roi(i,:) = avg_tooFast_roi(i,:)+shift;
%         end
%         for i = 1:size(avg_tooFast_roi,1);
%             hold on; subplot(2,3,6); plot(ts(1,:), avg_tooFast_roi(i,:), 'Color', colors(i,:));
%             axis tight;
%         end
%     elseif size(tooFast_roi,1) == 0
%         subplot(2,3,6); plot(ts(1,:), zeros(length(ts))); ylim([-0.01 0.01]);
%     else
%         subplot(2,3,6); bar(ts(1,:), mean(lick_trace_tooFast)/10); hold on
%         for i = 1:cluster.num_cluster  %baseline each curve so it passes through zero
%             shift = (-1)*avg_tooFast_roi(i,3);
%             avg_tooFast_roi(i,:) = avg_tooFast_roi(i,:)+shift;
%         end
%         for i = 1:size(avg_tooFast_roi,1);
%             hold on; subplot(2,3,6); errorbar(ts(1,:), avg_tooFast_roi(i,:), sm_tooFast(i,:), 'Color', colors(i,:));
%             axis tight;
%         end
%     end
%     xlabel('Time from release (ms)');
%     ylabel('dF/F');
%     title(['tooFast success n=' num2str(size(tooFast_roi,1))]);
%     axis tight;
    
    %set y scale to be equal for all 3 subplots----
%     YL = [];
%     for i =2:6
%         subplot(2,3,i); YL(i-1,:) = ylim;
%     end
%     
    %STANDARDIZE YLIMS, SAVE VARIABLES, REPORT Ns---------------------------------
%     subplot(2,3,2); ylim([min(YL(:,1)) max(YL(:,2))]);
%     subplot(2,3,3); ylim([min(YL(:,1)) max(YL(:,2))]);
%     subplot(2,3,4); ylim([min(YL(:,1)) max(YL(:,2))]); 
%     subplot(2,3,5); ylim([min(YL(:,1)) max(YL(:,2))]); 
%     subplot(2,3,6); ylim([min(YL(:,1)) max(YL(:,2))]);
     
%     success_roi = squeeze(success_roi);
%     fail_roi = squeeze(fail_roi);
%     fidget_roi = squeeze(fidget_roi);
%     destySucc = strcat(ANALYSIS_DIR, 'LeverSummaryFolder\', days{kk}, '_success');
%     destyFail = strcat(ANALYSIS_DIR, 'LeverSummaryFolder\', days{kk}, '_fail');
%     destyFidget = strcat(ANALYSIS_DIR, 'LeverSummaryFolder\', days{kk}, '_fidget');
%     destyFig = strcat(ANALYSIS_DIR, 'LeverFigureFolder\', days{kk}, '_fig');
%     save([destySucc], 'success_roi');
%     save([destyFail], 'fail_roi');
%     save([destyFidget], 'fidget_roi');
%     savefig([destyFig]);
%     
%     disp(['day/animal: ' num2str(ROI_name)])
%     disp(['# of successful trials = ' num2str(size(success_roi,1))])
%     disp(['# of failed trials = ' num2str(size(fail_roi,1))])
%     disp(['# of fidget trials = ' num2str(size(fidget_roi,1))])
%     disp(['# of tooFast_successes = ' num2str(size(tooFast_roi,1))])

end
subplot(2,1,1);
xlabel('Time from release (ms)');
ylabel('dF/F');
title('correct');
axis tight;

subplot(2,1,2);
xlabel('Time from release (ms)');
ylabel('dF/F');
title('early');
axis tight;

YL = [];
for i =1:2
    subplot(2,1,i); YL(i,:) = ylim;
end
subplot(2,1,1); ylim([min(YL(:,1)) max(YL(:,2))]);
subplot(2,1,2); ylim([min(YL(:,1)) max(YL(:,2))]);