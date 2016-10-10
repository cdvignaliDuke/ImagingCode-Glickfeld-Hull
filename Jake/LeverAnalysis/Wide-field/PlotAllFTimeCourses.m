%%plot all calclium traces on one graph
%PLOTTING GRAPHS 
clear
%close all 

%days = {'150717_img28', '151021_img29', '160131_img36', '160131_img35', '160315_img38', '160320_img41'}; %rand=1000 
%days = {'160129_img36', '151009_img30', '151011_img30', '160314_img38'};  %rand=4500
days = {'150716_img28', '150717_img28', '151021_img29', '151022_img29', '151009_img30', '151011_img30', '151211_img32', '151212_img32', '160129_img35', '160131_img35', '160129_img36','160131_img36', '160314_img38', '160315_img38', '160319_img41', '160320_img41', '160606_img46'}; %'150718_img27', '150719_img27',
colors = {'r', 'r', 'b', 'b', 'b', 'b', 'm', 'm', 'g', 'g', 'k', 'k', 'c', 'c', 'y', 'y', 'r', 'r', 'b', 'b', 'k'};
ROIcell = {[1], [1:5], [2], [2], [1:3], [1,3], [1:2], [1:2], [1:2], [1:2], [1:2], [1:2], [3:6], [2,3,5], [4], [1:2], [3:4]};

DATA_DIR =  'C:\Users\jake\TempData\';
FRAME_TIME_DIR = 'C:\Users\jake\TempData\';
BEHAVE_DIR = 'Z:\Data\WidefieldImaging\GCaMP\behavior\';
ANALYSIS_DIR ='Z:\Analysis\LeverAnalysis\';

% ---- do simple movie analysis.
func = @mean;
%func = @median;
%func = @std;
pre_frames = 5;
post_frames = 10;
%f1 = figure(1);
%f2 = figure(2); 
figure; hold on; %subplot(1,3,1);
title('cyan');
days_roi_matcher = struct();
ind_peak = [];
ind_peak_fail = [];
corr_roi_max = [];
fail_roi_max = [];
for kk=1:length(days)
    ROI_name  =  days{kk};
    currROI = ROIcell{kk};
    %currROI = days_roi_matcher.(strcat('i', days{kk}));
    bfile = dir([BEHAVE_DIR 'data-*i9' days{kk}(end-1:end) '-' days{kk}(1:6) '*' ]);
    behave_dest = [BEHAVE_DIR bfile.name];
    b_data = load(behave_dest);
    load(strcat(ANALYSIS_DIR, 'BxAndAnalysisOutputs\BxOutputs\', days{kk}, '_bx_outputs'));
    num_ROIs = cluster.num_cluster;
    ts = (-pre_frames:post_frames)*1000/round(sampling_rate);
    tot_frame = pre_frames + post_frames+1;
    
    % PLOT SUCCESSFUL TRIALS----------------------------------------------
    time_before = 500; % in ms, time before event w/o release
    time_after =1000; % in ms, time after event w/o press
    load(strcat(ANALYSIS_DIR, 'LeverSummaryFolder\', days{kk}, '_success'));
    ts = repmat(ts,[cluster.num_cluster 1]);
    avg_success_roi = squeeze(func(success_roi,1));
    if num_ROIs == 1
        avg_success_roi = avg_success_roi';
    end
    std_success = squeeze(std(squeeze(success_roi),1));
    sm_success = std_success./sqrt(size(success_roi,1));
    for i = 1:num_ROIs  %baseline each curve so it passes through zero
        shift = (-1)*avg_success_roi(i,3);
        avg_success_roi(i,:) = avg_success_roi(i,:)+shift;
    end
%     for i = currROI;
%         errorbar(ts(i,:), avg_success_roi(i,:), sm_success(i,:), 'g'); 
%         %subplot(1,3,1); hold on; errorbar(ts(i,:), avg_success_roi(i,:), sm_success(i,:), 'Color', colors{kk}); 
%     end
    
%     figure; hold on   %plots all trials from a single animal on the same plot
%     for iii = 1:size(success_roi,1)
%         plot(squeeze(success_roi(iii,ROIcell{kk}(1),:)));
%         title(days(kk));
%     end

    %PLOT FAILED TRIALS---------------------------------------------------
    load(strcat(ANALYSIS_DIR, 'LeverSummaryFolder\', days{kk}, '_fail'));
    avg_fail_roi = squeeze(func(fail_roi,1));
    if num_ROIs == 1
        avg_fail_roi = avg_fail_roi';
    end
    std_fail = squeeze(std(squeeze(fail_roi),1));
    sm_fail = std_fail./sqrt(size(fail_roi,1));
    for i = 1:num_ROIs  %baseline each curve so it passes through zero
        shift = (-1)*avg_fail_roi(i,3);
        avg_fail_roi(i,:) = avg_fail_roi(i,:)+shift;
    end
%     for i = currROI;
%         errorbar(ts(i,:), avg_fail_roi(i,:), sm_fail(i,:), 'r');
%         %subplot(1,3,2); hold on; errorbar(ts(i,:), avg_fail_roi(i,:), sm_fail(i,:), 'Color', colors{kk});
%     end
   
    %PLOT CUE TRIGGERED TRACES---------------------------------------------
%     load(strcat(ANALYSIS_DIR, 'LeverSummaryFolder\', days{kk}, '_cue'));
%     avg_cue_roi = squeeze(mean(cue_roi,1));
%     if num_ROIs == 1
%         avg_cue_roi = avg_cue_roi';
%     end
%     std_cue = squeeze(std(squeeze(cue_roi),1));
%     sm_cue = std_fail./sqrt(size(cue_roi,1));
%     for i = 1:num_ROIs  %baseline each curve so it passes through zero
%         shift = (-1)*avg_cue_roi(i,3);
%         avg_cue_roi(i,:) = avg_cue_roi(i,:)+shift;
%     end
%      for i = currROI;
%          subplot(1,3,3); hold on; errorbar(ts(i,:), avg_cue_roi(i,:), sm_cue(i,:), 'Color', colors{kk}); 
%      end
    
    %PLOT MAX SLOPE AND PEAK OF EACH PLOT. THEN COMPARE LATENCIES FOR CORR VS EARLY FOR EACH----- 
%     days{kk}
%     for i = ROIcell{kk};
%         for ii = 1:size(success_roi,1)
%             peak_frame = find(success_roi(ii,i,:)==max(success_roi(ii,i,6:11))); %finds the peak frame for each ROI from each trial. 
%             if length(peak_frame)>1
%                 peak_frame(peak_frame<6)=[];
%                 peak_frame(peak_frame>11)=[];
%                 peak_frame = peak_frame(1);
%             end
%             corr_roi_max = [corr_roi_max, peak_frame]; %stores the index of the frame num of peak df/f
%         end
%         avg_corr_roi_max = mean(corr_roi_max);
%         ind_peak = [ind_peak, avg_corr_roi_max]; % stores the avg frame num of the peak df/f for each ROI for each animal
%         for ii = 1:size(fail_roi,1)
%             peak_frame_fail = find(fail_roi(ii,i,:)==max(fail_roi(ii,i,6:11))); %finds the peak frame for each ROI from each trial.
%             if length(peak_frame_fail)>1
%                 peak_frame_fail(peak_frame_fail<6)=[];
%                 peak_frame_fail(peak_frame_fail>11)=[];
%                 peak_frame_fail = peak_frame_fail(1);
%             end
%             fail_roi_max = [fail_roi_max, peak_frame_fail]; %stores the index of the frame num of peak df/f
%         end
%         avg_fail_roi_max = mean(fail_roi_max);
%         ind_peak_fail = [ind_peak_fail, avg_fail_roi_max];
%         %peak_frame = find(avg_success_roi(i,:)==max(avg_success_roi(i,6:11))); %finds the peak frame for each ROI.
%         %ind_peak = [ind_peak, peak_frame];   %collects all the peak frames for correct trials for valid ROIs
%         %peak_frame_fail = find(avg_fail_roi(i,:)==max(avg_fail_roi(i,6:11)));
%         %ind_peak_fail = [ind_peak_fail, peak_frame_fail];
%         %vline(0,'k'); hold on;
%         %vline(peak_frame-(pre_frames+1), 'g'); hold on;   %subtract off (pre_frames+1) so the frame numbers will be aligned to lever release. +1 is for frame 0.
%         %vline(peak_frame_fail-(pre_frames+1), 'r'); hold on;
%         corr_fail_peak_latency = peak_frame - peak_frame_fail;
%     end

end
% vline(0,'k'); hold on;
% vline(ind_peak-(pre_frames+1), 'g'); hold on;   %subtract off (pre_frames+1) so the frame numbers will be aligned to lever release. +1 is for frame 0.
% vline(ind_peak_fail-(pre_frames+1), 'r'); hold on;
% xlim([-5 10]);
% disp(['latency to peak df/f correct - fail = ' num2str(ind_peak-ind_peak_fail)]);
% %title(['peak df/f for correct(green) and early(red) trials. avg diff=' num2str(corr_fail_peak_latency)]);
% xlabel('time from lever release');
% figure;
% %corr_peaks_plot = hist([ind_peak-(pre_frames+1)],ts(1,:)/100);
% corr_peaks_plot = hist([ind_peak-(pre_frames+1)],ts2);
% plot(ts2, corr_peaks_plot)
% hold on;
% %fail_peaks_plot = hist([ind_peak_fail-(pre_frames+1)],ts(1,:)/100);
% fail_peaks_plot = hist([ind_peak_fail-(pre_frames+1)],ts2);
% plot(ts2, fail_peaks_plot, 'r')
% xlabel('frames since lever release');
% ylabel('# of trials per 50ms bin');
% title('time to peak df/f for each ROI from each animal');
% sub_peaks_plot = (ind_peak-(pre_frames+1))-(ind_peak_fail-(pre_frames+1));
% figure; plot(ts2, sub_peaks_plot)
% title('time of avg peak df/f for corrects - same for earlies');
% xlabel('time of peak df/f for corrects (measure in 10hz frames)');



%figure(1);
%subplot(1,3,1);
xlabel('Time from release (ms)');
ylabel('dF/F');
%title('correct');
axis tight;
% 
% subplot(1,3,2);
% xlabel('Time from release (ms)');
% ylabel('dF/F');
% title('early');
% axis tight;
% 
% subplot(1,3,3);
% xlabel('Time from cue (ms)');
% ylabel('dF/F');
% title('cue aligned corrects');
% axis tight;
% 
% YL = [];
% for i =1:3
%     subplot(1,3,i); YL(i,:) = ylim;
% end
% subplot(1,3,1); ylim([min(YL(:,1)) max(YL(:,2))]);
% subplot(1,3,2); ylim([min(YL(:,1)) max(YL(:,2))]);
% subplot(1,3,3); ylim([min(YL(:,1)) max(YL(:,2))]);