%% script for taking correlation between df/f peak time and lever vs cue events
WF_plotting_lists_of_days;
bx_outputs_dir = ['Z:\Analysis\WF Lever Analysis\BxAndAnalysisOutputs\BxOutputs\'];
TC_dir = ['Z:\Analysis\WF Lever Analysis\LeverSummaryFolder\'];
step_size = 0.1;

colors = [1,0,0; 0,1,0; 0,0,1; 0.5,0.5,0.5; 1,0,1; 1,1,0; 0,1,1]; %sets up the color scheme for plotting multiple ROIs with errorbar 
plot_colors = {'r', 'b', 'g', 'k', 'c', 'm'};
colors_final_plot = {'k','k', 'b','b', 'g','g', 'c','c', 'm','m', 'r','r', 'y','y', 'k', 'b', 'r'}; 
plot_symbols = {'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 's', 's', 's'};
x_axis = [-4:3];
release_frame = 6;
release_frame_interp = 501;
analysis_window = [release_frame-1:release_frame+3];
analysis_window_cue = [release_frame+1:release_frame+5];
baseline_window = [1:3];
baseline_window_cue = [1:3];

corr_avg_lat_peak{1} = [];
fail_avg_lat_peak{1} = [];
cue_avg_lat_peak{1}  = [];
corr_std_lat_peak{1} = [];
fail_std_lat_peak{1} = [];
cue_std_lat_peak{1}  = [];

for session_num = 1:length(days);
    %load bx and corr trials TCs
    bx_out_dir  = [bx_outputs_dir days{session_num} '_bx_outputs'];
    load(bx_out_dir);
    load([TC_dir, days{session_num}, '_success']);
    load([TC_dir, days{session_num}, '_fail']);
    load([TC_dir, days{session_num}, '_tooFast']); %necessary?
    load([TC_dir, days{session_num}, '_cue']);
    lick_trace_corr = licking_data.lick_trace_succ;
    
    %remove non-LS ROIs
    fail_roi = fail_roi(:,[LS_ROIs{session_num}],:);
    success_roi = success_roi(:,[LS_ROIs{session_num}],:);
    cue_roi = cue_roi(:,[LS_ROIs{session_num}],:);
    
    %baseline TCs for non-interpolated data
    for ROI_num = 1:size(success_roi,2);
        for trial_num=1:size(success_roi,1)
            shift = squeeze(mean(success_roi(trial_num, ROI_num, baseline_window),3));
            success_roi(trial_num, ROI_num, :) = success_roi(trial_num, ROI_num, :)-shift;
            shift = squeeze(mean(cue_roi(trial_num, ROI_num, baseline_window_cue),3));
            cue_roi(trial_num, ROI_num, :) = cue_roi(trial_num, ROI_num, :)-shift;
        end
        for trial_num=1:size(fail_roi,1)
            shift = squeeze(mean(fail_roi(trial_num, ROI_num, baseline_window),3));
            fail_roi(trial_num, ROI_num, :) = fail_roi(trial_num, ROI_num, :)-shift;
        end
    end
    
    %find peak time and mags for either by ROI or for a single ROI or averaged across ROIs
    peak_mag_window = 0; %set to zero to get a single value. Otherwise it well get the mean of peak_ind-peak_mag_window : peak_ind+peak_mag_window
    [corr_tbyt_lat_peak, corr_tbyt_peak_mag] = get_peak_lat_and_mag(success_roi, analysis_window, release_frame, peak_mag_window);
    [fail_tbyt_lat_peak, fail_tbyt_peak_mag] = get_peak_lat_and_mag(fail_roi, analysis_window, release_frame, peak_mag_window);
    [cue_tbyt_lat_peak, cue_tbyt_peak_mag] = get_peak_lat_and_mag(cue_roi, analysis_window_cue, release_frame, peak_mag_window); %for cue_roi "release_frame" is actually the time of the cue
    
    %store lats and mags
    
    %Determine time of peak df/f for each ROI averaged across trials.
    corr_avg_lat_peak{session_num} = mean(corr_tbyt_lat_peak,2);
    fail_avg_lat_peak{session_num} = mean(fail_tbyt_lat_peak,2);
    cue_avg_lat_peak{session_num} = mean(cue_tbyt_lat_peak,2);
    
    %Determine std of peak df/f for corr vs cueCorr vs fail trials
    corr_std_lat_peak{session_num} = std(corr_tbyt_lat_peak, [],2);
    fail_std_lat_peak{session_num} = std(fail_tbyt_lat_peak, [],2);
    cue_std_lat_peak{session_num} = std(cue_tbyt_lat_peak, [],2);
    
    
end

%Plot the std values for each ROI. corr vs cue -------------------
figure;
subplot(1,2,1);
for session_num = 1:length(days)
    for ROI_num = 1:length(corr_std_lat_peak{session_num})
        if ~isempty(valid_LS_ROIs{session_num}) & ismember(ROI_num, valid_LS_ROIs{session_num});
           plot([1,2], [corr_std_lat_peak{session_num}(ROI_num), cue_std_lat_peak{session_num}(ROI_num)], strcat('-', plot_symbols{session_num}, 'g'));   hold on;
        else
           plot([1,2], [corr_std_lat_peak{session_num}(ROI_num), cue_std_lat_peak{session_num}(ROI_num)], strcat('-', plot_symbols{session_num}, 'r'));   hold on;
        end
        %plot([1,2], [corr_std_lat_peak{session_num}(ROI_num), cue_std_lat_peak{session_num}(ROI_num)], strcat('-', plot_symbols{session_num}, colors_final_plot{session_num}));   hold on;
    end
end
ylabel('st. Dev (100ms frames)');
xlabel('lever aligned        cue aligned');
xlim([0 3]); ylim([0 2]);
%Plot the std values for each ROI. corr vs early
subplot(1,2,2);
for session_num = 1:length(days)
    for ROI_num = 1:length(corr_std_lat_peak{session_num})
        if ~isempty(valid_LS_ROIs{session_num}) & ismember(ROI_num, valid_LS_ROIs{session_num});
          plot([1,2], [corr_std_lat_peak{session_num}(ROI_num), fail_std_lat_peak{session_num}(ROI_num)], strcat('-', plot_symbols{session_num}, 'g'));   hold on;
        else
           plot([1,2], [corr_std_lat_peak{session_num}(ROI_num), fail_std_lat_peak{session_num}(ROI_num)], strcat('-', plot_symbols{session_num}, 'r'));   hold on;
        end
        %plot([1,2], [corr_std_lat_peak{session_num}(ROI_num), fail_std_lat_peak{session_num}(ROI_num)], strcat('-', plot_symbols{session_num}, colors_final_plot{session_num}));   hold on;
    end
end
xlabel('correct        early');
xlim([0 3]); ylim([0 2]);
suptitle(['standard deviation of the latency to peak df/f']);

%plot the mean latency to peak df/f for each ROI: corr vs cue --------------------------------------------------
figure;
subplot(1,2,1);
for session_num = 1:length(days)
    for ROI_num = 1:length(corr_avg_lat_peak{session_num})
        if ~isempty(valid_LS_ROIs{session_num}) & ismember(ROI_num, valid_LS_ROIs{session_num});
            plot([1,2], [corr_avg_lat_peak{session_num}(ROI_num), cue_avg_lat_peak{session_num}(ROI_num)], strcat('-', plot_symbols{session_num}, 'g'));   hold on;
        else
            plot([1,2], [corr_avg_lat_peak{session_num}(ROI_num), cue_avg_lat_peak{session_num}(ROI_num)], strcat('-', plot_symbols{session_num}, 'r'));   hold on;
        end
        %plot([1,2], [corr_avg_lat_peak{session_num}(ROI_num), cue_avg_lat_peak{session_num}(ROI_num)], strcat('-', plot_symbols{session_num}, colors_final_plot{session_num}));   hold on;
    end
end
ylabel('latency to peak df/f (100ms frames)');
xlabel('lever aligned        cue aligned');
xlim([0 3]); ylim([-1 5]);
%correct vs early
subplot(1,2,2);
for session_num = 1:length(days)
    for ROI_num = 1:length(corr_avg_lat_peak{session_num})
        if ~isempty(valid_LS_ROIs{session_num}) & ismember(ROI_num, valid_LS_ROIs{session_num});
            plot([1,2], [corr_avg_lat_peak{session_num}(ROI_num), fail_avg_lat_peak{session_num}(ROI_num)], strcat('-', plot_symbols{session_num}, 'g'));   hold on;
        else
            plot([1,2], [corr_avg_lat_peak{session_num}(ROI_num), fail_avg_lat_peak{session_num}(ROI_num)], strcat('-', plot_symbols{session_num}, 'r'));   hold on;
        end
        %plot([1,2], [corr_avg_lat_peak{session_num}(ROI_num), fail_avg_lat_peak{session_num}(ROI_num)], strcat('-', plot_symbols{session_num}, colors_final_plot{session_num}));   hold on;
    end
end
ylabel('latency to peak df/f (100ms frames)');
xlabel('correct        early');
xlim([0 3]); ylim([-1 5]);
suptitle('latency to peak df/f');


 %uncomment if I decide to use interpolated data
%     %interpolate TC data and store in structure   %the new x-value for the lever release is 501
%     success_roi_interp = interpolate_TC_data(success_roi, step_size);
%     fail_roi_interp = interpolate_TC_data(fail_roi, step_size);
%     cue_roi_interp = interpolate_TC_data(cue_roi, step_size);
%     
%     %baseline TCs 
%     for trial_num = 1:size(success_roi_interp,2);  %dim1=time(10ms)  dim2=trial#
%         shift = mean(success_roi_interp([1:150], trial_num));
%         success_roi_interp(:, trial_num) = success_roi_interp(:,trial_num)-shift;
%     end
%     for trial_num = 1:size(fail_roi_interp,2);  %dim1=time(10ms)  dim2=trial#
%         shift = mean(fail_roi_interp([1:150], trial_num));
%         fail_roi_interp(:, trial_num) = fail_roi_interp(:,trial_num)-shift;
%     end
%     for trial_num = 1:size(cue_roi_interp,2);  %dim1=time(10ms)  dim2=trial#
%         shift = mean(cue_roi_interp([250:400], trial_num));
%         cue_roi_interp(:, trial_num) = cue_roi_interp(:,trial_num)-shift;
%     end