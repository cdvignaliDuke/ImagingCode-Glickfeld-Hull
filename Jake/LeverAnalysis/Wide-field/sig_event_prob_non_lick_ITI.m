%adopted from lick_triggered_averages
%script to isolate ITI periods without licks and then plot the probability
%of a significant event in that period. 
clear;

pre_frames = 8;
post_frames = 4;
bout_start_frame = pre_frames+1;
peak_window = [bout_start_frame-1:bout_start_frame+2];
min_lick_window = 2; %do not include frames from -min_lick_window:min_lick_window around a lick in the baseline calculation
non_lick_window = [1:13];

bx_source      = ['Z:\Data\WidefieldImaging\GCaMP\behavior\'];
bx_outputs_dir = ['Z:\Analysis\WF Lever Analysis\BxAndAnalysisOutputs\BxOutputs\'];  
baseline_output = ['Z:\Analysis\WF Lever Analysis\licking_investigation\baselines_ITI_no_licks'];
non_lick_output = ['Z:\Analysis\WF Lever Analysis\licking_investigation\lick_triggered_averages\non_lick_interval_TCs\'];
old_cd = cd; %save old cd so I can restore it later
WF_plotting_lists_of_days;

%datasets with good licking traces
days = {'151021_img29', '151022_img29', [], [], '151211_img32', '151212_img32', '160129_img35', '160131_img35', '160129_img36', '160131_img36', '160314_img38', '160315_img38', '160319_img41', '160320_img41', '160606_img46', '160722_img53', '160904_img55'}; %'150718_img27', '150719_img27', '150716_img28', '150717_img28', 
colors_roi = {'b', 'r', 'g', 'k', 'c', 'm'};
colors = {'k','k', 'b','b', 'g','g', 'c','c', 'm','m', 'r','r', 'y','y', 'k', 'b', 'r'}; 
plot_symbols = {'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 's', 's', 's'};

non_lick_mags = {}; 
percent_significant{1} = [];
percent_significant_ROI{1} = [];
num_intervals{1} = [];
num_trials_cell{1} = [];
for session_num = 1:length(days);
    if isempty(days{session_num})
        non_lick_mags{session_num} = [];
        continue
    end
    b_data = get_bx_data(bx_source, days{session_num});  %find the correct behavior file and loads it.
    bx_out_dir  = [bx_outputs_dir days{session_num} '_bx_outputs'];
    load(bx_out_dir);
    
    %assign variables
    imaging_start_MW_T = frame_info.imaging_start_MW_T;
    licksByFrame = licking_data.licksByFrame; 
    leverPress = cell2mat(b_data.holdStartsMs)-imaging_start_MW_T;
    holdDur = cell2mat(b_data.holdTimesMs);
    leverRelease = leverPress + holdDur;
    reqHoldTime = cell2mat(b_data.tTotalReqHoldTimeMs);
    changeOrientation = trial_outcome.change_orientation;
    
    %find and trim the trials which were not imaged. 
    trials_to_trim = [(1:frame_info.f_frame_trial_num), (frame_info.l_frame_trial_num:length(leverPress))];
    if ~isempty(trials_to_trim);
        leverPress(trials_to_trim) = [];
        holdDur(trials_to_trim) = [];
        leverRelease(trials_to_trim) = [];
        reqHoldTime(trials_to_trim) = [];
        changeOrientation(trials_to_trim) = [];
    end
    
    %collect all the lever events and index correct trials
    leverPressFrame = [];
    for iii = 1:length(leverPress)
        leverPressFrame = [leverPressFrame, frame_info.counter(leverPress(iii))];
    end
    leverReleaseFrame = [];
    for iii = 1:length(leverRelease)
        leverReleaseFrame = [leverReleaseFrame, frame_info.counter(leverRelease(iii))];
    end
    corrInx = zeros(1, length(holdDur));
    for iii = 1:length(holdDur)
        if holdDur(iii) > reqHoldTime(iii)
            corrInx(iii)=1;
        end
    end
    earlyInx = find(corrInx == 0);
    corrInx  = find(corrInx == 1);
    
    %find periods in the ITI without licking to get a baseline
    %identify ITI frames not near lever events
    ITI_bool = ones(1,length(tc_dfoverf));
    for lever_event = 1:length(leverReleaseFrame)
        ITI_bool(leverPressFrame(lever_event)-3:leverReleaseFrame(lever_event)+5) = 0;
    end
    ITI_inx_frames = find(ITI_bool);
    
    %sometimes there are missing frames at the end of licksbyframe. need to lengthen it and correct ITI_inx
    if length(ITI_bool) > length(licksByFrame)
        ITI_bool(length(licksByFrame)+1:end) = 0;
        licksByFrame(length(licksByFrame)+1:length(ITI_bool)) = 0;
    end
    
    %use ITI boolean and licksByFrame to get a boolean of ITI frames without licking
    non_lick_frames = ~licksByFrame;
    non_lick_ITI_bool = non_lick_frames.*ITI_bool;
    non_lick_ITI_inx = find(non_lick_ITI_bool);
    
    %look for periods of ITI without licking that match length of non_lick_window
    N_frames_min = length(non_lick_window); %required number of consecutive non-lick ITI frames 
    diff_non_lick_ITI_bool = diff(non_lick_ITI_inx)==1; %boolean indicating all datapoints in diff(non_lick_ITI_inx) followed by a squentially increasing value  (x)
    initial_value_in_sequence =find([false,diff_non_lick_ITI_bool]~=[diff_non_lick_ITI_bool,false]); %find the first value which does not match the value before it (f)
    valid_sequence_inx = find(initial_value_in_sequence(2:2:end)-initial_value_in_sequence(1:2:end-1)>=N_frames_min); %find which valuess in initial_value_in_sequence are the first in a sequence > N_frames_min (g)
    first_in_valid_sequence = non_lick_ITI_inx(initial_value_in_sequence(2*valid_sequence_inx-1)); 
    first_in_valid_sequence = first_in_valid_sequence(first_in_valid_sequence>20 & first_in_valid_sequence<size(tc_dfoverf,2)-20);
    if isempty(first_in_valid_sequence)
        continue
    end
    
    %plot the lick traces from first_in_valid_sequence-5:first_in_valid_sequence+15 to make sure algorithm worked properly
    all_lever_events = sort([leverPressFrame, leverReleaseFrame]);
    figure;
    for iii = 1:length(first_in_valid_sequence)
        plot([-5:15],  [ licksByFrame(first_in_valid_sequence(iii)-5:first_in_valid_sequence(iii)+15) ]); hold on;
        assert(isempty(intersect(all_lever_events, licksByFrame(first_in_valid_sequence(iii)-2:first_in_valid_sequence(iii)+N_frames_min+2)))); %checks to make sure there are no lever events near the identified intervals. 
    end
    title([days{session_num}, ' no-lick window=', num2str(N_frames_min)]);
    
    %remove non-LS ROIs
    tc_dfoverf = tc_dfoverf([LS_ROIs{session_num}],:);
    
    %use first_in_valid_sequence to collect frame intervals and measure df/f
    for iii = 1:length(first_in_valid_sequence) 
        if iii ==1 
            non_lick_TCs = tc_dfoverf(:, first_in_valid_sequence(iii):first_in_valid_sequence(iii)+non_lick_window(end)-1);
        else
            non_lick_TCs(:,:,iii) = tc_dfoverf(:, first_in_valid_sequence(iii):first_in_valid_sequence(iii)+non_lick_window(end)-1); %dim1=ROI dim2=frame# dim3=interval#
        end
    end
    non_lick_TCs = reshape(non_lick_TCs, size(non_lick_TCs,3),size(non_lick_TCs,1),size(non_lick_TCs,2)); %dim1=interval# dim2=ROI dim3=frame#
   
    %shift each TC so 0 is the mean of...
    shift = mean(non_lick_TCs(:,:,[1:6]),3);
    shift = repmat(shift, 1, 1, non_lick_window(end));
    non_lick_TCs = non_lick_TCs-shift;

    %get mean and sem across extracted intervals for each ROI
    non_lick_TC_sem = squeeze(std(non_lick_TCs,1)/sqrt(size(non_lick_TCs,1)));
    non_lick_TC_avg = squeeze(mean(non_lick_TCs)); %averaging across trials
    if length(LS_ROIs{session_num}) == 1
        non_lick_TC_avg = non_lick_TC_avg';
        non_lick_TC_sem = non_lick_TC_sem';
    end
    
    %find and store the peak values (averaged across trials)
    if size(non_lick_TCs,1) >= 4
        non_lick_mags{session_num} = max(non_lick_TC_avg(:,peak_window),[],2)';
    else
        non_lick_mags{session_num} = [];
    end
    
    %isolate the licks which occur in the ITI and use them to exclude frames from which to measure the baseline df/f
    ITI_licks = find(licksByFrame.*ITI_bool);
    for lick_num = 1:length(ITI_licks)
        if ITI_licks(lick_num)<3
            ITI_bool([1:4]) =0;
            continue
        end
        ITI_bool(ITI_licks(lick_num)-min_lick_window:ITI_licks(lick_num)+min_lick_window) = 0;
    end
    baseline_frames = tc_dfoverf(:, [find(ITI_bool)]); %non-LS ROIs already excluded before this
    baseline_frames_std = std(baseline_frames,[],2)';
    %save([baseline_output, '\', days{session_num}], 'baseline_frames', 'baseline_frames_std', 'min_lick_window');
        
    %determine what percent of trials exceed 2 stdevs first by ROI
    tbyt_peak_mag = max(non_lick_TCs(:,:,peak_window),[],3);
    this_session_sig = [];
    for ROI_num = 1:size(tbyt_peak_mag,2);
        this_session_sig = [this_session_sig, length(find(tbyt_peak_mag(:,ROI_num)>2*baseline_frames_std(ROI_num)))/size(tbyt_peak_mag,1)];
    end
    percent_significant_ROI{session_num} = this_session_sig;
        
    %then averaged across ROIs
    int_start_trig_roi_avg = squeeze(mean(non_lick_TCs,2));
    tbyt_peak_mag_avg = max(int_start_trig_roi_avg(:,peak_window),[],2);
    percent_significant{session_num} = length(find(tbyt_peak_mag_avg>2*mean(baseline_frames_std)))/length(tbyt_peak_mag_avg);
    num_intervals{session_num} = length(tbyt_peak_mag_avg);
    
    %identify which trials have a significant response
    sig_non_lick_resp_mat = zeros(size(tbyt_peak_mag));
    for ROI_num = 1:size(tbyt_peak_mag,2)
        sig_non_lick_resp_mat([find(tbyt_peak_mag(:,ROI_num)>2*baseline_frames_std(ROI_num))], ROI_num) = 1; %find all trials with significant response and log them in this matrix on an ROI basis
    end
    
    %make a color coded plot: trial by trial TCs of sig vs non-sig trials
    figure;
    for ROI_num = 1:size(tbyt_peak_mag,2)
        subplot(1,size(tbyt_peak_mag,2),ROI_num); hold on;
        this_roi_sig_frac = sum(sig_non_lick_resp_mat(:,ROI_num))/size(sig_non_lick_resp_mat,1);
        title(['ROI #', num2str(ROI_num), ' % of sig. trials=',num2str(this_roi_sig_frac)]);
        ylabel('df/f'); xlabel('frames relative to arbitrarly selected frame (no lick)');
        for trial_num = 1:size(tbyt_peak_mag,1)
            if sig_non_lick_resp_mat(trial_num, ROI_num) == 1
                plot([-pre_frames:post_frames], squeeze(non_lick_TCs(trial_num, ROI_num, :)), 'b');
            elseif sig_non_lick_resp_mat(trial_num, ROI_num) == 0
                plot([-pre_frames:post_frames], squeeze(non_lick_TCs(trial_num, ROI_num, :)), 'k');
            end
        end
        errorbar([-pre_frames:post_frames], non_lick_TC_avg(ROI_num,:), non_lick_TC_sem(ROI_num,:), 'm', 'LineWidth', 1.5);
    end
    suptitle([days(session_num), 'trial by trial lick triggered TCs: blue=significant black=non-sig (>2 St.Dev)']);
    %savefig([non_lick_output, days{session_num}]);
    
    %store the number of ITI licks in this session
    num_trials_cell{session_num} = size(sig_non_lick_resp_mat,1);
    
    %use first_in_valid_sequence to make a lick trace for ITI licks. Will feed
    %into another script to plot along with lever lick traces
    ITI_lick_trace = NaN(length(first_in_valid_sequence),16);
    for bout_num = 1:length(first_in_valid_sequence)
        ITI_lick_trace(bout_num, :) = licksByFrame(first_in_valid_sequence(bout_num)-5:first_in_valid_sequence(bout_num)+10);
    end
    assert(sum(sum(isnan(ITI_lick_trace)))==0);
    %save(['Z:\Analysis\WF Lever Analysis\licking_investigation\lick_triggered_averages\individual lick bout TCs\', days{session_num}, '_lick_trace'], 'ITI_lick_trace');
end

figure;
for session_num = 1:length(days)
    if isempty(percent_significant_ROI{session_num}) | num_trials_cell{session_num} < 10
        continue
    %elseif isempty(valid_LS_ROIs{session_num}) %only plot ROIs with a minimum corr:early ratio of 1.15
        %continue
    end
    plot(1, percent_significant{session_num}, strcat(plot_symbols{session_num}, colors{session_num})); hold on;
end
ylim([0 1]); ylabel('percent of extracted non-lick intervals with a significant positive df/f response');
title('trial peak df/f > 2 standard deviations above baseline. Averaged across ROIs');
%savefig([non_lick_output, 'percent_sig_trials_2_StDev']);
%save(non_lick_output, 'non_lick_mags');

%plot each ROI
norm_non_lick_prob_mat=[];
figure;
for session_num = 1:length(days)
    if isempty(percent_significant_ROI{session_num})  | num_trials_cell{session_num} < 10
        continue
%     elseif isempty(valid_LS_ROIs{session_num}) %only plot ROIs with a minimum corr:early ratio of 1.15
%         continue
    end
    norm_non_lick_prob_mat = [norm_non_lick_prob_mat , percent_significant_ROI{session_num}];
    plot(ones(1, length(percent_significant_ROI{session_num})), percent_significant_ROI{session_num}, strcat(plot_symbols{session_num}, colors{session_num})); hold on;
    %norm_LTA_prob_mat = [norm_LTA_prob_mat , percent_significant_ROI{session_num}([valid_LS_ROIs{session_num}])];
    %plot(ones(1, length(percent_significant_ROI{session_num}([valid_LS_ROIs{session_num}]))),   percent_significant_ROI{session_num}([valid_LS_ROIs{session_num}]),    strcat(plot_symbols{session_num}, colors{session_num})); hold on;
end
errorbar(1, mean(norm_non_lick_prob_mat), std(norm_non_lick_prob_mat)/sqrt(length(norm_non_lick_prob_mat)), 'o', 'MarkerFaceColor', 'k');
%plot(ones(1, length(plotting_sig_sessions)), plotting_sig_sessions, 'o');
xlim([0.5 1.5]);
ylim([0 1]); ylabel('percent of non-lick intervals with a significant positive df/f response');
title('trial peak df/f > 2 standard deviations above baseline.');
%savefig([non_lick_output, 'percent_sig_trials_2_StDev_ROI']);



