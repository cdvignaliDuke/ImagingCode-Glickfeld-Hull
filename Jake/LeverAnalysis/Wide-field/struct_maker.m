%make structure to hold all variables needed for comparing bx criteria to F
%responses
%function [all_days, LS_ROI] = struct_maker(fieldname, threshold)

bx_dir = 'Z:\Data\WidefieldImaging\GCaMP\behavior\';
bx_outputs_dir = 'Z:\Analysis\WF Lever Analysis\BxAndAnalysisOutputs\BxOutputs\';
TC_dir = 'Z:\Analysis\WF Lever Analysis\LeverSummaryFolder\';
struct_dir = 'Z:\Analysis\WF Lever Analysis\StructuresForPlotting\'; 
xls_dir = 'Z:\Data\WidefieldImaging\GCaMP\WF_exp_spreadsheet';
%all_days = {'151009_img30', '151011_img30', '151013_img30', '151015_img29','151019_img30','151021_img29','151022_img29','151026_img30','151210_img32','151211_img32','151212_img32','151222_img32','160129_img36','160129_img35','160131_img35','160131_img36','160205_img35','160207_img35','160207_img36','160208_img35','160208_img36','160209_img36','160228_img36','160314_img38','160315_img38','160319_img41','160320_img41','160512_img47','160515_img44','160515_img47','160515_img48','160516_img47','160606_img46','160704_img54','160722_img53', '160725_img53'}; 
all_days = {'150716_img28', '150717_img28', '151021_img29', '151022_img29', '151009_img30', '151011_img30', '151211_img32', '151212_img32', '160129_img35', '160131_img35', '160129_img36','160131_img36', '160314_img38', '160315_img38', '160319_img41', '160320_img41', '160606_img46'}; %not really all days but it was the days I had LS ROIs selected for already
ROIcell = {[1], [1:5], [2], [2], [1:3], [1,3], [1:2], [1:2], [1:2], [1:2], [1:2], [1:2], [3:4], [2,3,5], [1], [1:2], [3:4]};  %LS areas
did not analyze session 151021_img29

%load excel spreadsheet data
[~, ~, xl_spread] = xlsread('Z:\Data\WidefieldImaging\GCaMP\WF_exp_spreadsheet'); 
row_heading = cell(1,size(xl_spread,2));
%extract row headings
for i = 1:size(xl_spread,2)
    row_heading{1, i} = xl_spread{1, i};
end

xl_struct = cell2struct(xl_spread, row_heading, 2);
for ii = 4:length(all_days);
    subj_num = ['i9' num2str(all_days{ii}(end-1:end))]; %useful for naming the structure
    session_date = all_days{ii}(1:6); 
    struct_name = [subj_num, '_', session_date]; %unfortunately have to wait until saving the struct to name it this
    curr_struct = struct('corr_TCs', [], 'corr_licking_TCs', [], 'corr_peak_latency_release', [], 'corr_rise_times', [], ...
        'corr_magnitude', [], 'corr_onset_latency', [], 'corr_react_times', [], 'corr_hold_dur', [], ...
        'corr_cue_aligned_TCs', [], 'corr_cue_aligned_licking_TCs', [], 'corr_cue_aligned_peak_latency', [], 'corr_cue_aligned_rise_times', [], ...
        'corr_cue_aligned_magnitude', [], 'corr_cue_aligned_onset_latency', [], ...
        'early_TCs', [], 'early_licking_TCs', [], 'early_peak_latency_release', [], 'early_rise_times', [], 'early_magnitude', [], 'early_onset_latency', [], ...
        'early_hold_dur', [], ...
        'cue_TCs', [], 'cue_licking_TCs', [], 'cue_peak_latency', [], 'cue_rise_times', [], 'cue_magnitude', [], 'cue_onset_latency', [], 'cue_react_times', [], ...
        'cue_hold_dur', [], ...
        'percent_correct', [], 'training_day_num', [], 'peak_percent_corr', [], 'recent_percent_corr', [], 'corr_react_time_mean', [], 'corr_react_time_std', [], ...
        'cue_react_time_mean', [], 'cue_react_time_std', [], ...
        'corr_magnitude_mean', [], 'early_magnitude_mean', [], 'cue_magnitude_mean', [], 'corr_magnitude_std', [], 'early_magnitude_std', [], 'cue_magnitude_std', [], 'corr_early_dfoverf_ratio', [], ...
        'corr_rise_time_mean', [], 'early_rise_time_mean', [], 'cue_rise_time_mean', [], 'corr_onset_latency_mean', [], 'early_onset_latency_mean', [], 'cue_onset_latency_mean', [], 'corr_cue_aligned_onset_latency_mean', [],...
        'corr_onset_latency_std', [], 'early_onset_latency_std', [], 'cue_onset_latency_std', [], 'corr_cue_aligned_onset_latency_std', [], ...
        'corr_peak_latency_mean', [], 'early_peak_latency_mean', [],'cue_peak_latency_mean', [], 'corr_peak_latency_std', [], 'early_peak_latency_std', [], 'cue_peak_latency_std', []);
    
    %load some files
    b_file = dir([bx_dir, 'data-', subj_num, '-' session_date, '*']); %find and load behavior file
    load([bx_dir, b_file.name]); 
    bx_data = input; clear input; %bx file is stored as input so I rename it and clear input
    load([bx_outputs_dir, all_days{ii}, '_bx_outputs']);   %frame_info, lever, trial_outcome, cluster, lickTimes, ...
    
    %locate react times which correspond to the successful lever releases
    %used for TCs in trial_outcome
    %make sure counterValues are appropriate as are the first and last
    %trial nums
    old_f_frame_num = frame_info.f_frame_trial_num;
    if bx_data.counterValues{frame_info.f_frame_trial_num}(1)==2 & bx_data.counterValues{frame_info.f_frame_trial_num-1}(1)==1 & size(bx_data.counterValues{frame_info.f_frame_trial_num-1},2)==1;
        frame_info.f_frame_trial_num = frame_info.f_frame_trial_num-1; %there was one instance where the 1st frame occured in one trial and the second frame occured in the next trial
    end
    assert(bx_data.counterValues{frame_info.f_frame_trial_num}(1)==1); %making sure counter is working as it should so we locate the correct frame time
    if size(bx_data.counterValues{frame_info.f_frame_trial_num},2) >1
        assert(bx_data.counterValues{frame_info.f_frame_trial_num}(2)==2) ;
    else
        assert(bx_data.counterValues{frame_info.f_frame_trial_num+1}(1)==2);
    end
    
    %find trial start times
    first_frame_time = bx_data.counterTimesUs{frame_info.f_frame_trial_num}(1);
    first_frame_time = round(first_frame_time/1000); %converting to Ms
    trial_start_times = round(cell2mat(bx_data.tThisTrialStartTimeMs));
    trial_start_times = trial_start_times - double(first_frame_time); %store the trial start times zeroed to the time of the first frame
    
    %find the trial nums for correct trials 
    if ~isempty(trial_outcome.success_time)
         success_times = trial_outcome.success_time; %time of lever release on successful trials 
         success_times = success_times(find(success_times>0));
         success_times_trial_num = get_trial_num_from_time(trial_start_times, success_times);
    else
        success_times = [];
    end
    %use the trial nums for correct trials to get/store their reaction time data
    corr_react_times = cell2mat(bx_data.reactTimesMs(success_times_trial_num));
    curr_struct.corr_react_times = corr_react_times;  % store the reaction times for correct trial TCs in the structure
    curr_struct.corr_react_time_mean = round(mean(corr_react_times));
    curr_struct.corr_react_time_std = round(std(double(corr_react_times)));
    
   %find the trial nums for lapse and tooFast trials
    if ~isempty(trial_outcome.late_time)
        late_times = trial_outcome.late_time; %time of lever release on lapsed trials
        late_times = late_times(find(late_times>0));
        late_times_trial_num = get_trial_num_from_time(trial_start_times, late_times);
    else
        late_times_trial_num = [];
    end 
    if ~isempty(trial_outcome.tooFastCorrects)
        tooFast_times = trial_outcome.tooFastCorrects; %time of lever release on tooFast trials
        tooFast_times = tooFast_times(find(tooFast_times>0));
        tooFast_times_trial_num = get_trial_num_from_time(trial_start_times, tooFast_times);
    else
        tooFast_times_trial_num = [];
    end
    
    %need to find the first and last trials used to generate the _roi TCs 
    all_cue_change_times_imaged = round(trial_outcome.change_orientation);
    post_frames = 10; pre_frames = 5; 
    last_time=  find(frame_info.counter<size(tc_dfoverf,2)-post_frames, 1, 'last');
    all_cue_change_times_imaged = all_cue_change_times_imaged(all_cue_change_times_imaged < last_time);
    first_time = find(frame_info.counter>pre_frames+1, 1, 'first');
    all_cue_change_times_imaged = all_cue_change_times_imaged(all_cue_change_times_imaged >first_time);
    %use cue_change times to get their trial nums too. 
    all_cue_times_trial_num = get_trial_num_from_time(trial_start_times, all_cue_change_times_imaged);
    
    %use trial number of all cue trials to get reaction times and hold duration
    cue_react_times = cell2mat(bx_data.reactTimesMs(all_cue_times_trial_num));
    curr_struct.cue_react_times = cue_react_times;  % store the reaction times for all cue trial TCs in the structure
    curr_struct.cue_react_time_mean = round(mean(cue_react_times));
    curr_struct.cue_react_time_std = round(std(double(cue_react_times)));
    
    %find hold duration for cue trials 
    cue_hold_times = cell2mat(bx_data.holdTimesMs(all_cue_times_trial_num));
    curr_struct.cue_hold_times = cue_hold_times;  % store the hold times for all cue trial TCs in the structure
    curr_struct.cue_hold_time_mean = round(mean(cue_hold_times));
    curr_struct.cue_hold_time_std = round(std(double(cue_hold_times)));
    
    %get time course data
    load([TC_dir, all_days{ii}, '_success']);  % stored as success_roi
    load([TC_dir, all_days{ii}, '_fail']);     % stored as fail_roi
    load([TC_dir, all_days{ii}, '_cue']);      % stored as fail_roi
    success_roi = success_roi(:,ROIcell{ii},:); %remove non-LS areas
    fail_roi = fail_roi(:,ROIcell{ii},:);       %remove non-LS areas
    cue_roi = cue_roi(:,ROIcell{ii},:);         %remove non-LS areas
    
    %interpolate TC data and store in structure
    XVector = [1:.01:size(success_roi,3)];
    success_roi_interp = nan(length(XVector), size(success_roi,1), size(success_roi,2)); %dims: 1=T2 2=trial# 3=ROI#
    for i = 1:size(success_roi,2); %for each roi...
        succ_temp = squeeze(success_roi(:,i,:))';  %dims 1=Time 2=trial number
        succ_temp = interp1(succ_temp, XVector);
        success_roi_interp(:,:,i) = succ_temp;  %dims 1=Time 2=trialNumber  3=roi
    end
    fail_roi_interp = nan(length(XVector), size(fail_roi,1), size(fail_roi,2));
    for i = 1:size(fail_roi,2);
        fail_temp = squeeze(fail_roi(:,i,:))';  %dims 1=Time 2=trial number
        fail_temp = interp1(fail_temp, XVector);
        fail_roi_interp(:,:,i) = fail_temp;   %dim 3=roi
    end
    cue_roi_interp = nan(length(XVector), size(cue_roi,1), size(cue_roi,2));
    for i = 1:size(cue_roi,2);
        cue_temp = squeeze(cue_roi(:,i,:))';  %dims 1=Time 2=trial number
        cue_temp = interp1(cue_temp, XVector);
        cue_roi_interp(:,:,i) = cue_temp;   %dim 3=roi
    end 
    if size(ROIcell{ii},2) > 1
        success_roi_interp = mean(success_roi_interp,3); %average together the ROIs
        fail_roi_interp = mean(fail_roi_interp,3);
        cue_roi_interp = mean(cue_roi_interp,3);
    end
    
    %baseline TCs 
    for iii = 1:size(success_roi_interp,2);  %dim1=time(10ms)  dim2=trial#
        shift = mean(success_roi_interp([150:350], iii));
        success_roi_interp(:, iii) = success_roi_interp(:,iii)-shift;
    end
    for iii = 1:size(fail_roi_interp,2);  %dim1=time(10ms)  dim2=trial#
        shift = mean(fail_roi_interp([150:350], iii));
        fail_roi_interp(:, iii) = fail_roi_interp(:,iii)-shift;
    end
    for iii = 1:size(cue_roi_interp,2);  %dim1=time(10ms)  dim2=trial#
        shift = mean(cue_roi_interp([400:600], iii));
        cue_roi_interp(:, iii) = cue_roi_interp(:,iii)-shift;
    end
    curr_struct.corr_TCs = success_roi_interp; %store TCs
    curr_struct.early_TCs = fail_roi_interp;
    curr_struct.cue_TCs = cue_roi_interp;
    
    %use the trials nums/times in common between the corr condition and the cue condition to index cue_roi and create the corr_cue condition
    corr_cue_aligned_trial_num = success_times_trial_num; %the trial nums for corr_cue aligned and the successful conditions are the same 
    [trash, all_cue_ind_to_keep, trash2] = intersect(all_cue_times_trial_num, success_times_trial_num);   
    corr_cue_aligned_TCs = cue_roi_interp(:,all_cue_ind_to_keep); %only keep the cue trials which correspond to correct trials 
    clear trash trash2; 
    assert(size(corr_cue_aligned_TCs,2) == size(success_roi_interp,2));
    curr_struct.corr_cue_aligned_TCs = corr_cue_aligned_TCs; 
    
    %find %corr and other variables using the excel spreadsheet
    [row_ind_day, col_ind_day] = find(strcmp(xl_spread, all_days{ii})); %find location of the row corresponding to this imaging session
    [row_ind, col_ind] = find(strcmp(xl_spread, '% correct')); %find the column for %corr
    this_percent_corr = cell2mat(xl_spread(row_ind_day, col_ind));
    curr_struct.percent_correct = this_percent_corr;
    
    [row_ind, col_ind] = find(strcmp(xl_spread, 'Training day')); %find the column for training day
    this_training_day = cell2mat(xl_spread(row_ind_day, col_ind));
    curr_struct.training_day_num = this_training_day;
    
    [row_ind, col_ind] = find(strcmp(xl_spread, 'peak % correct')); %find the peak percent correct of any training day previous to the current imaging day
    this_peak_perc_corr = cell2mat(xl_spread(row_ind_day, col_ind));
    curr_struct.peak_percent_corr = this_peak_perc_corr;
    
    [row_ind, col_ind] = find(strcmp(xl_spread, 'avg %corr pre-imaging')); %
    this_recent_perc_corr = cell2mat(xl_spread(row_ind_day, col_ind));
    curr_struct.recent_percent_corr = this_recent_perc_corr;
    
    %get hold duration data which matches TC data
    curr_struct.corr_hold_dur  = trial_outcome.succ_hold_dur;
    curr_struct.early_hold_dur = trial_outcome.fail_hold_dur;
    
    %find tbyt peak times and values for correct trials     SHOULD ADAPT SO IT USES AN AVG VALUE IN A 100MS WINDOW AROUND THE PEAK TIME??
    corr_tbyt_lat_peak = [];
    corr_tbyt_peak_mag = [];
    for i = 1:size(success_roi_interp,2);
        temp_val = find(success_roi_interp(:,i)==max(success_roi_interp([510:1100],i)));  %finds the peak time of each trials TC
        if size(temp_val,1) > 1 %in case of finding multiple values it will select the first one
            temp_val = temp_val(find(temp_val >= 509,1, 'first'));
        end
        corr_tbyt_lat_peak = [corr_tbyt_lat_peak, temp_val]; %give peak time for each trial. Will occur in increments of 10 because the peak will always be on one of the actual frames. Not interpolated data points
        corr_tbyt_peak_mag = [corr_tbyt_peak_mag, success_roi_interp(temp_val,i)]; %collects all the trial by trial peak magnitudes
    end
    corr_tbyt_lat_peak = corr_tbyt_lat_peak-600;
    curr_struct.corr_peak_latency_release = corr_tbyt_lat_peak; %have to convert it to ms and zero to lever release
    curr_struct.corr_peak_latency_mean = mean(corr_tbyt_lat_peak);
    curr_struct.corr_peak_latency_std = std(corr_tbyt_lat_peak);
    curr_struct.corr_magnitude = corr_tbyt_peak_mag;
    curr_struct.corr_magnitude_mean = mean(corr_tbyt_peak_mag);
    curr_struct.corr_magnitude_std = std(corr_tbyt_peak_mag);
    
    %find tbyt peak times and values for early trials 
    early_tbyt_lat_peak = []; %SHOULD ADAPT SO IT USES AN AVG VALUE IN A 100MS WINDOW AROUND THE PEAK TIME??
    early_tbyt_peak_mag = [];
    for i = 1:size(fail_roi_interp,2);
        temp_val = find(fail_roi_interp(:,i)==max(fail_roi_interp([510:1100],i)));
        if size(temp_val,1) > 1
            temp_val = temp_val(find(temp_val >= 509,1, 'first'));
        end
        early_tbyt_lat_peak = [early_tbyt_lat_peak, temp_val]; %give peak time for each trial
        early_tbyt_peak_mag = [early_tbyt_peak_mag, fail_roi_interp(temp_val,i)]; %collects all the trial by trial peaks 
    end
    early_tbyt_lat_peak = early_tbyt_lat_peak-600;
    curr_struct.early_peak_latency_release = (early_tbyt_lat_peak); %have to convert it to ms and zero to lever release
    curr_struct.early_peak_latency_mean = mean(early_tbyt_lat_peak);
    curr_struct.early_peak_latency_std = std(early_tbyt_lat_peak);
    curr_struct.early_magnitude = early_tbyt_peak_mag;
    curr_struct.early_magnitude_mean = mean(early_tbyt_peak_mag);
    curr_struct.early_magnitude_std = std(early_tbyt_peak_mag);
    
    %calculate ratio of avg corr vs early df/f response magnitude
    curr_struct.corr_early_dfoverf_ratio =mean(corr_tbyt_peak_mag)/mean(early_tbyt_peak_mag);
    
    %find tbyt peak times and values for cue trials 
    cue_tbyt_lat_peak = []; %SHOULD ADAPT SO IT USES AN AVG VALUE IN A 100MS WINDOW AROUND THE PEAK TIME??
    cue_tbyt_peak_mag = [];
    for i = 1:size(cue_roi_interp,2);
        temp_val = find(cue_roi_interp(:,i)==max(cue_roi_interp([710:1300],i)));
        if size(temp_val,1) > 1
            temp_val = temp_val(find(temp_val >= 709,1, 'first'));
        end
        cue_tbyt_lat_peak = [cue_tbyt_lat_peak, temp_val]; %give peak time for each trial
        cue_tbyt_peak_mag = [cue_tbyt_peak_mag, cue_roi_interp(temp_val,i)]; %collects all the trial by trial peaks 
    end
    cue_tbyt_lat_peak = cue_tbyt_lat_peak-600;
    curr_struct.cue_peak_latency = (cue_tbyt_lat_peak); %have to convert it to ms and zero to cue
    curr_struct.cue_peak_latency_mean = mean(cue_tbyt_lat_peak);
    curr_struct.cue_peak_latency_std = std(cue_tbyt_lat_peak);
    curr_struct.cue_magnitude = cue_tbyt_peak_mag;
    curr_struct.cue_magnitude_mean = mean(cue_tbyt_peak_mag);
    curr_struct.cue_magnitude_std = std(cue_tbyt_peak_mag);
    
    %find tbyt peak times and values for corr_cue_aligned trials 
    corr_cue_tbyt_lat_peak = []; %SHOULD ADAPT SO IT USES AN AVG VALUE IN A 100MS WINDOW AROUND THE PEAK TIME??
    corr_cue_tbyt_peak_mag = [];
    for i = 1:size(corr_cue_aligned_TCs,2);
        temp_val = find(corr_cue_aligned_TCs(:,i)==max(corr_cue_aligned_TCs([710:1300],i)));
        if size(temp_val,1) > 1
            temp_val = temp_val(find(temp_val >= 709,1, 'first'));
        end
        corr_cue_tbyt_lat_peak = [corr_cue_tbyt_lat_peak, temp_val]; %give peak time for each trial
        corr_cue_tbyt_peak_mag = [corr_cue_tbyt_peak_mag, corr_cue_aligned_TCs(temp_val,i)]; %collects all the trial by trial peaks 
    end
    corr_cue_tbyt_lat_peak = corr_cue_tbyt_lat_peak-600;
    curr_struct.corr_cue_aligned_peak_latency = (corr_cue_tbyt_lat_peak); %have to convert it to ms and zero to cue
    curr_struct.corr_cue_aligned_peak_latency_mean = mean(corr_cue_tbyt_lat_peak);
    curr_struct.corr_cue_aligned_peak_latency_std = std(corr_cue_tbyt_lat_peak);
    curr_struct.corr_cue_aligned_magnitude = corr_cue_tbyt_peak_mag;
    curr_struct.corr_cue_aligned_magnitude_mean = mean(corr_cue_tbyt_peak_mag);
    curr_struct.corr_cue_aligned_magnitude_std = std(corr_cue_tbyt_peak_mag);
    
    %find the rise times and onset times
    rise_time_data = rise_time_finder(curr_struct);
    
    %store various variables
    curr_struct.corr_rise_times = rise_time_data.corr_rise_times;
    curr_struct.corr_onset_latency = rise_time_data.corr_onset_times-600;
    curr_struct.early_rise_times = rise_time_data.early_rise_times;
    curr_struct.early_onset_latency = rise_time_data.early_onset_times-600;
    curr_struct.cue_rise_times = rise_time_data.cue_rise_times;
    curr_struct.cue_onset_latency = rise_time_data.cue_onset_times-600;
    curr_struct.corr_cue_aligned_rise_times = rise_time_data.corr_cue_aligned_rise_times;
    curr_struct.corr_cue_aligned_onset_latency = rise_time_data.corr_cue_aligned_onset_times-600;
    
    curr_struct.corr_rise_time_mean = mean(curr_struct.corr_rise_times);
    curr_struct.early_rise_time_mean = mean(curr_struct.early_rise_times);
    curr_struct.cue_rise_time_mean = mean(curr_struct.cue_rise_times);
    
    curr_struct.corr_onset_latency_mean = mean(curr_struct.corr_onset_latency);
    curr_struct.early_onset_latency_mean = mean(curr_struct.early_onset_latency);
    curr_struct.cue_onset_latency_mean = mean(curr_struct.cue_onset_latency);
    curr_struct.corr_cue_aligned_onset_latency_mean = mean(curr_struct.corr_cue_aligned_onset_latency);
    
    curr_struct.corr_onset_latency_std = std(curr_struct.corr_onset_latency);
    curr_struct.early_onset_latency_std = std(curr_struct.early_onset_latency);
    curr_struct.cue_onset_latency_std = std(curr_struct.cue_onset_latency);
    curr_struct.corr_cue_aligned_onset_latency_std = std(curr_struct.corr_cue_aligned_onset_latency);
    
    save([struct_dir, all_days{ii}], 'curr_struct');
    disp(all_days(ii));
end
