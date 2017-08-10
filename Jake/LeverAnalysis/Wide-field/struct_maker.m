%make structure to hold all variables needed for comparing bx criteria to F
%responses
%function [all_days, LS_ROI] = struct_maker(fieldname, threshold)
clear
bx_dir = 'Z:\Data\WidefieldImaging\GCaMP\behavior\';
bx_outputs_dir = 'Z:\Analysis\WF Lever Analysis\BxAndAnalysisOutputs\BxOutputs\';
TC_dir = 'Z:\Analysis\WF Lever Analysis\LeverSummaryFolder\';
struct_dir = 'Z:\Analysis\WF Lever Analysis\StructuresForPlotting\'; 
xls_dir = 'Z:\Data\WidefieldImaging\GCaMP\WF_exp_spreadsheet';
%all_days = {'151009_img30', '151011_img30', '151013_img30', '151015_img29','151019_img30','151021_img29','151022_img29','151026_img30','151210_img32','151211_img32','151212_img32','151222_img32','160129_img36','160129_img35','160131_img35','160131_img36','160205_img35','160207_img35','160207_img36','160208_img35','160208_img36','160209_img36','160228_img36','160314_img38','160315_img38','160319_img41','160320_img41','160512_img47','160515_img44','160515_img47','160515_img48','160516_img47','160606_img46','160704_img54','160722_img53', '160725_img53'}; 
%old popultaion of sessions analyzed all_days = {'150716_img28', '150717_img28', '151021_img29', '151022_img29', '151009_img30', '151011_img30', '151211_img32', '151212_img32', '160129_img35', '160131_img35', '160129_img36','160131_img36', '160314_img38', '160315_img38', '160319_img41', '160320_img41', '160606_img46'}; %not really all days but it was the days I had LS ROIs selected for already
%ROIcell = {[1], [1:5], [2], [2], [1:3], [1,3], [1:2], [1:2], [1:2], [1:2], [1:2], [1:2], [3:4], [2,3,5], [1], [1:2], [3:4]};  %LS areas %did not analyze session 151021_img29

all_days = {'151021_img29', '151022_img29', '151009_img30', '151011_img30', '151211_img32', '151212_img32', '160129_img35', '160131_img35', '160129_img36', '160131_img36', '160314_img38', '160315_img38', '160319_img41', '160320_img41', '160606_img46', '160722_img53', }; %'150718_img27', '150719_img27', '150716_img28', '150717_img28', 
ROIcell = {[2], [2], [1:3], [1,3], [1:2], [1:2], [1:2], [1:2], [1:2], [1:2], [3:4], [2,3], [1], [1:2], [3:4], [1,2,3]}; % LS areas for all valid days of default bx

% all_days = {'160904_img55', '160905_img55', '160916_img61', '160918_img61', '160920_img61', '160921_img61', '161030_img62'};
% ROIcell = {[3:5], [3:5], [1:2], [1:2], [1:2], [1:2], [1:3]};

%load excel spreadsheet data
[~, ~, xl_spread] = xlsread('Z:\Data\WidefieldImaging\GCaMP\WF_exp_spreadsheet'); 
row_heading = cell(1,size(xl_spread,2));
%extract row headings
for i = 1:size(xl_spread,2)
    row_heading{1, i} = xl_spread{1, i};
end

%xl_struct = cell2struct(xl_spread, row_heading, 2);
for ii = 1%4:length(all_days);
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
        'corr_peak_latency_mean', [], 'early_peak_latency_mean', [],'cue_peak_latency_mean', [], 'corr_peak_latency_std', [], 'early_peak_latency_std', [], 'cue_peak_latency_std', [], ...
        'tooFast_peak_latency_release', [], 'tooFast_peak_latency_mean', [], 'tooFast_peak_latency_std', [], 'tooFast_magnitude', [], 'tooFast_magnitude_mean', [], 'tooFast_magnitude_std', []);
    
    %load some files
    b_file = dir([bx_dir, 'data-', subj_num, '-' session_date, '*']); %find and load behavior file
    load([bx_dir, b_file.name]); 
    bx_data = input; clear input; %bx file is stored as input so I rename it and clear input
    load([bx_outputs_dir, all_days{ii}, '_bx_outputs']);   %frame_info, lever, trial_outcome, cluster, lickTimes, ...
    
    %make sure counterValues are appropriate as are the first and last
    %trial nums
    old_f_frame_num = frame_info.f_frame_trial_num;
    if bx_data.counterValues{frame_info.f_frame_trial_num}(1)==2 & bx_data.counterValues{frame_info.f_frame_trial_num-1}(1)==1 & size(bx_data.counterValues{frame_info.f_frame_trial_num-1},2)==1;
        frame_info.f_frame_trial_num = frame_info.f_frame_trial_num-1; %there was one instance where the 1st frame occured in one trial and the second frame occured in the next trial
    end
    
    %load time of the first frame
    first_frame_time = frame_info.f_frame_MWorks_time;

    %find the trial nums for correct trials 
    success_times_trial_num = find(trial_outcome.corr_inx); % will include all trials, not just the imaged ones
    success_times_trial_num(success_times_trial_num <= frame_info.f_frame_trial_num) = [];
    success_times_trial_num(success_times_trial_num >= frame_info.l_frame_trial_num) = [];
    if length(success_times_trial_num) == size(trial_outcome.success_time,2)+1;  % sometimes WF_lever_plotting_TCs will chop off a trial from the success_roi because there were not enough frames after the lever release. When this happens it never edits trial_outcome.corr_inx
        success_times_trial_num(end) = [];  %
    end
    
    %use the trial nums for correct trials to get/store their reaction time data
    corr_react_times = cell2mat(bx_data.reactTimesMs(success_times_trial_num));
    curr_struct.corr_react_times = corr_react_times;  % store the reaction times for correct trial TCs in the structure
    curr_struct.corr_react_time_mean = round(mean(corr_react_times));
    curr_struct.corr_react_time_std = round(std(double(corr_react_times)));
    
   %find the trial nums for lapse and tooFast trials
    if ~isempty(trial_outcome.late_time)
        late_times_trial_num = find(trial_outcome.corr_inx); % will include all trials, not just the imaged ones
        late_times_trial_num(late_times_trial_num <= frame_info.f_frame_trial_num) = [];
        late_times_trial_num(late_times_trial_num >= frame_info.l_frame_trial_num) = [];
    else
        late_times_trial_num = [];
    end 
    if ~isempty(trial_outcome.tooFastCorrects)
        tooFast_times_trial_num = find(trial_outcome.corr_inx); % will include all trials, not just the imaged ones
        tooFast_times_trial_num(tooFast_times_trial_num <= frame_info.f_frame_trial_num) = [];
        tooFast_times_trial_num(tooFast_times_trial_num >= frame_info.l_frame_trial_num) = [];
    else
        tooFast_times_trial_num = [];
    end
    
    %find all the cue changes which were imaged 
    all_cue_times_trial_num = find(trial_outcome.corr_inx + trial_outcome.tooFast_inx);
    all_cue_times_trial_num(all_cue_times_trial_num <= frame_info.f_frame_trial_num) = [];
    all_cue_times_trial_num(all_cue_times_trial_num >= frame_info.l_frame_trial_num) = [];
    
    %use trial number of all cue trials to get reaction times 
    %cue_react_times = cell2mat(bx_data.reactTimesMs(all_cue_times_trial_num));
    cue_react_times = corr_react_times;
    curr_struct.cue_react_times = cue_react_times;  % store the reaction times for all cue trial TCs in the structure
    curr_struct.cue_react_time_mean = round(mean(cue_react_times));
    curr_struct.cue_react_time_std = round(std(double(cue_react_times)));
    
    %find hold duration for cue trials 
    cue_hold_times = cell2mat(bx_data.holdTimesMs(all_cue_times_trial_num));
    curr_struct.cue_hold_times = cue_hold_times;  % store the hold times for all cue trial TCs in the structure
    curr_struct.cue_hold_time_mean = round(mean(cue_hold_times));
    curr_struct.cue_hold_time_std = round(std(double(cue_hold_times)));
    curr_struct.cue_hold_dur = cue_hold_times; 
    
    %get time course data
    load([TC_dir, all_days{ii}, '_success']);    %stored as success_roi
    load([TC_dir, all_days{ii}, '_fail']);       %stored as fail_roi
    load([TC_dir, all_days{ii}, '_cue']);        %stored as fail_roi
    load([TC_dir, all_days{ii}, '_tooFast']);    %stored as tooFast
    success_roi = success_roi(:,ROIcell{ii},:); %remove non-LS areas
    fail_roi = fail_roi(:,ROIcell{ii},:);       %remove non-LS areas
    cue_roi = cue_roi(:,ROIcell{ii},:);         %remove non-LS areas
    tooFast_roi = tooFast_roi(:,ROIcell{ii},:); %remove non-LS areas
    
    %interpolate TC data and store in structure   %the new x-value for the lever release is 501
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
    if ~isempty(tooFast_roi)
        tooFast_roi_interp = nan(length(XVector), size(tooFast_roi,1), size(tooFast_roi,2)); %dims: 1=T2 2=trial# 3=ROI#
        for i = 1:size(tooFast_roi,2); %for each roi...
            tooFast_temp = squeeze(tooFast_roi(:,i,:))';  %dims 1=Time 2=trial number
            tooFast_temp = interp1(tooFast_temp, XVector);
            tooFast_roi_interp(:,:,i) = tooFast_temp;  %dims 1=Time 2=trialNumber  3=roi
        end
    end
    %if there is more than one LS ROI then average them together
    if size(ROIcell{ii},2) > 1
        success_roi_interp = mean(success_roi_interp,3); %average together the ROIs
        fail_roi_interp = mean(fail_roi_interp,3);
        cue_roi_interp = mean(cue_roi_interp,3);
        if ~isempty(tooFast_roi)
        tooFast_roi_interp = mean(tooFast_roi_interp,3);
        end
    end
    
    %baseline TCs 
    for iii = 1:size(success_roi_interp,2);  %dim1=time(10ms)  dim2=trial#
        shift = mean(success_roi_interp([1:150], iii));
        success_roi_interp(:, iii) = success_roi_interp(:,iii)-shift;
    end
    for iii = 1:size(fail_roi_interp,2);  %dim1=time(10ms)  dim2=trial#
        shift = mean(fail_roi_interp([1:150], iii));
        fail_roi_interp(:, iii) = fail_roi_interp(:,iii)-shift;
    end
    for iii = 1:size(cue_roi_interp,2);  %dim1=time(10ms)  dim2=trial#
        shift = mean(cue_roi_interp([250:400], iii));
        cue_roi_interp(:, iii) = cue_roi_interp(:,iii)-shift;
    end
    if ~isempty(tooFast_roi)
    for iii = 1:size(tooFast_roi_interp,2);  %dim1=time(10ms)  dim2=trial#
        shift = mean(tooFast_roi_interp([1:150], iii));
        tooFast_roi_interp(:, iii) = tooFast_roi_interp(:,iii)-shift;
    end
    end
    curr_struct.corr_TCs = success_roi_interp; %store TCs
    curr_struct.early_TCs = fail_roi_interp;
    curr_struct.cue_TCs = cue_roi_interp;
    if ~isempty(tooFast_roi)
    curr_struct.tooFast_TCs = tooFast_roi_interp;
    end
    
    %all cue trials correspond to successful trials. Thier indeces should be the same already. 
%     corr_cue_aligned_trial_num = success_times_trial_num; %the trial nums for corr_cue aligned and the successful conditions are the same 
%     [~, all_cue_ind_to_keep, ~] = intersect(all_cue_times_trial_num, corr_cue_aligned_trial_num);   
%     assert(length(all_cue_ind_to_keep)==size(cue_roi_interp,2))
    %corr_cue_aligned_TCs = cue_roi_interp(:,all_cue_ind_to_keep); %only keep the cue trials which correspond to correct trials 
    corr_cue_aligned_TCs = cue_roi_interp;
    assert(size(corr_cue_aligned_TCs,2) == size(success_roi_interp,2));
    curr_struct.corr_cue_aligned_TCs = corr_cue_aligned_TCs; 
    
    %find %corr and other variables using the excel spreadsheet
    [row_ind_day, col_ind_day] = find(strcmp(xl_spread, all_days{ii})); %find location of the row corresponding to this imaging session
    [row_ind, col_ind] = find(strcmp(xl_spread, 'percent correct')); %find the column for %corr
    this_percent_corr = cell2mat(xl_spread(row_ind_day, col_ind));
    curr_struct.percent_correct = this_percent_corr(1);
    
    [row_ind, col_ind] = find(strcmp(xl_spread, 'Training day')); %find the column for training day
    this_training_day = cell2mat(xl_spread(row_ind_day, col_ind));
    curr_struct.training_day_num = this_training_day(1);
    
    [row_ind, col_ind] = find(strcmp(xl_spread, 'peak percent correct')); %find the peak percent correct of any training day previous to the current imaging day
    this_peak_perc_corr = cell2mat(xl_spread(row_ind_day, col_ind));
    curr_struct.peak_percent_corr = this_peak_perc_corr(1);
    
    [row_ind, col_ind] = find(strcmp(xl_spread, 'avg percent corr pre-imaging')); %
    this_recent_perc_corr = cell2mat(xl_spread(row_ind_day, col_ind));
    curr_struct.recent_percent_corr = this_recent_perc_corr(1);
    
    %get hold duration data which matches TC data
    curr_struct.corr_hold_dur  = trial_outcome.succ_hold_dur;
    curr_struct.early_hold_dur = trial_outcome.fail_hold_dur;
    
    %find tbyt peak times and values for correct trials     SHOULD ADAPT SO IT USES AN AVG VALUE IN A 100MS WINDOW AROUND THE PEAK TIME??
    corr_tbyt_lat_peak = [];
    corr_tbyt_peak_mag = [];
    for i = 1:size(success_roi_interp,2);
        temp_val = find(success_roi_interp(:,i)==max(success_roi_interp([399:902],i)));  %finds the peak time of each trials TC
        if size(temp_val,1) > 1 %in case of finding multiple values it will select the first one
            temp_val = temp_val(find(temp_val >= 399,1, 'first'));
        end
        corr_tbyt_lat_peak = [corr_tbyt_lat_peak, temp_val]; %give peak time for each trial. Will occur in increments of 10 because the peak will always be on one of the actual frames. Not interpolated data points
        corr_tbyt_peak_mag = [corr_tbyt_peak_mag, success_roi_interp(temp_val,i)]; %collects all the trial by trial peak magnitudes
    end
    corr_tbyt_lat_peak = corr_tbyt_lat_peak-500;
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
        temp_val = find(fail_roi_interp(:,i)==max(fail_roi_interp([399:902],i)));
        if size(temp_val,1) > 1
            temp_val = temp_val(find(temp_val >= 399,1, 'first'));
        end
        early_tbyt_lat_peak = [early_tbyt_lat_peak, temp_val]; %give peak time for each trial
        early_tbyt_peak_mag = [early_tbyt_peak_mag, fail_roi_interp(temp_val,i)]; %collects all the trial by trial peaks 
    end
    early_tbyt_lat_peak = early_tbyt_lat_peak-500;
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
        temp_val = find(cue_roi_interp(:,i)==max(cue_roi_interp([599:1302],i)));
        if size(temp_val,1) > 1
            temp_val = temp_val(find(temp_val >= 599,1, 'first'));
        end
        cue_tbyt_lat_peak = [cue_tbyt_lat_peak, temp_val]; %give peak time for each trial
        cue_tbyt_peak_mag = [cue_tbyt_peak_mag, cue_roi_interp(temp_val,i)]; %collects all the trial by trial peaks 
    end
    cue_tbyt_lat_peak = cue_tbyt_lat_peak-500;
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
        temp_val = find(corr_cue_aligned_TCs(:,i)==max(corr_cue_aligned_TCs([599:1302],i)));
        if size(temp_val,1) > 1
            temp_val = temp_val(find(temp_val >= 599,1, 'first'));
        end
        corr_cue_tbyt_lat_peak = [corr_cue_tbyt_lat_peak, temp_val]; %give peak time for each trial
        corr_cue_tbyt_peak_mag = [corr_cue_tbyt_peak_mag, corr_cue_aligned_TCs(temp_val,i)]; %collects all the trial by trial peaks 
    end
    corr_cue_tbyt_lat_peak = corr_cue_tbyt_lat_peak-500;
    curr_struct.corr_cue_aligned_peak_latency = (corr_cue_tbyt_lat_peak); %have to convert it to ms and zero to cue
    curr_struct.corr_cue_aligned_peak_latency_mean = mean(corr_cue_tbyt_lat_peak);
    curr_struct.corr_cue_aligned_peak_latency_std = std(corr_cue_tbyt_lat_peak);
    curr_struct.corr_cue_aligned_magnitude = corr_cue_tbyt_peak_mag;
    curr_struct.corr_cue_aligned_magnitude_mean = mean(corr_cue_tbyt_peak_mag);
    curr_struct.corr_cue_aligned_magnitude_std = std(corr_cue_tbyt_peak_mag);
    
    %find tbyt peak times and values for tooFast trials     SHOULD ADAPT SO IT USES AN AVG VALUE IN A 100MS WINDOW AROUND THE PEAK TIME??
    if ~isempty(tooFast_roi)
    tooFast_tbyt_lat_peak = [];
    tooFast_tbyt_peak_mag = [];
    for i = 1:size(tooFast_roi_interp,2);
        temp_val = find(tooFast_roi_interp(:,i)==max(tooFast_roi_interp([399:902],i)));  %finds the peak time of each trials TC
        if size(temp_val,1) > 1 %in case of finding multiple values it will select the first one
            temp_val = temp_val(find(temp_val >= 399,1, 'first'));
        end
        tooFast_tbyt_lat_peak = [tooFast_tbyt_lat_peak, temp_val]; %give peak time for each trial. Will occur in increments of 10 because the peak will always be on one of the actual frames. Not interpolated data points
        tooFast_tbyt_peak_mag = [tooFast_tbyt_peak_mag, tooFast_roi_interp(temp_val,i)]; %collects all the trial by trial peak magnitudes
    end
    tooFast_tbyt_lat_peak = tooFast_tbyt_lat_peak-500;
    curr_struct.tooFast_peak_latency_release = tooFast_tbyt_lat_peak; %have to convert it to ms and zero to lever release
    curr_struct.tooFast_peak_latency_mean = mean(tooFast_tbyt_lat_peak);
    curr_struct.tooFast_peak_latency_std = std(tooFast_tbyt_lat_peak);
    curr_struct.tooFast_magnitude = tooFast_tbyt_peak_mag;
    curr_struct.tooFast_magnitude_mean = mean(tooFast_tbyt_peak_mag);
    curr_struct.tooFast_magnitude_std = std(tooFast_tbyt_peak_mag);
    end
    
    %find the rise times and onset times
    rise_time_data = rise_time_finder(curr_struct);
    
    %store the licking data
    curr_struct.corr_licking_TCs = licking_data.lick_trace_succ_10ms;
    curr_struct.corr_cue_aligned_licking_TCs = licking_data.lick_trace_cue_10ms;
    curr_struct.early_licking_TCs  = licking_data.lick_trace_fail_10ms;
    curr_struct.cue_licking_TCs  = licking_data.lick_trace_cue_10ms;
    
    %store various variables
    curr_struct.corr_rise_times = rise_time_data.corr_rise_times;
    curr_struct.corr_onset_latency = rise_time_data.corr_onset_times-500;
    curr_struct.early_rise_times = rise_time_data.early_rise_times;
    curr_struct.early_onset_latency = rise_time_data.early_onset_times-500;
    curr_struct.cue_rise_times = rise_time_data.cue_rise_times;
    curr_struct.cue_onset_latency = rise_time_data.cue_onset_times-500;
    curr_struct.corr_cue_aligned_rise_times = rise_time_data.corr_cue_aligned_rise_times;
    curr_struct.corr_cue_aligned_onset_latency = rise_time_data.corr_cue_aligned_onset_times-500;
    
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
