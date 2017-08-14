function rise_time_data = rise_time_finder(curr_struct);
corr_TCs  = curr_struct.corr_TCs;
early_TCs = curr_struct.early_TCs;
cue_TCs   = curr_struct.cue_TCs;
corr_cue_aligned_TCs = curr_struct.corr_cue_aligned_TCs;
corr_mag  = curr_struct.corr_magnitude;
early_mag = curr_struct.early_magnitude;
cue_mag   = curr_struct.cue_magnitude;
corr_cue_aligned_mag = curr_struct.corr_cue_aligned_magnitude;
corr_lat  = curr_struct.corr_peak_latency_release+500;
early_lat = curr_struct.early_peak_latency_release+500;
cue_lat   = curr_struct.cue_peak_latency+500;
corr_cue_aligned_lat = curr_struct.corr_cue_aligned_peak_latency+500;

rise_time_data = struct('corr_rise_times', [], 'corr_onset_times', [],'early_rise_times', [], 'early_onset_times', [], 'cue_rise_times', [], 'cue_onset_times', [],'corr_cue_aligned_rise_times', [], 'corr_cue_aligned_onset_times', []);

%find baseline values and indeces for each trial's F response
corr_TC_baselines = [];
corr_min_ind_vec = []; % this is a vector which will collect the indexed value of the min value calculated for each trial
for k = 1:size(corr_TCs,2)
    corr_min_ind = find(corr_TCs([corr_lat(k)-350:corr_lat(k)],k) == min(corr_TCs([corr_lat(k)-350:corr_lat(k)],k)),1,'last') + corr_lat(k)-350; %find the minimum values within a 350ms window before the peak
    if corr_min_ind < 151;
        corr_TC_baselines = [corr_TC_baselines, mean(corr_TCs(1:150,k))];
    else
        corr_TC_baselines = [corr_TC_baselines, mean(corr_TCs(corr_min_ind-150:corr_min_ind,k))]; %the baseline is the mean of the 150ms preceding the minimum 
    end
    corr_min_ind_vec = [corr_min_ind_vec, corr_min_ind];
end

early_TC_baselines = [];
early_min_ind_vec = [];
for k = 1:size(early_TCs,2)
    early_min_ind = find(early_TCs([early_lat(k)-350:early_lat(k)],k) == min(early_TCs([early_lat(k)-350:early_lat(k)],k)),1,'last') + early_lat(k)-350; %find the minimum values within a 350ms window before the peak
    if early_min_ind < 151;
        early_TC_baselines = [early_TC_baselines, mean(early_TCs(1:150,k))]; %if early_min_ind is less than 151 then the baseline is taken from the first 150ms
    else
        early_TC_baselines = [early_TC_baselines, mean(early_TCs(early_min_ind-150:early_min_ind,k))]; %the baseline is the mean of the 150ms preceding the minimum 
    end
    early_min_ind_vec = [early_min_ind_vec, early_min_ind];
end

cue_TC_baselines = [];
cue_min_ind_vec = []; 
for k = 1:size(cue_TCs,2)
    cue_min_ind = find(cue_TCs([cue_lat(k)-350:cue_lat(k)],k) == min(cue_TCs([cue_lat(k)-350:cue_lat(k)],k)),1,'last') + cue_lat(k)-350; %find the minimum values within a 350ms window before the peak
    cue_TC_baselines = [cue_TC_baselines, mean(cue_TCs(cue_min_ind-150:cue_min_ind,k))]; %the baseline is the mean of the 150ms preceding the minimum 
    cue_min_ind_vec = [cue_min_ind_vec, cue_min_ind];
end

corr_cue_TC_baselines = [];
corr_cue_min_ind_vec = []; 
for k = 1:size(corr_cue_aligned_TCs,2)
    corr_cue_min_ind = find(corr_cue_aligned_TCs([corr_cue_aligned_lat(k)-350:corr_cue_aligned_lat(k)],k) == min(corr_cue_aligned_TCs([corr_cue_aligned_lat(k)-350:corr_cue_aligned_lat(k)],k)),1,'last') + corr_cue_aligned_lat(k)-350; %find the minimum values within a 350ms window before the peak
    if corr_cue_min_ind < 151;
        corr_cue_TC_baselines = [corr_cue_TC_baselines, mean(corr_cue_aligned_TCs(1:150,k))]; 
    else
        corr_cue_TC_baselines = [corr_cue_TC_baselines, mean(corr_cue_aligned_TCs(corr_cue_min_ind-150:corr_cue_min_ind,k))]; %the baseline is the mean of the 150ms preceding the minimum 
    end
    corr_cue_min_ind_vec = [corr_cue_min_ind_vec, corr_cue_min_ind];
end

% find the 10% and 90% values and times. baseline=0% peak=100%
% need to take into account that some values could cross zero. 
% code needs to find the difference between the peak and baseline values.
% Then it will look the left of the peak time and log the 1st time the
% value hits 90% or less of the peak and the 1st time at which it hits 10%
% or less of the peak. 

%rise time data for correct trials 
corr_peak_to_base_diff = corr_mag-corr_TC_baselines; %find the difference between the peak magnitude and the magnitude of the baseline F
corr_ninety_val = (0.9*corr_peak_to_base_diff)+corr_TC_baselines;
corr_ten_val    = (0.1*corr_peak_to_base_diff)+corr_TC_baselines;
corr_ninety_time = [];
corr_ten_time = [];
for k = 1:size(corr_TCs,2);    
    corr_ninety_time = [corr_ninety_time, find(corr_TCs([1:corr_lat(k)],k)<=corr_ninety_val(k),1,'last')];  %find the last value from 1 to the time of the peak where the TC is less than the 90% value. 
    corr_ten_time    = [corr_ten_time, find(corr_TCs([1:corr_lat(k)],k)<=corr_ten_val(k),1,'last')];  %find the last value from 1 to the time of the peak where the TC is less than the 10% value. 
    assert(corr_ninety_time(k)<=corr_lat(k) | corr_ninety_time(k)>=corr_min_ind_vec(k));
    assert(corr_ten_time(k)<=corr_lat(k) | corr_ten_time(k)>=corr_min_ind_vec(k));
end
rise_time_data.corr_rise_times = corr_ninety_time-corr_ten_time; 
rise_time_data.corr_onset_times = corr_ten_time;

%Rise time data for early trials
early_ninety_time = early_mag-early_TC_baselines; %find the difference between the peak magnitude and the magnitude of the baseline F
early_ninety_val = (0.9*early_ninety_time)+early_TC_baselines;
early_ten_val    = (0.1*early_ninety_time)+early_TC_baselines;
early_ninety_time = [];
early_ten_time = [];
for k = 1:size(early_TCs,2);
    early_ninety_time = [early_ninety_time, find(early_TCs([1:early_lat(k)],k)<=early_ninety_val(k),1,'last')];  %find the last value from 1 to the time of the peak where the TC is less than the 90% value. 
    early_ten_time    = [early_ten_time, find(early_TCs([1:early_lat(k)],k)<=early_ten_val(k),1,'last')];  %find the last value from 1 to the time of the peak where the TC is less than the 10% value. 
    assert(early_ninety_time(k)<=early_lat(k) | early_ninety_time(k)>=early_min_ind_vec(k));
    assert(early_ten_time(k)<=early_lat(k) | early_ten_time(k)>=early_min_ind_vec(k));
end
rise_time_data.early_rise_times = early_ninety_time-early_ten_time; 
rise_time_data.early_onset_times = early_ten_time;

%rise time data for all cue
cue_ninety_time = cue_mag-cue_TC_baselines; %find the difference between the peak magnitude and the magnitude of the baseline F
cue_ninety_val = (0.9*cue_ninety_time)+cue_TC_baselines;
cue_ten_val    = (0.1*cue_ninety_time)+cue_TC_baselines;
cue_ninety_time = [];
cue_ten_time = [];
for k = 1:size(cue_TCs,2);
    cue_ninety_time = [cue_ninety_time, find(cue_TCs([1:cue_lat(k)],k)<=cue_ninety_val(k),1,'last')];  %find the last value from 1 to the time of the peak where the TC is less than the 90% value. 
    cue_ten_time    = [cue_ten_time, find(cue_TCs([1:cue_lat(k)],k)<=cue_ten_val(k),1,'last')];  %find the last value from 1 to the time of the peak where the TC is less than the 10% value. 
    assert(cue_ninety_time(k)<=cue_lat(k) | cue_ninety_time(k)>=cue_min_ind_vec(k));
    assert(cue_ten_time(k)<=cue_lat(k) | cue_ten_time(k)>=cue_min_ind_vec(k));
end
rise_time_data.cue_rise_times  = cue_ninety_time-cue_ten_time; 
rise_time_data.cue_onset_times = cue_ten_time;

%rise time data for all cue
corr_cue_ninety_time = corr_cue_aligned_mag-corr_cue_TC_baselines; %find the difference between the peak magnitude and the magnitude of the baseline F
corr_cue_ninety_val = (0.9*corr_cue_ninety_time)+corr_cue_TC_baselines;
corr_cue_ten_val    = (0.1*corr_cue_ninety_time)+corr_cue_TC_baselines;
corr_cue_ninety_time = [];
corr_cue_ten_time = [];
for k = 1:size(corr_cue_aligned_TCs,2);
    corr_cue_ninety_time = [corr_cue_ninety_time, find(corr_cue_aligned_TCs([1:corr_cue_aligned_lat(k)],k)<=corr_cue_ninety_val(k),1,'last')];  %find the last value from 1 to the time of the peak where the TC is less than the 90% value. 
    corr_cue_ten_time    = [corr_cue_ten_time, find(corr_cue_aligned_TCs([1:corr_cue_aligned_lat(k)],k)<=corr_cue_ten_val(k),1,'last')];  %find the last value from 1 to the time of the peak where the TC is less than the 10% value. 
    assert(corr_cue_ninety_time(k)<=corr_cue_aligned_lat(k) | corr_cue_ninety_time(k)>=corr_cue_min_ind_vec(k));
    assert(corr_cue_ten_time(k)<=corr_cue_aligned_lat(k) | corr_cue_ten_time(k)>=corr_cue_min_ind_vec(k));
end
rise_time_data.corr_cue_aligned_rise_times  = corr_cue_ninety_time-corr_cue_ten_time;
rise_time_data.corr_cue_aligned_onset_times = corr_cue_ten_time;
end