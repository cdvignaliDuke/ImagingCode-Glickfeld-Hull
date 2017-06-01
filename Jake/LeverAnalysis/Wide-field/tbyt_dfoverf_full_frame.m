function img_dfoverf = tbyt_dfoverf_full_frame(b_data, bx_out_dir, img);
%Uses the trial by trial baseline times to create a df/f timecourse

%load variables
load(bx_out_dir, 'lever', 'frame_info');

shift = 10; %this is the number of times which t_range is shifted. In lever trials the trial ends almost immediately after lever release. In CRP trials it ends 6s after the cue termination
baseline_timesMs = lever.baseline_timesMs;
first_baseline = find(~isnan(baseline_timesMs(1,:)),1, 'first');    %find the first trial / baseline_timesMs window that is not NaN
StartT = frame_info.imaging_start_MW_T; %time of imaging onset in MWorks time. 
img = double(img);

img_dfoverf = nan(size(img));
F_range = [];
for iT=frame_info.f_frame_trial_num+1: frame_info.l_frame_trial_num-1;    %only looks at the first and last fully imaged trials
    
    %finding F_range (F_range is the # of each frame which we will use to generate f. It uses an iti based f)
    if ~isnan(baseline_timesMs(1,iT));   %if there is a valid baseline interval then make that the new F
        F_range = frame_info.counter(baseline_timesMs(1,iT)):frame_info.counter(baseline_timesMs(2,iT));
    elseif isempty(F_range)   %if these are trials before there was ever a valid F_range then use the first valid F_range as this trials F_range also.
        F_range = frame_info.counter(baseline_timesMs(1,first_baseline)):frame_info.counter(baseline_timesMs(2,first_baseline));
    end    %if there was no valid F_range but there was previously a valid F_range then F_range will remain unchanged and the most recent one will be use.
    
    %Find t_range (the time interval encompasing the whole trial. Used to find the frames to which we will apply this df/f).
    if iT == frame_info.f_frame_trial_num+1; %if this is the first fully imaged trial then t_range includes all frames up to this point
        t_range = 1:(frame_info.counter(round(cell2mat(b_data.tThisTrialStartTimeMs(iT+1)))-StartT)-1);
    else
        t_range = frame_info.counter(round(cell2mat(b_data.tThisTrialStartTimeMs(iT))-StartT)):(frame_info.counter(round(cell2mat(b_data.tThisTrialStartTimeMs(iT+1))-StartT))-1);
    end
    if iT == frame_info.l_frame_trial_num-1; %if this is the last fully imaged trial
        t_range = (t_range(1)+shift):size(img,2); %then apply this F to all subsequent frames (with a frame shift)
    elseif iT == frame_info.f_frame_trial_num+1; %if this is the first fully imaged trial then apply this F to all frames up to this point with the shift
        t_range = 1:(t_range(end)+shift); 
    else
        t_range = t_range + shift;    %added this shift because we have a 2s analysis window post release but trial ends immediately after release. 
    end
    
    %calculate df/f
    F_avg= mean(img(:,F_range),2);
    t_df = bsxfun(@minus, img(:, t_range), F_avg);   %use bsxfun because sometimes F_avg is a vector
    t_dfoverf = bsxfun(@rdivide, t_df, F_avg);
    img_dfoverf(:,t_range) = t_dfoverf;
    iT
end