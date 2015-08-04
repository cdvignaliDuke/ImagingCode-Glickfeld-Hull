function [ frame_times frame_count sub_input ] = get_frame_time_by_counters(input, info)
% get the frame timing information from counter event times
 nframes = info.config.frames;
    % find start of collection
    ntrials = size(input.trialOutcomeCell,2);
    for i = 1:ntrials
        if ~isempty(input.counterValues{i});
            trial_start = i;
            break
        end
    end

    % find end of collection
    trial_end = 1;
    for i = 1:ntrials
        if find(input.counterValues{i}==nframes);
            trial_end = i;
            break
        end
    end
    % in case behavior ends before imaging
    if trial_end == 1;
        trial_end = ntrials;
    end
    
    % create structure with only imaged trials
    sub_input = trialChopper(input, [trial_start trial_end]);
    ntrials = size(sub_input.trialOutcomeCell,2);
    
    % find all frame times and counts
    frames_temp.time = [];
    frames_temp.count = [];
    for i = 1:ntrials
        frames_temp.time = [frames_temp.time round(sub_input.counterTimesUs{i}./1000)];
        frames_temp.count = [frames_temp.count sub_input.counterValues{i}];
    end

    % correct for repeated frames
    y = diff(frames_temp.count);
    rep_frame = find(y==0)+1;
    frames_temp.time_fix = frames_temp.time;
    frames_temp.time_fix(rep_frame) = [];
    frames_temp.count_fix = frames_temp.count;
    frames_temp.count_fix(rep_frame) = [];

    % correct for skipped frames
    y = diff(frames_temp.count_fix);
    x = diff(frames_temp.time_fix);
    miss_frame = find(y==2);
    avg_ifi = mean(x(1:100),2);
    for i = 1:length(miss_frame)
        frames_temp.count_fix = [frames_temp.count_fix(1:(miss_frame(i)+i-1)) (miss_frame(i)+i) frames_temp.count_fix(miss_frame(i)+i:end)];
        frames_temp.time_fix = [frames_temp.time_fix(1:(miss_frame(i)+i-1)) (frames_temp.time_fix(miss_frame(i)+i-1)+avg_ifi) frames_temp.time_fix(miss_frame(i)+i:end)];
    end

    %zero time stamp to start of first frame
    frame_times = frames_temp.time_fix-frames_temp.time_fix(1);
    frame_count = frames_temp.count_fix;
end

