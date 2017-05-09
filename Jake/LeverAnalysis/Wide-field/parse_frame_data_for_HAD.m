function [frame, lever] = parse_frame_data_for_HAD(b_data, frame_times_dir);
%Processes the behavior data for frames, lever info, and licking data so it can be used for TCs later. 
%frame.counter, NaN if undefined. In 1 ms resolution

%Obtain frame info from behavior file  %seperate this process from parse_bx and make them their own functions within the overview
[counterValues, counterTimesMs, f_frame_trial_num, l_frame_trial_num] = extract_counters_bx(b_data);
[counterValues, counterTimesMs] = counter_fixer(counterValues, counterTimesMs);

%store frame data
frame.counter_by_time = counterValues; 
frame.times = counterTimesMs; 
frame.counter = counter_calculator(counterValues, counterTimesMs);  %each slot represents a ms of time during imaging. Each value in a slot represents the frame# being collected at that time. 
first_frame =1;
last_frame = counterValues(end); 

%This code looks at frame.times while it is still in free floating MWtimes. It extracts the values needed to cut baseline_timesMs to remove
%all events without associated frames and align it to frame.counter after frame.counter has been cut to align with the camera start.
f_frame_MWorks_time  = double(frame.times(1)); %time of the first frame in free floating MWorks time
l_frame_MWorks_time  = double(frame.times(end));
imaging_start_MW_T = double(f_frame_MWorks_time - mode(diff(counterTimesMs))); %finds the IFI and subtracts from time of first frame. This finds the time that imaging first began. Frames times are the times of the end of each period of photon collection. 

%align bx frame times to the time of the first frame
assert(length(frame.times) == length(frame.counter_by_time));  %check to make sure there are the same number of frame counters and times

%load frame_times
load(frame_times_dir);

%looks at frame times from the camera
cam_IFI = mode(diff(frame_times));
frame_times = round(frame_times - frame_times(1) + cam_IFI);  %aligns camera frame times to the time of imgaging initiation
counter_by_frame_times(frame_times) = 1:length(frame_times);  %creates a vector with length = to num ms in session and interspersed values equal to frame times
for i=1:length(frame_times)-1   %fills in all of the slots representing the ms during a given frame with that frame's number
    if i == 1;  %frame time records the time at the END of the duration of each frame
        inx = 1:frame_times(i);
    else
        inx = frame_times(i-1)+1:(frame_times(i));
    end
    counter_by_frame_times(inx) = i;
end

%define and store values for the 'frame' or 'frame_info' set of variables 
inx = find(frame.counter_by_time>=first_frame & frame.counter_by_time<=last_frame);
frame.counter_by_time = frame.counter_by_time(inx) - first_frame+1;
frame.times = frame.times(inx) - imaging_start_MW_T+1;
frame.f_frame_MWorks_time = f_frame_MWorks_time;
frame.l_frame_MWorks_time = l_frame_MWorks_time;
frame.imaging_start_MW_T = imaging_start_MW_T;

%compare the frame times from the camera to those from the bx
min_l = min(length(counter_by_frame_times), length(frame.counter));  %One has been corrected for missing frames etc and the other has not
bad_inx = find(abs(counter_by_frame_times(1:min_l) - double(frame.counter(1:min_l))) >=3);
if(length(bad_inx)/min_l > 0.01)  %if there is a difference of greater than 3ms between the camera and bx frame times for more than 1% of trials...
    warning('something wrong with the camera synchronization! see parse_frame_data_for_HAD')
end

%store counter for the camera and trim both counters. 
frame.counter_cam = counter_by_frame_times;
frame.counter_cam = frame.counter_cam - first_frame + 1;  %so that the first frame to be anaylzed = 1. 
frame.counter = frame.counter - first_frame + 1; 

%Determine baseline times: Time windows in which to take an F for df/f
baseline_timesMs = find_baseline_times(b_data);

%adapt baseline times so it only includes imaged frames and is zeroed to MWtime of camera start. Aligns baseline_times wiht frame.counter
baseline_timesMs = adapt_baseline_times(baseline_timesMs, imaging_start_MW_T, l_frame_MWorks_time, f_frame_trial_num, l_frame_trial_num);

%store baseline_timesMs and last/first frame trial nums
lever.baseline_timesMs = baseline_timesMs;
frame.f_frame_trial_num = f_frame_trial_num;
frame.l_frame_trial_num = l_frame_trial_num;
return