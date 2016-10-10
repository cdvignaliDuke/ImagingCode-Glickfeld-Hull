function [lever, frame, trial_outcome, lickTimes] = parse_behavior_for_HAD(b_data, first_frame, last_frame, frame_times, holdT_min)
%Processes the behavior data for frames, lever info, licking data, and
%trial outcome so it can be used for TCs later. 

trial_outcome.ind_press_prerelease = [];
trial_start = round(double(cell2mat(b_data.tThisTrialStartTimeMs)));  %finds each trial's start time in free floating MWorks time
num_trials  = length(b_data.tThisTrialStartTimeMs); %gets vectors for # of trials. Includes unimaged trials. 

% --- find state of lever before first trial
prev_trial_state = [];
for i=1:num_trials
    if(isempty(b_data.leverValues{i})); %looks for the first trial with lever values. Should only be necessary for animals without solenoid. 
        continue;
    end
    prev_trial_state = double(~(b_data.leverValues{i}(1)));        %prev_trial_state becomes important later on but I have no idea what it is doing here?????
    break;
end
if(isempty(prev_trial_state))
    warning('no lever press');
    prev_trial_state = 0;
end

% --- stores the time of the beginning of the first trial in MWorks time. 
bx_start_MWorks_time  = trial_start(1);   

% -- Obtain frame info from behavior file
[counterValues, counterTimesMs, f_frame_trial_num, l_frame_trial_num] = extract_counters_bx(b_data);
if ~ismonotonic(counterValues) | ~ismonotonic(counterTimesMs);
    disp('error detected in counterTimes or counterValues: diagnosing');
    [counterValues, counterTimesMs] = counter_fixer(counterValues, counterTimesMs);
end
frame.counter_by_time = counterValues; 
frame.times = counterTimesMs; 
frame.counter = counter_calculator(counterValues, counterTimesMs);  %each slot represents a ms of time during imaging. Each value in a slot represents the frame# being collected at that time. 
last_frame = counterValues(end);

% -- for loop which goes over all trials and extracts lever info
%NEED TO WRITE IN A CONDITION FOR ANALOG LEVERS?
lever.state   = [];
lever.release = [];
lever.press   = [];
for i=1:num_trials-1
    t_begin = trial_start(i);                                      %times for trial start, end
    t_end = trial_start(i+1);
    t_len = t_end - t_begin;                                       %trial duration
    
    % ------ get the vector of lever presses
    lever_time = ceil(double(b_data.leverTimesUs{i})/1000 - t_begin)+1;                %extracts lever times (relative to trial start) and values for the current trial
    lever_value = b_data.leverValues{i};
    bool_press = zeros(t_len,1) + prev_trial_state;  %init with prev trial state       for each millisecond of the trial bool_press will record lever state
    pre_press_time = 1;
    for j=1:length(lever_time)                                      %1 to # of lever events in that trial
        if(j==1)
            if(lever_value(1) == prev_trial_state)
                warning(['check this perhaps a lever press/release was missed on trial ' num2str(i)]);
            end
        end
        bool_press(pre_press_time:lever_time(j)-1) = prev_trial_state;   %all values in bool press remain unchanged until the first lever event
        pre_press_time  = lever_time(j);
        prev_trial_state =  double(lever_value(j));
        if(j==length(lever_time))                    %last press
            bool_press(lever_time(j):end) = prev_trial_state;
        end
    end
    
    % ----- save a binary vector with state of the lever
    lever.state(end+1:end+length(bool_press)) = bool_press;
    % ----- save release and press times
    release = lever_value == 0;
    release_time = double(b_data.leverTimesUs{i}(release))/1000  -  bx_start_MWorks_time;
    lever.release(end+1:end+length(release_time)) = release_time';
    
    press = ~release;
    press_time = double(b_data.leverTimesUs{i}(press))/1000  -  bx_start_MWorks_time;
    lever.press(end+1:end+length(press_time)) = press_time';
end

%This code looks at frame.times while it is still in free floating MWtimes. It extracts the values needed to cut baseline_timesMs to remove
%all events without associated frames and align it to frame.counter after frame.counter has been cut to align with the camera start.
f_frame_MWorks_time  = double(frame.times(1)); %time of the first frame in free floating MWorks time
l_frame_MWorks_time  = double(frame.times(end));
imaging_start_MW_T = double(f_frame_MWorks_time - mode(diff(counterTimesMs))); %finds the IFI and subtracts from time of first frame. This finds the time that imaging first began. Frames times are the times of the end of each period of photon collection. 

%======================================================================================
%use time of first frame to align licking times to start of imaging
lickTimes=[];
if isfield(b_data, 'lickometerTimesUs');  %img24 and 25 use datasets which have no licking fields 
    for kk = 1:length(b_data.lickometerTimesUs);
        lickTimes = [lickTimes cell2mat(b_data.lickometerTimesUs(kk))/1000];
    end
    %lickTimes = lickTimes-f_frame_MWorks_time;
    lickTimes = double(lickTimes)-imaging_start_MW_T;
end
%=======================================================================

%align bx frame times to the time of the first frame
assert(length(frame.times) == length(frame.counter_by_time));  %check to make sure there are the same number of frame counters and times
%frame.times = frame.times-f_frame_MWorks_time;   %zeroes frame.times to the beginning of the first trial in MWorks time

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

%Collects and calculates various categories of trial outcome
%does not trim the unimaged trials
hold_start = double(cell2mat(b_data.holdStartsMs)) - imaging_start_MW_T;
hold_time  = double(cell2mat(b_data.holdTimesMs));   %duration of the lever hold on that trial
react_time = double(cell2mat(b_data.reactTimesMs));
req_hold   = double(cell2mat(b_data.tTotalReqHoldTimeMs));
rnd_hold   = double(cell2mat(b_data.tRandReqHoldTimeMs));
tot_req_hold = req_hold + rnd_hold;
release_time = hold_start + hold_time;
early_time = hold_time<req_hold & hold_time>=200;  %finds non-fidget early releases. Fidget < 350ms hold.   %LARGE BUG in early times. early inx values are all non zero and negative. 
has_reward = ~cellfun(@isempty, b_data.juiceTimesMsCell ); %using this as a proxy to isolate correct trials 
fidget_time = release_time(hold_time<350);
late_time = release_time(find(strcmp('ignore',b_data.trialOutcomeCell)));

%find all the correct trials which were too fast to be responding to visual cue
tooFastCorrects = zeros(size(b_data.reactTimesMs)); 
for i = 1:length(b_data.reactTimesMs)
    if react_time(i) < 176   %Removes tooFastCorrects from the successful trials count. 
        has_reward(i)=0;
    end
    if isequal(b_data.trialOutcomeCell{i}, 'success');  %will have to change this in order to run img24/25
        tooFastCorrects(i) = react_time(i)>0 & react_time(i)<176; %definition of and isolation of tooFastCorrects
    end
end

%remove all events that are not between first and last frame, adjust indices
%frame.counter = frame.counter(f_frame_MWorks_time:l_frame_MWorks_time);  %trims frame.counter and lever.state to exclude unimaged timepoints
%lever.state = lever.state(f_frame_MWorks_time:l_frame_MWorks_time); %1 if lever pressed 0 if not, in ms resolution

inx = find(frame.counter_by_time>=first_frame & frame.counter_by_time<=last_frame);
frame.counter_by_time = frame.counter_by_time(inx) - first_frame+1;
frame.times = frame.times(inx) - imaging_start_MW_T+1;
frame.f_frame_MWorks_time = f_frame_MWorks_time;
frame.imaging_start_MW_T = imaging_start_MW_T;

%============================================================
%PERHAPS I CAN REMOVE LEVER EVENTS BASED ON TRIAL NUM INSTEAD OF EVENT TIME
%remove unimaged frames.
lever.press = remove_events_and_adjust_inx(lever.press, imaging_start_MW_T, l_frame_MWorks_time);  %this is a nested function at the end of parse_bx
lever.release = remove_events_and_adjust_inx(lever.release, imaging_start_MW_T, l_frame_MWorks_time);
%======================================================

%compare the frame times from the camera to those from the bx
min_l = min(length(counter_by_frame_times), length(frame.counter));  %One has been corrected for missing frames etc and the other has not
bad_inx = find(abs(counter_by_frame_times(1:min_l) - double(frame.counter(1:min_l))) >=3);
if(length(bad_inx)/min_l > 0.01)  %if there is a difference of greater than 3ms between the camera and bx frame times for more than 1% of trials...
    warning('something wrong with the camera synchronization!')
end

%store counter for the camera and trim both counters. 
frame.counter_cam = counter_by_frame_times;
frame.counter_cam = frame.counter_cam - first_frame + 1;  %so that the first frame to be anaylzed = 1. 
frame.counter = frame.counter - first_frame + 1; 

%Determine baseline times: Time windows in which to take an F for df/f
baseline_timesMs = find_baseline_times(b_data, trial_outcome, holdT_min);

%remove events from baseline_times which do not have associated frames. JH
for t = 1:length(baseline_timesMs)
    if baseline_timesMs(1,t)<imaging_start_MW_T;
        baseline_timesMs(:,t)=NaN;
    elseif baseline_timesMs(2,t)>l_frame_MWorks_time;
        baseline_timesMs(:,t)=NaN;
    end
end

%convert 1st and last trials' baseline_times to NaNs  
baseline_timesMs(:,f_frame_trial_num)=NaN;    %the first and last trials with camera pulses will be incompletely imaged. Therefore they are unusable baseline_times
baseline_timesMs(:,l_frame_trial_num)=NaN;

%subtract the MWtime of camera start in order to align baseline_timesMS to frame.counter
baseline_timesMs = baseline_timesMs - imaging_start_MW_T;
lever.baseline_timesMs = baseline_timesMs;
frame.f_frame_trial_num = f_frame_trial_num;
frame.l_frame_trial_num = l_frame_trial_num;

%calculate trial outcome for all categories for all trials. Put into indeces. Values = time of lever release zeroed relative to f frame time MW
early_inx = (early_time.*release_time);  %right now early time is just a bunch of 1s and 0s indexing the early trials. 
corr_inx = (has_reward.*release_time);
tooFast_inx = (tooFastCorrects.*release_time);
late_inx = ismember(release_time, late_time).*release_time;
fidget_inx = ismember(release_time, fidget_time).*release_time;

%store the indexed vectors
trial_outcome.early_inx = early_inx;
trial_outcome.corr_inx = corr_inx;
trial_outcome.tooFast_inx = tooFast_inx;
trial_outcome.late_inx = late_inx;
trial_outcome.fidget_inx = fidget_inx;

%remove unimaged frames from hold duration variables and store them in trial_outcome
has_reward([1:f_frame_trial_num, l_frame_trial_num:end])=0;
early_time([1:f_frame_trial_num, l_frame_trial_num:end])=0;
trial_outcome.succ_hold_dur = hold_time(has_reward);
trial_outcome.fail_hold_dur = hold_time(early_time);

%store the various categories of trial outcome in the for of the lever times for only those categories agnostic of trial num
trial_outcome.success_time = release_time(has_reward);
trial_outcome.tooFastCorrects = release_time(1,find(tooFastCorrects)); %isolates correct trials with reaction times >0 but <175ms. Stores them as a matrix of their release times.
trial_outcome.early_time = release_time(early_time);                   %altered to exclude fidgets 2/27/16
trial_outcome.fidget = fidget_time;
is_ignore =@(x)isequal(x, 'ignore');
trial_outcome.late_time  = late_time; 
trial_outcome.change_orientation = hold_start + req_hold;  %zeroed relative to first frame time in MWorks time
trial_outcome.change_orientation(react_time < 0) = NaN;               %if respond before change then no change in orientation
if b_data.doLever == 0;                             %CONTROL trials 
    trial_outcome.early_time = release_time(hold_time<req_hold);
    trial_outcome.fidget = [];
end

%subtract off the time of the first frame from the lever times
trial_outcome.late_time    = trial_outcome.late_time - imaging_start_MW_T;

%explanation of certain variabels
%lever.state 1 if lever down 0 if up, in 1ms increments
%frame.counter, NaN if undefined. In 1 ms resolution
%trial_outcome.success_time   lever release for successes 
%trial_outcome.early_time     time of lever release in early trials 
%trial_outcome.late_time      time of lever solenoid coming up.
return; 

function res = remove_events_and_adjust_inx(vec, imaging_start_MW_T, l_frame_MWorks_time)  %takes a given vector of events and removes the events which did not occur during imaging. Also subtracts off the time of the first frame
res = vec(vec >=imaging_start_MW_T & vec <= l_frame_MWorks_time)  - imaging_start_MW_T;