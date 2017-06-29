function [lever, frame, trial_outcome, licking_data] = cleanBehav(b_data, ifi, holdT_min)
%Processes the behavior data for frames, lever info, licking data, and
%trial outcome so it can be used for TCs later.

trial_outcome.ind_press_prerelease = [];

% --- Obtain frame info from behavior file

[counterValues, counterTimesMs] = counter_fixer(cell2mat(cellfun(@int64,b_data.counterValues,'UniformOutput',0)), cell2mat(cellfun(@int64,b_data.counterTimesUs,'UniformOutput',0))/1000);

%store frame data
frame.counter_by_time = counterValues;
frame.times = counterTimesMs;
frame.counter = counter_calculator(counterValues, counterTimesMs);  %each slot represents a ms of time during imaging. Each value in a slot represents the frame# being collected at that time.
last_frame = counterValues(end);

%This code looks at frame.times while it is still in free floating MWtimes. It extracts the values needed to cut baseline_timesMs to remove
%all events without associated frames and align it to frame.counter after frame.counter has been cut to align with the camera start.
f_frame_MWorks_time  = double(frame.times(1)); %time of the first frame in free floating MWorks time
l_frame_MWorks_time  = double(frame.times(end));
imaging_start_MW_T = f_frame_MWorks_time; %finds the IFI and subtracts from time of first frame. This finds the time that imaging first began. Frames times are the times of the end of each period of photon collection.

%======================================================================================
%use time of first frame to align licking times to start of imaging
lickTimes=[];
if isfield(b_data, 'lickometerTimesUs');  %img24 and 25 use datasets which have no licking fields
    for kk = 1:length(b_data.lickometerTimesUs);
        lickTimes = [lickTimes b_data.lickometerTimesUs{kk}/1000];
    end
    %lickTimes = lickTimes-f_frame_MWorks_time;
    lickTimes = double(lickTimes);
    frameTimes = frame.times; %frame times is monotonic by the time it gets here.
    
    licksByFrame = [];
%     lickStart = min(find((lickTimes(1) - frameTimes<0))) - 1;
    for i = 1:length(frameTimes)-1;
        
        licksThisFrame = sum(lickTimes>frameTimes(i) & lickTimes<=frameTimes(i+1));
        licksByFrame = [licksByFrame licksThisFrame];

    end
    licking_data = [];
    licking_data.lickTimes = lickTimes-imaging_start_MW_T;
    licking_data.licksByFrame = licksByFrame;
end
%=======================================================================


if b_data.doLeverSolenoidAllTrials == 1
    cLeverDown = cell2mat(b_data.cLeverDown);
    cLeverDown = cLeverDown - counterValues(1) + 1;
    cLeverUp = cell2mat(b_data.cLeverUp);
    cLeverUp = cLeverUp - counterValues(1) + 1;
    cLeverDown(cLeverDown == 0) = 1;
    cLeverUp(cLeverUp == 0) = 1;
    lever.press = counterTimesMs(cLeverDown);
    lever.release   = counterTimesMs(cLeverUp);
    
    lever.press = remove_events_and_adjust_inx(lever.press, imaging_start_MW_T, l_frame_MWorks_time);  %this is a nested function at the end of parse_bx
    lever.release = remove_events_and_adjust_inx(lever.release, imaging_start_MW_T, l_frame_MWorks_time);

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
    fidget_time = release_time(hold_time<200);
    late_time = release_time(find(strcmp('ignore',b_data.trialOutcomeCell)));
    
    %find all the correct trials which were too fast to be responding to visual cue
    tooFastCorrects = zeros(size(b_data.reactTimesMs));
    for i = 1:length(b_data.reactTimesMs)
        if react_time(i) < 200   %Removes tooFastCorrects from the successful trials count.
            has_reward(i)=0;
        end
        if isequal(b_data.trialOutcomeCell{i}, 'success');  %will have to change this in order to run img24/25
            tooFastCorrects(i) = react_time(i)>0 & react_time(i)<200; %definition of and isolation of tooFastCorrects
        end
    end
    
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
    % has_reward([1:f_frame_trial_num, l_frame_trial_num:end])=0;
    % early_time([1:f_frame_trial_num, l_frame_trial_num:end])=0;
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
    trial_outcome.late_time    = trial_outcome.late_time;
    
else
    cTargetOn = cell2mat(b_data.cTargetOn);
    cTargetOn = cTargetOn - counterValues(1) + 1;
    cLeverUp = cell2mat(b_data.cLeverUp);
    cLeverUp = cLeverUp - counterValues(1) + 1;
    cTargetOn(cTargetOn == 0) = 1;
    cLeverUp(cLeverUp == 0) = 1;
    lever.press = counterTimesMs(cTargetOn);
    lever.release   = counterTimesMs(cLeverUp);
    
    lever.press = remove_events_and_adjust_inx(lever.press, imaging_start_MW_T, l_frame_MWorks_time);  %this is a nested function at the end of parse_bx
    lever.release = remove_events_and_adjust_inx(lever.release, imaging_start_MW_T, l_frame_MWorks_time);
    
    
    omitRewardIndx = cell2mat(b_data.tRewardOmissionTrial);
    unexpRewardIndx = cell2mat(b_data.tDoNoStimulusChange);
    normalRewardIndx = ones(size(omitRewardIndx));
    normalRewardIndx(omitRewardIndx | unexpRewardIndx) = 0;
    
    trial_outcome.normalReward = lever.press(logical(normalRewardIndx));
    trial_outcome.omitReward = lever.press(logical(omitRewardIndx));
    trial_outcome.unexpReward = lever.press(logical(unexpRewardIndx));
    
end
% inx = find(frame.counter_by_time>=first_frame & frame.counter_by_time<=last_frame);
% frame.counter_by_time = frame.counter_by_time(inx) - first_frame+1;
% frame.times = frame.times(inx) - imaging_start_MW_T+1;
frame.f_frame_MWorks_time = f_frame_MWorks_time;
frame.imaging_start_MW_T = imaging_start_MW_T;

%store counter for the camera and trim both counters.
% frame.counter_cam = counter_by_frame_times;
% frame.counter_cam = frame.counter_cam;  %so that the first frame to be anaylzed = 1.
frame.counter = frame.counter;

%Determine baseline times: Time windows in which to take an F for df/f
[baseline_timesMs,  trial_outcome] = find_baseline_times_2P(b_data, trial_outcome, holdT_min);

%remove events from baseline_times which do not have associated frames. JH
for t = 1:length(baseline_timesMs)
    if baseline_timesMs(1,t)<imaging_start_MW_T;
        baseline_timesMs(:,t)=NaN;
    elseif baseline_timesMs(2,t)>l_frame_MWorks_time;
        baseline_timesMs(:,t)=NaN;
    end
end

% %convert 1st and last trials' baseline_times to NaNs
% baseline_timesMs(:,f_frame_trial_num)=NaN;    %the first and last trials with camera pulses will be incompletely imaged. Therefore they are unusable baseline_times
% baseline_timesMs(:,l_frame_trial_num)=NaN;

%subtract the MWtime of camera start in order to align baseline_timesMS to frame.counter
baseline_timesMs = baseline_timesMs - imaging_start_MW_T;
lever.baseline_timesMs = baseline_timesMs;
% frame.f_frame_trial_num = f_frame_trial_num;
% frame.l_frame_trial_num = l_frame_trial_num;
return;

function res = remove_events_and_adjust_inx(vec, imaging_start_MW_T, l_frame_MWorks_time)  %takes a given vector of events and removes the events which did not occur during imaging. Also subtracts off the time of the first frame
res = vec(vec >=imaging_start_MW_T & vec <= l_frame_MWorks_time)  - imaging_start_MW_T;