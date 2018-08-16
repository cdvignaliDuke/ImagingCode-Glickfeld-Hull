function [lever, frame, trial_outcome, licking_data] = cleanBehav_CRP(b_data, ifi, holdT_min)
%Processes the behavior data for frames, lever info, licking data, and
%trial outcome so it can be used for TCs later.

trial_outcome.ind_press_prerelease = [];

% --- Obtain frame info from behavior file

[counterValues, counterTimesMs, shift_info] = counter_fixer_2P_CRP(cell2mat(cellfun(@int64,b_data.counterValues,'UniformOutput',0)), cell2mat(cellfun(@int64,b_data.counterTimesUs,'UniformOutput',0))/1000, b_data.saveTime(1:6), b_data.subjectNum);

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
doLick = sum(cell2mat(cellfun(@int64,b_data.lickometerTimesUs,'UniformOutput',0)));
licking_data = [];

if isfield(b_data, 'lickometerTimesUs') && doLick ~= 0  %img24 and 25 use datasets which have no licking fields
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
    
    licking_data.lickTimes = lickTimes-imaging_start_MW_T;
    licking_data.licksByFrame = licksByFrame;
    
    % find lick bout during ITI
    % lick bout criteria: at least 3 licks and each within 6 frames
    % and no lick during 300ms before lick bout start
    licking_data.bout_onset = []; licking_data.bout_rate= []; licking_data.bout_onsetTime = [];
    licking_data.single_lickFrame = []; licking_data.single_lickTime = [];
    licking_data.single_lickFrameAll = []; licking_data.single_lickTimeAll = [];
    licking_data.lickRate = [];
    buffer_lick = round(500/double(ifi));
    buffer_press = round(500/double(ifi));
    licking_data.buffer_lick = buffer_lick;
    for i = 2:length(b_data.cTrialStart)-1
        
        lick_itiTrial = licksByFrame(b_data.cTrialStart{i} : b_data.cLeverDown{i} - buffer_press);
        lick_itiTrial2 = licksByFrame(b_data.cTrialStart{i} : b_data.cTrialStart{i} + b_data.tItiWaitFrames{i});
        lick_idx = find(lick_itiTrial);
        lick_idx2 = find(lick_itiTrial2);
        
        for f = 1:length(lick_idx2)
            single_idx = (lick_idx2(f) + b_data.cTrialStart{i}-1);
            % to be removed from spontaneous events
            licking_data.single_lickFrameAll = [licking_data.single_lickFrameAll single_idx];
            licking_data.single_lickTimeAll = [licking_data.single_lickTimeAll counterTimesMs(single_idx)];
            
            if sum(licksByFrame(single_idx - buffer_lick : single_idx -1)) < 1 && sum(licksByFrame(single_idx + 1: single_idx + buffer_lick)) < 1
                licking_data.single_lickFrame = [licking_data.single_lickFrame single_idx];
                licking_data.single_lickTime = [licking_data.single_lickTime counterTimesMs(single_idx)];
            end
        end
        
        if length(lick_idx) >= 3
            inter_lick = diff(lick_idx);
            %i
            for li = 1:length(inter_lick)
                % li
                if length(inter_lick(li:end)) >= 3 && mean(inter_lick(li: li+2 )) <= round(300/double(ifi)) % each lick is less than 300ms apart
                    % and at least three licks 
                    onset_idx = lick_idx(li) + b_data.cTrialStart{i} -1;
                    if sum(licksByFrame(onset_idx - buffer_lick : onset_idx -1)) < 1 % no licking 500ms before the onset
                        licking_data.bout_onset = [licking_data.bout_onset onset_idx]; % match index with trigger_licks
                        licking_data.bout_onsetTime = [licking_data.bout_onsetTime counterTimesMs(onset_idx)];
                        licking_data.bout_rate = [licking_data.bout_rate 1/mean(inter_lick(li: li+2))/double(ifi)*1000];
%                         licking_data.bout = cat(licking_data.bout, 3, lick_idx-lick_idx(1) + onset_idx);
                        break;
                    end
                    
                end
            end
        end
    end
    licking_data.bout_onsetTime = remove_events_and_adjust_inx(licking_data.bout_onsetTime, imaging_start_MW_T, l_frame_MWorks_time); 
    licking_data.single_lickTime = remove_events_and_adjust_inx(licking_data.single_lickTime, imaging_start_MW_T, l_frame_MWorks_time); 
    licking_data.single_lickTimeAll = remove_events_and_adjust_inx(licking_data.single_lickTimeAll, imaging_start_MW_T, l_frame_MWorks_time);
else
    licking_data = [];
end

%========================================================

    cTargetOn = cell2mat(b_data.cTargetOn); % cue onset frame
    cTargetOn = cTargetOn - counterValues(1) + 1;
    cLeverUp = cell2mat(b_data.cLeverUp);  % Reward frame
    cLeverUp = cLeverUp - counterValues(1) + 1;
    cTrialEnd = cell2mat(b_data.cTrialEnd);
    cTrialStart = cell2mat(b_data.cTrialStart);
    if ~isempty(shift_info)
        mod_vals_ind = find(cTargetOn >= shift_info.first_bad_frame);
        cTargetOn(mod_vals_ind) = cTargetOn(mod_vals_ind) - shift_info.frame_num_shift; 
        mod_vals_ind = find(cLeverUp >= shift_info.first_bad_frame);
        cLeverUp(mod_vals_ind) = cLeverUp(mod_vals_ind) - shift_info.frame_num_shift; 
        mod_vals_ind = find(cTrialEnd >= shift_info.first_bad_frame);
        cTrialEnd(mod_vals_ind) = cTrialEnd(mod_vals_ind) - shift_info.frame_num_shift; 
        mod_vals_ind = find(cTrialStart >= shift_info.first_bad_frame);
        cTrialStart(mod_vals_ind) = cTrialStart(mod_vals_ind) - shift_info.frame_num_shift; 
    end
   
    lever.press = counterTimesMs(cTargetOn); % cue onset counter times
    lever.release   = counterTimesMs(cLeverUp); % reward counter times
    
    lever.press = remove_events_and_adjust_inx(lever.press, imaging_start_MW_T, l_frame_MWorks_time);  %this is a nested function at the end of parse_bx
    lever.release = remove_events_and_adjust_inx(lever.release, imaging_start_MW_T, l_frame_MWorks_time);
    
    eIdx = cellfun(@isempty, b_data.tRewardOmissionTrial);  %fix empty cell 
    temp = b_data.tRewardOmissionTrial;
    temp(eIdx) = {int64(0)};
    omitRewardIndx = cell2mat(temp);
    
    unexpRewardIndx = cell2mat(b_data.tDoNoStimulusChange);
    normalRewardIndx = ones(size(omitRewardIndx));
    normalRewardIndx(omitRewardIndx | unexpRewardIndx) = 0;
    
    % get reward timestamp
    trial_outcome.normalRewardCue = lever.press(logical(normalRewardIndx));
    trial_outcome.omitRewardCue = lever.press(logical(omitRewardIndx));
    trial_outcome.unexpRewardCue = lever.press(logical(unexpRewardIndx));
    
    delayTime = double(cell2mat(b_data.tRewardDelayDurationMs));
    
    trial_outcome.normalReward = lever.release(logical(normalRewardIndx)) + delayTime(logical(normalRewardIndx));
    trial_outcome.omitReward = lever.release(logical(omitRewardIndx)) + delayTime(logical(omitRewardIndx));
    trial_outcome.unexpReward = lever.release(logical(unexpRewardIndx)) + delayTime(logical(unexpRewardIndx));
    
    trial_outcome.trialLen = cTrialEnd - cTrialStart %========================================
    trial_outcome.omitRewardIndx = omitRewardIndx;

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
[baseline_timesMs] = find_baseline_times_2P_CRP(b_data);

%remove events from baseline_times which do not have associated frames. JHl
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