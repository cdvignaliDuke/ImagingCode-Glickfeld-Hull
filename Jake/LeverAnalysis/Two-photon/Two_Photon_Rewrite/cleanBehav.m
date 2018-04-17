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
%=======================================================================


if b_data.doLeverSolenoidAllTrials == 1 || b_data.doLever == 1
    cTrialStart = cell2mat(b_data.cTrialStart);
    cTrialStart = cTrialStart - counterValues(1) + 1;
    cTrialEnd = cell2mat(b_data.cTrialEnd);
    cTrialEnd = cTrialEnd - counterValues(1) + 1;
    if isfield(b_data, 'tItiWaitFrames')
        
        cItiReward = cTrialStart + int64(cell2mat(b_data.tItiWaitFrames)/2);
    else
        cItiReward = cTrialStart + int64(cell2mat(b_data.tItiWaitTimeMs)/(2*double(ifi)));
    end
    
    cItiReward(cItiReward > cTrialEnd) = [];
    
    cItiReward = unique(cItiReward);
    cLeverDown = cell2mat(b_data.cLeverDown);
    cLeverDown = cLeverDown - counterValues(1) + 1;
    cLeverUp = cell2mat(b_data.cLeverUp);
    cLeverUp = cLeverUp - counterValues(1) + 1;
    cLeverDown = unique(cLeverDown);
    [cLeverUp, icc] = unique(cLeverUp);
    cLeverDown(cLeverDown == 0) = 1;
    cLeverUp(cLeverUp == 0) = 1;
    lever.press = counterTimesMs(cLeverDown);
    lever.release   = counterTimesMs(cLeverUp);
    lever.itiReward = counterTimesMs(cItiReward);
    
    cueOn = cellfun(@isempty, b_data.cTargetOn);
    cTargetOn = b_data.cTargetOn;
    cTargetOn(cueOn) = {int64(1)};
    cTargetOn = cell2mat(cTargetOn);
    
    cTargetOn = cTargetOn - counterValues(1) + 1;
    
    cTargetOn(cTargetOn <=0) = 1;
    
    lever.cue = counterTimesMs(cTargetOn);

    lever.press = remove_events_and_adjust_inx(lever.press, imaging_start_MW_T, l_frame_MWorks_time);  %this is a nested function at the end of parse_bx
    lever.release = remove_events_and_adjust_inx(lever.release, imaging_start_MW_T, l_frame_MWorks_time);
    lever.cue = remove_events_and_adjust_inx(lever.cue, imaging_start_MW_T, l_frame_MWorks_time);
    lever.itiReward = remove_events_and_adjust_inx(lever.itiReward, imaging_start_MW_T, l_frame_MWorks_time);
    
    reward_trial = ~cellfun(@isempty, regexp(b_data.trialOutcomeCell, 'success') );
    
%     if length(reward_trial) ~= lever.release
        reward_trial_reduce = reward_trial(icc);
        lever.cue_release = lever.release(reward_trial_reduce);
        lever.cue = lever.cue(icc);
%     else
%         lever.cue_release = lever.release(~cellfun(@isempty, b_data.juiceTimesMsCell ));
%     end
    
    lever.cue2release = lever.release - lever.cue;
    lever.cue2release = lever.cue2release(reward_trial_reduce);
    
    if ~isempty(licking_data)
        for ll = 1:length(cLeverUp) 
            if cLeverUp(ll) + round(500 / double(ifi)) > length(licksByFrame)
                licking_data.lickRate = [licking_data.lickRate 1/mean(licksByFrame(cLeverUp(ll):end))/double(ifi)*1000];
            else
                licking_data.lickRate = [licking_data.lickRate sum(licksByFrame(cLeverUp(ll) : cLeverUp(ll) + round(500 / double(ifi))))/0.5];
            end
        end
        licking_data.lickRate(isnan(licking_data.lickRate)) = 0;
    end
    
    
    
%     lever.cue_release = remove_events_and_adjust_inx(lever.cue_release, imaging_start_MW_T, l_frame_MWorks_time);
    
    %Collects and calculates various categories of trial outcome
    %does not trim the unimaged trials
    hold_start = double(cell2mat(b_data.holdStartsMs)) - imaging_start_MW_T;
    hold_time  = double(cell2mat(b_data.holdTimesMs));   %duration of the lever hold on that trial
    react_time = double(cell2mat(b_data.reactTimesMs));
    req_hold   = double(cell2mat(b_data.tTotalReqHoldTimeMs));
    rnd_hold   = double(cell2mat(b_data.tRandReqHoldTimeMs));
    stimOnMS = double(cell2mat(b_data.stimOnUs)/1000) - imaging_start_MW_T;
    
    release_time = hold_start + hold_time;
    early_time = hold_time<req_hold & hold_time>=200;  %finds non-fidget early releases. Fidget < 350ms hold.   %LARGE BUG in early times. early inx values are all non zero and negative.
    early_early_time = hold_time<=1000 & hold_time>=200 & hold_time<req_hold;%hold_time<=1000 & hold_time>=200 & hold_time<req_hold;
    late_early_time = hold_time<req_hold & hold_time>=3500;
    early_correct = hold_time>req_hold & hold_time<= 1000;
    late_correct = hold_time>req_hold & hold_time >=3500;
    missed = cellfun(@(s) strcmp(b_data.trialOutcomeCell, s), {'ignore'}, 'UniformOutput', 0);
    missed = missed{1};
    early_correct = early_correct & (~missed);
    late_correct = late_correct & (~missed);
    
    has_reward = ~cellfun(@isempty, b_data.juiceTimesMsCell ); %using this as a proxy to isolate correct trials
    if isfield(b_data, 'itiRewardUnexpectPercent') && b_data.itiRewardUnexpectPercent > 0 
        has_reward = cellfun(@(s) strcmp(b_data.trialOutcomeCell, s), {'success'}, 'UniformOutput', 0);
        has_reward = has_reward{1};
    end
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
    
    early_tooFast = tooFastCorrects & early_correct;
    late_tooFast = tooFastCorrects & late_correct;
    %calculate trial outcome for all categories for all trials. Put into indeces. Values = time of lever release zeroed relative to f frame time MW
    early_inx = (early_time.*release_time);  %right now early time is just a bunch of 1s and 0s indexing the early trials.
    corr_inx = (has_reward.*release_time);
    tooFast_inx = (tooFastCorrects.*release_time);
    late_inx = ismember(release_time, late_time).*release_time;
    fidget_inx = ismember(release_time, fidget_time).*release_time;
    
    %store the indexed vectors
    trial_outcome.hasReward = has_reward;
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
    trial_outcome.tf_hold_dur = hold_time(logical(tooFastCorrects));
    trial_outcome.RT_cue = release_time - stimOnMS;
    trial_outcome.req_holdTime = req_hold(icc);
    trial_outcome.req_holdTime = trial_outcome.req_holdTime(reward_trial_reduce);
    
    %store the various categories of trial outcome in the for of the lever times for only those categories agnostic of trial num
    trial_outcome.success_time = release_time(has_reward);
    trial_outcome.tooFastCorrects = release_time(1,find(tooFastCorrects)); %isolates correct trials with reaction times >0 but <175ms. Stores them as a matrix of their release times.
    trial_outcome.early_time = release_time(early_time);                   %altered to exclude fidgets 2/27/16
    trial_outcome.early_early_time = release_time(early_early_time);
    trial_outcome.late_early_time = release_time(late_early_time);
    trial_outcome.early_correct_time = release_time(early_correct);
    trial_outcome.late_correct_time = release_time(late_correct);
    trial_outcome.early_tooFast_time = release_time(early_tooFast);
    trial_outcome.late_tooFast_time = release_time(late_tooFast);
    trial_outcome.fidget = fidget_time;
    trial_outcome.success_ptime = trial_outcome.success_time - hold_time(has_reward);
    trial_outcome.success_ptimeEndF = round(trial_outcome.succ_hold_dur/double(ifi));
    trial_outcome.early_ptime = trial_outcome.early_time - hold_time(early_time);
    trial_outcome.early_ptimeEndF = round(trial_outcome.fail_hold_dur/double(ifi));
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
    
    if ~isempty(licking_data)
        licking_data.lickRateS = licking_data.lickRate(has_reward);
        licking_data.lickRateF = licking_data.lickRate(early_time);
    end
    
    if isfield(b_data, 'rewardOmissionPercent') && b_data.rewardOmissionPercent ~= 0
        eIdx = cellfun(@isempty, b_data.tRewardOmissionTrial);  %fix empty cell
        temp = b_data.tRewardOmissionTrial;
        temp(eIdx) = {int64(0)};
        omitRewardIndx = cell2mat(temp);
        trial_outcome.omitReward = lever.release(logical(omitRewardIndx));
        trial_outcome.omitRewardIndx = omitRewardIndx;
    end
    
    if isfield(b_data,'itiRewardUnexpectPercent') && b_data.itiRewardUnexpectPercent ~= 0
        unexpRewardIndx = cell2mat(b_data.tItiUnexpectedRewardTrial);
        trial_outcome.ItiunexpReward = lever.itiReward(logical(unexpRewardIndx));
    end
    
else
    cTargetOn = cell2mat(b_data.cTargetOn); % cue onset frame
    cTargetOn = cTargetOn - counterValues(1) + 1;
    cLeverUp = cell2mat(b_data.cLeverUp);  % Reward frame
    cLeverUp = cLeverUp - counterValues(1) + 1;
   
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
    
    trial_outcome.trialLen = cell2mat(b_data.cTrialEnd) - cell2mat(b_data.cTrialStart);
    trial_outcome.omitRewardIndx = omitRewardIndx;
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