function [trial_outcome, lickTimes] = parse_behavior_for_HAD(b_data, frame_info);
%Processes the behavior data for trial_outcome and lickTimes so it can be used for TCs later. More documentation at the end of the function

%use time of first frame to align licking times to start of imaging
imaging_start_MW_T = frame_info.imaging_start_MW_T;
lickTimes=[];
if isfield(b_data, 'lickometerTimesUs');  %img24 and 25 use datasets which have no licking fields 
    for kk = 1:length(b_data.lickometerTimesUs);
        lickTimes = [lickTimes cell2mat(b_data.lickometerTimesUs(kk))/1000];
    end
    lickTimes = double(lickTimes)-imaging_start_MW_T;
end

%Collects and calculates various categories of trial outcome. Does not trim the unimaged trials
hold_start = double(cell2mat(b_data.holdStartsMs)) - imaging_start_MW_T;
hold_time  = double(cell2mat(b_data.holdTimesMs));   %duration of the lever hold on that trial
react_time = double(cell2mat(b_data.reactTimesMs));
req_hold   = double(cell2mat(b_data.tTotalReqHoldTimeMs));
release_time = hold_start + hold_time;
early_time = hold_time<req_hold & hold_time>=200;  %finds non-fidget early releases. Fidget < 350ms hold.   %LARGE BUG in early times. early inx values are all non zero and negative. 
has_reward = ~cellfun(@isempty, b_data.juiceTimesMsCell ); %using this as a proxy to isolate correct trials 
fidget_time = release_time(hold_time<200);
late_time = release_time(find(strcmp('ignore',b_data.trialOutcomeCell)));
unexp_rew = double(cell2mat(b_data.tDoNoStimulusChange));

%find all the correct trials which were too fast to be responding to visual cue
tooFastCorrects = zeros(size(b_data.reactTimesMs)); 
if isfield(b_data, 'doFakeMouseSuccessOnly');
    1
    fake_mouse = b_data.doFakeMouseSuccessOnly;
else
    fake_mouse = 0;
end
for i = 1:length(b_data.reactTimesMs)
    if react_time(i) < 176 & fake_mouse ==0;   %Removes tooFastCorrects from the successful trials count. But not for CRP trials. 
        has_reward(i)=0;
    end
    if isequal(b_data.trialOutcomeCell{i}, 'success'); 
        tooFastCorrects(i) = react_time(i)>0 & react_time(i)<176; %definition of and isolation of tooFastCorrects
    end
end

%calculate trial outcome for all categories for all trials. Put into indeces. Values = time of lever release zeroed relative to f frame time MW
early_inx = (early_time.*release_time);  %right now early time is just a bunch of 1s and 0s indexing the early trials. 
corr_inx = (has_reward.*release_time); 
corr_inx(find(unexp_rew)) = 0;
tooFast_inx = (tooFastCorrects.*release_time);
late_inx = ismember(release_time, late_time).*release_time;
fidget_inx = ismember(release_time, fidget_time).*release_time;
reward_omission_inx = release_time;
reward_omission_inx(find(has_reward)) = 0;
unexp_inx = (unexp_rew.*release_time);


%store the indexed vectors
trial_outcome = [];
if fake_mouse == 0; %condition set for regular lever trials and no-lever controls
    trial_outcome.early_inx = early_inx;
    trial_outcome.corr_inx = corr_inx;
    trial_outcome.tooFast_inx = tooFast_inx;
    trial_outcome.late_inx = late_inx;
    trial_outcome.fidget_inx = fidget_inx;
elseif fake_mouse == 1; %condition set for cue-reward pairing trials
    trial_outcome.rewarded_trials_inx = corr_inx;
    trial_outcome.rew_omission_inx = reward_omission_inx;
    trial_outcome.unexp_inx = unexp_inx;
end

%remove unimaged frames from hold duration variables and store them in trial_outcome
has_reward([1:frame_info.f_frame_trial_num, frame_info.l_frame_trial_num:end])=0;
early_time([1:frame_info.f_frame_trial_num, frame_info.l_frame_trial_num:end])=0;
trial_outcome.succ_hold_dur = hold_time(has_reward);
trial_outcome.fail_hold_dur = hold_time(early_time);

%store the various categories of trial outcome in the for of the lever times for only those categories agnostic of trial num
if b_data.doLever== 1; %normal lever sessions
    trial_outcome.success_time = release_time(has_reward);
    trial_outcome.tooFastCorrects = release_time(1,find(tooFastCorrects)); %isolates correct trials with reaction times >0 but <175ms. Stores them as a matrix of their release times.
    trial_outcome.early_time = release_time(early_time);                   %altered to exclude fidgets 2/27/16
    trial_outcome.fidget = fidget_time;
    trial_outcome.late_time = late_time;               %if respond before change then no change in orientation
elseif fake_mouse ==1; %if this is a cue reward pairing session
    trial_outcome.rewarded_trial_time = release_time(find(corr_inx));
    trial_outcome.rew_omission_trial_time = release_time(~has_reward);
    trial_outcome.unexpected_rew_time = release_time(find(unexp_rew));
elseif b_data.doLever == 0;                             %no-lever control session
    trial_outcome.early_time = release_time(hold_time<req_hold);
    trial_outcome.fidget = [];
end
trial_outcome.change_orientation = hold_start + req_hold;  %zeroed relative to first frame time in MWorks time
trial_outcome.change_orientation(react_time < 0) = NaN;

%explanation of certain variabels
%trial_outcome.success_time   lever release for successes 
%trial_outcome.early_time     time of lever release in early trials 
%trial_outcome.late_time      time of lever solenoid coming up.
return; 
