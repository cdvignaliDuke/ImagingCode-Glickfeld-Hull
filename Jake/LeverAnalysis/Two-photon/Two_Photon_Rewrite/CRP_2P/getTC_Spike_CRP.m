%code to detect calcium transients
% 1- Find periods with no lever activity and detect "spontaneous" events
% 2- Find periods around release events and detect
% 3- Find periods around press events and detect
% note: % buffer is a longer window and data_start, data_end is the real range, in order to plot 
% time course of all event types, they all need to have the same data_start
% and data_end to be well aligned to each other.

%load data and behavior info
load([dest, '_dFOverF_TC.mat'])
load([dest '_cue_movies.mat']) %get ifi
load([dest '_cue_movies_lick.mat'])
load([dest 'parse_behavior.mat'])
if isfield(frame_info, 'laser_on_ind_conserv')
    laser_on_ind_conserv = frame_info.laser_on_ind_conserv;
else
    load([dest, 'img_reg.mat'], 'laser_on_ind_conserv');
end
nIC = size(data_tc,1);
%data_tc = tc_avg;

%% 1- Find periods with no lever activity and detect "spontaneous" events
%extract spontaneous windows

%define some variables related to frame rate and number
Sampeling_rate = 1000/ifi;
pre_buffer = ceil(200/double(ifi));  %old value for pre_buffer=1000ms. Changed 181016
post_buffer = ceil(4000/double(ifi));

%find the first trial to use
empty_ind = find(~cellfun(@isempty, input.counterValues));
if input.counterValues{empty_ind(1)}(1) == 0
    start_i  =  empty_ind(2) + 1;
else
    start_i  =  empty_ind(1) + 1;
end

%find the last trial to use
trial_len = cell2mat(cellfun(@length, input.counterValues, 'UniformOutput', 0));
if strfind(input.savedDataName, 'i077-180403')
    trial_len(77) = []; %this trial was incompletely imaged. Was cut out of behavior variables in cleanBehav_CRP
end
valid_trial = find(trial_len > 2);
stop_i = min(valid_trial(end) - 1, length(lever.press));

%find cue onset timestamps +buffers   (lever.press = cue onset   lever.release is cue off)
fr_lever = [];   %fr_lever will be a series of frame numbers corresponding to lever events +&- the buffers
for ipress = start_i:stop_i    %this could be problematic due to unremoved NaNs
    fr_lever = [fr_lever frame_info.counter(round(lever.press(ipress)))-pre_buffer:frame_info.counter(round(lever.press(ipress))) + post_buffer ]; %used to be press-pre : press-post
end

%only take unique frames and frames which were imaged. 
fr_lever = unique(fr_lever);
fr_lever((fr_lever<1)) = [];
fr_lever((fr_lever>size(data_tc,1))) = [];
% fr_lever(find(fr_lever>length(frame_times))) = [];  %only looks at values that correspond to actual frames.
 
%takes a subset of frames which correspond to iti.
data_tc_spont = data_tc'; %data_tc_spont dim1=frame dim2=cell
laser_off_ind = 1:size(data_tc,2);
laser_off_ind(laser_on_ind_conserv) = [];
data_tc_spont(laser_off_ind,:) = NaN;
data_tc_spont(fr_lever,:) = NaN;
use_frames_ind = find(~isnan(sum(data_tc_spont,2)));

%identify frames on the borders between splices
use_fr_diff = diff(use_frames_ind);
border_fr = find(use_fr_diff>1);   %gives the indeces of use_frame_ind which are border frames. SHould also correspond to frame nums

%take the derivative of the raw tc
data_tc_spont = data_tc_spont([use_frames_ind],:);
events_diff_x = diff(data_tc_spont, [], 1);

%find all events to get rate
events = [];
events_ind = {};
events_rate = zeros(1,nIC);

%Search each derivative of cell's raw F tc for Ca events
opt.thresh = 2.1; % standard deviation of spike threshold
opt.normalization = 1;
opt.deconvtau = 0;
opt.dt = 1;
for ic = 1:nIC
    %     events_ind{ic} = find(events_diff_x(:,ic)>thresh(ic));
    [~, events_ind_temp, ~] = CellsortFindspikes(events_diff_x(:,ic), opt.thresh, opt.dt, opt.deconvtau, opt.normalization);
    %remove any false events on the borders of splices
    if ~isempty(intersect(events_ind_temp, border_fr));
        [~,cut_fr,~] = intersect(events_ind_temp, border_fr);
        events_ind_temp(cut_fr) = [];
    end
    events_ind{ic} = events_ind_temp;
    events_ind{ic}(find(diff(events_ind{ic})==1)+1)=[];  %if you have events in two consecutive frames, this removes the second event.
    events_rate(1,ic) = length(events_ind{ic})./((double(ifi)./1000)*size(data_tc_spont,1));
end

% %find all spontaneous events to get waveform 
min_iei = round(650./double(ifi)); %remove all events that occur within 650 ms
data_start = round(300./double(ifi));
data_end = round(600./double(ifi));    %looks at frames 300ms before and 600ms after the spont events
events_ind_good = {};
for ic = 1:nIC
    double_events = find(diff(events_ind{ic})<min_iei);
    events_ind_good{ic} = events_ind{ic};
    events_ind_good{ic}([double_events (double_events+1)])=[];
    ind_out = [];
    events(ic).f_chunk = [];
    for ii = 1: length(events_ind_good{ic})
        %event window may not include frames before 1st frame or after last frame. 
        %may not include border frames between splices
        if and(events_ind_good{ic}(ii)-data_start>0, events_ind_good{ic}(ii)+data_end<size(data_tc_spont,1)) &  isempty(intersect( [events_ind_good{ic}(ii)-data_start:events_ind_good{ic}(ii)+ round(350./double(ifi))], border_fr))         
            events(ic).f_chunk = [events(ic).f_chunk; data_tc_spont(events_ind_good{ic}(ii)-data_start:events_ind_good{ic}(ii)+data_end,ic)'];
        else
            ind_out = [ind_out ii];
        end
    end   %event.f_chunk will consist of frame ranges before and after a spont event that are sufficently far away from the begining or end of imaging. dim1=event#  dim2=frame#
    events_ind_good{ic}(ind_out)=[];
    events(ic).f_base = mean(events(ic).f_chunk(:,data_start-round(50./double(ifi)):data_start+1),2); %13:15),2); %
    events(ic).df_chunk = bsxfun(@minus, events(ic).f_chunk, events(ic).f_base);
    events(ic).dfoverf_chunk = bsxfun(@rdivide, events(ic).df_chunk, events(ic).f_base);
end

save([dest_sub '_spont_events.mat'], 'events', 'data_tc_spont', 'fr_lever', 'opt', 'events_rate', 'data_start', 'data_end');

%% 2- Find periods around cue/reward

%correct for vis stim artifact in three datasets
if strcmp(session, '180502_img084') | strcmp(session, '180403_img077') | strcmp(session, '180404_img077')
    if strcmp(session, '180502_img084')
        artifact_ind = [2:5];
    elseif strcmp(session, '180403_img077')
        artifact_ind = [2:4];
    elseif strcmp(session, '180404_img077')
        artifact_ind = [3:6];
    end
    for trial_ind = 1:length(lever.press)
        for cell_num = 1:nIC
            data_tc(cell_num,    [frame_info.counter([lever.press(trial_ind)])+artifact_ind]            ) = ...
                linspace(  data_tc(cell_num,    frame_info.counter(lever.press(trial_ind))+artifact_ind(1)-1   ), ...
                data_tc(cell_num,     frame_info.counter(lever.press(trial_ind))+artifact_ind(end)+1     ), ...
                length(artifact_ind));    
        end
    end
end

%Define some input variables
pre_buffer_rew = ceil(2000/double(ifi));
post_buffer_rew = ceil(2000/double(ifi));
NRdata_start_rew = ceil(300/double(ifi));
NRdata_end_rew = ceil(600/double(ifi));
% =================need to rewrite so it doesnt look at small snippets?
%Find Normal Reward events
[NR_tc, normalR] = extractEvent(data_tc, frame_info, trial_outcome.normalReward, NR_lick_info, opt, pre_buffer_rew, post_buffer_rew, NRdata_start_rew, NRdata_end_rew, min_iei);

%Find Reward omission events
if ~isempty(trial_outcome.omitReward)
    [OR_tc, omitR] = extractEvent(data_tc, frame_info, trial_outcome.omitReward, OR_lick_info, opt, pre_buffer_rew, post_buffer_rew, NRdata_start_rew, NRdata_end_rew, min_iei);
else
    OR_tc = nan(size(NR_tc));
    omitR = [];
end

%Find Unexpected Reward events
if ~isempty(trial_outcome.unexpReward)
    [UR_tc, unexpR] = extractEvent(data_tc, frame_info, trial_outcome.unexpReward, UR_lick_info, opt, pre_buffer_rew, post_buffer_rew, NRdata_start_rew, NRdata_end_rew, min_iei);
else
    UR_tc = nan(size(NR_tc));
    unexpR = [];
end

% find Cue for NR and OR
pre_buffer_cue = ceil(1000/double(ifi));
post_buffer_cue = ceil(2000/double(ifi));
Cuedata_start = ceil(300/double(ifi));  %widnow around individual Ca events which will be extracted. 
Cuedata_end = ceil(600/double(ifi));

% time courses aligned to the cue
[NRCue_tc, normalCue] = extractEvent(data_tc, frame_info, trial_outcome.normalRewardCue, NR_lick_info, opt, pre_buffer_cue, post_buffer_cue, Cuedata_start, Cuedata_end, min_iei);
if ~isempty(trial_outcome.omitRewardCue)
    [ORCue_tc, omitCue] = extractEvent(data_tc, frame_info, trial_outcome.omitRewardCue, OR_lick_info, opt, pre_buffer_cue, post_buffer_cue, Cuedata_start, Cuedata_end, min_iei);
else
    ORCue_tc = nan(size(NRCue_tc));
    omitCue = [];
end
if ~isempty(trial_outcome.unexpRewardCue)
    [URCue_tc, unexpCue] = extractEvent(data_tc, frame_info, trial_outcome.unexpRewardCue, UR_lick_info, opt, pre_buffer_cue, post_buffer_cue, Cuedata_start, Cuedata_end, min_iei);
else
    URCue_tc = nan(size(NRCue_tc));
    unexpCue = [];
end

save([dest_sub '_evoked_events.mat'], 'NR_tc', 'OR_tc', 'UR_tc', 'NRCue_tc', 'ORCue_tc', 'URCue_tc', 'normalR', 'omitR', 'unexpR', 'normalCue', 'omitCue', 'unexpCue', 'min_iei', 'pre_buffer_cue', 'post_buffer_cue',...
    'NRdata_start_rew', 'NRdata_end_rew', 'Cuedata_start', 'Cuedata_end', 'post_buffer_rew', 'pre_buffer_rew', 'opt')

