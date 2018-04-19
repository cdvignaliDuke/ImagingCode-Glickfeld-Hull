%code to detect calcium transients
% 1- Find periods with no lever activity and detect "spontaneous" events
% 2- Find periods around release events and detect
% 3- Find periods around press events and detect
% note: % buffer is a longer window and data_start, data_end is the real range, in order to plot 
% time course of all event types, they all need to have the same data_start
% and data_end to be well aligned to each other.

%load data and behavior info
load([dest, '_dFOverF_TC.mat'])
data_tc = tc_avg;
nIC = size(data_tc,2);
load([dest '_cue_movies.mat']) %get ifi
load([dest '_cue_movies_lick.mat'])
load([dest 'parse_behavior.mat'])
%% 1- Find periods with no lever activity and detect "spontaneous" events
%extract spontaneous windows

Sampeling_rate = 1000/ifi;
pre_buffer = ceil(1000/double(ifi));  
post_buffer = ceil(2000/double(ifi));

empty_ind = find(~cellfun(@isempty, input.counterValues));
if input.counterValues{empty_ind(1)}(1) == 0
    
    start_i  =  empty_ind(2) + 1;
else
    
    start_i  =  empty_ind(1) + 1;
end

trial_len = cell2mat(cellfun(@length, input.counterValues, 'UniformOutput', 0));
valid_trial = find(trial_len > 2);
stop_i = min(valid_trial(end) - 1, length(lever.press));

fr_lever = [];   %fr_lever will be a series of frame numbers corresponding to lever events +&- the buffers
for ipress = start_i:stop_i    %this could be problematic due to unremoved NaNs

    fr_lever = [fr_lever frame_info.counter(round(lever.press(ipress)))-pre_buffer:frame_info.counter(round(lever.press(ipress))) - post_buffer ];
end
% for irelease = start_i:stop_i    %this could be problematic due to unremoved NaNs
%     fr_lever = [fr_lever frame_info.counter(round(lever.release(irelease)))-pre_buffer:frame_info.counter(round(lever.release(irelease)))+post_buffer];
% end
% 
fr_lever = unique(fr_lever);
fr_lever((fr_lever<1)) = [];
fr_lever((fr_lever>size(data_tc,1))) = [];
% % fr_lever(find(fr_lever>length(frame_times))) = [];  %only looks at values that correspond to actual frames.
% 
% 
data_tc_spont = data_tc;
data_tc_spont(fr_lever,:) = [];

%find all events to get rate
events = [];
events_diff_x = zeros(size(data_tc_spont,1)-1,nIC);
% events_diff_x2 = zeros(size(data_tc_spont,1)-2,nIC);
events_ind = {};
events_rate = zeros(1,nIC);

for ic = 1:nIC
    events_diff_x(:,ic) = diff(data_tc_spont(:,ic),[],1);
%     events_diff_x1(:,ic) = diff(data_tc_spont(:,ic),[],1);
%     events_diff_x2(:,ic) = diff(events_diff_x1(:,ic),[],1);
end


thresh = 2.1; % standard deviation of spike threshold

normalization = 1;
deconvtau = 0;
dt = 1;

for ic = 1:nIC
    %     events_ind{ic} = find(events_diff_x(:,ic)>thresh(ic));
    [~, events_ind{ic}, ~] = CellsortFindspikes(events_diff_x(:,ic), thresh, dt, deconvtau, normalization);
    events_ind{ic}(find(diff(events_ind{ic})==1)+1)=[];
    events_rate(1,ic) = length(events_ind{ic})./((double(ifi)./1000)*size(data_tc_spont,1));
end

% %find all spontaneous events to get waveform
min_iei = round(650./double(ifi)); %remove all events that occur within 650 ms
data_start = round(300./double(ifi));
data_end = round(600./double(ifi));    %looks at frames 300ms before and 700ms after the spont events
events_ind_good = {};
for ic = 1:nIC
    double_events = find(diff(events_ind{ic})<min_iei);
    events_ind_good{ic} = events_ind{ic};
    events_ind_good{ic}([double_events (double_events+1)])=[];
    ind_out = [];
    events(ic).f_chunk = [];
    for ii = 1: length(events_ind_good{ic})
        if and(events_ind_good{ic}(ii)-data_start>0, events_ind_good{ic}(ii)+data_end<size(data_tc_spont,1))
            events(ic).f_chunk = [events(ic).f_chunk; data_tc_spont(events_ind_good{ic}(ii)-data_start:events_ind_good{ic}(ii)+data_end,ic)'];
        else
            ind_out = [ind_out ii];
        end
    end   %event.f_chunk will consist of frame ranges before and after a spont event that are sufficently far away from the begining or end of imaging
    events_ind_good{ic}(ind_out)=[];
    events(ic).f_base = mean(events(ic).f_chunk(:,data_start-round(50./double(ifi)):data_start+1),2); %13:15),2); %
    events(ic).df_chunk = bsxfun(@minus, events(ic).f_chunk, events(ic).f_base);
    events(ic).dfoverf_chunk = bsxfun(@rdivide, events(ic).df_chunk, events(ic).f_base);
end

save([dest_sub '_spont_events.mat'], 'events', 'data_tc_spont', 'fr_lever', 'thresh', 'events_rate', 'data_start', 'data_end');

%% 2- Find periods around release events and detect
%find Normal Reward and Reward omission
opt.thresh = 2.1;

opt.normalization = 1;
opt.deconvtau = 0;
opt.dt = 1;
pre_buffer = ceil(2000/double(ifi));
post_buffer = ceil(2000/double(ifi));
NRdata_start = ceil(300/double(ifi));
NRdata_end = ceil(600/double(ifi));

[NR_tc, normalR] = extractEvent(data_tc, frame_info, trial_outcome.normalReward, NR_lick_info, opt, pre_buffer, post_buffer, NRdata_start, NRdata_end, min_iei);

if ~isempty(trial_outcome.omitReward)
    [OR_tc, omitR] = extractEvent(data_tc, frame_info, trial_outcome.omitReward, OR_lick_info, opt, pre_buffer, post_buffer, NRdata_start, NRdata_end, min_iei);
else
    OR_tc = nan(size(NR_tc));
    omitR = [];
end

if ~isempty(trial_outcome.unexpReward)
    [UR_tc, unexpR] = extractEvent(data_tc, frame_info, trial_outcome.unexpReward, UR_lick_info, opt, pre_buffer, post_buffer, NRdata_start, NRdata_end, min_iei);
else
    UR_tc = nan(size(NR_tc));
    unexpR = [];
end

% find Cue for NR and OR
pre_buffer = ceil(1000/double(ifi));
post_buffer = ceil(2000/double(ifi));
Cuedata_start = ceil(300/double(ifi));
Cuedata_end = ceil(600/double(ifi));

% time course extracted here not used
[NRCue_tc, normalCue] = extractEvent(data_tc, frame_info, trial_outcome.normalRewardCue, NR_lick_info, opt, pre_buffer, post_buffer, Cuedata_start, Cuedata_end, min_iei);

if ~isempty(trial_outcome.omitRewardCue)
    [ORCue_tc, omitCue] = extractEvent(data_tc, frame_info, trial_outcome.omitRewardCue, OR_lick_info, opt, pre_buffer, post_buffer, Cuedata_start, Cuedata_end, min_iei);
else
    ORCue_tc = nan(size(NRCue_tc));
    omitCue = [];
end

if ~isempty(trial_outcome.unexpRewardCue)
    [URCue_tc, unexpCue] = extractEvent(data_tc, frame_info, trial_outcome.unexpRewardCue, UR_lick_info, opt, pre_buffer, post_buffer, Cuedata_start, Cuedata_end, min_iei);
else
    URCue_tc = nan(size(NRCue_tc));
    unexpCue = [];
end

save([dest_sub '_evoked_events.mat'], 'NR_tc', 'OR_tc', 'UR_tc', 'NRCue_tc', 'ORCue_tc', 'URCue_tc', 'normalR', 'omitR', 'unexpR', 'normalCue', 'omitCue', 'unexpCue', 'min_iei', 'thresh', 'pre_buffer', 'post_buffer',...
    'NRdata_start', 'NRdata_end', 'Cuedata_start', 'Cuedata_end')

