%code to detect calcium transients
% 1- Find periods with no lever activity and detect "spontaneous" events
% 2- Find periods around release events and detect
% 3- Find periods around press events and detect

%load data and behavior info
% load([dest, '_ROI_TCs.mat'])
load([dest, '_dFOverF_TC.mat'])
data_tc = tc_avg;
nIC = size(data_tc,2);
% load([dest '_frame_times.mat'])
load([dest 'parse_behavior.mat'])
%% 1- Find periods with no lever activity and detect "spontaneous" events
%extract spontaneous windows
if exist('frame_times') == 1
    ifi = (frame_times(end)-frame_times(1))/length(frame_times);
else
    ifi = mode(diff(cell2mat(cellfun(@int64,input.counterTimesUs,'UniformOutput',0))))/1000;
end
Sampeling_rate = 1000/ifi;
pre_buffer = ceil(300/double(ifi));
post_buffer = ceil(1000/double(ifi));

if isfield(input, 'tItiWaitFrames')
    
    itiFrames = cell2mat(input.tItiWaitFrames);
else
    itiFrames = cell2mat(input.tItiWaitTimeMs)/double(ifi);
end
    
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
fr_iti = [];
for ipress = start_i:stop_i    %this could be problematic due to unremoved NaNs
    
    fr_lever = [fr_lever frame_info.counter(round(lever.press(ipress)))-pre_buffer:frame_info.counter(round(lever.press(ipress)))+post_buffer];
    fr_iti = [fr_iti; frame_info.counter(round(lever.press(ipress)))-itiFrames(ipress) frame_info.counter(round(lever.press(ipress)))];
end
for irelease = start_i:stop_i    %this could be problematic due to unremoved NaNs
    fr_lever = [fr_lever frame_info.counter(round(lever.release(irelease)))-pre_buffer:frame_info.counter(round(lever.release(irelease)))+post_buffer];
end

fr_lever = unique(fr_lever);
fr_lever(find(fr_lever<1)) = [];
fr_lever(find(fr_lever>size(data_tc,1))) = [];
% fr_lever(find(fr_lever>length(frame_times))) = [];  %only looks at values that correspond to actual frames.
fr_iti(find(fr_iti<1),:) = [];
fr_iti(find(fr_iti(:,2)>size(data_tc,1)),:) = [];
% fr_lick = [];
if exist('lick_data', 'var') && ~isempty(lick_data)
    for ilick = 1:length(lick_data.single_lickTimeAll)    %this could be problematic due to unremoved NaNs
        %         fr_lick = [fr_lick frame_info.counter(round(lick_data.bout_onsetTime(ilick)))-pre_buffer:frame_info.counter(round(lick_data.bout_onsetTime(ilick)))+post_buffer];
        fr_lever = [fr_lever frame_info.counter(round(lick_data.single_lickTimeAll(ilick)))-pre_buffer:frame_info.counter(round(lick_data.single_lickTimeAll(ilick)))+post_buffer];
    end
end

data_tc_spont = data_tc;
%  data_tc_spont = tc_dfoverf';
data_tc_spont(fr_lever,:) = [];  %removes frames surrounding lever events



% % [row, col] = find(data_tc_spont==0);
% % data_tc_spont(row,:) = [];

% meanF = mean(data_tc_spont, 2);
%
% df_f = bsxfun(@minus, data_tc_spont, meanF);
%
% data_tc_spont = bsxfun(@rdivide, df_f, meanF);


%find all events to get rate
events = [];
events_diff_x1 = zeros(size(data_tc_spont,1)-1,nIC);
events_diff_x2 = zeros(size(data_tc_spont,1)-2,nIC);
events_ind = {};
events_rate = zeros(1,nIC);

for ic = 1:nIC
    %     events_diff_x(:,ic) = diff(data_tc_spont(:,ic),[],1);
    events_diff_x1(:,ic) = diff(data_tc_spont(:,ic),[],1);
    events_diff_x2(:,ic) = diff(events_diff_x1(:,ic),[],1);
end

% disp('manually select thresholds for each TC')
% for ii = size(events_diff_x,2):-1:1;
%     figure; plot(events_diff_x(:,ii));
%     title(['diffTC for TC #' num2str(ii)]);
% end
%pause

%%%%% 2P
if sub < 5
    thresh = 1.5;
elseif sub == 5 || sub == 6 || sub == 7 || sub ==12 || sub == 20
    thresh = 3; % two standard deviation
elseif sub >= 21 && sub <=26
    thresh = 2.1;
elseif sub == 29 || sub == 30
    thresh = 2.1;
else
    thresh = 2.5;
end

%%%% WF
% thresh = 1.5;
opt.thresh = thresh;

opt.normalization = 1;
opt.deconvtau = 0;
opt.dt = 1;

normalization = 1;
deconvtau = 0;
dt = 1;
iti_acor_all = [];
for ii = 1:size(fr_iti,1)
    iti_tc = data_tc(fr_iti(ii,1):fr_iti(ii,2),:);
    for ic = 1:nIC
        diff_iti(:,ic) = diff(iti_tc(:,ic), [], 1);
    end
    iti_events = 0*diff_iti;
    for ic = 1:nIC
        [~, iti_ind, ~] = CellsortFindspikes(diff_iti(:,ic), 2, dt, deconvtau, normalization);
        iti_ind(find(diff(iti_ind)==1)+1)=[];
        iti_events(iti_ind,ic) = 1;
        iti_acor(ic,:) = autocorr(iti_events(:,ic), round(2000/double(ifi)))';
    end
    iti_acor_all = cat(3, iti_acor_all, iti_acor);
end
iti_cor_mean = squeeze(nanmean(iti_acor_all,1));


events_inds_x1 = 0*events_diff_x1;
for ic = 1:nIC
    %     events_ind{ic} = find(events_diff_x(:,ic)>thresh(ic));
    [~, events_ind{ic}, ~] = CellsortFindspikes(events_diff_x1(:,ic), thresh, dt, deconvtau, normalization);
    events_ind{ic}(find(diff(events_ind{ic})==1)+1)=[];
    events_inds_x1(events_ind{ic},ic) = 1;
    events_rate(1,ic) = length(events_ind{ic})./((double(ifi)./1000)*size(data_tc_spont,1));
end

%find all spontaneous events to get waveform
min_iei = round(650./double(ifi)); %remove all events that occur within 650 ms
data_start = round(300./double(ifi));
data_end = round(700./double(ifi));    %looks at frames 300ms before and 700ms after the spont events
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

save([dest_sub '_spont_events.mat'], 'events', 'data_tc_spont', 'events_inds_x1', 'iti_cor_mean', 'fr_lever', 'thresh', 'events_rate', 'data_start', 'data_end');

% %% 2- Find periods around release events and detect
% %find successes, earlies
% pre_buffer = ceil(1000/double(ifi));
% post_buffer = ceil(1000/double(ifi));
% success_fr = zeros(length(trial_outcome.success_time),pre_buffer+post_buffer+1);
% fail_fr = zeros(length(trial_outcome.early_time),pre_buffer+post_buffer+1);
% success_tc = [];
% fail_tc = [];
% neg_ind = find(trial_outcome.success_time<0);
% if isempty(neg_ind)
%     start_is = 1;
% else
%     start_is = neg_ind(end) + 1;
% end
% for is = start_is:length(trial_outcome.success_time)    %this could be problematic due to unremoved NaNs
%     
%     if trial_outcome.success_time(is) < length(frame_info.counter) && frame_info.counter(trial_outcome.success_time(is))-pre_buffer > 0
%         success_fr(is,:) = [frame_info.counter(trial_outcome.success_time(is))-pre_buffer:frame_info.counter(trial_outcome.success_time(is))+post_buffer];
%         if and(min(success_fr(is,:),[],2)>0, max(success_fr(is,:),[],2)<size(data_tc,1))
%             success_tc = cat(3, success_tc, data_tc(success_fr(is,:),:));
%         end
%     elseif trial_outcome.success_time(is) >= length(frame_info.counter)
%         success_fr(is:end,:) = [];
%         break
%     end
% end
% if ~isempty(neg_ind)
%     success_fr(1:start_is-1 ,:)=[];
% end
% [xx, yy] = find(success_fr == 0);
% success_fr(unique(xx),:) = [];
% 
% neg_ind = find(trial_outcome.early_time<0);
% if isempty(neg_ind)
%     start_ie = 1;
% else
%     start_ie = neg_ind(end) + 1;
% end
% for ie = start_ie:length(trial_outcome.early_time)    %this could be problematic due to unremoved NaNs
%     if trial_outcome.early_time(ie) < length(frame_info.counter) && frame_info.counter(trial_outcome.early_time(ie))-pre_buffer > 0
%         fail_fr(ie,:) = [frame_info.counter(trial_outcome.early_time(ie))-pre_buffer:frame_info.counter(trial_outcome.early_time(ie))+post_buffer];
%         if and(min(fail_fr(ie,:),[],2)>0, max(fail_fr(ie,:),[],2)<size(data_tc,1))
%             fail_tc = cat(3, fail_tc, data_tc(fail_fr(ie,:),:));
%         end
%     elseif trial_outcome.early_time(ie) >= length(frame_info.counter)
%         fail_fr(ie:end,:) = [];
%     end
% end
% if ~isempty(neg_ind)
%     fail_fr(1:start_ie-1,:)=[];
% end
% 
% [xx, yy] = find(fail_fr == 0);
% fail_fr(unique(xx),:) = [];
% 
% %find all events during success trials
% success = [];
% win_release_start = round(33./double(ifi));
% win_release_end = round(200./double(ifi));
% for ic = 1:nIC
%     success_diff = diff(squeeze(success_tc(:,ic,:)),[],1);
%     success(ic).f_chunk = [];
%     success(ic).ind = {};
%     success(ic).good_ind = {};
%     for is = 1:size(success_diff,2)
%         [~,ind,~] = CellsortFindspikes(success_diff(:,is), thresh, dt, deconvtau, normalization);
%         %         ind = find(success_diff(:,is)>thresh(ic));
%         ind(find(diff(ind)==1)+1)=[];
%         double_events = find(diff(ind)<min_iei);
%         ind_good = ind;
%         ind_good([double_events (double_events+1)])=[];
%         success(ic).ind{is} = ind;
%         success(ic).good_ind{is} = ind_good;
%         if intersect(find(ind>pre_buffer-win_release_start), find(ind<pre_buffer+win_release_end))
%             wind = success_fr(is,pre_buffer-win_release_start:pre_buffer+win_release_end);
%             if intersect(find(ind_good>pre_buffer-win_release_start), find(ind_good<pre_buffer+win_release_end))
%                 success(ic).event_ind(:,is) =ind_good(intersect(find(ind_good>pre_buffer-win_release_start), find(ind_good<pre_buffer+win_release_end)));
%                 first_event = success_fr(is,success(ic).event_ind(:,is));
%                 success(ic).good_event(:,is) = 1;
%                 success(ic).event(:,is) = 1;
%             else
%                 success(ic).event_ind(:,is) = min(ind(intersect(find(ind>pre_buffer-win_release_start), find(ind<pre_buffer+win_release_end))),[],1);
%                 first_event = success_fr(is,success(ic).event_ind(:,is));
%                 success(ic).good_event(:,is) = 0;
%                 success(ic).event(:,is) = 1;
%             end
%             success(ic).f_chunk = [success(ic).f_chunk; data_tc(first_event-data_start:first_event+data_end,ic)']; %data_tc(success_fr(is,:),ic)']; %
%         else
%             success(ic).event_ind(:,is) = NaN;
%             wind = success_fr(is,pre_buffer-win_release_start:pre_buffer+win_release_end);
%             first_event = wind(1);
%             success(ic).good_event(:,is) = NaN;
%             success(ic).event(:,is) = 0;
%             success(ic).f_chunk = [success(ic).f_chunk; data_tc(first_event-data_start:first_event+data_end,ic)']; %data_tc(success_fr(is,:),ic)']; %
%         end
%     end
%     success(ic).f_base = mean(success(ic).f_chunk(:,data_start-round(50./double(ifi)):data_start+1),2); %13:15),2); %
%     success(ic).df_chunk = bsxfun(@minus, success(ic).f_chunk, success(ic).f_base);
%     success(ic).dfoverf_chunk = bsxfun(@rdivide, success(ic).df_chunk, success(ic).f_base);
% end
% 
% %find all events during early trials
% fail = [];
% for ic = 1:nIC
%     fail_diff = diff(squeeze(fail_tc(:,ic,:)),[],1);
%     fail(ic).f_chunk = [];
%     fail(ic).ind = {};
%     fail(ic).good_ind = {};
%     for ie = 1:size(fail_diff,2)
%         [~,ind,~] = CellsortFindspikes(fail_diff(:,ie), thresh, dt, deconvtau, normalization);
%         %         ind = find(fail_diff(:,ie)>thresh(ic));
%         ind(find(diff(ind)==1)+1)=[];
%         double_events = find(diff(ind)<min_iei);
%         ind_good = ind;
%         ind_good([double_events (double_events+1)])=[];
%         fail(ic).ind{ie} = ind;
%         fail(ic).good_ind{ie} = ind_good;
%         if intersect(find(ind>pre_buffer-win_release_start), find(ind<pre_buffer+win_release_end))
%             wind = fail_fr(ie,pre_buffer-win_release_start:pre_buffer+win_release_end);
%             if intersect(find(ind_good>pre_buffer-win_release_start), find(ind_good<pre_buffer+win_release_end))
%                 fail(ic).event_ind(:,ie) =ind_good(intersect(find(ind_good>pre_buffer-win_release_start), find(ind_good<pre_buffer+win_release_end)));
%                 first_event = fail_fr(ie,fail(ic).event_ind(:,ie));
%                 fail(ic).good_event(:,ie) = 1;
%                 fail(ic).event(:,ie) = 1;
%             else
%                 fail(ic).event_ind(:,ie) = min(ind(intersect(find(ind>pre_buffer-win_release_start), find(ind<pre_buffer+win_release_end))),[],1);
%                 first_event = fail_fr(ie,fail(ic).event_ind(:,ie));
%                 fail(ic).good_event(:,ie) = 0;
%                 fail(ic).event(:,ie) = 1;
%             end
%             fail(ic).f_chunk = [fail(ic).f_chunk; data_tc(first_event-data_start:first_event+data_end,ic)']; %data_tc(fail_fr(ie,:),ic)']; %
%         else
%             fail(ic).event_ind(:,ie) = NaN;
%             wind = fail_fr(ie,pre_buffer-win_release_start:pre_buffer+win_release_end);
%             first_event = wind(1);
%             fail(ic).good_event(:,ie) = NaN;
%             fail(ic).event(:,ie) = 0;
%             fail(ic).f_chunk = [fail(ic).f_chunk; data_tc(first_event-data_start:first_event+data_end,ic)']; %data_tc(fail_fr(ie,:),ic)']; %
%         end
%     end
%     fail(ic).f_base = mean(fail(ic).f_chunk(:,data_start-round(50./double(ifi)):data_start+1),2); %13:15),2); %
%     fail(ic).df_chunk = bsxfun(@minus, fail(ic).f_chunk, fail(ic).f_base);
%     fail(ic).dfoverf_chunk = bsxfun(@rdivide, fail(ic).df_chunk, fail(ic).f_base);
% end
% 
% 
% %% 3- Find periods around press events and detect
% %find presses
% pre_buffer = ceil(1000/double(ifi));
% post_buffer = ceil(1000/double(ifi));
% press_fr = [];  %press fr will equal a series of frame numbers around each lever press
% press_tc = [];
% ind_long = intersect(find(cell2mat(input.holdTimesMs)>=500),find(~isnan(trial_outcome.ind_press_prerelease)));
% for ip = 1:length(ind_long)
%     iT = ind_long(ip);
%     if isempty(cell2mat(cellfun(@double,input.leverTimesUs, 'UniformOutput', 0)) );
%         leverTimes = round((cell2mat(input.holdStartsMs(iT))*1000-input.counterTimesUs{start_i - 1}(1))./1000);
%     else
%         
%         leverTimes = round((cell2mat(input.leverTimesUs(iT))-input.counterTimesUs{start_i - 1}(1))./1000);
%     end
%     if leverTimes(trial_outcome.ind_press_prerelease(iT)) > 0
%         if leverTimes(trial_outcome.ind_press_prerelease(iT)) < length(frame_info.counter) && frame_info.counter(leverTimes(trial_outcome.ind_press_prerelease(iT)))-pre_buffer > 0
%             fr_range = [frame_info.counter(leverTimes(trial_outcome.ind_press_prerelease(iT)))-pre_buffer:frame_info.counter(leverTimes(trial_outcome.ind_press_prerelease(iT)))+post_buffer];
%             press_fr = [press_fr; fr_range];
%             if and(min(fr_range,[],2)>0, max(fr_range,[],2)<size(data_tc,1))
%                 press_tc = cat(3, press_tc, data_tc(fr_range,:));
%             end
%         end
%     end
% end
% 
% %find all events around the press
% press = [];
% win_press_start = round(300./double(ifi));  %why this window?
% win_press_end = round(100./double(ifi));
% for ic = 1:nIC
%     press_diff = diff(squeeze(press_tc(:,ic,:)),[],1);   %calcium traces for time around the lever press. Lever press happens in the center.
%     press(ic).f_chunk = [];
%     press(ic).ind = {};
%     press(ic).good_ind = {};
%     for ip = 1:size(press_diff,2)
%         [~,ind,~] = CellsortFindspikes(press_diff(:,ip), thresh, dt, deconvtau, normalization);
%         %         ind = find(press_diff(:,ip)>thresh(ic));  %finds all calcium events around a lever press that are > threshold for that neuron
%         ind(find(diff(ind)==1)+1)=[];     %events must be separated by at least one frame.
%         double_events = find(diff(ind)<min_iei);
%         ind_good = ind;
%         ind_good([double_events (double_events+1)])=[];
%         press(ic).ind{ip} = ind;
%         press(ic).good_ind{ip} = ind_good;
%         if intersect(find(ind>pre_buffer-win_press_start), find(ind<pre_buffer+win_press_end))  %if there are events within the window
%             wind = press_fr(ip,pre_buffer-win_press_start:pre_buffer+win_press_end);
%             if intersect(find(ind_good>pre_buffer-win_press_start), find(ind_good<pre_buffer+win_press_end))  %if ind_good has an event within the window
%                 ii = ind_good(intersect(find(ind_good>pre_buffer-win_press_start), find(ind_good<pre_buffer+win_press_end))); % ii = events within the window
%                 if size(ii,1)>1
%                     [minn, ind] = min(abs(ii-pre_buffer),[],1);
%                     ii = ii(ind);
%                 end
%                 press(ic).event_ind(:,ip) = ii;
%                 first_event = press_fr(ip,press(ic).event_ind(:,ip));
%                 press(ic).good_event(:,ip) = 1;
%                 press(ic).event(:,ip) = 1;
%             else
%                 ii = ind(intersect(find(ind>pre_buffer-win_press_start), find(ind<pre_buffer+win_press_end)));
%                 if size(ii,1)>1
%                     [nmin, imin] = min(abs(ii-pre_buffer),[],1);
%                     ii = ii(imin);
%                 end
%                 press(ic).event_ind(:,ip) = ii;
%                 first_event = press_fr(ip,press(ic).event_ind(:,ip));
%                 press(ic).good_event(:,ip) = 0;
%                 press(ic).event(:,ip) = 1;
%             end
%             if and(ip == 1, (first_event-data_start)<1)
%                 press(ic).f_chunk = NaN(1,[data_start+data_end+1]);
%                 continue
%             end
%             press(ic).f_chunk = [press(ic).f_chunk; data_tc(first_event-data_start:first_event+data_end,ic)']; %data_tc(press_fr(ip,:),ic)']; %
%         else   %if there is no event within the window then produce a random fake event (not counting fake events as real, using them as place holders only)
%             press(ic).event_ind(:,ip) = NaN;
%             %             wind = press_fr(ip,pre_buffer-win_press_start:pre_buffer+win_press_end);
%             %             if (wind(1)-data_start)<1
%             %                 wind = wind(find((wind-data_start)==1):length(wind));
%             %             end
%             %             first_event = wind(randsample(length(wind),1));
%             wind = press_fr(ip,pre_buffer-win_press_start:pre_buffer+win_release_end);
%             first_event = wind(1);
%             press(ic).good_event(:,ip) = NaN;
%             press(ic).event(:,ip) = 0;
%             press(ic).f_chunk = [press(ic).f_chunk; data_tc([first_event-data_start]:[first_event+data_end],ic)']; %data_tc(press_fr(ip,:),ic)']; %
%         end
%     end
%     press(ic).f_base = mean(press(ic).f_chunk(:,data_start-round(50./double(ifi)):data_start+1),2); %13:15),2); %
%     press(ic).df_chunk = bsxfun(@minus, press(ic).f_chunk, press(ic).f_base);
%     press(ic).dfoverf_chunk = bsxfun(@rdivide, press(ic).df_chunk, press(ic).f_base);
% end
% 
% [press_success_tc, press_success] = extractEvent_2P(data_tc, trial_outcome.success_ptime, frame_info, pre_buffer, ceil(3000/double(ifi)), opt, double(ifi), data_start,data_end);
% [press_fail_tc, press_fail] = extractEvent_2P(data_tc, trial_outcome.early_ptime, frame_info, pre_buffer, ceil(3000/double(ifi)), opt, double(ifi), data_start,data_end);
% 
% %% 3- Find periods around lick ITI events and detect
% %find licks
% lick = [];
% if exist('lick_data', 'var') && ~isempty(lick_data)
%     pre_buffer = ceil(1000/double(ifi));
%     post_buffer = ceil(1000/double(ifi));
%     lick_tc = [];
%     lick_fr = zeros(length(lick_data.bout_onsetTime),pre_buffer+post_buffer+1);
%     for is = 1:length(lick_data.bout_onsetTime)    %this could be problematic due to unremoved NaNs
%         if lick_data.bout_onsetTime(is) < length(frame_info.counter)
%             lick_fr(is,:) = [frame_info.counter(lick_data.bout_onsetTime(is))-pre_buffer:frame_info.counter(lick_data.bout_onsetTime(is))+post_buffer];
%             if and(min(lick_fr(is,:),[],2)>0, max(lick_fr(is,:),[],2)<size(data_tc,1))
%                 lick_tc = cat(3, lick_tc, data_tc(lick_fr(is,:),:));
%             end
%         else
%             lick_fr(is:end,:) = [];
%             break
%         end
%     end
%     
%     single_lick_tc = [];
%     single_lick_fr = zeros(length(lick_data.single_lickTime),pre_buffer+post_buffer+1);
%     for is = 1:length(lick_data.single_lickTime)    %this could be problematic due to unremoved NaNs
%         if lick_data.single_lickTime(is) < length(frame_info.counter)
%             single_lick_fr(is,:) = [frame_info.counter(lick_data.single_lickTime(is))-pre_buffer:frame_info.counter(lick_data.single_lickTime(is))+post_buffer];
%             if and(min(single_lick_fr(is,:),[],2)>0, max(single_lick_fr(is,:),[],2)<size(data_tc,1))
%                 single_lick_tc = cat(3, single_lick_tc, data_tc(single_lick_fr(is,:),:));
%             end
%         else
%             single_lick_fr(is:end,:) = [];
%             break
%         end
%     end
%     
%     %find all events around the lick
%     % thresh = 2.5;
%     win_lick_start = round(100./double(ifi));  %why this window?
%     win_lick_end = round(100./double(ifi));
%     for ic = 1:nIC
%         lick_diff = diff(squeeze(lick_tc(:,ic,:)),[],1);   %calcium traces for time around the lever press. Lever press happens in the center.
%         lick(ic).f_chunk = [];
%         lick(ic).ind = {};
%         lick(ic).good_ind = {};
%         for ip = 1:size(lick_diff,2)
%             [~,ind,~] = CellsortFindspikes(lick_diff(:,ip), thresh, dt, deconvtau, normalization);
%             %         ind = find(press_diff(:,ip)>thresh(ic));  %finds all calcium events around a lever press that are > threshold for that neuron
%             ind(find(diff(ind)==1)+1)=[];     %events must be separated by at least one frame.
%             double_events = find(diff(ind)<min_iei);
%             ind_good = ind;
%             ind_good([double_events (double_events+1)])=[];
%             lick(ic).ind{ip} = ind;
%             lick(ic).good_ind{ip} = ind_good;
%             if intersect(find(ind>pre_buffer-win_lick_start), find(ind<pre_buffer+win_lick_end))  %if there are events within the window
%                 wind = lick_fr(ip,pre_buffer-win_lick_start:pre_buffer+win_lick_end);
%                 if intersect(find(ind_good>pre_buffer-win_lick_start), find(ind_good<pre_buffer+win_lick_end))  %if ind_good has an event within the window
%                     ii = ind_good(intersect(find(ind_good>pre_buffer-win_lick_start), find(ind_good<pre_buffer+win_lick_end))); % ii = events within the window
%                     if size(ii,1)>1
%                         [minn, ind] = min(abs(ii-pre_buffer),[],1);
%                         ii = ii(ind);
%                     end
%                     lick(ic).event_ind(:,ip) = ii;
%                     first_event = lick_fr(ip,lick(ic).event_ind(:,ip));
%                     lick(ic).good_event(:,ip) = 1;
%                     lick(ic).event(:,ip) = 1;
%                 else
%                     ii = ind(intersect(find(ind>pre_buffer-win_lick_start), find(ind<pre_buffer+win_lick_end)));
%                     if size(ii,1)>1
%                         [nmin, imin] = min(abs(ii-pre_buffer),[],1);
%                         ii = ii(imin);
%                     end
%                     lick(ic).event_ind(:,ip) = ii;
%                     first_event = lick_fr(ip,lick(ic).event_ind(:,ip));
%                     lick(ic).good_event(:,ip) = 0;
%                     lick(ic).event(:,ip) = 1;
%                 end
%                 if and(ip == 1, (first_event-data_start)<1)
%                     lick(ic).f_chunk = NaN(1,[data_start+data_end+1]);
%                     continue
%                 end
%                 lick(ic).f_chunk = [lick(ic).f_chunk; data_tc(first_event-data_start:first_event+data_end,ic)']; %data_tc(press_fr(ip,:),ic)']; %
%             else   %if there is no event within the window then produce a random fake event (not counting fake events as real, using them as place holders only)
%                 lick(ic).event_ind(:,ip) = NaN;
%                 %             wind = press_fr(ip,pre_buffer-win_press_start:pre_buffer+win_press_end);
%                 %             if (wind(1)-data_start)<1
%                 %                 wind = wind(find((wind-data_start)==1):length(wind));
%                 %             end
%                 %             first_event = wind(randsample(length(wind),1));
%                 wind = lick_fr(ip,pre_buffer-win_lick_start:pre_buffer+win_lick_end);
%                 first_event = wind(1);
%                 lick(ic).good_event(:,ip) = NaN;
%                 lick(ic).event(:,ip) = 0;
%                 lick(ic).f_chunk = [lick(ic).f_chunk; data_tc([first_event-data_start]:[first_event+data_end],ic)']; %data_tc(press_fr(ip,:),ic)']; %
%             end
%         end
%         lick(ic).f_base = mean(lick(ic).f_chunk(:,data_start-round(50./double(ifi)):data_start+1),2); %13:15),2); %
%         lick(ic).df_chunk = bsxfun(@minus, lick(ic).f_chunk, lick(ic).f_base);
%         lick(ic).dfoverf_chunk = bsxfun(@rdivide, lick(ic).df_chunk, lick(ic).f_base);
%     end
%     
%     for ic = 1:nIC
%         single_lick_diff = diff(squeeze(single_lick_tc(:,ic,:)),[],1);   %calcium traces for time around the lever press. Lever press happens in the center.
%         single_lick(ic).f_chunk = [];
%         single_lick(ic).ind = {};
%         single_lick(ic).good_ind = {};
%         for ip = 1:size(single_lick_diff,2)
%             [~,ind,~] = CellsortFindspikes(single_lick_diff(:,ip), thresh, dt, deconvtau, normalization);
%             %         ind = find(press_diff(:,ip)>thresh(ic));  %finds all calcium events around a lever press that are > threshold for that neuron
%             ind(find(diff(ind)==1)+1)=[];     %events must be separated by at least one frame.
%             double_events = find(diff(ind)<min_iei);
%             ind_good = ind;
%             ind_good([double_events (double_events+1)])=[];
%             single_lick(ic).ind{ip} = ind;
%             single_lick(ic).good_ind{ip} = ind_good;
%             if intersect(find(ind>pre_buffer-win_lick_start), find(ind<pre_buffer+win_lick_end))  %if there are events within the window
%                 wind = single_lick_fr(ip,pre_buffer-win_lick_start:pre_buffer+win_lick_end);
%                 if intersect(find(ind_good>pre_buffer-win_lick_start), find(ind_good<pre_buffer+win_lick_end))  %if ind_good has an event within the window
%                     ii = ind_good(intersect(find(ind_good>pre_buffer-win_lick_start), find(ind_good<pre_buffer+win_lick_end))); % ii = events within the window
%                     if size(ii,1)>1
%                         [minn, ind] = min(abs(ii-pre_buffer),[],1);
%                         ii = ii(ind);
%                     end
%                     single_lick(ic).event_ind(:,ip) = ii;
%                     first_event = single_lick_fr(ip,single_lick(ic).event_ind(:,ip));
%                     single_lick(ic).good_event(:,ip) = 1;
%                     single_lick(ic).event(:,ip) = 1;
%                 else
%                     ii = ind(intersect(find(ind>pre_buffer-win_lick_start), find(ind<pre_buffer+win_lick_end)));
%                     if size(ii,1)>1
%                         [nmin, imin] = min(abs(ii-pre_buffer),[],1);
%                         ii = ii(imin);
%                     end
%                     single_lick(ic).event_ind(:,ip) = ii;
%                     first_event = single_lick_fr(ip,single_lick(ic).event_ind(:,ip));
%                     single_lick(ic).good_event(:,ip) = 0;
%                     single_lick(ic).event(:,ip) = 1;
%                 end
%                 if and(ip == 1, (first_event-data_start)<1)
%                     single_lick(ic).f_chunk = NaN(1,[data_start+data_end+1]);
%                     continue
%                 end
%                 single_lick(ic).f_chunk = [single_lick(ic).f_chunk; data_tc(first_event-data_start:first_event+data_end,ic)']; %data_tc(press_fr(ip,:),ic)']; %
%             else   %if there is no event within the window then produce a random fake event (not counting fake events as real, using them as place holders only)
%                 single_lick(ic).event_ind(:,ip) = NaN;
%                 %             wind = press_fr(ip,pre_buffer-win_press_start:pre_buffer+win_press_end);
%                 %             if (wind(1)-data_start)<1
%                 %                 wind = wind(find((wind-data_start)==1):length(wind));
%                 %             end
%                 %             first_event = wind(randsample(length(wind),1));
%                 wind = single_lick_fr(ip,pre_buffer-win_lick_start:pre_buffer+win_lick_end);
%                 first_event = wind(1);
%                 single_lick(ic).good_event(:,ip) = NaN;
%                 single_lick(ic).event(:,ip) = 0;
%                 single_lick(ic).f_chunk = [single_lick(ic).f_chunk; data_tc([first_event-data_start]:[first_event+data_end],ic)']; %data_tc(press_fr(ip,:),ic)']; %
%             end
%         end
%         single_lick(ic).f_base = mean(single_lick(ic).f_chunk(:,data_start-round(50./double(ifi)):data_start+1),2); %13:15),2); %
%         single_lick(ic).df_chunk = bsxfun(@minus, single_lick(ic).f_chunk, single_lick(ic).f_base);
%         single_lick(ic).dfoverf_chunk = bsxfun(@rdivide, single_lick(ic).df_chunk, single_lick(ic).f_base);
%     end
%     
% end
% 
% itiR = [];
% if isfield(trial_outcome, 'ItiunexpReward')
%     pre_buffer = ceil(1000/double(ifi));
%     post_buffer = ceil(1000/double(ifi));
%     itiR_tc = [];
%     itiR_fr = zeros(length(trial_outcome.ItiunexpReward),pre_buffer+post_buffer+1);
%     for is = 1:length(trial_outcome.ItiunexpReward)    %this could be problematic due to unremoved NaNs
%         if trial_outcome.ItiunexpReward(is) < length(frame_info.counter)
%             itiR_fr(is,:) = [frame_info.counter(trial_outcome.ItiunexpReward(is))-pre_buffer:frame_info.counter(trial_outcome.ItiunexpReward(is))+post_buffer];
%             if and(min(itiR_fr(is,:),[],2)>0, max(itiR_fr(is,:),[],2)<size(data_tc,1))
%                 itiR_tc = cat(3, itiR_tc, data_tc(itiR_fr(is,:),:));
%             end
%         else
%             itiR_fr(is:end,:) = [];
%             break
%         end
%     end
%     
%     %find all events around the iti unexpected reward
%     % thresh = 2.5;
%    
%     for ic = 1:nIC
%         itiR_diff = diff(squeeze(itiR_tc(:,ic,:)),[],1);
%         itiR(ic).f_chunk = [];
%         itiR(ic).ind = {};
%         itiR(ic).good_ind = {};
%         for is = 1:size(itiR_diff,2)
%             [~,ind,~] = CellsortFindspikes(itiR_diff(:,is), thresh, dt, deconvtau, normalization);
%             %         ind = find(success_diff(:,is)>thresh(ic));
%             ind(find(diff(ind)==1)+1)=[];
%             double_events = find(diff(ind)<min_iei);
%             ind_good = ind;
%             ind_good([double_events (double_events+1)])=[];
%             itiR(ic).ind{is} = ind;
%             itiR(ic).good_ind{is} = ind_good;
%             if intersect(find(ind>pre_buffer-win_release_start), find(ind<pre_buffer+win_release_end))
%                 wind = itiR_fr(is,pre_buffer-win_release_start:pre_buffer+win_release_end);
%                 if intersect(find(ind_good>pre_buffer-win_release_start), find(ind_good<pre_buffer+win_release_end))
%                     itiR(ic).event_ind(:,is) =ind_good(intersect(find(ind_good>pre_buffer-win_release_start), find(ind_good<pre_buffer+win_release_end)));
%                     first_event = itiR_fr(is,itiR(ic).event_ind(:,is));
%                     itiR(ic).good_event(:,is) = 1;
%                     itiR(ic).event(:,is) = 1;
%                 else
%                     itiR(ic).event_ind(:,is) = min(ind(intersect(find(ind>pre_buffer-win_release_start), find(ind<pre_buffer+win_release_end))),[],1);
%                     first_event = itiR_fr(is,itiR(ic).event_ind(:,is));
%                     itiR(ic).good_event(:,is) = 0;
%                     itiR(ic).event(:,is) = 1;
%                 end
%                 itiR(ic).f_chunk = [itiR(ic).f_chunk; data_tc(first_event-data_start:first_event+data_end,ic)']; %data_tc(success_fr(is,:),ic)']; %
%             else
%                 itiR(ic).event_ind(:,is) = NaN;
%                 wind = itiR_fr(is,pre_buffer-win_release_start:pre_buffer+win_release_end);
%                 first_event = wind(1);
%                 itiR(ic).good_event(:,is) = NaN;
%                 itiR(ic).event(:,is) = 0;
%                 itiR(ic).f_chunk = [itiR(ic).f_chunk; data_tc(first_event-data_start:first_event+data_end,ic)']; %data_tc(success_fr(is,:),ic)']; %
%             end
%         end
%         itiR(ic).f_base = mean(itiR(ic).f_chunk(:,data_start-round(50./double(ifi)):data_start+1),2); %13:15),2); %
%         itiR(ic).df_chunk = bsxfun(@minus, itiR(ic).f_chunk, itiR(ic).f_base);
%         itiR(ic).dfoverf_chunk = bsxfun(@rdivide, itiR(ic).df_chunk, itiR(ic).f_base);
%     end
%     
% else
%     itiR_tc = nan(size(success_tc));
% end
% 
% omitR = [];
% if isfield(trial_outcome, 'omitReward')
%     pre_buffer = ceil(1000/double(ifi));
%     post_buffer = ceil(1000/double(ifi));
%     omitR_tc = [];
%     omitR_fr = zeros(length(trial_outcome.omitReward),pre_buffer+post_buffer+1);
%     for is = 1:length(trial_outcome.omitReward)    %this could be problematic due to unremoved NaNs
%         if trial_outcome.omitReward(is) < length(frame_info.counter)
%             omitR_fr(is,:) = [frame_info.counter(trial_outcome.omitReward(is))-pre_buffer:frame_info.counter(trial_outcome.omitReward(is))+post_buffer];
%             if and(min(omitR_fr(is,:),[],2)>0, max(omitR_fr(is,:),[],2)<size(data_tc,1))
%                 omitR_tc = cat(3, omitR_tc, data_tc(omitR_fr(is,:),:));
%             end
%         else
%             omitR_fr(is:end,:) = [];
%             break
%         end
%     end
%     
%     %find all events around the iti unexpected reward
%     % thresh = 2.5;
%    
%     for ic = 1:nIC
%         omitR_diff = diff(squeeze(omitR_tc(:,ic,:)),[],1);
%         omitR(ic).f_chunk = [];
%         omitR(ic).ind = {};
%         omitR(ic).good_ind = {};
%         for is = 1:size(omitR_diff,2)
%             [~,ind,~] = CellsortFindspikes(omitR_diff(:,is), thresh, dt, deconvtau, normalization);
%             %         ind = find(success_diff(:,is)>thresh(ic));
%             ind(find(diff(ind)==1)+1)=[];
%             double_events = find(diff(ind)<min_iei);
%             ind_good = ind;
%             ind_good([double_events (double_events+1)])=[];
%             omitR(ic).ind{is} = ind;
%             omitR(ic).good_ind{is} = ind_good;
%             if intersect(find(ind>pre_buffer-win_release_start), find(ind<pre_buffer+win_release_end))
%                 wind = omitR_fr(is,pre_buffer-win_release_start:pre_buffer+win_release_end);
%                 if intersect(find(ind_good>pre_buffer-win_release_start), find(ind_good<pre_buffer+win_release_end))
%                     omitR(ic).event_ind(:,is) =ind_good(intersect(find(ind_good>pre_buffer-win_release_start), find(ind_good<pre_buffer+win_release_end)));
%                     first_event = omitR_fr(is,itiR(ic).event_ind(:,is));
%                     omitR(ic).good_event(:,is) = 1;
%                     omitR(ic).event(:,is) = 1;
%                 else
%                     omitR(ic).event_ind(:,is) = min(ind(intersect(find(ind>pre_buffer-win_release_start), find(ind<pre_buffer+win_release_end))),[],1);
%                     first_event = omitR_fr(is,itiR(ic).event_ind(:,is));
%                     omitR(ic).good_event(:,is) = 0;
%                     omitR(ic).event(:,is) = 1;
%                 end
%                 omitR(ic).f_chunk = [omitR(ic).f_chunk; data_tc(first_event-data_start:first_event+data_end,ic)']; %data_tc(success_fr(is,:),ic)']; %
%             else
%                 omitR(ic).event_ind(:,is) = NaN;
%                 wind = omitR_fr(is,pre_buffer-win_release_start:pre_buffer+win_release_end);
%                 first_event = wind(1);
%                 omitR(ic).good_event(:,is) = NaN;
%                 omitR(ic).event(:,is) = 0;
%                 omitR(ic).f_chunk = [omitR(ic).f_chunk; data_tc(first_event-data_start:first_event+data_end,ic)']; %data_tc(success_fr(is,:),ic)']; %
%             end
%         end
%         omitR(ic).f_base = mean(omitR(ic).f_chunk(:,data_start-round(50./double(ifi)):data_start+1),2); %13:15),2); %
%         omitR(ic).df_chunk = bsxfun(@minus, omitR(ic).f_chunk, omitR(ic).f_base);
%         omitR(ic).dfoverf_chunk = bsxfun(@rdivide, omitR(ic).df_chunk, omitR(ic).f_base);
%     end
%     
% else
%     omitR_tc = nan(size(success_tc));
% end
% 
% if ~isempty(lick_data)
%     save([dest_sub '_evoked_events.mat'], 'success_tc', 'fail_tc', 'press_tc', 'success', 'fail', 'press', 'lick', 'single_lick', ...
%         'win_release_end', 'win_press_end', 'win_press_start', 'win_lick_start', 'win_lick_end','min_iei', 'thresh', 'pre_buffer', 'post_buffer', ...
%         'itiR_tc', 'itiR', 'omitR_tc', 'omitR')
% else
%     save([dest_sub '_evoked_events.mat'], 'success_tc', 'fail_tc', 'press_tc', 'success', 'fail', 'press', ...
%         'win_release_end', 'win_press_end', 'win_press_start', 'min_iei', 'thresh', 'pre_buffer', 'post_buffer', 'press_success', 'press_fail',...
%         'press_success_tc', 'press_fail_tc')
% end
% 
