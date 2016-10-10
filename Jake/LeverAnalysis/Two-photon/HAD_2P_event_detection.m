%code to detect calcium transients
% 1- Find periods with no lever activity and detect "spontaneous" events
% 2- Find periods around release events and detect
% 3- Find periods around press events and detect

%load data and behavior info
load([dest, '_ROI_TCs.mat'])
load([dest, '_dFOverF_TC.mat'])
nIC = size(data_tc,2);
load([dest '_frame_times.mat'])
load([dest '_parse_behavior.mat'])
%% 1- Find periods with no lever activity and detect "spontaneous" events
%extract spontaneous windows
ifi = (frame_times(end)-frame_times(1))/length(frame_times);
Sampeling_rate = 1000/ifi;
pre_buffer = ceil(300/double(ifi));
post_buffer = ceil(1000/double(ifi));
fr_lever = [];   %fr_lever will be a series of frame numbers corresponding to lever events +&- the buffers
for ipress = 1:length(lever.press)    %this could be problematic due to unremoved NaNs
    fr_lever = [fr_lever frame_info.counter(round(lever.press(ipress)))-pre_buffer:frame_info.counter(round(lever.press(ipress)))+post_buffer];
end 
for irelease = 1:length(lever.release)    %this could be problematic due to unremoved NaNs
    fr_lever = [fr_lever frame_info.counter(round(lever.release(irelease)))-pre_buffer:frame_info.counter(round(lever.release(irelease)))+post_buffer];
end 

fr_lever = unique(fr_lever);
fr_lever(find(fr_lever<1)) = [];
fr_lever(find(fr_lever>size(data_tc,1))) = [];
fr_lever(find(fr_lever>length(frame_times))) = [];  %only looks at values that correspond to actual frames. 

data_tc_spont = data_tc; 
data_tc_spont(fr_lever,:) = [];  %removes frames surrounding lever events 

%find all events to get rate
events = [];
events_diff_x = zeros(size(data_tc_spont,1)-1,nIC);
events_ind = {};
events_rate = zeros(1,nIC);

for ic = 1:nIC
    events_diff_x(:,ic) = diff(data_tc_spont(:,ic),[],1);
end

disp('manually select thresholds for each TC')
for ii = size(events_diff_x,2):-1:1; 
    figure; plot(events_diff_x(:,ii));
    title(['diffTC for TC #' num2str(ii)]);
end
%pause
%thresh = [1100, 600, 300, 450, 350, 500, 500, 250, 500, 300, 400, 400, 500, 500, 500, 500, 500, 700, 500, 500, 400, 350]; % manually selected thresholds for 160319_img38
%thresh = [500, 450, 250, 400, 450, 300, 200, 200, 300, 350, 300, 300, 300, 250, 300, 300, 300, 250]; %160220_img38
%thresh= [250, 400, 400, 300, 400, 300, 400, 200, 400, 250, 600, 250, 300, 500, 400, 500, 450, 750, 1800, 250, 600, 250, 250, 700]; %160324_img41
%thresh = [500, 500, 600, 600, 600, 500, 500, 500, 600, 600, 1000, 1000, 1000, 400, 500, 750, 600]; %160203_img36
%thresh = [350, 200, 350, 350, 300, 400, 500, 600, 500, 400, 500, 350, 400, 400, 500]; %160202_img36
%thresh = [700, 750, 700, 700, 1000, 1100, 1200, 1000, 1100, 550, 1000]; %151215_img32
%thresh = [500, 1000, 1000, 600, 300, 400, 450, 700, 1000, 900, 800, 700, 750, 800]; %160609_img46
%thresh = [500, 400, 900, 750, 800, 500, 1000, 750, 500, 400, 400, 700, 500, 750, 700, 800, 600, 600, 500, 450, 500, 450, 600, 750, 900, 500, 500]; %160607_img46
%thresh = [750, 600, 600, 2500, 1000, 1200, 900, 750, 900, 1000, 2000]; %img53_160723 Not usable since TCs will not be in the same order. Find saved variables and use them to make this set usable again. 
thresh = [1500, 750, 600, 700, 650, 1000, 750, 500, 450, 1000, 700, 600, 750, 800, 1000, 1100]; %img53_160724

for ic = 1:nIC
    events_ind{ic} = find(events_diff_x(:,ic)>thresh(ic));
    events_ind{ic}(find(diff(events_ind{ic})==1)+1)=[];
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

save([dest_sub '_spont_events.mat'], 'events', 'data_tc_spont', 'fr_lever', 'thresh', 'events_rate', 'data_start', 'data_end');

%% 2- Find periods around release events and detect
%find successes, earlies
pre_buffer = ceil(1000/double(ifi));
post_buffer = ceil(1000/double(ifi));
success_fr = zeros(length(trial_outcome.success_time),pre_buffer+post_buffer+1);
fail_fr = zeros(length(trial_outcome.early_time),pre_buffer+post_buffer+1);
success_tc = [];
fail_tc = [];
for is = 1:length(trial_outcome.success_time)    %this could be problematic due to unremoved NaNs
    success_fr(is,:) = [frame_info.counter(trial_outcome.success_time(is))-pre_buffer:frame_info.counter(trial_outcome.success_time(is))+post_buffer];
    if and(min(success_fr(is,:),[],2)>0, max(success_fr(is,:),[],2)<size(data_tc,1)) 
        success_tc = cat(3, success_tc, data_tc(success_fr(is,:),:));
    end
end 
for ie = 1:length(trial_outcome.early_time)    %this could be problematic due to unremoved NaNs
    fail_fr(ie,:) = [frame_info.counter(trial_outcome.early_time(ie))-pre_buffer:frame_info.counter(trial_outcome.early_time(ie))+post_buffer];
    if and(min(fail_fr(ie,:),[],2)>0, max(fail_fr(ie,:),[],2)<size(data_tc,1)) 
        fail_tc = cat(3, fail_tc, data_tc(fail_fr(ie,:),:));
    end
end 


%find all events during success trials
success = [];
win_release_end = round(200./double(ifi));
for ic = 1:nIC
    success_diff = diff(squeeze(success_tc(:,ic,:)),[],1);
    success(ic).f_chunk = [];
    success(ic).ind = {};
    success(ic).good_ind = {};
    for is = 1:size(success_diff,2)
        ind = find(success_diff(:,is)>thresh(ic));
        ind(find(diff(ind)==1)+1)=[];
        double_events = find(diff(ind)<min_iei);
        ind_good = ind;
        ind_good([double_events (double_events+1)])=[];
        success(ic).ind{is} = ind;
        success(ic).good_ind{is} = ind_good;
        if intersect(find(ind>pre_buffer-1), find(ind<pre_buffer+win_release_end))
            wind = success_fr(is,pre_buffer-1:pre_buffer+win_release_end);
            if intersect(find(ind_good>pre_buffer-1), find(ind_good<pre_buffer+win_release_end))
                success(ic).event_ind(:,is) =ind_good(intersect(find(ind_good>pre_buffer-1), find(ind_good<pre_buffer+win_release_end)));
                first_event = success_fr(is,success(ic).event_ind(:,is));
                success(ic).good_event(:,is) = 1;
                success(ic).event(:,is) = 1;
            else
                success(ic).event_ind(:,is) = min(ind(intersect(find(ind>pre_buffer-1), find(ind<pre_buffer+win_release_end))),[],1);
                first_event = success_fr(is,success(ic).event_ind(:,is));
                success(ic).good_event(:,is) = 0;
                success(ic).event(:,is) = 1;
            end
            success(ic).f_chunk = [success(ic).f_chunk; data_tc(first_event-data_start:first_event+data_end,ic)']; %data_tc(success_fr(is,:),ic)']; %
        else
            success(ic).event_ind(:,is) = NaN;
            wind = success_fr(is,pre_buffer-1:pre_buffer+1);
            first_event = wind(1);
            success(ic).good_event(:,is) = NaN;
            success(ic).event(:,is) = 0;
            success(ic).f_chunk = [success(ic).f_chunk; data_tc(first_event-data_start:first_event+data_end,ic)']; %data_tc(success_fr(is,:),ic)']; %
        end
    end
    success(ic).f_base = mean(success(ic).f_chunk(:,data_start-round(50./double(ifi)):data_start+1),2); %13:15),2); %
    success(ic).df_chunk = bsxfun(@minus, success(ic).f_chunk, success(ic).f_base);
    success(ic).dfoverf_chunk = bsxfun(@rdivide, success(ic).df_chunk, success(ic).f_base);
end

%find all events during early trials
fail = [];
for ic = 1:nIC
    fail_diff = diff(squeeze(fail_tc(:,ic,:)),[],1);
    fail(ic).f_chunk = [];
    fail(ic).ind = {};
    fail(ic).good_ind = {};
    for ie = 1:size(fail_diff,2)
        ind = find(fail_diff(:,ie)>thresh(ic));
        ind(find(diff(ind)==1)+1)=[];
        double_events = find(diff(ind)<min_iei);
        ind_good = ind;
        ind_good([double_events (double_events+1)])=[];
        fail(ic).ind{ie} = ind;
        fail(ic).good_ind{ie} = ind_good;
        if intersect(find(ind>pre_buffer-1), find(ind<pre_buffer+win_release_end))
            wind = fail_fr(ie,pre_buffer-1:pre_buffer+win_release_end);
            if intersect(find(ind_good>pre_buffer-1), find(ind_good<pre_buffer+win_release_end))
                fail(ic).event_ind(:,ie) =ind_good(intersect(find(ind_good>pre_buffer-1), find(ind_good<pre_buffer+win_release_end)));
                first_event = fail_fr(ie,fail(ic).event_ind(:,ie));
                fail(ic).good_event(:,ie) = 1;
                fail(ic).event(:,ie) = 1;
            else
                fail(ic).event_ind(:,ie) = min(ind(intersect(find(ind>pre_buffer-1), find(ind<pre_buffer+win_release_end))),[],1);
                first_event = fail_fr(ie,fail(ic).event_ind(:,ie));
                fail(ic).good_event(:,ie) = 0;
                fail(ic).event(:,ie) = 1;
            end
            fail(ic).f_chunk = [fail(ic).f_chunk; data_tc(first_event-data_start:first_event+data_end,ic)']; %data_tc(fail_fr(ie,:),ic)']; %
        else
            fail(ic).event_ind(:,ie) = NaN;
            wind = fail_fr(ie,pre_buffer-1:pre_buffer+win_release_end);
            first_event = wind(1);
            fail(ic).good_event(:,ie) = NaN;
            fail(ic).event(:,ie) = 0;
            fail(ic).f_chunk = [fail(ic).f_chunk; data_tc(first_event-data_start:first_event+data_end,ic)']; %data_tc(fail_fr(ie,:),ic)']; %
        end
    end
    fail(ic).f_base = mean(fail(ic).f_chunk(:,data_start-round(50./double(ifi)):data_start+1),2); %13:15),2); %
    fail(ic).df_chunk = bsxfun(@minus, fail(ic).f_chunk, fail(ic).f_base);
    fail(ic).dfoverf_chunk = bsxfun(@rdivide, fail(ic).df_chunk, fail(ic).f_base);
end

%% 3- Find periods around press events and detect
%find presses
pre_buffer = ceil(1000/double(ifi));
post_buffer = ceil(1000/double(ifi));
press_fr = [];  %press fr will equal a series of frame numbers around each lever press
press_tc = [];
ind_long = intersect(find(cell2mat(input.holdTimesMs)>=500),find(~isnan(trial_outcome.ind_press_prerelease)));
for ip = 1:length(ind_long)
    iT = ind_long(ip);
    leverTimes = round((cell2mat(input.leverTimesUs(iT))-input.counterTimesUs{1}(1))./1000);
    if leverTimes(trial_outcome.ind_press_prerelease(iT)) > 0
        fr_range = [frame_info.counter(leverTimes(trial_outcome.ind_press_prerelease(iT)))-pre_buffer:frame_info.counter(leverTimes(trial_outcome.ind_press_prerelease(iT)))+post_buffer];
        press_fr = [press_fr; fr_range];
        if and(min(fr_range,[],2)>0, max(fr_range,[],2)<size(data_tc,1)) 
           press_tc = cat(3, press_tc, data_tc(fr_range,:));
        end
    end
end 

%find all events around the press
press = [];
win_press_start = round(300./double(ifi));  %why this window?
win_press_end = round(100./double(ifi));
for ic = 1:nIC
    press_diff = diff(squeeze(press_tc(:,ic,:)),[],1);   %calcium traces for time around the lever press. Lever press happens in the center.
    press(ic).f_chunk = [];
    press(ic).ind = {};
    press(ic).good_ind = {};
    for ip = 1:size(press_diff,2)
        ind = find(press_diff(:,ip)>thresh(ic));  %finds all calcium events around a lever press that are > threshold for that neuron
        ind(find(diff(ind)==1)+1)=[];     %events must be separated by at least one frame. 
        double_events = find(diff(ind)<min_iei);
        ind_good = ind;
        ind_good([double_events (double_events+1)])=[];
        press(ic).ind{ip} = ind;
        press(ic).good_ind{ip} = ind_good;
        if intersect(find(ind>pre_buffer-win_press_start), find(ind<pre_buffer+win_press_end))  %if there are events within the window
            wind = press_fr(ip,pre_buffer-win_press_start:pre_buffer+win_press_end);
            if intersect(find(ind_good>pre_buffer-win_press_start), find(ind_good<pre_buffer+win_press_end))  %if ind_good has an event within the window
                ii = ind_good(intersect(find(ind_good>pre_buffer-win_press_start), find(ind_good<pre_buffer+win_press_end))); % ii = events within the window
                if size(ii,1)>1
                    [minn, ind] = min(abs(ii-pre_buffer),[],1);
                    ii = ii(ind);
                end
                press(ic).event_ind(:,ip) = ii;                 
                first_event = press_fr(ip,press(ic).event_ind(:,ip)); 
                press(ic).good_event(:,ip) = 1;
                press(ic).event(:,ip) = 1;
                else
                ii = ind(intersect(find(ind>pre_buffer-win_press_start), find(ind<pre_buffer+win_press_end)));
                if size(ii,1)>1
                    [nmin, imin] = min(abs(ii-pre_buffer),[],1);
                    ii = ii(imin);
                end
                press(ic).event_ind(:,ip) = ii;
                first_event = press_fr(ip,press(ic).event_ind(:,ip));
                press(ic).good_event(:,ip) = 0;
                press(ic).event(:,ip) = 1;
            end
            if and(ip == 1, (first_event-data_start)<1)
                    press(ic).f_chunk = NaN(1,[data_start+data_end+1]);
                    continue
            end
            press(ic).f_chunk = [press(ic).f_chunk; data_tc(first_event-data_start:first_event+data_end,ic)']; %data_tc(press_fr(ip,:),ic)']; %
        else   %if there is no event within the window then produce a random fake event (not counting fake events as real, using them as place holders only)
            press(ic).event_ind(:,ip) = NaN;
            wind = press_fr(ip,pre_buffer-win_press_start:pre_buffer+win_press_end);
            if (wind(1)-data_start)<1
                wind = wind(find((wind-data_start)==1):length(wind));
            end
            first_event = wind(randsample(length(wind),1));
            press(ic).good_event(:,ip) = NaN;
            press(ic).event(:,ip) = 0;
            press(ic).f_chunk = [press(ic).f_chunk; data_tc([first_event-data_start]:[first_event+data_end],ic)']; %data_tc(press_fr(ip,:),ic)']; %
        end
    end
    press(ic).f_base = mean(press(ic).f_chunk(:,data_start-round(50./double(ifi)):data_start+1),2); %13:15),2); %
    press(ic).df_chunk = bsxfun(@minus, press(ic).f_chunk, press(ic).f_base);
    press(ic).dfoverf_chunk = bsxfun(@rdivide, press(ic).df_chunk, press(ic).f_base);
end

save([dest_sub '_evoked_events.mat'], 'success_tc', 'fail_tc', 'press_tc', 'success', 'fail', 'press', 'win_release_end', 'win_press_end', 'win_press_start', 'min_iei', 'thresh', 'pre_buffer', 'post_buffer')

