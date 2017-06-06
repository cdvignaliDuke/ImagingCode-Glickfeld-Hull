function [trigger_movie, use_times, trigger_licks, trigger_licks_10ms] = trigger_movie_by_event(movie, frame_info, ...
    event_times, pre_frames, post_frames, licking_data)
% ---- cut the movie around event times
% movie - the movie in 2-dimantional matrix 
% frame_info: output of parse behavior syncs frame number to time  in behavior

% --- use only events that are in the times of the movie, leave space for last event
last_time=  find(frame_info.counter<size(movie,2)-post_frames, 1, 'last');
use_event_times = event_times(event_times < last_time);
first_time = find(frame_info.counter>pre_frames+1, 1, 'first');
% first_time = 758175; %for img38_160320 remove first 25 trials
use_event_times = use_event_times(use_event_times >first_time);
trigger_movie =nan([length(use_event_times) size(movie,1), pre_frames + post_frames+1]);
use_times =[];

% if there is licking_data 
if ~isempty(licking_data)
    lickTimes = licking_data.lickTimes;
    licksByFrame = licking_data.licksByFrame;
    trigger_licks = nan([length(use_event_times), pre_frames + post_frames+1]);
    trigger_licks_10ms = nan([length(use_event_times), pre_frames*10 + post_frames*10+1]);
end

for i=1:length(use_event_times)
    c_time = use_event_times(i);
    % inefficient - but simple in the future use indexes
    % find the corresponding frame
    frame_no = frame_info.counter(c_time);
    if(isnan(frame_no))
        continue; % discard
    end    
    use_times(end+1) = c_time;
    
    trigger_movie(i,:,:) = movie(:, frame_no-pre_frames:frame_no+post_frames);
    
    trigger_licks(i,:) = licksByFrame(frame_no-pre_frames:frame_no+post_frames);
    
    %gather licks by 10ms bins
    frame_time = frame_info.times(frame_no);
    time_range = (frame_time-pre_frames*10):10:(frame_time+post_frames*10);
    licks_this_TC = [];
    for kk = 2:length(time_range); %find all the # of licks which occur in each time bin. 
        licks_this_bin = lickTimes(find(lickTimes>=time_range(kk-1) & lickTimes<time_range(kk)));
        licks_this_TC  = [licks_this_TC, licks_this_bin];
    end
    trigger_licks_10ms(i,:);
end

%Attempted to access frame_info.counter(629.805); index must be a positive integer or logical.

%Error in trigger_movie_by_event (line 18)
 %   frame_no = frame_info.counter(c_time);

%Error in HAD_cmp_success_fail (line 159)
 %   press_movie = trigger_movie_by_event(img, frame_info, ...
 