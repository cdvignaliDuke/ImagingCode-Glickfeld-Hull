function [trigger_movie, use_times, trigger_licks, trigger_licks_10ms] = trigger_movie_by_event(movie, frame_info, licking_data, event_times, pre_frames, post_frames)
% ---- cut the movie around event times
% movie - the movie in 2-dimantional matrix 
% frame_info: output of parse behavior syncs frame number to time in behavior
 
% --- use only events that are in the times of the movie, leave space for last event
lickTimes = licking_data.lickTimes;  %MWorks time of each recorded lick. 
licksByFrame = licking_data.licksByFrame; %# of licks which happened in each frame. Binned my MWorks counter times. 

%find frame times for the first and last frames that fall within boundaries imposed by size of the analysis window
last_time =  find(frame_info.counter<size(movie,2)-post_frames, 1, 'last');
first_time = find(frame_info.counter>pre_frames+1, 1, 'first');

%exclude event times which fall outside these boundaries
use_event_times = event_times(event_times < last_time);
use_event_times = use_event_times(use_event_times >first_time);

%allocate memory
trigger_movie =nan([length(use_event_times) size(movie,1), pre_frames + post_frames+1]);  %dims: 
trigger_licks = nan([length(use_event_times), pre_frames + post_frames+1]);
trigger_licks_10ms = nan([length(use_event_times), pre_frames*10 + post_frames*10+1]);
use_times =[];

%main forloop which uses event times to extract specific frames and lick traces
for i=1:length(use_event_times);
    c_time = use_event_times(i);
    
    %checks to make sure the event time has frames associated with it
    frame_no = frame_info.counter(c_time); %find the frame num of the event time
    if(isnan(frame_no));
        continue; % discard
    end
    use_times(end+1) = c_time; %keeps track of all the event times which were actually used
    
    %use specific frame #s to extract intervals of frames and lick trace
    trigger_movie(i,:,:) = movie(:, frame_no-pre_frames:frame_no+post_frames);
    if frame_no+post_frames > length(licksByFrame)
        ii=i;
        continue;
    end
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
if exist('ii')
    trigger_licks = trigger_licks(1:ii-1,:);
end
end
 