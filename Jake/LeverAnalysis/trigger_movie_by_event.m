function [trigger_movie, use_times] = trigger_movie_by_event(movie, frame_info, event_times  , pre_frames, post_frames)
% extract sequences of frames corresponding to trial types. Used in frame_by_frame_ROI_fig
 
%use only events that are in the times of the movie, leave space for last event
last_time=  find(frame_info.counter<size(movie,2)-post_frames, 1, 'last');
first_time = find(frame_info.counter>pre_frames+1, 1, 'first');
use_event_times = event_times(event_times >first_time & event_times < last_time);

%allocate memory
trigger_movie =nan([length(use_event_times) size(movie,1), pre_frames + post_frames+1]);
use_times =[];

%convert times to frames and extract a series of frames for each trial of the selected outcome type
for i=1:length(use_event_times)
    c_time = use_event_times(i);
    frame_no = frame_info.counter(c_time); %use event time to find frame#
    if(isnan(frame_no))
        continue; % discard
    end    
    use_times(end+1) = c_time;
    trigger_movie(i,:,:) = movie(:, frame_no-pre_frames:frame_no+post_frames); %trigger movie dim1=trial   dim2=linear pixels   dim3=frame num within that trial
end