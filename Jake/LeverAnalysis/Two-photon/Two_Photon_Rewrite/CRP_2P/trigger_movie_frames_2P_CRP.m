function [trigger_movie, remove_event_idx, hold_dur, use_times] = trigger_movie_by_event_2P(movie, frame_info, ...
    event_times, pre_frames, post_frames, licking_data, hold_dur, do_alignLick)

% ---- cut the movie around event times
% movie - the movie in 2-dimantional matrix
% frame_info: output of parse behavior syncs frame number to time  in behavior

% --- use only events that are in the times of the movie, leave space for last event
last_time=  find(frame_info.counter<size(movie,2)-post_frames, 1, 'last');
use_event_times = event_times(event_times < last_time);
first_time = find(frame_info.counter>pre_frames+1, 1, 'first');
% first_time = 758175; %for img38_160320 remove first 25 trials
use_event_times = use_event_times(use_event_times >first_time);
if ~isempty(hold_dur)
    hold_dur = hold_dur(event_times < last_time & event_times > first_time);
end
trigger_movie =nan([length(use_event_times) size(movie,1), pre_frames + post_frames+1]);
use_times =[]; 

% if there is licking_data
if ~isempty(licking_data)
    licksByFrame = licking_data.licksByFrame;
end

%for each event time extract the frames and licking info
for event_num=1:length(use_event_times)
    
    %checks to make sure this is a valid event time and also logs the event time
    c_time = use_event_times(event_num);
    frame_no = frame_info.counter(c_time); % find the corresponding frame
    if(isnan(frame_no))
        continue; % discard
    end
    use_times(end+1) = c_time;
    
    %Either align TCs to the first lick or to the reward/cue
    if ~do_alignLick %input variable 
        trigger_movie(event_num,:,:) = movie(:, frame_no-pre_frames:frame_no+post_frames);
    else
        temp_lick = licksByFrame(frame_no:frame_no+post_frames);
        if ~isempty( find(temp_lick, 1, 'first') )
            frame_no = find(temp_lick, 1, 'first') + frame_no - 1;
        end
        trigger_movie(event_num,:,:) = movie(:, frame_no-pre_frames:frame_no+post_frames);
    end
end
remove_event_idx = find(~ismember(event_times, use_event_times));


