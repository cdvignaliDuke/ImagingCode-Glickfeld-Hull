function [trigger_piezo, piezo_info] = trigger_movie_piezo_2P_CRP(piezo_data, piezo_frame_ind, frame_info, use_event_times, trigger_vars);
% movement analysis

%allocae structure
 piezo_info.bout_onset = []; piezo_info.cue_onset = []; piezo_info.reward_onset = []; piezo_info.reward_late_onset = [];
 piezo_info.lickTrial = []; piezo_info.lickTrial_1000 = [];    piezo_info.no_lick_cue_to_500=[]; piezo_info.reward_late_onset_500 = [];
 piezo_info.reward_late_onset_250 = []; piezo_info.reward_prox = [];
 %define variables
 pre_frames = trigger_vars.pre_cue_frames;
 post_frames = trigger_vars.post_cue_frames;
 cue_2_rew_delay_ms = trigger_vars.cue_2_rew_delay_ms;
 trigger_piezo = nan([length(use_event_times), pre_frames + post_frames+1]);
 
 %exxtract movement trace from each trial
 for event_num=1:length(use_event_times)  
     
    %get licking TC for the interval around event time
    trigger_piezo(event_num,:) = licksByFrame(frame_no-pre_frames:frame_no+post_frames);
     
     
 end
 
 
 
