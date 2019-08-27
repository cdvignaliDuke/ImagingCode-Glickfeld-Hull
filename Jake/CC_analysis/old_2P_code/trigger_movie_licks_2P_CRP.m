 function [trigger_licks, trigger_licks_10ms, lick] = trigger_movie_by_event_2P(frame_info, licking_data, use_event_times, ...
     pre_frames, post_frames, do_lickAna, cue_2_rew_delay_ms)
 
 % licking analysis
 lick.bout_onset = []; lick.cue_onset = []; lick.reward_onset = []; lick.reward_late_onset = [];
 lick.lickTrial = []; lick.lickTrial_1000 = [];    lick.no_lick_cue_to_500=[]; lick.reward_late_onset_500 = [];
 lick.reward_late_onset_250 = []; lick.reward_prox = [];
 
 % if there is licking_data
 if ~isempty(licking_data)
     lickTimes = licking_data.lickTimes;
     licksByFrame = licking_data.licksByFrame;
     trigger_licks = nan([length(use_event_times), pre_frames + post_frames+1]);
     trigger_licks_10ms = nan([length(use_event_times), pre_frames*10 + post_frames*10+1]);
 else
     trigger_licks = [];
     trigger_licks_10ms = [];
 end

for event_num=1:length(use_event_times)  
    
    %get licking TC for the interval around event time
    trigger_licks(event_num,:) = licksByFrame(frame_no-pre_frames:frame_no+post_frames);
    
    %Check if licking analysis is indicated -----------------------------
    if do_lickAna ~= 1
        continue
    end
    
    % find lick bout onset
    lick_reward = trigger_licks(event_num, pre_frames +1 :end); %index licks relative to cue
    lick_idx = find(lick_reward);
    if length(lick_idx) >= 4  %must be a minimum of four licks in the whole trace
        inter_lick = diff(lick_idx);
        for li = 1:length(inter_lick)
            if inter_lick(li) <= floor(200/frame_info.ifi) && length(inter_lick(li:end)) >= 3 && max(inter_lick(li: li+2 )) <= floor(200/frame_info.ifi) % each lick is less than 6 frames apart and at least three inter licks in bewteen
                lick.bout_onset(event_num) = lick_idx(li) + pre_frames + 1; % index relative to cue onset frame. So idx of 1 = first frame after cue onset
                break;
            elseif li == length(inter_lick)
                lick.bout_onset(event_num) = NaN;
            end
        end
    else
        lick.bout_onset(event_num) = NaN;
    end
    
    %find trials without any licks from cue : cue+500ms
    if sum(trigger_licks(event_num, pre_frames+1:pre_frames + 1 + floor(500/frame_info.ifi))) < 1
        lick.no_lick_cue_to_500(event_num) = 1;
    else
        lick.no_lick_cue_to_500(event_num) = 0;
    end
    
    % find first lick between cue and reward  
    lick_cue = trigger_licks(event_num, pre_frames + 1: pre_frames + floor(cue_2_rew_delay_ms/frame_info.ifi)+1);
    if sum(lick_cue) > 0
        lick.cue_onset(event_num) =  find(lick_cue,1,'first') + pre_frames;
    else
        lick.cue_onset(event_num) = NaN;
    end
    
    % find lick onset 150ms around reward
    lick_reward = trigger_licks(event_num, pre_frames +1 + floor(cue_2_rew_delay_ms/frame_info.ifi)- floor(200/frame_info.ifi): pre_frames + 1+ floor(cue_2_rew_delay_ms/frame_info.ifi) + floor(200/frame_info.ifi));
    if sum(lick_reward) > 0
        lick.reward_prox(event_num) = find(lick_reward,1,'first') + pre_frames +1 +floor(cue_2_rew_delay_ms/frame_info.ifi)- floor(200/frame_info.ifi);
    else
        lick.reward_prox(event_num) = NaN;
    end
    
    % find lick onset 250ms after reward
    lick_reward_late = trigger_licks(event_num, pre_frames + floor(250/frame_info.ifi) + 1: end);
    %must be no licking from rew-100:rew+250ms
    if sum(trigger_licks(event_num, pre_frames +1 + floor(cue_2_rew_delay_ms/frame_info.ifi)- floor(100/frame_info.ifi) : pre_frames +1 + floor(cue_2_rew_delay_ms/frame_info.ifi)- floor(250/frame_info.ifi)))==0
        if sum(lick_reward_late) > 1
            lick.reward_late_onset_250(event_num) = find(lick_reward_late,1,'first') + pre_frames + floor(250/frame_info.ifi) + 1;
        else
            lick.reward_late_onset_250(event_num) = NaN;
        end
    else
        lick.reward_late_onset_250(event_num) = NaN;
    end
    
    % find lick onset 500ms after reward
    lick_reward_late = trigger_licks(event_num, pre_frames + floor(500/frame_info.ifi) + 1: end);
    %must be no licking from rew-100:rew+250ms
    if sum(trigger_licks(event_num, pre_frames +1 + floor(cue_2_rew_delay_ms/frame_info.ifi)- floor(100/frame_info.ifi) : pre_frames +1 + floor(cue_2_rew_delay_ms/frame_info.ifi)- floor(500/frame_info.ifi)))==0
        if sum(lick_reward_late) > 1
            lick.reward_late_onset_500(event_num) = find(lick_reward_late,1,'first') + pre_frames + floor(500/frame_info.ifi) + 1;
        else
            lick.reward_late_onset_500(event_num) = NaN;
        end
    else
        lick.reward_late_onset_500(event_num) = NaN;
    end
    
    % find trial without lick bout  =============== just look for the NaNs  in the lick bout onset variable
%     if sum(trigger_licks(event_num, pre_frames - floor(100/frame_info.ifi):pre_frames + floor(500/frame_info.ifi))) > 1 % no more than 1 lick 100ms before
%         % cue and 100ms before reward
%         lick_cuereward = trigger_licks(event_num, pre_frames - floor(100/frame_info.ifi):pre_frames + floor(500/frame_info.ifi));
%         lick_idx = find(lick_cuereward);
%         if length(lick_idx) >= 4
%             inter_lick = diff(lick_idx);
%             %i
%             for li = 1:length(inter_lick)
%                 % li
%                 if inter_lick(li) <= 6 && length(inter_lick(li:end)) >= 3 && mean(inter_lick(li: li+2 )) <=6 % each lick is less than 6 frames apart
%                     % and at least three inter licks in bewteen
%                     lick.lickTrial = [lick.lickTrial 1];
%                     break;
%                 end
%                 
%                 if length(inter_lick(li:end)) < 3  % if less than 3 licks in the end, count as no lick trial
%                     lick.lickTrial = [lick.lickTrial 0];
%                     break;
%                 end
%                 
%             end
%         else
%             lick.lickTrial = [lick.lickTrial 0];
%         end
%     else
%         lick.lickTrial = [lick.lickTrial 0];
%     end
    
    % find trial without lick bout for 1000ms case
    if sum(trigger_licks(event_num, pre_frames - floor(100/frame_info.ifi):pre_frames + floor(1600/frame_info.ifi))) > 1 % no more than 1 lick 100ms before
        % cue and 100ms before reward
        lick_cuereward = trigger_licks(event_num, pre_frames - floor(100/frame_info.ifi):pre_frames + floor(1600/frame_info.ifi));
        lick_idx = find(lick_cuereward);
        if length(lick_idx) >= 4
            inter_lick = diff(lick_idx);
            %i
            for li = 1:length(inter_lick)
                % li
                if inter_lick(li) <= 6 && length(inter_lick(li:end)) >= 3 && mean(inter_lick(li: li+2 )) <=6 % each lick is less than 6 frames apart
                    %and at least three inter licks in bewteen
                    lick.lickTrial_1000 = [lick.lickTrial_1000 1];
                    break;
                end
                
                if length(inter_lick(li:end)) < 3  % if less than 3 licks in the end, count as no lick trial
                    lick.lickTrial_1000 = [lick.lickTrial_1000 0];
                    break;
                end
                
            end
        else
            lick.lickTrial_1000 = [lick.lickTrial_1000 0];
        end
    else
        lick.lickTrial_1000 = [lick.lickTrial_1000 0];
    end
    
    %gather licks by 10ms bins
    frame_time = frame_info.times(frame_no);
    time_range = (frame_time-pre_frames*10):10:(frame_time+post_frames*10);
    licks_this_TC = [];
    for kk = 2:length(time_range); %find all the # of licks which occur in each time bin.
        licks_this_bin = lickTimes(find(lickTimes>=time_range(kk-1) & lickTimes<time_range(kk)));
        licks_this_TC  = [licks_this_TC, licks_this_bin];
    end
    trigger_licks_10ms(event_num,:);
end