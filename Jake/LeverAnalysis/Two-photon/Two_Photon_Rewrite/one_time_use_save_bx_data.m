% this script was use in Oct 2017 because in the main analysis loop the
% licking_data was no longer being saved. So I reran the bx analysis and
% resaved all the related variables. 
% May have caused a 500ms misalignment between licking data and the F TCs??
clear;
file_info_CRP; clear days89 days90 days91 days92 days93 days94 days_1 days_post days_UR days_1000 days_1000_post comp_500ms irun;
behav_dir = 'Z:\Data\2P_imaging\behavior\';
for sub = [9, 22, 23, 31]%[3, 6, 10, 25, 27, 32]
    for rID = 1
        if sub ==9 
            rID=2;
        end
        sub
        out_dir  = fullfile('Z:', 'Analysis','Cue_reward_pairing_analysis','2P',[dates{sub}, '_', runID{rID}, '_', mouseID{sub}],'\');
        dest = out_dir;
        mouse = mouseID{sub};
        this_date = dates{sub};
        
        b_data = get_bx_data(behav_dir, [dates{sub}, '_', mouseID{sub}]);
        
        %% 1. find frame and lever times
        if exist('frame_times') == 1
            ifi = (frame_times(end)-frame_times(1))/length(frame_times);
        else
            ifi = mode(diff(cell2mat(cellfun(@int64,b_data.counterTimesUs,'UniformOutput',0))))/1000;
        end
        
        Sampeling_rate = 1000/ifi;
        
        % ---- parse behavior
        holdT_min  = 500000;  %us
        [lever, frame_info, trial_outcome, lick_data] = cleanBehav(b_data, ifi, holdT_min);
        frame_info.ifi = ifi;
        
        data_dest = [dest 'parse_behavior.mat'];
        save(data_dest, 'lever', 'frame_info', 'trial_outcome', 'lick_data', 'Sampeling_rate', 'holdT_min', 'ifi');
    end
end