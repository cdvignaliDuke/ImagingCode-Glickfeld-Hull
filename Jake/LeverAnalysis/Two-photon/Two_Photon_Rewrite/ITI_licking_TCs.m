%script for using ITI lick bouts to ID licking related neurons and ID them
%will us dfoverf_TC

%this code was originally designed to isolate lick bouts in the ITI and
%plot df/f licking from those bouts. Was adapted to plot liick triggered
%avgs from the ITI

clear;
file_info_CRP; clear days89 days90 days91 days92 days93 days94 days_1 days_post days_UR days_1000 days_1000_post comp_500ms irun;
dir_base = 'Z:\Analysis\Cue_reward_pairing_analysis\2P\';
bx_dir_base = 'Z:\Data\2P_imaging\behavior\';
post_lick_frames = round(1500/33.33); %ms of df/f and licking to be plotted after the first lick of the lick bout
pre_lick_frames = round(1500/33.33); 
post_reward_buffer = 5000;
max_diff =6; %maximum number of difference in frame # between licks in a bout
min_licks_per_bout = 3;
pre_bout_int = 22;

for ii = [3, 6, 10, 25, 27, 32,   2, 5, 9, 22, 23, 31]  %     2, 5, 9, 22, 23, 31]%  % %post learning          []; day 1
    %assign an output directory for the results
    out_dir = [dir_base, 'ITI licking\', dates{ii}, '_', mouseID{ii}, '\'];
    if ~exist(out_dir);
        mkdir(out_dir);
    end
    
    %check to see if the file exists and select currect runID
    for rID = 1:3;
        data_dir = [dir_base, dates{ii}, '_', runID{rID}, '_', mouseID{ii}, '\'];
        TC_fn = dir([data_dir, 'Results.mat']);
        if ~isempty(TC_fn);
            break
        end
    end
    
    %load TCs  bx data and licking data
    load([data_dir, TC_fn.name], 'tc_avg');
    load([data_dir, 'parse_behavior'], 'lick_data', 'lever', 'frame_info');
    b_data = get_bx_data(bx_dir_base, [dates{ii}, '_', mouseID{ii}]);
    rew_delay = b_data.RewardDelayDurationMs;
    cue_on = frame_info.counter(lever.press); %converts cue on/off times to framesfa
    cue_off = frame_info.counter(lever.release);
    ifi = mean(diff(frame_info.times)); assert(ifi>32 && ifi<34);
    licks_by_frame = lick_data.licksByFrame;
    
    %generate index of frames which could be used to ITI licks
    valid_frame_inx = ones(1,length(tc_avg));
    assert( length(cue_on) == length(cue_off) );
    assert( cue_on(1) < cue_off(1) );
    for iii = 1:length(cue_on)
        invalid_frames = [ (cue_on(iii)-post_lick_frames) : cue_off(iii)+round((rew_delay+post_reward_buffer)/ifi) ];
        if invalid_frames(end) > length(valid_frame_inx);  %if the invalid frame inx exceeds matrix dims then just go to the end of the movie
            invalid_frames = [ (cue_on(iii)-post_lick_frames) : length(valid_frame_inx) ];
        end
        valid_frame_inx(invalid_frames) = 0;
        valid_frame_inx([1:pre_bout_int]) = 0;
        valid_frame_inx([end-(1+min_licks_per_bout*max_diff):end]) = 0;
    end
    
    %identify frames which correspond to the beginning of a lick bout
    bout_start_inx = [];
    if length(licks_by_frame) > length(valid_frame_inx)
        licks_by_frame = licks_by_frame(1:length(valid_frame_inx));
    end
    for iii = find(valid_frame_inx(1:length(licks_by_frame)));
        if licks_by_frame(iii)==0 %if there is no lick in this frame then skip to the next
            continue
        elseif ~isempty(bout_start_inx) && bout_start_inx(end) >= iii-max_diff; %prevents counting the second lick in a bout as a bout onset
            continue
        elseif sum( licks_by_frame([(iii-pre_bout_int):(iii-1)]) ) > 0 %if there is a lick in the time before the current frame
            continue
        %elseif sum( licks_by_frame( [(iii+1) : (iii+1+min_licks_per_bout*max_diff)] )) < min_licks_per_bout;  %minimum # of licks per bout 
            %continue
        %elseif max(diff(find( licks_by_frame([iii:(iii+1+min_licks_per_bout*max_diff)])  ))) > max_diff %maximum difference between lick frames in order to avoid confounds due one lick counting as multiple
            %continue 
        elseif sum( licks_by_frame([(iii+2):(iii+pre_bout_int)]) ) > 0 %if there is a lick in the time AFTER the current frame
            continue
        end
        bout_start_inx = [bout_start_inx, iii];
    end
    
    if isempty(bout_start_inx)
        disp(['no ITI lick bouts for ' dates{ii}, ' ', mouseID{ii}]);
        continue
    end
    
    %save bout start inx, parameters for determining valid bouts and valid frames
    %save([out_dir, 'ITI_lick_parameters_'], 'bout_start_inx', 'max_diff', 'min_licks_per_bout', 'post_lick_frames', 'pre_lick_frames', 'post_reward_buffer');
   
    %pull out TCs from each neuron aligned to the lick bout starts
    use_ev_ITI_bout = frame_info.times([bout_start_inx]) -frame_info.imaging_start_MW_T +10; %added the extra 10ms so that the reported ev times will be well within their own fram_num's duration instead of on the border. 
    [ITI_bout_movie, ~, lick_trace_ITI_bout, ~] = trigger_movie_by_event(tc_avg', frame_info, ...
        use_ev_ITI_bout, pre_lick_frames, post_lick_frames, lick_data);
    
    %calculate df/f for each trial and each neuron
    base_f_window = [pre_lick_frames-20:pre_lick_frames-10];
    base_f = mean( ITI_bout_movie(:,:,[base_f_window]) ,3);
    base_f = repmat(base_f,1,1,size(ITI_bout_movie,3));
    ITI_bout_movie_dfof = (ITI_bout_movie-base_f)./base_f; %dim1=trial# dim2=cell# dim3=frame#
    
    %determine which cells significantly respond to licking
    [LR_h, LR_p, LR_resp_cells, LR_resp_avg, LR_resp_sem, LR_base, LR_resp] = findRespCell2(ITI_bout_movie_dfof, pre_lick_frames, ifi);
    LR_non_resp_cells = find(LR_h==0);
    
    %plot avg df/f and licking for non responsive neurons
    x_axis = ([1:size(ITI_bout_movie_dfof,3)]-pre_lick_frames-1)*round(ifi);
%     figure; subplot(1,2,1); 
%     suptitle([' avg df/f and licking ' , dates{ii}, ' ', mouseID{ii}, ' ']);
%     bar( x_axis, mean(lick_trace_ITI_bout,1)*30/100 );  hold on;
%     plot( x_axis, squeeze(mean(mean(ITI_bout_movie_dfof(:,LR_non_resp_cells,:),1),2)) ,'g');
%     errorbar( x_axis, squeeze(mean(mean(ITI_bout_movie_dfof(:,LR_non_resp_cells,:),1),2)), std(squeeze(mean(ITI_bout_movie_dfof(:,LR_non_resp_cells,:),1)),[],1)/sqrt(size(ITI_bout_movie_dfof(:,LR_non_resp_cells,:),2)),'g');
%     title(['All cells:  ',  'lick bouts=', num2str(size(ITI_bout_movie_dfof,1)), ' neurons=', num2str(size(ITI_bout_movie_dfof(:,LR_non_resp_cells,:),2))]);
%     xlabel('time relative to bout start (ms)');
%     ylabel('df/f  and  lick rate (Hz/100)');
%     xlim([-1500 1500]); 
%     ylim([-0.02 max(squeeze(mean(mean(ITI_bout_movie_dfof,1),2)))*1.3]);
    
    %plot df/f and licking for lick responsive neurons
%     subplot(1,2,2); 
%     if ~isempty(LR_resp_cells)
%     bar( x_axis, mean(lick_trace_ITI_bout,1)*30/100 );  hold on;
%     plot( x_axis, squeeze(mean(mean(ITI_bout_movie_dfof(:,[LR_resp_cells],:),1),2)) ,'g');
%     errorbar( x_axis, squeeze(mean(mean(ITI_bout_movie_dfof(:,[LR_resp_cells],:),1),2)), std(squeeze(mean(ITI_bout_movie_dfof(:,[LR_resp_cells],:),1)),[],1)/sqrt(length(LR_resp_cells)),'g');
%     title(['Lick Responsive cells: ', 'lick bouts=', num2str(size(ITI_bout_movie_dfof,1)), ' neurons=', num2str(length(LR_resp_cells))]);
%     xlabel('time relative to bout start (ms)');
%     ylabel('df/f  and  lick rate (Hz/100)');
%     xlim([-1500 1500]); ylim([-0.02 max(squeeze(mean(mean(ITI_bout_movie_dfof(:,[LR_resp_cells],:),1),2)))*1.3]);
%     end
    %save figs and data
    %savefig([out_dir, 'avg_TCs_no_prelick_exclusion']);
    %save([out_dir, 'ITI_licking_outputs'], 'ITI_bout_movie', 'ITI_bout_movie_dfof', 'lick_trace_ITI_bout', 'LR_h', 'LR_p', 'LR_resp_cells', 'LR_resp_avg', 'LR_resp_sem', 'LR_base', 'LR_resp', 'bout_start_inx');
    
%     %plot all trials for a subset of neurons 
%     cell_subset = round(linspace(1,size(ITI_bout_movie_dfof,2),5));
%     for cell_num = 1:length(cell_subset)
%         shift = 0;
%         figure; hold on;
%         for bout_num = 1:size(ITI_bout_movie_dfof,1)
%             plot(x_axis, squeeze(ITI_bout_movie_dfof(bout_num, cell_num, :)) + shift);
%             licks_this_TC = lick_trace_ITI_bout(bout_num,:);
%             plot(x_axis(find(licks_this_TC)), repmat(shift,1,length(find(licks_this_TC))), 'c*');   %==================need to plot a dot at y value shift and x value where the lick occurred. 
%             shift = shift + 0.2;
%         end
%         xlabel('time relative to bout start (ms)');
%         ylabel('df/f  and  lick rate (Hz/100)');
%         title(['df/f and licking ' , dates{ii}, ' ', mouseID{ii}, ' ' 'cell number ' num2str(cell_subset(cell_num)), ' all ITI lick bouts(',num2str(bout_num), ') offset']);
%         savefig([out_dir, 'cell_', num2str(cell_num), 'all_ITI_bouts']);
%     end

    %plot all cells averaged across all licks for lick aligned TCs with no prelick exclusion criteria
    bout_num = size(ITI_bout_movie_dfof,1);
    shift=0;
    figure; 
    bar( x_axis, mean(lick_trace_ITI_bout,1)*30/100 );  hold on;
    for cell_num = 1:size(ITI_bout_movie_dfof,2)
        this_cell_TC = squeeze(mean(ITI_bout_movie_dfof(:, cell_num, :),1));
        plot(x_axis, this_cell_TC + shift, 'g');
        shift = shift - 0.03;
    end
    plot( x_axis, squeeze(mean(mean(ITI_bout_movie_dfof,1),2)), 'r');
    errorbar( x_axis, squeeze(mean(mean(ITI_bout_movie_dfof,1),2)), std(squeeze(mean(ITI_bout_movie_dfof,1)),[],1)/sqrt(size(ITI_bout_movie_dfof,2)),'r');
    xlabel('time from lick (ms)');
    ylabel('df/f and lick rate');
    title([dates{ii}, ' ', mouseID{ii}, ' ', 'df/f traces for each cell in a FoV avg across all ITI licks n=', num2str(bout_num), '. Red = grand mean']);
    savefig([out_dir, 'cell_', num2str(cell_num), 'isolated_ITI_licks_single_cells']);
end


%% Ashley and Lindsey ran some experiments on the 2P for frame/bx alignment. So it appears there is about a 1 
% frame delay between the command HAD_2P_frames sends to deliver the cue and when 
% the cue actually changes. I will need to correct for that. However, they
% seem to think that the frame pulse 1s before the 2nd pulse actually does
% correspond to an imaged frame. They said that the two extra pulses come
% from frames remaining in the buffer which were never saved. So the last
% two pulses do not have associated frames. 




