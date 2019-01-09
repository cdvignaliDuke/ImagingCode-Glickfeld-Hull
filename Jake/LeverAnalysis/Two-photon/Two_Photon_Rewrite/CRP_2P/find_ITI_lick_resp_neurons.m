% 
% 1) identify ITI frames from original 6 CRP datasets
%      -get trials start/end frames and cue/rew frames
%      -make sure and log an index
% 2) load those ITI frames from .sbx file
% 3) register these frames to the img_ref from the trial frame analysis
% 4) apply mask from main analysis to get TCs
% 5) Take a df/f (use baseline_times df/f values?)
% 6) extract movies cued to single licks and lick bouts
% 7) Determine which neurons are responsive



clear
%file_info_CRP2;
file_info_CRP_all;
usFacs = 100;
behav_dir = 'Z:\Data\2P_imaging\behavior\';
crp_dir = 'Z:\Analysis\Cue_reward_pairing_analysis\2P\PCA_ICA_recalibration';
lick_dir_base = 'Z:\Analysis\Cue_reward_pairing_analysis\2P\PCA_ICA_recalibration\ITI_licking\'
doRerun = 0; useWarp = 0; doPCA = 0; useGPU = 0; checkImg = 0;
for session_subset = [1:2]
    if session_subset ==1
        stable_int = stable_int_1;
        all_nPCA = all_nPCA_1;
        all_PCuse = all_PCuse_1;
        all_PCuse2 = all_PCuse2_1;
        days_subset = days_1;
        nPCA_to_use = nPCA_to_use_1;
    elseif session_subset ==2
        stable_int = stable_int_post;
        all_nPCA = all_nPCA_post;
        all_PCuse = all_PCuse_post;
        all_PCuse2 = all_PCuse2_post;
        days_subset = days_post;
        nPCA_to_use = nPCA_to_use_post;
    elseif session_subset ==3
        stable_int = stable_int_UR;
        all_nPCA = all_nPCA_UR;
        all_PCuse = all_PCuse_UR;
        all_PCuse2 = all_PCuse2_UR;
        days_subset = days_UR;
        nPCA_to_use = nPCA_to_use_UR;
    end
    
    for sub = [1:6]%
        %% collect session/mouse information
        session = days_subset{sub};
        session_date = days_subset{sub}(1:6);
        if session(end-2) == 'g'
            mouse_num = ['9', session(end-1:end)];
            mouse_ID = ['img', session(end-1:end)];
        elseif session(end-2) =='0'
            mouse_num = session(end-2:end);
            mouse_ID = ['img', session(end-2:end)];
        end
        session
        
        %set pathnames
        data_dir = fullfile('Z:\Data\2P_imaging',session, mouse_ID, '\');
        for rID = 1:2
            if exist([data_dir, mouse_ID, '_000_', runID{rID}, '.sbx'], 'file') == 2
                break
            end
        end
        if strcmp(session, '170426_img91') | strcmp(session, '170522_img89')
            rID=2;
        end
        out_dir =  fullfile(crp_dir, [session_date, '_', runID{rID}, '_', mouse_ID], '\');
        data_dir = fullfile('Z:\Data\2P_imaging',session, mouse_ID, '\');
        lick_out_dir = fullfile(lick_dir_base, session, '\');
        if ~exist(lick_out_dir)
            mkdir(lick_out_dir);
        end
        
        %% 1) Load parse_behavior data and determine which frames to load
        %         load([out_dir, 'parse_behavior.mat']);
        %         load([out_dir, 'img_reg.mat'], 'img_ref', 'img_nframes');
        %         subMat = dir([behav_dir, '*', mouse_num, '-', session_date, '*']);
        %         load([behav_dir, subMat.name]);
        %
        %         % Identify ITI frame numbers
        %         ifi = mode(diff(frame_info.times));
        %         buffer_frames = ceil(500/ifi);
        %         postRewardFrames = ceil(input.postRewardMs/ifi);
        %         fixedReqHoldTimeFrames = ceil(input.fixedReqHoldTimeMs/ifi);
        %         RewardDelayDurationFrames = ceil(input.RewardDelayDurationMs/ifi);
        %         cue_frames = frame_info.counter(lever.press);
        %         trial_frames = [];
        %         for trial_num = 1:length(lever.press)
        %             trial_frames = [trial_frames, cue_frames(trial_num)-buffer_frames : cue_frames(trial_num)+postRewardFrames+buffer_frames];
        %         end
        %         ITI_frames_ind = ones(1,img_nframes);
        %         ITI_frames_ind([trial_frames]) = 0;
        %         ITI_frame_nums = find(ITI_frames_ind);
        %         assert(size(ITI_frame_nums,2) == size(unique(ITI_frame_nums),2));
        
        %%  Motion registration or check for exisint motion reg outputs
        
        %         % 2)  load sbx file
        %         config_fn = dir(fullfile(data_dir,['*' runID{rID} '.mat']));
        %         [img, ~, ~] = loadFile(data_dir, runID{rID}, [], []);
        %         assert(img_nframes == size(img,3));
        %
        %         % 3)  motion register
        %         [reg_out, img_reg] = stackRegister(img(:,:,ITI_frame_nums), img_ref);
        %         figure; plot(squeeze(mean(mean(img_reg,2),1)));
        %         title([session, ' mean raw F for all motion registered iti frames']);
        %         save([lick_out_dir, 'ITI_img_reg.mat'], 'img_reg', 'reg_out', 'img_ref', 'img_nframes', 'ITI_frame_nums', '-v7.3');
        
        %% 4) apply existing pixel masks to registered movies and get TCs
        
        %load variables
%         nPCA = nPCA_to_use{sub};
%         load([lick_out_dir, 'ITI_img_reg.mat']);
%         load([out_dir, session, '_nPCA_', num2str(nPCA), '_post_process_outputs.mat']);
%         
%         %get TCs
%         nmask = size(mask_3D_buffed,3);
%         tc_avg = getTC(img_reg, mask_3D_buffed, nmask);
%         save([lick_out_dir, session_date, '_', mouse_ID, '_nPCA_', num2str(nPCA), '_Results.mat'], 'tc_avg');

        %% 5) Take a df/f of time courses
        
%         %load variables
%         nPCA = nPCA_to_use{sub};
%         load([lick_out_dir, 'ITI_img_reg.mat'], 'ITI_frame_nums');
%         load([lick_out_dir, session_date, '_', mouse_ID, '_nPCA_', num2str(nPCA), '_Results.mat']);
%         
%         %use diff in frame nums to find trial borders. 
%         assert(length(unique(diff(ITI_frame_nums))) == 2);
%         border_inds = find(diff(ITI_frame_nums) > 1);
%         ITI_start_fnum = [1, (border_inds+1)]; %indeces of tc_avg which correspond to the beginning and end of the ITI
%         ITI_end_fnum = [border_inds, length(ITI_frame_nums)];
%         
%         %take an average f from the entire ITI period. 
%         num_ITIs = length(ITI_start_fnum);
%         num_cells = size(tc_avg,2);
%         avg_f = NaN(num_cells, num_ITIs);
%         for iti_num = 1:num_ITIs
%             %ITI_std = std(tc_avg([ITI_start_fnum(iti_num):ITI_end_fnum(iti_num)],:),1);
%             %ITI_mean = mean(tc_avg([ITI_start_fnum(iti_num):ITI_end_fnum(iti_num)],:),1);
%             avg_f(:,iti_num) = mean(tc_avg([ITI_start_fnum(iti_num):ITI_end_fnum(iti_num)],:),1);
%         end
%         assert(isempty(find(isnan(avg_f))));
%         
%         %calculate df/f for each frame for each cell
%         dfoverf_tc = NaN(size(tc_avg));
%         avg_f = avg_f';
%         for iti_num = 1:num_ITIs
%             int_length = length([ITI_start_fnum(iti_num):ITI_end_fnum(iti_num)]);
%             df_int = tc_avg([ITI_start_fnum(iti_num):ITI_end_fnum(iti_num)],:) - repmat(avg_f(iti_num,:),int_length,1);   % 1)frames,cells   2) 1,cells ->  frames,cells
%             dfoverf_tc([ITI_start_fnum(iti_num):ITI_end_fnum(iti_num)],:) = df_int./repmat(avg_f(iti_num,:),int_length,1);
%         end
%         assert(isempty(find(isnan(dfoverf_tc))));
%         save([lick_out_dir, 'dfoverf_tc'], 'dfoverf_tc', 'border_inds', 'ITI_start_fnum', 'ITI_end_fnum');
%         
        %% 6) extract movies cued to single licks and lick bouts
        
        %load and define some variables
        load([out_dir, 'parse_behavior.mat']);
        load([lick_out_dir, 'dfoverf_tc']);
        load([lick_out_dir, 'ITI_img_reg.mat'], 'ITI_frame_nums');
        ITI_frame_nums = ITI_frame_nums(ITI_frame_nums < length(lick_data.licksByFrame));
        licksByFrame = lick_data.licksByFrame([ITI_frame_nums]); %trim licksByFrame so that it has only ITI frames
        lickFrames = find(licksByFrame);
        num_ITIs = length(ITI_start_fnum);
        num_cells = size(dfoverf_tc,2);
        ifi = frame_info.ifi;
        pre_lick_fr = round(750/ifi);
        post_lick_fr = round(750/ifi);
        lick_buffer = round(500/ifi);
        lick_min_for_bout = 4;
        inter_lick_int = round(200/ifi);
        lick_window = (lick_min_for_bout-1)*inter_lick_int;
        
        %Look through iti frame nums of licksByFrame to find ITI licks and bouts               %make sure licks are buffer_distance away from border_frames
        single_lick_frame = [];
        bout_frame = [];
        for iti_num = 1:num_ITIs
            %find licks during this ITI which far enough away from the beginning or end of ITI to have buffers
            this_ITI_licks = find(lickFrames > ITI_start_fnum(iti_num)+pre_lick_fr &  lickFrames < ITI_end_fnum(iti_num)-post_lick_fr);
            %this_ITI_licks = lickFrames(this_ITI_licks);
            
            %find single licks with no other licks within 500ms
            for this_lick = this_ITI_licks
                if this_lick==1 
                    continue
                end
                if lickFrames(this_lick)-lickFrames(this_lick-1) > lick_buffer
                    if  this_lick == length(lickFrames) | lickFrames(this_lick+1)-lickFrames(this_lick) > lick_buffer    
                        single_lick_frame = [single_lick_frame, lickFrames(this_lick)];
                    end
                end
            end
            
            %find lick bouts   (4 licks within 200ms of each other)
            for this_lick = this_ITI_licks
                if this_lick==1 
                    continue
                end
                if lickFrames(this_lick)-lickFrames(this_lick-1) > lick_buffer
                    licks_this_window = licksByFrame(lickFrames(this_lick): lickFrames(this_lick)+lick_window);
                    if  sum(licks_this_window) >= lick_min_for_bout  %there are a total of four licks (including this lick) within the following 600ms
                        bout_frame = [bout_frame, lickFrames(this_lick)];
                    end
                end
            end
        end
        
        %save licking info
        ITI_lick_info.pre_lick_fr = pre_lick_fr;
        ITI_lick_info.post_lick_fr = post_lick_fr;
        ITI_lick_info.lick_buffer = lick_buffer;
        ITI_lick_info.single_lick_frame = single_lick_frame;
        ITI_lick_info.bout_frame = bout_frame;
        ITI_lick_info.lick_min_for_bout = lick_min_for_bout;
        ITI_lick_info.inter_lick_int = inter_lick_int;
        %save([lick_out_dir, 'ITI_lick_info.mat'], 'ITI_lick_info');

        %% 7) Ectract Ca and licking traces. Determine significance
        
        %single licks - extract licking and df/f traces using frame nums
        if length(single_lick_frame) > 9
            single_lick_ca_traces = NaN(length([-pre_lick_fr:post_lick_fr]), size(dfoverf_tc,2), length(single_lick_frame));  % dim1=frames   dim2=cells   dim3=lick number
            single_lick_lick_traces = NaN(length([-pre_lick_fr:post_lick_fr]), length(single_lick_frame));  % dim1=frames    dim2=lick number
            for this_lick = 1:length(single_lick_frame)
                single_lick_ca_traces(:,:,this_lick) = dfoverf_tc([single_lick_frame(this_lick)-pre_lick_fr : single_lick_frame(this_lick)+post_lick_fr],:);
                single_lick_lick_traces(:,this_lick) = licksByFrame([single_lick_frame(this_lick)-pre_lick_fr : single_lick_frame(this_lick)+post_lick_fr]);
            end
            assert(isempty(find(isnan(single_lick_ca_traces))));
            assert(isempty(find(isnan(single_lick_lick_traces))));
            
            %ttest for significance
            single_lick_ca_traces = permute(single_lick_ca_traces, [3,2,1]);
            [ITI_lick_h, ITI_lick_p, ITI_lick_resp_cells, ITI_lick_resp_avg, ITI_lick_resp_sem, ITI_lick_base, ITI_lick_resp] = findRespCell_ITI(single_lick_ca_traces, pre_lick_fr, ifi);   %dim1=trial# dim2=cell#  dim3=frames
        else
            ITI_lick_h=[]; ITI_lick_p=[]; ITI_lick_resp_cells=[NaN]; ITI_lick_resp_avg=[]; ITI_lick_resp_sem=[];  ITI_lick_base=[]; ITI_lick_resp=[];
        end
        
        %lick bouts - extract licking and df/f traces using frame nums
        if length(bout_frame) > 9
            bout_ca_traces = NaN(length([-pre_lick_fr:post_lick_fr]), size(dfoverf_tc,2), length(bout_frame));  % dim1=frames   dim2=cells   dim3=lick number
            bout_lick_traces = NaN(length([-pre_lick_fr:post_lick_fr]), length(bout_frame));  % dim1=frames     dim2=lick number
            for this_lick = 1:length(bout_frame)
                bout_ca_traces(:,:,this_lick) = dfoverf_tc([bout_frame(this_lick)-pre_lick_fr : bout_frame(this_lick)+post_lick_fr],:);
                bout_lick_traces(:,this_lick) = licksByFrame([bout_frame(this_lick)-pre_lick_fr : bout_frame(this_lick)+post_lick_fr]);
            end
            assert(isempty(find(isnan(bout_ca_traces))));
            assert(isempty(find(isnan(bout_lick_traces))));
            
            %ttest for significance
            bout_ca_traces = permute(bout_ca_traces, [3,2,1]);
            [ITI_bout_h, ITI_bout_p, ITI_bout_resp_cells, ITI_bout_resp_avg, ITI_bout_resp_sem, ITI_bout_base, ITI_bout_resp] = findRespCell_ITI(bout_ca_traces, pre_lick_fr, ifi);
        else
            ITI_bout_h=[]; ITI_bout_p=[]; ITI_bout_resp_cells=[NaN]; ITI_bout_resp_avg=[]; ITI_bout_resp_sem=[]; ITI_bout_base=[]; ITI_bout_resp=[];
        end
        
        %save lick responsive neurons in format useable in get_TC_events etc
        save([lick_out_dir, 'lick_resp_cells.mat'], 'ITI_lick_h', 'ITI_lick_p', 'ITI_lick_resp_cells', 'ITI_lick_resp_avg', 'ITI_lick_resp_sem', 'ITI_lick_base', 'ITI_lick_resp', ...
            'ITI_bout_h', 'ITI_bout_p', 'ITI_bout_resp_cells', 'ITI_bout_resp_avg', 'ITI_bout_resp_sem', 'ITI_bout_base', 'ITI_bout_resp');
        
        %plot mean lick trace and mean lick aligned ca trace 
        x_axis = [-pre_lick_fr:post_lick_fr]*ifi;
        color_cell = { 'r',  'g',  'b', 'k', 'm', 'c'};
        color_cell = repmat(color_cell,1,50);
        shift_window = [pre_lick_fr+1-round(750/33):pre_lick_fr+1-round(450/33)];
        figure; 
        if length(single_lick_frame) > 9
            %SINGLE LICK RESP
            subplot(2,2,1);
            bar( x_axis, squeeze(mean(single_lick_lick_traces,2))' ); hold on;
            errorbar(x_axis, squeeze(mean(single_lick_lick_traces,2))', std(single_lick_lick_traces,[],2)'/sqrt(size(single_lick_lick_traces,1)), 'g', 'LineStyle', 'none');
            for cell_num= ITI_lick_resp_cells
                shift_mag = mean( squeeze(mean(single_lick_ca_traces(:,cell_num,[shift_window]),1)) );
                plot(x_axis, squeeze(mean(single_lick_ca_traces(:,cell_num,:),1))-shift_mag, color_cell{cell_num});
                %errorbar(x_axis, squeeze(mean(single_lick_ca_traces(:,cell_num,:),1)),   squeeze(std(single_lick_ca_traces(:,cell_num,:),[],1))/sqrt(size(single_lick_ca_traces,1)),   'k'); % dim1=lick number  dim2=cells   dim3=frames
            end
            xlabel('time from lick (ms)');
            ylabel('lick rate or df/f');
            title(['Responsive: single lick: n=', num2str(size(single_lick_ca_traces,1)), ' licks. n=', num2str(length(ITI_lick_resp_cells)),' cells']);
            ylim([-0.1 0.15])
            
            %SINGLE LICK NON-RESP
            subplot(2,2,2);
            bar( x_axis, squeeze(mean(single_lick_lick_traces,2))' ); hold on;
            errorbar(x_axis, squeeze(mean(single_lick_lick_traces,2))', std(single_lick_lick_traces,[],2)'/sqrt(size(single_lick_lick_traces,1)), 'b', 'LineStyle', 'none');
            ITI_lick_non_resp_cells = 1:size(single_lick_ca_traces,2);
            ITI_lick_non_resp_cells([ITI_lick_resp_cells])=[];
            for cell_num= ITI_lick_non_resp_cells
                shift_mag = mean( squeeze(mean(single_lick_ca_traces(:,cell_num,[shift_window]),1)) );
                plot(x_axis, squeeze(mean(single_lick_ca_traces(:,cell_num,:),1))-shift_mag, color_cell{cell_num});
                %errorbar(x_axis, squeeze(mean(single_lick_ca_traces(:,cell_num,:),1)),   squeeze(std(single_lick_ca_traces(:,cell_num,:),[],1))/sqrt(size(single_lick_ca_traces,1)),   'k'); % dim1=lick number  dim2=cells   dim3=frames
            end
            xlabel('time from lick (ms)');
            ylabel('lick rate or df/f');
            title(['Non-responsive: single lick:  n=', num2str(size(single_lick_ca_traces,1)), ' licks. n=', num2str(length(ITI_lick_non_resp_cells)),' cells']);
            ylim([-0.1 0.15])
        end
        
        %BOUT RESP
        if length(bout_frame) > 9
            subplot(2,2,3);
            bar(x_axis, squeeze(mean(bout_lick_traces,2))); hold on;
            errorbar(x_axis, squeeze(mean(bout_lick_traces,2)), std(bout_lick_traces,[],2)'/sqrt(size(bout_lick_traces,2)), 'b', 'LineStyle', 'none');
            for cell_num=ITI_bout_resp_cells;
                shift_mag = mean( squeeze(mean(bout_ca_traces(:,cell_num,[shift_window]),1)) );
                plot(x_axis, squeeze(mean(bout_ca_traces(:,cell_num,:),1))-shift_mag, color_cell{cell_num});
                %errorbar(x_axis, squeeze(mean(mean(bout_ca_traces,2),1))', squeeze(std(mean(bout_ca_traces,1),[],2))'/sqrt(size(bout_ca_traces,2)),'g');
            end
            xlabel('time from lick (ms)');
            ylabel('lick rate or df/f');
            title(['Responsive: lick bout: n=', num2str(size(bout_ca_traces,1)), ' licks. n=', num2str(length(ITI_bout_resp_cells)),' cells']);
            ylim([-0.1 0.15]);
            
            %BOUT NON-RESP
            subplot(2,2,4);
            bar(x_axis, squeeze(mean(bout_lick_traces,2))); hold on;
            errorbar(x_axis, squeeze(mean(bout_lick_traces,2)), std(bout_lick_traces,[],2)'/sqrt(size(bout_lick_traces,2)), 'b', 'LineStyle', 'none');
            ITI_bout_non_resp_cells = 1:size(bout_ca_traces,2);
            ITI_bout_non_resp_cells([ITI_bout_resp_cells])=[];
            for cell_num=ITI_bout_non_resp_cells;
                shift_mag = mean( squeeze(mean(bout_ca_traces(:,cell_num,[shift_window]),1)) );
                plot(x_axis, squeeze(mean(bout_ca_traces(:,cell_num,:),1))-shift_mag, color_cell{cell_num});
                %errorbar(x_axis, squeeze(mean(mean(bout_ca_traces,2),1))', squeeze(std(mean(bout_ca_traces,1),[],2))'/sqrt(size(bout_ca_traces,2)),'g');
            end
            xlabel('time from lick (ms)');
            ylabel('lick rate or df/f');
            title(['Non-resp: lick bout: n=', num2str(size(bout_ca_traces,1)), ' licks. n=', num2str(length(ITI_bout_non_resp_cells)),' cells']);
            ylim([-0.1 0.15]);
        end
        suptitle([session_date, ' ',  mouse_ID])
        if length(bout_frame) > 9 | length(single_lick_frame) > 9
            savefig([lick_out_dir, 'lick aligned dfoverf mean traces.fig']);
        end
    end
end









