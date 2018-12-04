%script for plotting population summaries for the CRP data. 

clear
file_info_CRP_all;
out_base = fullfile('Z:\Analysis\Cue_reward_pairing_analysis\2P\PCA_ICA_recalibration\Summary_folder\');
crp_dir = 'Z:\Analysis\Cue_reward_pairing_analysis\2P\PCA_ICA_recalibration\';
ITI_dir = 'Z:\Analysis\Cue_reward_pairing_analysis\2P\PCA_ICA_recalibration\ITI_licking\';
%load([out_base 'cell_count.mat']);
normalR_nevents_D1 = 0; spont_nevents_D1 = 0;
sess_subset =  [1:15];  % [1:14]%1:6%  [1,2,3,4,5,6,7,8,9,10,11,13,15]; %1:size(days_1,2) %

%load and group variables for DAY 1
pre_buffer_all = [];
post_buffer_all = [];
for id = sess_subset
    %determine session info
    session = days_1{id};
    session_date = days_1{id}(1:6);
    if session(end-2) == 'g'
        mouse_num = ['9', session(end-1:end)];
        mouse_ID = ['img', session(end-1:end)];
    elseif session(end-2) =='0'
        mouse_num = session(end-2:end);
        mouse_ID = ['img', session(end-2:end)];
    end
    
    %set pathnames
    for rID = 1:2
        if exist([crp_dir, session_date, '_', runID{rID}, '_', mouse_ID], 'file') == 7
            break
        end
    end
    if strcmp(session, '170426_img91') | strcmp(session, '170522_img89')
        rID=2;
    end
    dest_sub  = fullfile(crp_dir, [session_date, '_', runID{rID}, '_', mouse_ID],'\');
    dest_sub_spikes  = fullfile(dest_sub, 'spike_outputs', '\');
    iti_lick_dir = fullfile(ITI_dir, session, '\');
    
    %load variables
    load([dest_sub '_cell_categories.mat']);
    load([dest_sub_spikes '_event_summary.mat']);
    load([dest_sub_spikes '_spont_events.mat']);
    load([dest_sub_spikes '_event_hist.mat']);
    load([dest_sub_spikes '_evoked_events.mat'], 'pre_buffer_cue', 'post_buffer_cue');
    load([dest_sub '_cue_movies.mat']);
    %load([iti_lick_dir, 'lick_resp_cells.mat']);
    %load([dest_sub '_cell_TCs.mat']);
    %             load([dest_sub 'parse_behavior']);
    %             load([dest_sub '_event_peaks.mat']);
    %             load([dest_sub_spikes '_norm2spont.mat']);
    %             load([out_base, 'cell_count.mat']);
    
    
    pre_buffer_all = [pre_buffer_all, pre_buffer_cue];
    post_buffer_all = [post_buffer_all, post_buffer_cue];
    ncells(id) = nCells;
    
    %Define lick responsive cells
    RS_cells_D1{id} = allresp_cells; %
%     if isnan(ITI_bout_resp_cells);   %    ITI_lick_resp_cells
%         RS_cells_D1{id} = [];
%     else 
%         lick_resp_cells = logical(zeros(1,nCells)); %  lick_resp_cells
%         lick_resp_cells([ITI_bout_resp_cells]) = 1;
%         non_lick_resp_cells = ~lick_resp_cells;
%         RS_cells_D1{id} = allresp_cells; %      |||    allresp_cells_pos  cue_cells_pos  rew_cells_pos  allresp_cells_neg  rew_cells_neg  cue_cells_neg
%     end
    
    %             TC_length(id) = size(all_success_event,2);
    %             TC_ifi(id) = ifi;
    
    %Group event hists for all NR/OR trials. all/RS. Aligned to cue. 
    normalR_hist_D1{id} = all_NR_hist.*(1000/double(ifi));
    omitR_hist_D1{id} = all_OR_hist.*(1000/double(ifi));
    normalR_hist_RS_D1{id} = all_NR_hist(:,RS_cells_D1{id}).*(1000/double(ifi));
    omitR_hist_RS_D1{id} = all_OR_hist(:,RS_cells_D1{id}).*(1000/double(ifi));
    
    %Group event hists for NO LICK NR/OR trials. all/RS
    normalR_hist_nolick_D1{id} = all_NR_nolick_hist.*(1000/double(ifi));
    omitR_hist_nolick_D1{id} = all_OR_nolick_hist.*(1000/double(ifi));
    normalR_hist_nolick_RS_D1{id} = all_NR_nolick_hist(:,RS_cells_D1{id}).*(1000/double(ifi));
    omitR_hist_nolick_RS_D1{id} = all_OR_nolick_hist(:,RS_cells_D1{id}).*(1000/double(ifi));
    
    NR_ind_D1{id} = (NR_ind_all-1).*double(ifi);
    OR_ind_D1{id} = (OR_ind_all-1).*double(ifi);
    NR_RS_ind_D1{id} = (NR_ind_RS-1).*double(ifi);
    OR_RS_ind_D1{id} = (OR_ind_RS-1).*double(ifi);
    
    NR_D1_plot{id} = ones(size(NR_ind_all));
    OR_D1_plot{id} = ones(size(OR_ind_all));
    NR_RS_D1_plot{id} = ones(size(NR_ind_RS));
    OR_RS_D1_plot{id} = ones(size(OR_ind_RS));
    
    %             normalR_norm_all_D1{id} = normalR_norm_avg;
    %             omitR_norm_all_D1{id} = omitR_norm_avg;
    %             spont_norm_all_D1{id} = spont_norm_avg;
    %             normalCue_norm_all_D1{id} = normalCue_norm_avg;
    %             omitCue_norm_all_D1{id} = omitCue_norm_avg;
    
    %             normalR_norm_RS_D1{id} = normalR_norm_avg(:,RS_cells_D1{id});
    %             omitR_norm_RS_D1{id} = omitR_norm_avg(:,RS_cells_D1{id});
    %             spont_norm_RS_D1{id} = spont_norm_avg(:,RS_cells_D1{id});
    %             normalCue_norm_RS_D1{id} = normalCue_norm_avg(:,RS_cells_D1{id});
    %             omitCue_norm_RS_D1{id} = omitCue_norm_avg(:,RS_cells_D1{id});
    %
    %             normalR_norm1_all_D1{id} = normalR_norm1_avg;
    %             omitR_norm1_all_D1{id} = omitR_norm1_avg;
    %             spont_norm1_all_D1{id} = spont_norm1_avg;
    %             normalCue_norm1_all_D1{id} = normalCue_norm1_avg;
    %             omitCue_norm1_all_D1{id} = omitCue_norm1_avg;
    
    %             normalR_norm1_RS_D1{id} = normalR_norm1_avg(:,RS_cells_D1{id});
    %             omitR_norm1_RS_D1{id} = omitR_norm1_avg(:,RS_cells_D1{id});
    %             spont_norm1_RS_D1{id} = spont_norm1_avg(:,RS_cells_D1{id});
    %             normalCue_norm1_RS_D1{id} = normalCue_norm1_avg(:,RS_cells_D1{id});
    %             omitCue_norm1_RS_D1{id} = omitCue_norm1_avg(:,RS_cells_D1{id});
    %
    rate_D1{id} = events_rate;
    
    normalR_nevents_D1 = normalR_nevents_D1 + sum(normalR_nevents);
    spont_nevents_D1 = spont_nevents_D1 + sum(spont_nevents);
    
    %             ts{id} = [-data_start:data_end].*double(ifi);
    %             data_start_frame{id} = data_start;
    %             data_end_frame{id} = data_end;
    %             sz1 = size(all_success_hist,1);
    %             th{id} = [1-(sz1/2):(sz1/2)].*double(ifi);
end
assert(min(pre_buffer_all)==max(pre_buffer_all));

normalCue_nevents_P1 = 0;  spont_nevents_P1 =0;
for id = sess_subset %  [1:14]%  [1,2,3,4,5,6,7,8,9,10,11,13,15]; % 1:size(days_post,2) %1:6 %
    mouse = days_post{id};

    %determine session info
    session = days_post{id};
    session_date = days_post{id}(1:6);
    if session(end-2) == 'g'
        mouse_num = ['9', session(end-1:end)];
        mouse_ID = ['img', session(end-1:end)];
    elseif session(end-2) =='0'
        mouse_num = session(end-2:end);
        mouse_ID = ['img', session(end-2:end)];
    end
    
    %set pathnames
    for rID = 1:2
        if exist([crp_dir, session_date, '_', runID{rID}, '_', mouse_ID], 'file') == 7
            break
        end
    end
    if strcmp(session, '170426_img91') | strcmp(session, '170522_img89')
        rID=2;
    end
    dest_sub  = fullfile(crp_dir, [session_date, '_', runID{rID}, '_', mouse_ID],'\');
    dest_sub_spikes  = fullfile(dest_sub, 'spike_outputs', '\');
    iti_lick_dir = fullfile(ITI_dir, session, '\');
    
    %load variables
    load([dest_sub '_cell_categories.mat']);
    load([dest_sub_spikes '_spont_events.mat']);
    load([dest_sub_spikes '_event_hist.mat']);
    load([dest_sub '_cue_movies.mat']);
    %load([iti_lick_dir, 'lick_resp_cells.mat']);
%     load([dest_sub '_norm2spont.mat']);
%     load([dest_sub '_cell_TCs.mat']);
    
    %Define responsive cells
     RS_cells_P1{id} = allresp_cells;
%     if isnan(ITI_bout_resp_cells);   %   ITI_lick_resp_cells
%         RS_cells_P1{id} = [];
%     else
%         lick_resp_cells = logical(zeros(1,nCells));
%         lick_resp_cells([ITI_bout_resp_cells]) = 1;  %  ITI_bout_resp_cells  lick_resp_cells
%         non_lick_resp_cells = ~lick_resp_cells;
%         RS_cells_P1{id} = allresp_cells;   %     |||   allresp_cells_pos   cue_cells_pos  rew_cells_pos  allresp_cells_neg  rew_cells_neg  cue_cells_neg    lick_resp_cells
%     end
    
    normalR_hist_P1{id} = all_NR_hist.*(1000/double(ifi));
    omitR_hist_P1{id} = all_OR_hist.*(1000/double(ifi));
    normalR_hist_RS_P1{id} = all_NR_hist(:,RS_cells_P1{id}).*(1000/double(ifi));
    omitR_hist_RS_P1{id} = all_OR_hist(:,RS_cells_P1{id}).*(1000/double(ifi));
    
    %             normalR_hist_rewPos_P1{id} = all_NR_hist(:,NR_Rew_resp_pos_cells).*(1000/double(ifi));
    %             normalR_hist_rewNeg_P1{id} = all_NR_hist(:,NR_Rew_resp_neg_cells).*(1000/double(ifi));
    %             omitR_hist_rewPos_P1{id} = all_OR_hist(:,OR_Rew_resp_pos_cells).*(1000/double(ifi));
    %             omitR_hist_rewNeg_P1{id} = all_OR_hist(:,OR_Rew_resp_neg_cells).*(1000/double(ifi));
    %             normalR_hist_ORNeg_P1{id} = all_NR_hist(:,OR_Rew_resp_neg_cells).*(1000/double(ifi));
    
    normalR_hist_nolick_P1{id} = all_NR_nolick_hist.*(1000/double(ifi));
    omitR_hist_nolick_P1{id} = all_OR_nolick_hist.*(1000/double(ifi));
    normalR_hist_nolick_RS_P1{id} = all_NR_nolick_hist(:,RS_cells_P1{id}).*(1000/double(ifi));
    omitR_hist_nolick_RS_P1{id} = all_OR_nolick_hist(:,RS_cells_P1{id}).*(1000/double(ifi));
    
    %             normalR_hist_nolick_rewPos_P1{id} = all_NR_nolick_hist(:,NR_Rew_resp_pos_cells).*(1000/double(ifi));
    %             normalR_hist_nolick_rewNeg_P1{id} = all_NR_nolick_hist(:,NR_Rew_resp_neg_cells).*(1000/double(ifi));
    %             omitR_hist_nolick_rewPos_P1{id} = all_OR_nolick_hist(:,OR_Rew_resp_pos_cells).*(1000/double(ifi));
    %             omitR_hist_nolick_rewNeg_P1{id} = all_OR_nolick_hist(:,OR_Rew_resp_neg_cells).*(1000/double(ifi));
    %             normalR_hist_nolick_ORNeg_P1{id} = all_NR_nolick_hist(:,OR_Rew_resp_neg_cells).*(1000/double(ifi));
    
    omitR_hist_nolick_NegRate_P1{id} = all_OR_hist_negR.*(1000/double(ifi));
    
    NR_ind_P1{id} = (NR_ind_all-1).*double(ifi);
    OR_ind_P1{id} = (OR_ind_all-1).*double(ifi);
    NR_RS_ind_P1{id} = (NR_ind_RS-1).*double(ifi);
    OR_RS_ind_P1{id} = (OR_ind_RS-1).*double(ifi);
    
    NR_P1_plot{id} = 2*ones(size(NR_ind_all));
    OR_P1_plot{id} = 2*ones(size(OR_ind_all));
    NR_RS_P1_plot{id} = 2*ones(size(NR_ind_RS));
    OR_RS_P1_plot{id} = 2*ones(size(OR_ind_RS));
    
    %             normalR_norm_all_P1{id} = normalR_norm_avg;
    %             omitR_norm_all_P1{id} = omitR_norm_avg;
    %             spont_norm_all_P1{id} = spont_norm_avg;
    %             normalCue_norm_all_P1{id} = normalCue_norm_avg;
    %             omitCue_norm_all_P1{id} = omitCue_norm_avg;
    %
    %             normalR_norm_RS_P1{id} = normalR_norm_avg(:,RS_cells_P1{id});
    %             omitR_norm_RS_P1{id} = omitR_norm_avg(:,RS_cells_P1{id});
    %             spont_norm_RS_P1{id} = spont_norm_avg(:,RS_cells_P1{id});
    %             normalCue_norm_RS_P1{id} = normalCue_norm_avg(:,RS_cells_P1{id});
    %             omitCue_norm_RS_P1{id} = omitCue_norm_avg(:,RS_cells_P1{id});
    %
    %             normalR_norm1_all_P1{id} = normalR_norm1_avg;
    %             omitR_norm1_all_P1{id} = omitR_norm1_avg;
    %             spont_norm1_all_P1{id} = spont_norm1_avg;
    %             normalCue_norm1_all_P1{id} = normalCue_norm1_avg;
    %             omitCue_norm1_all_P1{id} = omitCue_norm1_avg;
    %
    %             normalR_norm1_RS_P1{id} = normalR_norm1_avg(:,RS_cells_P1{id});
    %             omitR_norm1_RS_P1{id} = omitR_norm1_avg(:,RS_cells_P1{id});
    %             spont_norm1_RS_P1{id} = spont_norm1_avg(:,RS_cells_P1{id});
    %             normalCue_norm1_RS_P1{id} = normalCue_norm1_avg(:,RS_cells_P1{id});
    %             omitCue_norm1_RS_P1{id} = omitCue_norm1_avg(:,RS_cells_P1{id});
    
    rate_P1{id} = events_rate;
    
    %             normalCue_nevents_P1 = normalCue_nevents_P1 + sum(normalCue_nevents);
    %             spont_nevents_P1 = spont_nevents_P1 + sum(spont_nevents);
end

%Unexpected Reward Condition
for id = 1:size(days_UR,2)
    %determine session info
    mouse = days_UR{id};
    session = days_UR{id};
    session_date = days_UR{id}(1:6);
    if session(end-2) == 'g'
        mouse_num = ['9', session(end-1:end)];
        mouse_ID = ['img', session(end-1:end)];
    elseif session(end-2) =='0'
        mouse_num = session(end-2:end);
        mouse_ID = ['img', session(end-2:end)];
    end
    
    %set pathnames
    for rID = 1:2
        if exist([crp_dir, session_date, '_', runID{rID}, '_', mouse_ID], 'file') == 7
            break
        end
    end
    if strcmp(session, '170426_img91') | strcmp(session, '170522_img89')
        rID=2;
    end
    dest_sub  = fullfile(crp_dir, [session_date, '_', runID{rID}, '_', mouse_ID],'\');
    dest_sub_spikes  = fullfile(dest_sub, 'spike_outputs', '\');
    
    %load data for this session
    load([dest_sub '_cell_categories.mat']);
    load([dest_sub_spikes '_spont_events.mat']);
    load([dest_sub_spikes '_event_hist.mat']);
    load([dest_sub '_cue_movies.mat']);
    %load([dest_sub '_cell_TCs.mat']);
    %load([dest_sub '_spont_events.mat']);
    
    %fill cells with data
    RS_cells_PN1{id} = allresp_cells;   %   allresp_cells
    %all trials
    normalR_hist_PN1{id} = all_NR_hist.*(1000/double(ifi));
    unexpR_hist_PN1{id} = all_UR_hist.*(1000/double(ifi));
    normalR_hist_RS_PN1{id} = all_NR_hist(:,RS_cells_PN1{id}).*(1000/double(ifi));
    unexpR_hist_RS_PN1{id} = all_UR_hist(:,RS_cells_PN1{id}).*(1000/double(ifi));
    %no lick condition
    normalR_hist_nolick_PN1{id} = all_NR_nolick_hist.*(1000/double(ifi));
    unexpR_hist_nolick_PN1{id} = all_UR_nolick_hist.*(1000/double(ifi));
    normalR_hist_nolick_RS_PN1{id} = all_NR_nolick_hist(:,RS_cells_PN1{id}).*(1000/double(ifi));
    unexpR_hist_nolick_RS_PN1{id} = all_UR_nolick_hist(:,RS_cells_PN1{id}).*(1000/double(ifi));
    %calculate rate
    rate_PN1{id} = events_rate;
    
end

%1000ms delay
for id = 1:size(days_1000,2)
    %determine session info
    mouse = days_1000{id};
    session = days_1000{id};
    session_date = days_1000{id}(1:6);
    if session(end-2) == 'g'
        mouse_num = ['9', session(end-1:end)];
        mouse_ID = ['img', session(end-1:end)];
    elseif session(end-2) =='0'
        mouse_num = session(end-2:end);
        mouse_ID = ['img', session(end-2:end)];
    end
    
    %set pathnames
    for rID = 1:2
        if exist([crp_dir, session_date, '_', runID{rID}, '_', mouse_ID], 'file') == 7
            break
        end
    end
    dest_sub  = fullfile(crp_dir, [session_date, '_', runID{rID}, '_', mouse_ID],'\');
    dest_sub_spikes  = fullfile(dest_sub, 'spike_outputs', '\');
    
    
    %load data for this session
    load([dest_sub '_cell_categories.mat']);
    load([dest_sub_spikes '_spont_events.mat']);
    load([dest_sub_spikes '_event_hist.mat']);
    load([dest_sub '_cue_movies.mat']);
    %load([dest_sub '_cell_TCs.mat']);
    %load([dest_sub '_spont_events.mat']);
    
    %fill cells with data
    RS_cells_PN2{id} = allresp_cells;
    %all trials
    normalR_hist_PN2{id} = all_NR_hist.*(1000/double(ifi));
    omitR_hist_PN2{id} = all_OR_hist.*(1000/double(ifi));
    normalR_hist_RS_PN2{id} = all_NR_hist(:,RS_cells_PN2{id}).*(1000/double(ifi));
    omitR_hist_RS_PN2{id} = all_OR_hist(:,RS_cells_PN2{id}).*(1000/double(ifi));
    
    %             normalR_hist_rewPos_PN2{id} = all_NR_hist(:,NR_Rew_resp_pos_cells).*(1000/double(ifi));
    %             normalR_hist_rewNeg_PN2{id} = all_NR_hist(:,NR_Rew_resp_neg_cells).*(1000/double(ifi));
    %             omitR_hist_rewPos_PN2{id} = all_OR_hist(:,OR_Rew_resp_pos_cells).*(1000/double(ifi));
    %             omitR_hist_rewNeg_PN2{id} = all_OR_hist(:,OR_Rew_resp_neg_cells).*(1000/double(ifi));
    %             normalR_hist_ORNeg_PN2{id} = all_NR_hist(:,OR_Rew_resp_neg_cells).*(1000/double(ifi));
    %
    
    %no lick condition
    omitR_hist_nolick_NegRate_PN2{id} = all_OR_hist_negR.*(1000/double(ifi));
    normalR_hist_nolick_PN2{id} = all_NR_nolick_hist.*(1000/double(ifi));
    omitR_hist_nolick_PN2{id} = all_OR_nolick_hist.*(1000/double(ifi));
    normalR_hist_nolick_RS_PN2{id} = all_NR_nolick_hist(:,RS_cells_PN2{id}).*(1000/double(ifi));
    omitR_hist_nolick_RS_PN2{id} = all_OR_nolick_hist(:,RS_cells_PN2{id}).*(1000/double(ifi));
    
    %             normalR_hist_nolick_rewPos_PN2{id} = all_NR_nolick_hist(:,NR_Rew_resp_pos_cells).*(1000/double(ifi));
    %             normalR_hist_nolick_rewNeg_PN2{id} = all_NR_nolick_hist(:,NR_Rew_resp_neg_cells).*(1000/double(ifi));
    %             omitR_hist_nolick_rewPos_PN2{id} = all_OR_nolick_hist(:,OR_Rew_resp_pos_cells).*(1000/double(ifi));
    %             omitR_hist_nolick_rewNeg_PN2{id} = all_OR_nolick_hist(:,OR_Rew_resp_neg_cells).*(1000/double(ifi));
    %             normalR_hist_nolick_ORNeg_PN2{id} = all_NR_nolick_hist(:,OR_Rew_resp_neg_cells).*(1000/double(ifi));
    
    rate_PN2{id} = events_rate;
    cue_rew_cells{id} = (cue_cells & rew_cells);
    cue_cells_fr{id} = sum(cue_rew_cells{id})/sum(cue_cells);
    rew_cells_fr{id} = sum(cue_rew_cells{id})/sum(rew_cells);
    plot_ind{id} = 1;
end

%% plotting
% ['peach' 'army green' 'yellow' 'purple' 'black']
col_mat = [ 0.9  0.9  0;
    1  0  1;
    0  1  1;
    0.5  0  0;
    0  1  0;
    0  0  1;
    1  0.6  1;
    0  0  0;
    1  0.8 0.4
    0  0.5 0.7
    0.5 0.4 0];

col_mat_s = 0.5*ones(size(col_mat));
fig=figure;scatter_plot(days_1_mouse, plot_ind, cue_cells_fr, col_mat_s)
ylim([0,1]);
set(gca, 'XTick', 1, 'XTickLabels', 'Cue Responsive');
ylabel('fraction of cells responsive for reward and cue');
saveas(fig, [out_base 'Summary_100ms_rewardcue_cell_1.fig']);
print([out_base 'Summary_100ms_rewardcue_cell_1.eps'], '-depsc');
print([out_base 'Summary_100ms_rewardcue_cell_1.pdf'], '-dpdf');

fig=figure;scatter_plot(days_1_mouse, plot_ind, rew_cells_fr, col_mat_s)
ylim([0,1]);
set(gca, 'XTick', 1, 'XTickLabels', 'Reward Responsive');
ylabel('fraction of cells responsive for reward and cue');
saveas(fig, [out_base 'Summary_100ms_rewardcue_cell_2.fig']);
print([out_base 'Summary_100ms_rewardcue_cell_2.eps'], '-depsc');
print([out_base 'Summary_100ms_rewardcue_cell_2.pdf'], '-dpdf');

%summary of event rate
% fig=figure;
rate_all = [];
% cid = 1;
for id = 1:size(days_1_mouse,2)
%     if id > 1 && ~strcmp(mouseID{id},mouseID{id-1})
%             
%             cid = cid + 1;
%              
%     end
%         plot(id*ones(size(rate{id})), rate{id}, 'o', 'color', col_mat(cid,:));
    
    rate_all = [rate_all rate_D1{id}];
    hold on
end

for id = 1:size(days_1_mouse,2)
%     if id > 1 && ~strcmp(mouseID{id},mouseID{id-1})
%             
%             cid = cid + 1;
%              
%     end
%         plot(id*ones(size(rate{id})), rate{id}, 'o', 'color', col_mat(cid,:));
    
    rate_all = [rate_all rate_P1{id}];
    hold on
end
% ylim([0 2])
% xlim([0.5 size(mouseID,2)+.5])
% ylabel('Spike rate (Hz)')
% xlabel('session')
% title(['Avg rate: ' num2str(chop(mean(rate_all,2),3)) ' +/-' num2str(chop(std(rate_all,[],2)./(sqrt(size(rate_all,2))),2)) ' n = ' num2str(size(rate_all,2)) ' cells'])
% saveas(fig, [out_base 'Summary_event_rate.fig']);
% print([out_base 'Summary_event_rate.eps'], '-depsc');
% print([out_base 'Summary_event_rate.pdf'], '-dpdf');

fig = figure;
histogram(rate_all, 0:0.1:1.4)
xlabel('Firing rate (HZ)')
ylabel('#Cell')
title(['total cell=', num2str(length(rate_all)), ' mean firing rate is ', num2str(mean(rate_all)), ' and Standard error is ', num2str(std(rate_all)/sqrt(length(rate_all)))])
% title(['Avg rate: ' num2str(chop(mean(rate_all,2),3)) ' +/-' num2str(chop(std(rate_all,[],2)./(sqrt(size(rate_all,2))),2)) ' n = ' num2str(size(rate_all,2)) ' cells'])
saveas(fig, [out_base 'Summary_event_rate_hist.fig']);
print([out_base 'Summary_event_rate_hist.eps'], '-depsc');
print([out_base 'Summary_event_rate_hist.pdf'], '-dpdf');


%summary of average events
% fig=figure;
% for id = 1:size(mouseID,2)
%     press_events = press_event_TC{id};
%     success_events = success_event_TC{id};
%     fail_events = fail_event_TC{id};
%     cell_use = find(~isnan(mean([press_events success_events fail_events],2)));
%     subplot(4,5,id)
%     errorbar(ts{id}, nanmean(press_events(cell_use,:),1), std(press_events(cell_use,:),[],1)./sqrt(length(cell_use)), 'c');
%     hold on
%     errorbar(ts{id}, nanmean(success_events(cell_use,:),1), std(success_events(cell_use,:),[],1)./sqrt(length(cell_use)), 'k');
%     hold on
%     errorbar(ts{id}, nanmean(fail_events(cell_use,:),1), std(fail_events(cell_use,:),[],1)./sqrt(length(cell_use)), 'r');
%     title([date{id} ' ' mouseID{id} '- ' num2str(length(cell_use)) ' cells'])
% end
% xlabel('Time (ms)')
% ylabel('dF/F')
% supertitle(['Good events- All cells- Black: success; Red: early; Cyan: press'])
% saveas(fig, [out_base 'Summary_allcells_event_TC.fig']);
% print([out_base 'Summary_allcells_event_TC.eps'], '-depsc');
% print([out_base 'Summary_allcells_event_TC.pdf'], '-dpdf');
% fig=figure;
% for id = 1:size(mouseID,2)
%     press_events = press_event_TC_RS{id};
%     success_events = success_event_TC_RS{id};
%     fail_events = fail_event_TC_RS{id};
%     cell_use = find(~isnan(mean([press_events success_events fail_events],2)));
%     subplot(4,5,id)
%     errorbar(ts{id}, nanmean(press_events(cell_use,:),1), std(press_events(cell_use,:),[],1)./sqrt(length(cell_use)), 'c');
%     hold on
%     errorbar(ts{id}, nanmean(success_events(cell_use,:),1), std(success_events(cell_use,:),[],1)./sqrt(length(cell_use)), 'k');
%     hold on
%     errorbar(ts{id}, nanmean(fail_events,1), std(fail_events(cell_use,:),[],1)./sqrt(length(cell_use)), 'r');
%     title([date{id} ' ' mouseID{id} '- ' num2str(length(cell_use)) ' cells'])
% end
% xlabel('Time (ms)')
% ylabel('dF/F')
% supertitle(['Good events- Resp cells- Black: success; Red: early; Cyan: press'])
% saveas(fig, [out_base 'Summary_respcells_event_TC.fig']);
% print([out_base 'Summary_respcells_event_TC.eps'], '-depsc');
% print([out_base 'Summary_respcells_event_TC.pdf'], '-dpdf');


% %% commented to fix acquisition rates
% %average across days
% fig=figure;
% ra = [];
% pa = [];
% rr = [];
% pr = [];
% sa = [];
% fa = [];
% sr = [];
% fr = [];
% cid = 1;
% for id = 1:size(mouseID,2)
%     ra = [ra release_event_peak{id}];
%     pa = [pa press_event_peak{id}];
%     rr = [rr release_event_peak_RS{id}];
%     pr = [pr press_event_peak_RS{id}];
%     sa = [sa success_event_peak{id}];
%     fa = [fa fail_event_peak{id}];
%     sr = [sr success_event_peak_RS{id}];
%     fr = [fr fail_event_peak_RS{id}];
%     if id > 1 && ~strcmp(mouseID{id},mouseID{id-1})
%             
%             cid = cid + 1;
%              
%     end
%     subplot(2,2,1)
%     scatter(release_event_peak{id}, press_event_peak{id}, 'MarkerEdgeColor', col_mat(cid,:));
%     hold on
%     subplot(2,2,2)
%     scatter(success_event_peak{id}, fail_event_peak{id}, 'MarkerEdgeColor', col_mat(cid,:));
%     hold on
%     subplot(2,2,3)
%     scatter(release_event_peak_RS{id}, press_event_peak_RS{id}, 'MarkerEdgeColor', col_mat(cid,:));
%     hold on
%     subplot(2,2,4)
%     scatter(success_event_peak_RS{id}, fail_event_peak_RS{id}, 'MarkerEdgeColor', col_mat(cid,:));
%     hold on
%     
% end
% x = 0:.1:1;
% y = x;
% subplot(2,2,1)
% plot(x,y,'-k')
% xlim([0 1])
% ylim([0 1])
% xlabel('Release event peak dF/F')
% ylabel('Press event peak dF/F')
% [h_rpa p_rpa]= ttest(ra,pa);
% title(['All cells- n = ' num2str(total_cells) '; p = ' num2str(chop(p_rpa,2))])
% subplot(2,2,2)
% plot(x,y,'-k')
% xlim([0 1])
% ylim([0 1])
% xlabel('Success event peak dF/F')
% ylabel('Fail event peak dF/F') 
% [h_sfa p_sfa]= ttest(sa,fa);
% title(['All cells- n = ' num2str(total_cells) '; p = ' num2str(chop(p_sfa,2))])
% subplot(2,2,3)
% plot(x,y,'-k')
% xlim([0 1])
% ylim([0 1])
% xlabel('Release event peak dF/F')
% ylabel('Press event peak dF/F')
% [h_rpr p_rpr]= ttest(rr,pr);
% title(['Responsive cells- n = ' num2str(total_resp) '; p = ' num2str(chop(p_rpr,2))])
% subplot(2,2,4)
% plot(x,y,'-k')
% xlim([0 1])
% ylim([0 1])
% xlabel('Success event peak dF/F')
% ylabel('Fail event peak dF/F') 
% [h_sfr p_sfr]= ttest(sr,fr);
% title(['Responsive cells- n = ' num2str(total_resp) '; p = ' num2str(chop(p_sfr,2))])
% suptitle('Peak event amplitude')
% saveas(fig, [out_base 'Summary_event_amp_scatter.fig']);
% print([out_base 'Summary_event_amp_scatter.eps'], '-depsc');
% print([out_base 'Summary_event_amp_scatter.pdf'], '-dpdf');
% 
% %average event amplitude by mouse
% fig=figure;
% subplot(2,2,1)
% scatter_plot(mouseID, release_event_peak, press_event_peak, col_mat);
% 
% subplot(2,2,2)
% scatter_plot(mouseID, success_event_peak, fail_event_peak, col_mat);
% 
% subplot(2,2,3)
% scatter_plot(mouseID, release_event_peak_RS, press_event_peak_RS, col_mat);
% 
% subplot(2,2,4)
% scatter_plot(mouseID, success_event_peak_RS, fail_event_peak_RS, col_mat);

% for id = 1:size(mouseID,2)
%     subplot(2,2,1)
%     rp_peak = [release_event_peak{id}; press_event_peak{id}];
%     rp_use = find(~isnan(mean(rp_peak,1)));
%     errorbarxy(mean(rp_peak(1,rp_use),2), nanmean(rp_peak(2,rp_use),2), std(rp_peak(1,rp_use),[],2)./sqrt(length(rp_use)), std(rp_peak(2,rp_use),[],2)./sqrt(length(rp_use)),{['o' col_mat(id,:)], col_mat(id,:),col_mat(id,:)});
%     hold on
%     scatter(mean(rp_peak(1,rp_use),2), nanmean(rp_peak(2,rp_use),2),col_mat(id,:));
%     subplot(2,2,2)
%     sf_peak = [success_event_peak{id}; fail_event_peak{id}];
%     sf_use = find(~isnan(mean(sf_peak,1)));
%     errorbarxy(mean(sf_peak(1,sf_use),2), nanmean(sf_peak(2,sf_use),2),std(sf_peak(1,sf_use),[],2)./sqrt(length(sf_use)), std(sf_peak(2,sf_use),[],2)./sqrt(length(sf_use)),{['o' col_mat(id,:)], col_mat(id,:),col_mat(id,:)});
%     hold on
%     scatter(mean(sf_peak(1,sf_use),2), nanmean(sf_peak(2,sf_use),2),col_mat(id,:));
%     subplot(2,2,3)
%     rp_peak_RS = [release_event_peak_RS{id}; press_event_peak_RS{id}];
%     rp_use_RS = find(~isnan(mean(rp_peak_RS,1)));
%     errorbarxy(mean(rp_peak_RS(1,rp_use_RS),2), nanmean(rp_peak_RS(2,rp_use_RS),2), std(rp_peak_RS(1,rp_use_RS),[],2)./sqrt(length(rp_use_RS)), std(rp_peak_RS(2,rp_use_RS),[],2)./sqrt(length(rp_use_RS)),{['o' col_mat(id,:)], col_mat(id,:),col_mat(id,:)});
%     hold on
%     scatter(mean(rp_peak_RS(1,rp_use_RS),2), nanmean(rp_peak_RS(2,rp_use_RS),2), col_mat(id,:));
%     subplot(2,2,4)
%     sf_peak_RS = [success_event_peak_RS{id}; fail_event_peak_RS{id}];
%     sf_use_RS = find(~isnan(mean(sf_peak_RS,1)));
%     errorbarxy(mean(sf_peak_RS(1,sf_use_RS),2), nanmean(sf_peak_RS(2,sf_use_RS),2), std(sf_peak_RS(1,sf_use_RS),[],2)./sqrt(length(sf_use_RS)), std(sf_peak_RS(2,sf_use_RS),[],2)./sqrt(length(sf_use_RS)),{['o' col_mat(id,:)], col_mat(id,:),col_mat(id,:)});
%     hold on
%     scatter(mean(sf_peak_RS(1,sf_use_RS),2), nanmean(sf_peak_RS(2,sf_use_RS),2), col_mat(id,:));
% end
% x = 0:.1:1;
% y = x;
% subplot(2,2,1)
% plot(x,y,'-k')
% xlim([0 1])
% ylim([0 1])
% xlabel('Release event peak dF/F')
% ylabel('Press event peak dF/F')
% title(['All cells'])
% subplot(2,2,2)
% plot(x,y,'-k')
% xlim([0 1])
% ylim([0 1])
% xlabel('Success event peak dF/F')
% ylabel('Fail event peak dF/F') 
% title(['All cells'])
% subplot(2,2,3)
% plot(x,y,'-k')
% xlim([0 1])
% ylim([0 1])
% xlabel('Release event peak dF/F')
% ylabel('Press event peak dF/F')
% title(['Responsive cells'])
% subplot(2,2,4)
% plot(x,y,'-k')
% xlim([0 1])
% ylim([0 1])
% xlabel('Success event peak dF/F')
% ylabel('Fail event peak dF/F') 
% title(['Responsive cells'])
% supertitle('Peak event amplitude')
% saveas(fig, [out_base 'Summary_event_amp_avg_scatter.fig']);
% print([out_base 'Summary_event_amp_avg_scatter.eps'], '-depsc');
% print([out_base 'Summary_event_amp_avg_scatter.pdf'], '-dpdf');
% 
% %summary of PSTH- average
% fig=figure;
% for id = 1:size(mouseID,2)
%     subplot(4,5,id)
%     shadedErrorBar(th{id}, nanmean(success_hist{id},2), std(success_hist{id},[],2)./sqrt(size(success_hist{id},2)), 'k');
%     hold on;
%     shadedErrorBar(th{id}, nanmean(fail_hist{id},2), std(fail_hist{id},[],2)./sqrt(size(fail_hist{id},2)), 'r');
%     title([date{id} ' ' mouseID{id}])
% end
% xlabel('Time (ms)')
% ylabel('Firing rate (Hz)')
% supertitle(['PSTH- all cells- Black: success; Red: Failure'])
% saveas(fig, [out_base 'Summary_PSTH_allcells.fig']);
% print([out_base 'Summary_PSTH_allcells.eps'], '-depsc');
% print([out_base 'Summary_PSTH_allcells.pdf'], '-dpdf');

% fig=figure;
% for id = 1:size(mouseID,2)
%     subplot(4,5,id)
%     shadedErrorBar(th{id}, nanmean(success_hist_RS{id},2), std(success_hist_RS{id},[],2)./sqrt(size(success_hist_RS{id},2)), 'k');
%     hold on;
%     shadedErrorBar(th{id}, nanmean(fail_hist_RS{id},2), std(fail_hist_RS{id},[],2)./sqrt(size(fail_hist_RS{id},2)), 'r');
%     title([date{id} ' ' mouseID{id}])
% end
% xlabel('Time (ms)')
% ylabel('Firing rate (Hz)')
% supertitle(['PSTH- responsive cells- Black: success; Red: Failure'])
% saveas(fig, [out_base 'Summary_PSTH_respcells.fig']);
% print([out_base 'Summary_PSTH_respcells.eps'], '-depsc');
% print([out_base 'Summary_PSTH_respcells.pdf'], '-dpdf');


%% needs to be updated for different acquisition rates
% NR_hist_D1 = cellfun(@transpose, normalR_hist_D1, 'UniformOutput', 0);
% OR_hist_D1 = cellfun(@transpose, omitR_hist_D1, 'UniformOutput', 0);
% NR_hist_RS_D1 = cellfun(@transpose, normalR_hist_RS_D1, 'UniformOutput', 0);
% OR_hist_RS_D1 = cellfun(@transpose, omitR_hist_RS_D1, 'UniformOutput', 0);
% 
% NR_hist_P1 = cellfun(@transpose, normalR_hist_P1, 'UniformOutput', 0);
% OR_hist_P1 = cellfun(@transpose, omitR_hist_P1, 'UniformOutput', 0);
% NR_hist_RS_P1 = cellfun(@transpose, normalR_hist_RS_P1, 'UniformOutput', 0);
% OR_hist_RS_P1 = cellfun(@transpose, omitR_hist_RS_P1, 'UniformOutput', 0);
% frame_size = cell2mat(cellfun(@size, success_hist, 'UniformOutput', 0));
% max_frame = max(frame_size(2:2:end));
% 
% success_hist_all = interp_frame(success_hist, max_frame);
% fail_hist_all    = interp_frame(fail_hist, max_frame);
% success_hist_RS_all   = interp_frame(success_hist_RS, max_frame);12
% fail_hist_RS_all = interp_frame(fail_hist_RS, max_frame);

tth = [-pre_buffer_cue(1):post_buffer_cue].*double(min(ifi));
% tth = [1-(sz1/2):(sz1/2)].*double(min((TC_ifi)));
NR_hist_D1 = cell2mat(normalR_hist_D1)';
OR_hist_D1 = cell2mat(omitR_hist_D1)';
NR_hist_NL_D1 = cell2mat(normalR_hist_nolick_D1)';
OR_hist_NL_D1 = cell2mat(omitR_hist_nolick_D1)';
NR_hist_RS_D1 = cell2mat(normalR_hist_RS_D1)';
OR_hist_RS_D1 = cell2mat(omitR_hist_RS_D1)';
NR_hist_NL_RS_D1 = cell2mat(normalR_hist_nolick_RS_D1)';
OR_hist_NL_RS_D1 = cell2mat(omitR_hist_nolick_RS_D1)';

NR_hist_P1 = cell2mat(normalR_hist_P1)';
OR_hist_P1 = cell2mat(omitR_hist_P1)';
NR_hist_NL_P1 = cell2mat(normalR_hist_nolick_P1)';
OR_hist_NL_P1 = cell2mat(omitR_hist_nolick_P1)';
NR_hist_RS_P1 = cell2mat(normalR_hist_RS_P1)';
OR_hist_RS_P1 = cell2mat(omitR_hist_RS_P1)';
% normalR_hist_nolick_RS_P1{3} = [];
% omitR_hist_nolick_RS_P1{3}=[];
NR_hist_NL_RS_P1 = cell2mat(normalR_hist_nolick_RS_P1)';
OR_hist_NL_RS_P1 = cell2mat(omitR_hist_nolick_RS_P1)';
OR_hist_NL_NegCell_P1 = cell2mat(omitR_hist_nolick_NegRate_P1)';

% NR_hist_RewPos_P1 = cell2mat(normalR_hist_rewPos_P1)';
% NR_hist_RewNeg_P1 = cell2mat(normalR_hist_rewNeg_P1)';
% OR_hist_RewPos_P1 = cell2mat(omitR_hist_rewPos_P1)';
% OR_hist_RewNeg_P1 = cell2mat(omitR_hist_rewNeg_P1)';
% NR_hist_ORNeg_P1 = cell2mat(normalR_hist_ORNeg_P1)';
% 
% NR_hist_NL_RewPos_P1 = cell2mat(normalR_hist_nolick_rewPos_P1)';
% NR_hist_NL_RewNeg_P1 = cell2mat(normalR_hist_nolick_rewNeg_P1)';
% OR_hist_NL_RewPos_P1 = cell2mat(omitR_hist_nolick_rewPos_P1)';
% OR_hist_NL_RewNeg_P1 = cell2mat(omitR_hist_nolick_rewNeg_P1)';
% NR_hist_NL_ORNeg_P1 = cell2mat(normalR_hist_nolick_ORNeg_P1)';

NR_hist_PN1 = cell2mat(normalR_hist_PN1)';
UR_hist_PN1 = cell2mat(unexpR_hist_PN1)';
NR_hist_NL_PN1 = cell2mat(normalR_hist_nolick_PN1)';
UR_hist_NL_PN1 = cell2mat(unexpR_hist_nolick_PN1)';
NR_hist_RS_PN1 = cell2mat(normalR_hist_RS_PN1)';
UR_hist_RS_PN1 = cell2mat(unexpR_hist_RS_PN1)';
NR_hist_NL_RS_PN1 = cell2mat(normalR_hist_nolick_RS_PN1)';
UR_hist_NL_RS_PN1 = cell2mat(unexpR_hist_nolick_RS_PN1)';

NR_hist_PN2 = cell2mat(normalR_hist_PN2)';
OR_hist_PN2 = cell2mat(omitR_hist_PN2)';
NR_hist_NL_PN2 = cell2mat(normalR_hist_nolick_PN2)';
OR_hist_NL_PN2 = cell2mat(omitR_hist_nolick_PN2)';
NR_hist_RS_PN2 = cell2mat(normalR_hist_RS_PN2)';
OR_hist_RS_PN2 = cell2mat(omitR_hist_RS_PN2)';
% normalR_hist_nolick_RS_PN2{3}=[];
% omitR_hist_nolick_RS_PN2{3}=[];
NR_hist_NL_RS_PN2 = cell2mat(normalR_hist_nolick_RS_PN2)';
OR_hist_NL_RS_PN2 = cell2mat(omitR_hist_nolick_RS_PN2)';
OR_hist_NL_NegCell_PN2 = cell2mat(omitR_hist_nolick_NegRate_PN2)';

NR_hist_RewPos_PN2 = cell2mat(normalR_hist_rewPos_PN2)';
NR_hist_RewNeg_PN2 = cell2mat(normalR_hist_rewNeg_PN2)';
OR_hist_RewPos_PN2 = cell2mat(omitR_hist_rewPos_PN2)';
OR_hist_RewNeg_PN2 = cell2mat(omitR_hist_rewNeg_PN2)';
NR_hist_ORNeg_PN2 = cell2mat(normalR_hist_ORNeg_PN2)';

NR_hist_NL_RewPos_PN2 = cell2mat(normalR_hist_nolick_rewPos_PN2)';
NR_hist_NL_RewNeg_PN2 = cell2mat(normalR_hist_nolick_rewNeg_PN2)';
OR_hist_NL_RewPos_PN2 = cell2mat(omitR_hist_nolick_rewPos_PN2)';
OR_hist_NL_RewNeg_PN2 = cell2mat(omitR_hist_nolick_rewNeg_PN2)';
NR_hist_NL_ORNeg_PN2 = cell2mat(normalR_hist_nolick_ORNeg_PN2)';

%tth = [1-(size(NR_hist_D1,2)/2):(size(NR_hist_D1,2)/2)].*double(min((ifi)));

fig = figure;
subplot(1,2,1)
shadedErrorBar(tth, nanmean(OR_hist_NL_NegCell_P1,1), nanstd(OR_hist_NL_NegCell_P1,[],1)./sqrt(size(OR_hist_NL_NegCell_P1,1)), 'r')
xlim([-1500 1500])
ylim([0 3])
hold on; vline(600, '--k')
axis square
xlabel('Time (ms)')
ylabel('Firing rate (Hz)')
title(['Post Learning DayN+1 for Cells with Decreased Spike Rate, nCell = ', num2str(size(OR_hist_NL_NegCell_P1,1))]);

subplot(1,2,2)
shadedErrorBar(tth, nanmean(OR_hist_NL_NegCell_PN2,1), nanstd(OR_hist_NL_NegCell_PN2,[],1)./sqrt(size(OR_hist_NL_NegCell_PN2,1)), 'r')
xlim([-1500 1500])
ylim([0 3])
hold on; vline(1100, '--k')
axis square
xlabel('Time (ms)')
ylabel('Firing rate (Hz)')
title(['Post Learning DayN+3 for Cells with Decreased Spike Rate, nCell = ', num2str(size(OR_hist_NL_NegCell_PN2,1))]);
saveas(fig, [out_base 'Summary_PSTH_avgexpt_30Hz_OR.fig']);
print([out_base 'Summary_PSTH_avgexpt_30Hz_OR.eps'], '-depsc');
print([out_base 'Summary_PSTH_avgexpt_30Hz_OR.pdf'], '-dpdf');

% PSTHS ======================================================================================
fig=figure('rend', 'painters', 'pos', [50 150 1200 550]);
subplot(2,2,1)
shadedErrorBar(tth, nanmean(NR_hist_D1,1), nanstd(NR_hist_D1,[],1)./sqrt(size(NR_hist_D1,1)), 'k'); hold on;
shadedErrorBar(tth, nanmean(OR_hist_D1,1), nanstd(OR_hist_D1,[],1)./sqrt(size(OR_hist_D1,1)), 'r');
xlim([-1000 1950]);     ylim([0 1.5]); 
vline(0,'c');       vline(600, 'b'); 
xlabel('Time (ms)');       ylabel('Firing rate (Hz)');
[~,cell_countNR] = cellfun(@size, normalR_hist_D1);         [~,cell_countOR] = cellfun(@size, omitR_hist_D1); 
title(['Day1: ', num2str(size(NR_hist_D1,1)), 'neurons. NR: ', num2str(length(cell_countNR>1)), 'sessions.  OR: ', num2str(length(cell_countOR>1)), ' sessions']);

subplot(2,2,2)
shadedErrorBar(tth, nanmean(NR_hist_P1,1), nanstd(NR_hist_P1,[],1)./sqrt(size(NR_hist_P1,1)), 'k');
hold on;
shadedErrorBar(tth, nanmean(OR_hist_P1,1), nanstd(OR_hist_P1,[],1)./sqrt(size(OR_hist_P1,1)), 'r');
xlim([-1000 1950]);
ylim([0 1.5]); vline(0,'c'); vline(600, 'b'); 
xlabel('Time (ms)');
ylabel('Firing rate (Hz)');
title(['PL n=', num2str(size(NR_hist_P1,1)), ' neurons']);
[~,cell_countNR] = cellfun(@size, normalR_hist_P1);         [~,cell_countOR] = cellfun(@size, omitR_hist_P1); 
title(['PL: ', num2str(size(NR_hist_P1,1)), 'neurons. NR: ', num2str(length(cell_countNR>1)), 'sessions.  OR: ', num2str(length(cell_countOR>1)), ' sessions']);

subplot(2,2,3)
shadedErrorBar(tth, nanmean(NR_hist_PN1,1), nanstd(NR_hist_PN1,[],1)./sqrt(size(NR_hist_PN1,1)), 'k');    hold on;
shadedErrorBar(tth, nanmean(UR_hist_PN1,1), nanstd(UR_hist_PN1,[],1)./sqrt(size(UR_hist_PN1,1)), 'g');
xlim([-1000 1950]);    ylim([0 1.5]); 
vline(0, '--c');      vline(600,'b');
xlabel('Time (ms)')
ylabel('Firing rate (Hz)')
[~,cell_countNR] = cellfun(@size, normalR_hist_PN1);         [~,cell_countOR] = cellfun(@size, unexpR_hist_PN1); 
title(['Day1: ', num2str(size(NR_hist_PN1,1)), 'neurons. NR: ', num2str(length(cell_countNR>1)), 'sessions.  OR: ', num2str(length(cell_countOR>1)), ' sessions']);

subplot(2,2,4)
shadedErrorBar(tth, nanmean(NR_hist_PN2,1), nanstd(NR_hist_PN2,[],1)./sqrt(size(NR_hist_PN2,1)), 'k');
hold on;
shadedErrorBar(tth, nanmean(OR_hist_PN2,1), nanstd(OR_hist_PN2,[],1)./sqrt(size(OR_hist_PN2,1)), 'r');
xlim([-1000 1950]);
ylim([0 1.5]); vline(0, 'c'); vline(600,'--b'); vline(1100, 'b');
xlabel('Time (ms)')
ylabel('Firing rate (Hz)')
[~,cell_countNR] = cellfun(@size, normalR_hist_PN2);         [~,cell_countOR] = cellfun(@size, omitR_hist_PN2); 
title(['1000ms: ', num2str(size(NR_hist_PN2,1)), 'neurons. NR: ', num2str(length(cell_countNR>1)), 'sessions.  OR: ', num2str(length(cell_countOR>1)), ' sessions']);

suptitle(['All cells: All trials'] );
saveas(fig, [out_base 'Summary_PSTH_avgexpt_30Hz.fig']);
print([out_base 'Summary_PSTH_avgexpt_30Hz.eps'], '-depsc');
print([out_base 'Summary_PSTH_avgexpt_30Hz.pdf'], '-dpdf');


fig=figure('rend', 'painters', 'pos', [50 150 1200 550]);
subplot(2,2,1)
shadedErrorBar(tth, nanmean(NR_hist_NL_D1,1), nanstd(NR_hist_NL_D1,[],1)./sqrt(size(NR_hist_NL_D1,1)), 'k');   hold on;
shadedErrorBar(tth, nanmean(OR_hist_NL_D1,1), nanstd(OR_hist_NL_D1,[],1)./sqrt(size(OR_hist_NL_D1,1)), 'r');
xlim([-1000 1950]);
ylim([0 2]); vline(0,'c'); vline(600, 'b');  vline(500, 'k');
xlabel('Time (ms)');    ylabel('Firing rate (Hz)')
[~,cell_countNR] = cellfun(@size, normalR_hist_nolick_D1);         [~,cell_countOR] = cellfun(@size, omitR_hist_nolick_D1); 
title(['Day1: ', num2str(size(NR_hist_NL_D1,1)), 'neurons. NR: ', num2str(length(cell_countNR>1)), 'sessions.  OR: ', num2str(length(cell_countOR>1)), ' sessions']);

subplot(2,2,2)
shadedErrorBar(tth, nanmean(NR_hist_NL_P1,1), nanstd(NR_hist_NL_P1,[],1)./sqrt(size(NR_hist_NL_P1,1)), 'k');    hold on;
shadedErrorBar(tth, nanmean(OR_hist_NL_P1,1), nanstd(OR_hist_NL_P1,[],1)./sqrt(size(OR_hist_NL_P1,1)), 'r');
xlim([-1000 1950]);    ylim([0 2]); 
vline(0,'c'); vline(600, 'b');  vline(500, 'k');
xlabel('Time (ms)');    ylabel('Firing rate (Hz)');
[~,cell_countNR] = cellfun(@size, normalR_hist_nolick_P1);         [~,cell_countOR] = cellfun(@size, omitR_hist_nolick_P1); 
title(['PL: ', num2str(size(NR_hist_NL_P1,1)), 'neurons. NR: ', num2str(length(cell_countNR>1)), 'sessions.  OR: ', num2str(length(cell_countOR>1)), ' sessions']);

subplot(2,2,3)
shadedErrorBar(tth, nanmean(NR_hist_NL_PN1,1), nanstd(NR_hist_NL_PN1,[],1)./sqrt(size(NR_hist_NL_PN1,1)), 'k'); hold on;
shadedErrorBar(tth, nanmean(UR_hist_NL_PN1,1), nanstd(UR_hist_NL_PN1,[],1)./sqrt(size(UR_hist_NL_PN1,1)), 'g');
xlim([-1000 1950]);
ylim([0 2]); vline(0,'--c'); vline(600, 'b');  vline(500, 'k');
xlabel('Time (ms)')
ylabel('Firing rate (Hz)')
title(['UR n=', num2str(size(UR_hist_NL_PN1,1)), ' neurons']);
[~,cell_countNR] = cellfun(@size, normalR_hist_nolick_PN1);         [~,cell_countOR] = cellfun(@size, unexpR_hist_nolick_PN1); 
title(['UR: ', num2str(size(NR_hist_NL_PN1,1)), 'neurons. NR: ', num2str(length(cell_countNR>1)), 'sessions.  OR: ', num2str(length(cell_countOR>1)), ' sessions']);

subplot(2,2,4)
shadedErrorBar(tth, nanmean(NR_hist_NL_PN2,1), nanstd(NR_hist_NL_PN2,[],1)./sqrt(size(NR_hist_NL_PN2,1)), 'k'); hold on;
shadedErrorBar(tth, nanmean(OR_hist_NL_PN2,1), nanstd(OR_hist_NL_PN2,[],1)./sqrt(size(OR_hist_NL_PN2,1)), 'r');
xlim([-1000 1950]);
ylim([0 2]); vline(0, 'c'); vline(600, '--b'); vline(1100, 'b');
xlabel('Time (ms)')
ylabel('Firing rate (Hz)')
title(['1000ms delay n=', num2str(size(OR_hist_NL_PN2,1)), ' neurons']);
[~,cell_countNR] = cellfun(@size, normalR_hist_nolick_PN2);         [~,cell_countOR] = cellfun(@size, omitR_hist_nolick_PN2); 
title(['1000ms: ', num2str(size(NR_hist_NL_PN2,1)), 'neurons. NR: ', num2str(length(cell_countNR>1)), 'sessions.  OR: ', num2str(length(cell_countOR>1)), ' sessions']);

suptitle(['All cells: No lick condition: n=', num2str(length(sess_subset)), ' animals'] );
saveas(fig, [out_base 'Summary_PSTH_avgexpt_30Hz_nolick.fig']);
print([out_base 'Summary_PSTH_avgexpt_30Hz_nolick.eps'], '-depsc');
print([out_base 'Summary_PSTH_avgexpt_30Hz_nolick.pdf'], '-dpdf');

%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
fig=figure('rend', 'painters', 'pos', [50 150 1200 550]);
subplot(2,2,1); hold on;
shadedErrorBar(tth, nanmean(OR_hist_RS_D1,1), nanstd(OR_hist_RS_D1,[],1)./sqrt(size(OR_hist_RS_D1,1)), 'r');
shadedErrorBar(tth, nanmean(NR_hist_RS_D1,1), nanstd(NR_hist_RS_D1,[],1)./sqrt(size(NR_hist_RS_D1,1)), 'k');
xlim([-1000 1950])
ylim([0 2]); vline(600, 'b'); vline(0, 'c');
xlabel('Time (ms)')
ylabel('Firing rate (Hz)')
[~,cell_countNR] = cellfun(@size, normalR_hist_RS_D1);         [~,cell_countOR] = cellfun(@size, omitR_hist_RS_D1); 
title(['Day1: ', num2str(size(NR_hist_RS_D1,1)), 'neurons. NR: ', num2str(length(cell_countNR>1)), 'sessions.  OR: ', num2str(length(cell_countOR>1)), ' sessions']);

subplot(2,2,2)
hold on;
shadedErrorBar(tth, nanmean(OR_hist_RS_P1,1), nanstd(OR_hist_RS_P1,[],1)./sqrt(size(OR_hist_RS_P1,1)), 'r');
shadedErrorBar(tth, nanmean(NR_hist_RS_P1,1), nanstd(NR_hist_RS_P1,[],1)./sqrt(size(NR_hist_RS_P1,1)), 'k');
xlim([-1000 1950]);
ylim([0 2]); vline(600, 'b'); vline(0, 'c');
xlabel('Time (ms)')
ylabel('Firing rate (Hz)')
[~,cell_countNR] = cellfun(@size, normalR_hist_RS_P1);         [~,cell_countOR] = cellfun(@size, omitR_hist_RS_P1); 
title(['PL: ', num2str(size(NR_hist_RS_P1,1)), 'neurons. NR: ', num2str(length(cell_countNR>1)), 'sessions.  OR: ', num2str(length(cell_countOR>1)), ' sessions']);

subplot(2,2,3); hold on;
shadedErrorBar(tth, nanmean(UR_hist_RS_PN1,1), nanstd(UR_hist_RS_PN1,[],1)./sqrt(size(UR_hist_RS_PN1,1)), 'g');
shadedErrorBar(tth, nanmean(NR_hist_RS_PN1,1), nanstd(NR_hist_RS_PN1,[],1)./sqrt(size(NR_hist_RS_PN1,1)), 'k');
xlim([-1000 1950]);
ylim([0 2]); vline(0, '--c'); vline(600, 'b');
xlabel('Time (ms)');           ylabel('Firing rate (Hz)');
[~,cell_countNR] = cellfun(@size, normalR_hist_RS_PN1);         [~,cell_countOR] = cellfun(@size, unexpR_hist_RS_PN1); 
title(['UR: ', num2str(size(NR_hist_RS_PN1,1)), 'neurons. NR: ', num2str(length(cell_countNR>1)), 'sessions.  OR: ', num2str(length(cell_countOR>1)), ' sessions']);

subplot(2,2,4); hold on;
shadedErrorBar(tth, nanmean(OR_hist_RS_PN2,1), nanstd(OR_hist_RS_PN2,[],1)./sqrt(size(OR_hist_RS_PN2,1)), 'r');
shadedErrorBar(tth, nanmean(NR_hist_RS_PN2,1), nanstd(NR_hist_RS_PN2,[],1)./sqrt(size(NR_hist_RS_PN2,1)), 'k');
xlim([-1000 1950])
ylim([0 2]); vline(0, 'c'); vline(600, '--b'); vline(1000, 'b');
xlabel('Time (ms)');       ylabel('Firing rate (Hz)')
[~,cell_countNR] = cellfun(@size, normalR_hist_RS_PN2);         [~,cell_countOR] = cellfun(@size, omitR_hist_RS_PN2); 
title(['1000ms: ', num2str(size(NR_hist_RS_PN2,1)), 'neurons. NR: ', num2str(length(cell_countNR>1)), 'sessions.  OR: ', num2str(length(cell_countOR>1)), ' sessions']);

suptitle(['All resp neurons: all trials'] );
suptitle(['bout resp neurons: all trials: n=', num2str(length(sess_subset)), ' animals'] );
saveas(fig, [out_base 'Summary_PSTH_avgexpt_30Hz_resp.fig']);
print([out_base 'Summary_PSTH_avgexpt_30Hz_resp.eps'], '-depsc');
print([out_base 'Summary_PSTH_avgexpt_30Hz_resp.pdf'], '-dpdf');


fig=figure('rend', 'painters', 'pos', [50 150 1200 550]);
subplot(2,2,1);  hold on;
shadedErrorBar(tth, nanmean(OR_hist_NL_RS_D1,1), nanstd(OR_hist_NL_RS_D1,[],1)./sqrt(size(OR_hist_NL_RS_D1,1)), 'r');
shadedErrorBar(tth, nanmean(NR_hist_NL_RS_D1,1), nanstd(NR_hist_NL_RS_D1,[],1)./sqrt(size(NR_hist_NL_RS_D1,1)), 'k');
xlim([-1000 1950]);  ylim([0 2])
vline(600, 'b'); vline(500, 'k'); vline(0, 'c'); vline(600, 'b')
xlabel('Time (ms)'); ylabel('Firing rate (Hz)');
[~,cell_countNR] = cellfun(@size, normalR_hist_RS_PN2);         [~,cell_countOR] = cellfun(@size, omitR_hist_RS_PN2); 
title(['NR: ', num2str(size(NR_hist_NL_RS_D1,1)), 'neurons. NR: ', num2str(length(cell_countNR>1)), 'sessions.  OR: ', num2str(length(cell_countOR>1)), ' sessions']);

subplot(2,2,2); hold on;
shadedErrorBar(tth, nanmean(OR_hist_NL_RS_P1,1), nanstd(OR_hist_NL_RS_P1,[],1)./sqrt(size(OR_hist_NL_RS_P1,1)), 'r');
shadedErrorBar(tth, nanmean(NR_hist_NL_RS_P1,1), nanstd(NR_hist_NL_RS_P1,[],1)./sqrt(size(NR_hist_NL_RS_P1,1)), 'k');
xlim([-1000 1950]);    ylim([0 2]);
vline(600, 'b'); vline(500, 'k'); vline(0, 'c'); vline(600, 'b')
xlabel('Time (ms)');     ylabel('Firing rate (Hz)')
[~,cell_countNR] = cellfun(@size, normalR_hist_nolick_RS_P1);         [~,cell_countOR] = cellfun(@size, omitR_hist_nolick_RS_P1); 
title(['PL: ', num2str(size(NR_hist_NL_RS_P1,1)), 'neurons. NR: ', num2str(length(cell_countNR>1)), 'sessions.  OR: ', num2str(length(cell_countOR>1)), ' sessions']);

subplot(2,2,3);   hold on;
shadedErrorBar(tth, nanmean(UR_hist_NL_RS_PN1,1), nanstd(UR_hist_NL_RS_PN1,[],1)./sqrt(size(UR_hist_NL_RS_PN1,1)), 'g');
shadedErrorBar(tth, nanmean(NR_hist_NL_RS_PN1,1), nanstd(NR_hist_NL_RS_PN1,[],1)./sqrt(size(NR_hist_NL_RS_PN1,1)), 'k');
xlim([-1000 1950]);
ylim([0 2]);  vline(500, 'k'); vline(0, '--c'); vline(600, 'b');
xlabel('Time (ms)');       ylabel('Firing rate (Hz)')
[~,cell_countNR] = cellfun(@size, normalR_hist_nolick_RS_PN1);         [~,cell_countOR] = cellfun(@size, unexpR_hist_nolick_RS_PN1); 
title(['UR: ', num2str(size(NR_hist_NL_RS_PN1,1)), 'neurons. NR: ', num2str(length(cell_countNR>1)), 'sessions.  OR: ', num2str(length(cell_countOR>1)), ' sessions']);

subplot(2,2,4); hold on;
shadedErrorBar(tth, nanmean(OR_hist_NL_RS_PN2,1), nanstd(OR_hist_NL_RS_PN2,[],1)./sqrt(size(OR_hist_NL_RS_PN2,1)), 'r');
shadedErrorBar(tth, nanmean(NR_hist_NL_RS_PN2,1), nanstd(NR_hist_NL_RS_PN2,[],1)./sqrt(size(NR_hist_NL_RS_PN2,1)), 'k');
xlim([-1000 1950]);
ylim([0 2]); vline(500, 'k'); vline(0, 'c'); vline(600, '--b'); vline(1100,'b');
xlabel('Time (ms)')
ylabel('Firing rate (Hz)')
[~,cell_countNR] = cellfun(@size, normalR_hist_nolick_RS_PN2);         [~,cell_countOR] = cellfun(@size, omitR_hist_nolick_RS_PN2); 
title(['1000ms: ', num2str(size(NR_hist_NL_RS_PN2,1)), 'neurons. NR: ', num2str(length(cell_countNR>1)), 'sessions.  OR: ', num2str(length(cell_countOR>1)), ' sessions']);

suptitle(['bout resp neurons: no lick condition: n=', num2str(length(sess_subset)), ' animals'] );
saveas(fig, [out_base 'Summary_PSTH_avgexpt_30Hz_resp_nolick.fig']);
print([out_base 'Summary_PSTH_avgexpt_30Hz_resp_nolick.eps'], '-depsc');
print([out_base 'Summary_PSTH_avgexpt_30Hz_resp_nolick.pdf'], '-dpdf');

%===========================================================================================================

fig = figure('rend', 'painters', 'pos', [50 150 1200 550]);
subplot(2,2,1)
hold on;
shadedErrorBar(tth, nanmean(NR_hist_NL_RewNeg_P1,1), nanstd(NR_hist_NL_RewNeg_P1,[],1)./sqrt(size(NR_hist_NL_RewNeg_P1,1)), 'r');
shadedErrorBar(tth, nanmean(NR_hist_NL_RewPos_P1,1), nanstd(NR_hist_NL_RewPos_P1,[],1)./sqrt(size(NR_hist_NL_RewPos_P1,1)), 'g');
shadedErrorBar(tth, nanmean(NR_hist_NL_ORNeg_P1,1), nanstd(NR_hist_NL_ORNeg_P1,[],1)./sqrt(size(NR_hist_NL_ORNeg_P1,1)), 'm');
xlim([-1500 1500])
ylim([0 2.5])
xlabel('Time (ms)')
ylabel('Firing rate (Hz)')
title(['DayN Normal Reward P=', num2str(size(NR_hist_NL_RewPos_P1,1)), ' N=', num2str(size(NR_hist_NL_RewNeg_P1,1))]);

subplot(2,2,2)
hold on;
shadedErrorBar(tth, nanmean(OR_hist_NL_RewNeg_P1,1), nanstd(OR_hist_NL_RewNeg_P1,[],1)./sqrt(size(OR_hist_NL_RewNeg_P1,1)), 'r');
shadedErrorBar(tth, nanmean(OR_hist_NL_RewPos_P1,1), nanstd(OR_hist_NL_RewPos_P1,[],1)./sqrt(size(OR_hist_NL_RewPos_P1,1)), 'g');
xlim([-1500 1500])
ylim([0 2.5])
xlabel('Time (ms)')
ylabel('Firing rate (Hz)')
title(['DayN Reward Omission P=', num2str(size(OR_hist_NL_RewPos_P1,1)), ' N=', num2str(size(OR_hist_NL_RewNeg_P1,1))]);

subplot(2,2,3)
hold on;
shadedErrorBar(tth, nanmean(NR_hist_NL_RewNeg_PN2,1), nanstd(NR_hist_NL_RewNeg_PN2,[],1)./sqrt(size(NR_hist_NL_RewNeg_PN2,1)), 'r');
shadedErrorBar(tth, nanmean(NR_hist_NL_RewPos_PN2,1), nanstd(NR_hist_NL_RewPos_PN2,[],1)./sqrt(size(NR_hist_NL_RewPos_PN2,1)), 'g');
shadedErrorBar(tth, nanmean(NR_hist_NL_ORNeg_PN2,1), nanstd(NR_hist_NL_ORNeg_PN2,[],1)./sqrt(size(NR_hist_NL_ORNeg_PN2,1)), 'm');
xlim([-1500 1500])
ylim([0 2.5])
xlabel('Time (ms)')
ylabel('Firing rate (Hz)')
title(['1000ms Delay Normal Reward P=', num2str(size(NR_hist_NL_RewPos_PN2,1)), ' N=', num2str(size(NR_hist_NL_RewNeg_PN2,1))]);

subplot(2,2,4)
hold on;
shadedErrorBar(tth, nanmean(OR_hist_NL_RewNeg_PN2,1), nanstd(OR_hist_NL_RewNeg_PN2,[],1)./sqrt(size(OR_hist_NL_RewNeg_PN2,1)), 'r');
shadedErrorBar(tth, nanmean(OR_hist_NL_RewPos_PN2,1), nanstd(OR_hist_NL_RewPos_PN2,[],1)./sqrt(size(OR_hist_NL_RewPos_PN2,1)), 'g');
xlim([-1500 1500])
ylim([0 2.5])
xlabel('Time (ms)')
ylabel('Firing rate (Hz)')
title(['1000ms Delay Reward Omission P=', num2str(size(OR_hist_NL_RewPos_PN2,1)), ' N=', num2str(size(OR_hist_NL_RewNeg_PN2,1))]);

supertitle(['PSTH- responsive cells without lickbout- Green: Positive; Red: Negative; Magneta: Omission Negative'])
saveas(fig, [out_base 'Summary_PSTH_avgexpt_30Hz_nolick_PosNeg.fig']);
print([out_base 'Summary_PSTH_avgexpt_30Hz_nolick_PosNeg.eps'], '-depsc');
print([out_base 'Summary_PSTH_avgexpt_30Hz_nolick_PosNeg.pdf'], '-dpdf');

fig = figure('rend', 'painters', 'pos', [50 150 1200 550]);
subplot(2,2,1)
hold on;
shadedErrorBar(tth, nanmean(NR_hist_RewNeg_P1,1), nanstd(NR_hist_RewNeg_P1,[],1)./sqrt(size(NR_hist_RewNeg_P1,1)), 'r');
shadedErrorBar(tth, nanmean(NR_hist_RewPos_P1,1), nanstd(NR_hist_RewPos_P1,[],1)./sqrt(size(NR_hist_RewPos_P1,1)), 'g');
shadedErrorBar(tth, nanmean(NR_hist_ORNeg_P1,1), nanstd(NR_hist_ORNeg_P1,[],1)./sqrt(size(NR_hist_ORNeg_P1,1)), 'm');
xlim([-1500 1500])
ylim([0 2.5])
xlabel('Time (ms)')
ylabel('Firing rate (Hz)')
title(['DayN Normal Reward P=', num2str(size(NR_hist_RewPos_P1,1)), ' N=', num2str(size(NR_hist_RewNeg_P1,1))]);

subplot(2,2,2)
hold on;
shadedErrorBar(tth, nanmean(OR_hist_RewNeg_P1,1), nanstd(OR_hist_RewNeg_P1,[],1)./sqrt(size(OR_hist_RewNeg_P1,1)), 'r');
shadedErrorBar(tth, nanmean(OR_hist_RewPos_P1,1), nanstd(OR_hist_RewPos_P1,[],1)./sqrt(size(OR_hist_RewPos_P1,1)), 'g');
xlim([-1500 1500])
ylim([0 2.5])
xlabel('Time (ms)')
ylabel('Firing rate (Hz)')
title(['DayN Reward Omission P=', num2str(size(OR_hist_RewPos_P1,1)), ' N=', num2str(size(OR_hist_RewNeg_P1,1))]);

subplot(2,2,3)
hold on;
shadedErrorBar(tth, nanmean(NR_hist_RewNeg_PN2,1), nanstd(NR_hist_RewNeg_PN2,[],1)./sqrt(size(NR_hist_RewNeg_PN2,1)), 'r');
shadedErrorBar(tth, nanmean(NR_hist_RewPos_PN2,1), nanstd(NR_hist_RewPos_PN2,[],1)./sqrt(size(NR_hist_RewPos_PN2,1)), 'g');
shadedErrorBar(tth, nanmean(NR_hist_ORNeg_PN2,1), nanstd(NR_hist_ORNeg_PN2,[],1)./sqrt(size(NR_hist_ORNeg_PN2,1)), 'm');
xlim([-1500 1500])
ylim([0 2.5])
xlabel('Time (ms)')
ylabel('Firing rate (Hz)')
title(['1000ms Delay Normal Reward P=', num2str(size(NR_hist_RewPos_PN2,1)), ' N=', num2str(size(NR_hist_RewNeg_PN2,1))]);

subplot(2,2,4)
hold on;
shadedErrorBar(tth, nanmean(OR_hist_RewNeg_PN2,1), nanstd(OR_hist_RewNeg_PN2,[],1)./sqrt(size(OR_hist_RewNeg_PN2,1)), 'r');
shadedErrorBar(tth, nanmean(OR_hist_RewPos_PN2,1), nanstd(OR_hist_RewPos_PN2,[],1)./sqrt(size(OR_hist_RewPos_PN2,1)), 'g');
xlim([-1500 1500])
ylim([0 2.5])
xlabel('Time (ms)')
ylabel('Firing rate (Hz)')
title(['1000ms Delay Reward Omission P=', num2str(size(OR_hist_RewPos_PN2,1)), ' N=', num2str(size(OR_hist_RewNeg_PN2,1))]);

supertitle(['PSTH- responsive cells without lickbout- Green: Positive; Red: Negative; Magneta: Omission Negative'])
saveas(fig, [out_base 'Summary_PSTH_avgexpt_30Hz_withlick_PosNeg.fig']);
print([out_base 'Summary_PSTH_avgexpt_30Hz_withlick_PosNeg.eps'], '-depsc');
print([out_base 'Summary_PSTH_avgexpt_30Hz_withlick_PosNeg.pdf'], '-dpdf');
%scatter of peak rate
% fig=figure; 
% all_sr = [];
% all_si = [];
% RS_sr = [];
% RS_si = [];
% all_fr = [];
% all_fi = [];
% RS_fr = [];
% RS_fi = [];
% cid = 1;
% for id = 1:size(mouseID,2)
%     if id > 1 && ~strcmp(mouseID{id},mouseID{id-1})
%             
%             cid = cid + 1;
%              
%     end
%     subplot(2,2,1)
%     scatter(success_rate{id}, fail_rate{id},'MarkerEdgeColor',  col_mat(cid,:))
%     hold on
%     subplot(2,2,2)
%     scatter(success_RS_rate{id}, fail_RS_rate{id},'MarkerEdgeColor',  col_mat(cid,:))
%     hold on
%     subplot(2,2,3)
%     scatter(success_ind{id}, fail_ind{id}, 'MarkerEdgeColor', col_mat(cid,:))
%     hold on
%     subplot(2,2,4)
%     scatter(success_RS_ind{id}, fail_RS_ind{id},'MarkerEdgeColor',  col_mat(cid,:))
%     hold on
%     all_sr = [all_sr success_rate{id}];
%     all_si = [all_si success_ind{id}];
%     RS_sr = [RS_sr success_RS_rate{id}];
%     RS_si = [RS_si success_RS_ind{id}];
%     all_fr = [all_fr fail_rate{id}];
%     all_fi = [all_fi fail_ind{id}];
%     RS_fr = [RS_fr fail_RS_rate{id}];
%     RS_fi = [RS_fi fail_RS_ind{id}];
% end
% x = 0:.1:10;
% y = x;
% subplot(2,2,1)
% plot(x,y,'k')
% xlabel('Success rate (Hz)')
% ylabel('Fail rate (Hz)')
% xlim([-.5 8])
% ylim([-.5 8])
% [h_ra p_ra] = ttest(all_sr, all_fr);
% title(['All cells- n = ' num2str(total_cells) '; p = ' num2str(chop(p_ra,2))])
% subplot(2,2,2)
% plot(x,y,'k')
% xlim([-.5 8])
% ylim([-.5 8])
% xlabel('Success rate (Hz)')
% ylabel('Fail rate (Hz)')
% [h_rr p_rr] = ttest(RS_sr, RS_fr);
% title(['Responsive cells- n = ' num2str(total_resp) '; p = ' num2str(chop(p_rr,2))])
% subplot(2,2,3)
% xlabel('Success latency (ms)')
% ylabel('Fail latency (ms)')
% [h_la p_la] = ttest(all_si, all_fi);
% title(['All cells- n = ' num2str(total_cells) '; p = ' num2str(chop(p_la,2))])
% subplot(2,2,4)
% plot(x,y,'k')
% xlabel('Success latency (ms)')
% ylabel('Fail latency (ms)')
% [h_lr p_lr] = ttest(RS_si, RS_fi);
% title(['Responsive cells- n = ' num2str(total_resp) '; p = ' num2str(chop(p_lr,2))])
% supertitle(['Event rate and latency on releases- all cells'])
% saveas(fig, [out_base 'Summary_rate_latency_scatter.fig']);
% print([out_base 'Summary_rate_latency_scatter.eps'], '-depsc');
% print([out_base 'Summary_rate_latency_scatter.pdf'], '-dpdf');

%% plot latency 
% fig=figure;
% 
% subplot(1,2,1); axis([0 3 0 800]); axis square;
% scatter_plot(days_1, NR_D1_plot, NR_ind_D1, col_mat);
% ylabel('Latency of Response (ms)');
% title('Normal Reward Peak Onset Delay');
% 
% hold on;
% scatter_plot(days_post, NR_P1_plot, NR_ind_P1, col_mat);
% 
% % subplot(2,2,2); axis([0 3 0 800])
% % scatter_plot(days_1, NR_RS_D1_plot, NR_RS_ind_D1, col_mat);
% % 
% % hold on;
% % scatter_plot(days_post, NR_RS_P1_plot, NR_RS_ind_P1, col_mat);
% 
% subplot(1,2,2); axis([0 3 0 800]); axis square;
% scatter_plot(days_1, OR_D1_plot, OR_ind_D1, col_mat);
% ylabel('Latency of Response (ms)');
% title('Omission Reward Peak Onset Delay');
% 
% hold on;
% scatter_plot(days_post, OR_P1_plot, OR_ind_P1, col_mat);
% 
% supertitle(['for Day1- nCell = ', num2str(total_cellD1), ' and Post Day1- nCell = ', num2str(total_cellP1)]);
% saveas(fig, [out_base 'Summary_resp_latency.fig']);
% print([out_base 'Summary_resp_latency.eps'], '-depsc');
% print([out_base 'Summary_resp_latency.pdf'], '-dpdf');
%%
% % subplot(2,2,4); axis([0 3 0 800])
% % scatter_plot(days_1, OR_RS_D1_plot, OR_RS_ind_D1, col_mat);
% % 
% % hold on;
% % scatter_plot(days_post, OR_RS_P1_plot, OR_RS_ind_P1, col_mat);
% 
% 
% x = 0:.1:10;
% y = x;
% subplot(2,2,1)
% plot(x,y,'k')
% xlim([0 6])
% ylim([0 6])
% xlabel('Success rate (Hz)')
% ylabel('Fail rate (Hz)')
% title(['All cells- rate'])
% subplot(2,2,2)
% plot(x,y,'k')
% xlim([0 6])
% ylim([0 6])
% xlabel('Success rate (Hz)')
% ylabel('Fail rate (Hz)')
% title(['Responsive cells- rate'])
% subplot(2,2,3)
% x = 0:1:150;
% y = x;
% plot(x,y,'k')
% xlabel('Success latency (ms)')
% ylabel('Fail latency (ms)')
% title(['All cells- latency'])
% subplot(2,2,4)
% plot(x,y,'k')
% xlabel('Success latency (ms)')
% ylabel('Fail latency (ms)')
% title(['Responsive cells- latency'])
% supertitle(['Event rate and latency on releases- average within expts'])
% saveas(fig, [out_base 'Summary_rate_latency_scatter_avg.fig']);
% print([out_base 'Summary_rate_latency_scatter_avg.eps'], '-depsc');
% print([out_base 'Summary_rate_latency_scatter_avg.pdf'], '-dpdf');

% fig=figure;
% subplot(1,2,1)
% scatter_plot(mouseID, success_syn_all_RS, fail_syn_all_RS, col_mat);
% hold on
% x = 0:1:150;
% y = x;
% plot(x,y,'k')
% xlabel('Success standard deviation of latency')
% ylabel('Fail standard deviation of latency')
% title(['All cells- standard deviation of latency across cells'])
% 
% subplot(1,2,2)
% scatter_plot(mouseID, success_syn_c_RS, fail_syn_c_RS, col_mat);
% hold on
% x = 0:1:150;
% y = x;
% plot(x,y,'k')
% xlabel('Success standard deviation of latency')
% ylabel('Fail standard deviation of latency')
% title(['All cells- standard deviation of latency across events'])
% supertitle('Syncrony');
% saveas(fig, [out_base 'Summary_std_latency_scatter_avg_abs.fig']);
% print([out_base 'Summary_std_latency_scatter_avg_abs.eps'], '-depsc');
% print([out_base 'Summary_std_latency_scatter_avg_abs.pdf'], '-dpdf');
% 
% %summary of average event waveform relative to spont
% fig=figure;
for id = 1:size(days_1,2)
%     subplot(4,5,id)
    normalR_norm_D1 = bsxfun(@rdivide, normalR_norm_all_D1{id}, max(spont_norm_all_D1{id},[],1));
    omitR_norm_D1 = bsxfun(@rdivide,omitR_norm_all_D1{id}, max(spont_norm_all_D1{id},[],1));
%     spont_norm = bsxfun(@rdivide,spont_norm_all{id}, max(spont_norm_all{id},[],1));
    normalCue_norm_D1 = bsxfun(@rdivide,normalCue_norm_all_D1{id}, max(spont_norm_all_D1{id},[],1));
    omitCue_norm_D1 = bsxfun(@rdivide,omitCue_norm_all_D1{id}, max(spont_norm_all_D1{id},[],1));
    normalR_norm_peak_D1{id} = max(normalR_norm_D1,[],1);
    omitR_norm_peak_D1{id} = max(omitR_norm_D1,[],1);
%     spont_norm_peak{id} = max(spont_norm,[],1);
    normalCue_norm_peak_D1{id} = max(normalCue_norm_D1,[],1);
    omitCue_norm_peak_D1{id} = max(omitCue_norm_D1,[],1);
    
    
    normalR_norm_P1 = bsxfun(@rdivide, normalR_norm_all_P1{id}, max(spont_norm_all_P1{id},[],1));
    omitR_norm_P1 = bsxfun(@rdivide,omitR_norm_all_P1{id}, max(spont_norm_all_P1{id},[],1));
%     spont_norm = bsxfun(@rdivide,spont_norm_all{id}, max(spont_norm_all{id},[],1));
    normalCue_norm_P1 = bsxfun(@rdivide,normalCue_norm_all_P1{id}, max(spont_norm_all_P1{id},[],1));
    omitCue_norm_P1 = bsxfun(@rdivide,omitCue_norm_all_P1{id}, max(spont_norm_all_P1{id},[],1));
    normalR_norm_peak_P1{id} = max(normalR_norm_P1,[],1);
    omitR_norm_peak_P1{id} = max(omitR_norm_P1,[],1);
%     spont_norm_peak{id} = max(spont_norm,[],1);
    normalCue_norm_peak_P1{id} = max(normalCue_norm_P1,[],1);
    omitCue_norm_peak_P1{id} = max(omitCue_norm_P1,[],1);
%     
%     cell_use = find(~isnan(mean([release_norm; press_norm; spont_norm],1)));
%     errorbar(ts{id}, nanmean(release_norm(:,cell_use),2), std(release_norm(:,cell_use),[],2)./sqrt(length(cell_use)), 'k')
%     hold on;
%     errorbar(ts{id}, nanmean(press_norm(:,cell_use),2), std(press_norm(:,cell_use),[],2)./sqrt(length(cell_use)), 'c')
%     hold on;
%     errorbar(ts{id}, nanmean(spont_norm(:,cell_use),2), std(spont_norm(:,cell_use),[],2)./sqrt(length(cell_use)), 'b')
%     title([date{id} ' ' mouseID{id} '-'  num2str(length(cell_use)) ' cells'])   
%     xlabel('Time (ms)')
%     ylabel('Norm dF/F')
end
% supertitle('All cells- Event waveform normalized to spont')
% saveas(fig, [out_base 'Summary_avgevent_norm2spont.fig']);
% print([out_base 'Summary_avgevent_norm2spont.eps'], '-depsc');
% print([out_base 'Summary_avgevent_norm2spont.pdf'], '-dpdf');
% 
% fig=figure;
for id = 1:size(days_1,2)
%     subplot(4,5,id)
    normScale = max(spont_norm_RS_D1{id},[],1);
%      normScale = 1;
    normalR_renorm_RS_D1 = bsxfun(@rdivide,normalR_norm_RS_D1{id}, normScale);
    omitR_renorm_RS_D1 = bsxfun(@rdivide,omitR_norm_RS_D1{id}, normScale);
%     spont_renorm_RS = bsxfun(@rdivide,spont_norm_RS{id}, normScale);
    normalCue_renorm_RS_D1 = bsxfun(@rdivide,normalCue_norm_RS_D1{id}, normScale);
    omitCue_renorm_RS_D1 = bsxfun(@rdivide,omitCue_norm_RS_D1{id}, normScale);
    normalR_norm_peak_RS_D1{id} = max(normalR_renorm_RS_D1,[],1);
    omitR_norm_peak_RS_D1{id} = max(omitR_renorm_RS_D1,[],1);
%     spont_norm_peak_RS{id} = max(spont_renorm_RS,[],1);
    normalCue_norm_peak_RS_D1{id} = max(normalCue_renorm_RS_D1,[],1);
    omitCue_norm_peak_RS_D1{id} = max(omitCue_renorm_RS_D1,[],1);
    
    
    normScale = max(spont_norm_RS_P1{id},[],1);
%      normScale = 1;
    normalR_renorm_RS_P1 = bsxfun(@rdivide,normalR_norm_RS_P1{id}, normScale);
    omitR_renorm_RS_P1 = bsxfun(@rdivide,omitR_norm_RS_P1{id}, normScale);
%     spont_renorm_RS = bsxfun(@rdivide,spont_norm_RS{id}, normScale);
    normalCue_renorm_RS_P1 = bsxfun(@rdivide,normalCue_norm_RS_P1{id}, normScale);
    omitCue_renorm_RS_P1 = bsxfun(@rdivide,omitCue_norm_RS_P1{id}, normScale);
    normalR_norm_peak_RS_P1{id} = max(normalR_renorm_RS_P1,[],1);
    omitR_norm_peak_RS_P1{id} = max(omitR_renorm_RS_P1,[],1);
%     spont_norm_peak_RS{id} = max(spont_renorm_RS,[],1);
    normalCue_norm_peak_RS_P1{id} = max(normalCue_renorm_RS_P1,[],1);
    omitCue_norm_peak_RS_P1{id} = max(omitCue_renorm_RS_P1,[],1);
    
%     cell_use = find(~isnan(mean([release_renorm_RS; press_renorm_RS; spont_renorm_RS],1)));
%     errorbar(ts{id}, nanmean(release_renorm_RS(:,cell_use),2), std(release_renorm_RS(:,cell_use),[],2)./sqrt(length(cell_use)), 'k')
%     hold on;
%     errorbar(ts{id}, nanmean(press_renorm_RS(:,cell_use),2), std(press_renorm_RS(:,cell_use),[],2)./sqrt(length(cell_use)), 'c')
%     hold on;
%     errorbar(ts{id}, nanmean(spont_renorm_RS(:,cell_use),2), std(spont_renorm_RS(:,cell_use),[],2)./sqrt(length(cell_use)), 'b')
%     title([date{id} ' ' mouseID{id} '-'  num2str(length(cell_use)) ' cells'])   
%     xlabel('Time (ms)')
%     ylabel('Norm dF/F')
end
% supertitle('Responsive cells- Event waveform normalized to spont')
% saveas(fig, [out_base 'Summary_avgevent_norm2spont.fig']);
% print([out_base 'Summary_avgevent_norm2spont.eps'], '-depsc');
% print([out_base 'Summary_avgevent_norm2spont.pdf'], '-dpdf');
% 
% %peak scatters
% r_all = [];
% s_all = [];  % success
% f_all = [];  % fail
% sp_all = []; % spont
% p_all = [];
% r_RS = [];
% s_RS = [];
% f_RS = [];
% sp_RS = [];
% p_RS = [];
% fig=figure;
% cid = 1;
% for id = 1:size(mouseID,2)
%     if id > 1 && ~strcmp(mouseID{id},mouseID{id-1})
%             
%             cid = cid + 1;
%              
%     end
%     r_all = [r_all release_norm_peak{id}];
%     p_all = [p_all press_norm_peak{id}];
%     r_RS = [r_RS release_norm_peak_RS{id}];
%     p_RS = [p_RS press_norm_peak_RS{id}];   
%     subplot(1,2,1)
%     scatter(release_norm_peak{id}, press_norm_peak{id}, 'MarkerEdgeColor', col_mat(cid,:))
%     hold on;
%     subplot(1,2,2)
%     scatter(release_norm_peak_RS{id}, press_norm_peak_RS{id},'MarkerEdgeColor',  col_mat(cid,:))
%     hold on;
% end
% subplot(1,2,1)
% xlim([0 2.5])
% ylim([0 2])
% hline(1, '--k')
% hold on;
% vline(1, '--k')
% use_cells_all = find(~isnan(mean([p_all; r_all],1)));
% [h_rall p_rall] = ttest(r_all(1,use_cells_all),1);
% [h_pall p_pall] = ttest(p_all(1,use_cells_all),1);
% xlabel('Release amplitude')
% ylabel('Press amplitude')
% ncells = sum(~isnan(mean([p_all; r_all],1)),2);
% title(['All cells- n = ' num2str(total_cells) '; Release p = ' num2str(chop(p_rall,2)) '; Press p = ' num2str(chop(p_pall,2))])
% subplot(1,2,2)
% xlim([0 2.5])
% ylim([0 2])
% hline(1, '--k')
% hold on;
% vline(1, '--k')
% use_cells_RS = find(~isnan(mean([p_RS; r_RS],1)),2);
% [h_rRS p_rRS] = ttest(r_RS(1,use_cells_RS),1);
% [h_pRS p_pRS] = ttest(p_RS(1,use_cells_RS),1);
% xlabel('Release amplitude')
% ylabel('Press amplitude')
% ncells = sum(~isnan(mean([p_RS; r_RS],1)),2);
% title(['Responsive cells- n = ' num2str(total_resp) '; Release p = ' num2str(chop(p_rRS,2)) '; Press p = ' num2str(chop(p_pRS,2))])
% supertitle('Event waveform normalized to spont')
% saveas(fig, [out_base 'Summary_pr_scatter_norm2spont.fig']);
% print([out_base 'Summary_pr_scatter_norm2spont.eps'], '-depsc');
% print([out_base 'Summary_pr_scatter_norm2spont.pdf'], '-dpdf');
% 
% 
% %peak scatters for spont, fail, success and press
% NR_all = [];  % 
% OR_all = [];  % 
% NRCue_all = []; %
% NR_RS = [];
% OR_RS = [];
% ORCue_RS = [];

% x = [0:.01:1.5];
% y = x;
% fig=figure;
% col_mat_s = [0.5, 0.5, 0.5];
% cid = 1;
% for id = 1:size(mouseID,2)
%     if id > 1 && ~strcmp(mouseID{id},mouseID{id-1})
%             
%             cid = cid + 1;
%              
%     end
%     NR_all = [NR_all normalR_norm_peak{id}];
%     OR_all = [OR_all omitR_norm_peak{id}];
%     NRCue_all = [NRCue_all normalCue_norm_peak{id}];
%     ORCue_all = [ORCue_all omitCue_norm_peak{id}];
% %     sp_all = [sp_all spont_norm_peak{id}];
%     NR_RS = [NR_RS normalR_norm_peak_RS{id}];
%     OR_RS = [OR_RS omitR_norm_peak_RS{id}];
%     NRCue_RS = [NRCue_RS normalCue_norm_peak_RS{id}];
%     ORCue_RS = [ORCue_RS omitCue_norm_peak_RS{id}];
% %     sp_RS = [sp_RS spont_norm_peak_RS{id}];
%     
% %     subplot(1,3,1)
% %     scatter(press_norm_peak_RS{id}, spont_norm_peak_RS{id}, 4, 'MarkerEdgeColor', col_mat_s)
% %     hold on;
% %     subplot(1,3,2)
% %     scatter(success_norm_peak_RS{id}, spont_norm_peak_RS{id}, 4, 'MarkerEdgeColor',  col_mat_s)
% %     hold on;
% %     subplot(1,3,3)
% %     scatter(fail_norm_peak_RS{id}, spont_norm_peak_RS{id}, 4, 'MarkerEdgeColor',  col_mat_s)
% %     hold on;
% end
% subplot(1,3,1)
% hold on; plot(x,y,'-k');
% xlim([0 1.5])
% ylim([0 1.5])
% use_cells_RS = find(~isnan(mean([p_RS; sp_RS],1)));
% [h_pspRS p_pspRS] = ttest(sp_RS(1,use_cells_RS),p_RS(1,use_cells_RS));
% % [h_pRS p_pRS] = ttest(p_RS(1,use_cells_RS),1);
% xlabel('Press amplitude')
% ylabel('Spontaneous amplitude')
% ncells = sum(~isnan(mean([p_RS; sp_RS],1)),2);
% title(['Responsive cells- n = ' num2str(total_resp) '; p = ' num2str(chop(p_pspRS,2))])
% 
% subplot(1,3,2)
% hold on; plot(x,y,'-k');
% xlim([0 1.5])
% ylim([0 1.5])
% use_cells_RS = find(~isnan(mean([s_RS; sp_RS],1)),2);
% [h_sRS p_sspRS] = ttest(s_RS(1,use_cells_RS),sp_RS(1,use_cells_RS));
% % [h_spRS p_spRS] = ttest(sp_RS(1,use_cells_RS),1);
% xlabel('Correct amplitude')
% ylabel('Spontaneous amplitude')
% ncells = sum(~isnan(mean([s_RS; sp_RS],1)),2);
% title(['Responsive cells- n = ' num2str(total_resp) '; p = ' num2str(chop(p_sspRS,2))])
% 
% subplot(1,3,3)
% hold on; plot(x,y,'-k');
% xlim([0 1.5])
% ylim([0 1.5])
% use_cells_RS = find(~isnan(mean([f_RS; sp_RS],1)),2);
% [h_fRS p_fspRS] = ttest(f_RS(1,use_cells_RS),sp_RS(1,use_cells_RS));
% % [h_spRS p_spRS] = ttest(sp_RS(1,use_cells_RS),1);
% xlabel('Early amplitude')
% ylabel('Spontaneous amplitude')
% ncells = sum(~isnan(mean([f_RS; sp_RS],1)),2);
% title(['Responsive cells- n = ' num2str(total_resp) '; p = ' num2str(chop(p_fspRS,2))])
% 
% supertitle('Peak of event waveform')
% % saveas(fig, [out_base 'Summary_event_peak_scatter.fig']);
% % print([out_base 'Summary_event_peak_scatter.eps'], '-depsc');
% % print([out_base 'Summary_event_peak_scatter.pdf'], '-dpdf');
% 
col_mat_s = 0.5*ones(10,3);
fig = figure;
cid =1;
subplot(1,2,1)
scatter_plot(days_1_mouse, normalR_norm_peak_D1, normalCue_norm_peak_P1,col_mat);
axis square

subplot(1,2,2)
scatter_plot(days_1_mouse, normalR_norm_peak_RS_D1, normalCue_norm_peak_RS_P1, col_mat);
axis square
% for id = 1:size(mouseID,2)
%     if id > 1 && ~strcmp(mouseID{id},mouseID{id-1})
%             
%             cid = cid + 1;
%              
%     end
% scatter(normalR_norm_peak_P1{id}, normalCue_norm_peak_D1{id}, 4, 'MarkerEdgeColor',  col_mat_s)
% hold on;
% end
% 
% cid =1;
% subplot(1,2,2)
% for id = 1:size(mouseID,2)
%     if id > 1 && ~strcmp(mouseID{id},mouseID{id-1})
%             
%             cid = cid + 1;
%              
%     end
% scatter(normalR_norm_peak_RS_P1{id}, normalCue_norm_peak_RS_D1{id}, 4, 'MarkerEdgeColor',  col_mat_s)
% hold on;
% end

ax1 = subplot(1,2,1);
xlim([0.1 10])
ylim([0.1 10])
hline(1, '--k')
hold on;
vline(1, '--k')
set(ax1, 'XScale', 'log', 'YScale', 'log');
xlabel('Day1 Reward amplitude')
ylabel('Day5 Cue amplitude')
title('All cells');

ax2 = subplot(1,2,2);
xlim([0.1 10])
ylim([0.1 10])
hline(1, '--k')
hold on;
vline(1, '--k')
set(ax2, 'XScale', 'log', 'YScale', 'log');
xlabel('Day1 Reward amplitude')
ylabel('Day5 Cue amplitude')
title('Responsive cells');
% use_cells_RS = find(~isnan(mean([s_RS; f_RS],1)),2);
% [h_sfRS p_sfRS] = ttest(s_RS(1,use_cells_RS),f_RS(1,use_cells_RS));
% title(['Responsive cells- n = ' num2str(total_resp) '; p = ' num2str(chop(p_sfRS,2))])
saveas(fig, [out_base 'Summary_event_peak_scatter_CE_norm.fig']);
print([out_base 'Summary_event_peak_scatter_CE_norm.eps'], '-depsc');
print([out_base 'Summary_event_peak_scatter_CE_norm.pdf'], '-dpdf');
% 
% %average across mice
% fig=figure;
% cid = 1;
% for id = 1:size(mouseID,2) 
%     if id > 1 && ~strcmp(mouseID{id},mouseID{id-1})
%             
%             cid = cid + 1;
%              
%     end
%     subplot(1,2,1)
%     temp_release = release_norm_peak{id};
%     temp_press = press_norm_peak{id};
%     use_cells_all = find(~isnan(mean([temp_release; temp_press],1)));
%     temp_release_RS = release_norm_peak_RS{id};
%     temp_press_RS = press_norm_peak_RS{id};
%     use_cells_RS = find(~isnan(mean([temp_release_RS; temp_press_RS],1)));
%     errorbarxy(mean(temp_release(:,use_cells_all),2), mean(temp_press(:,use_cells_all),2), std(temp_release(:,use_cells_all),[],2)./sqrt(length(use_cells_all)), std(temp_press(:,use_cells_all),[],2)./sqrt(length(use_cells_all)),{'o', col_mat(cid,:), col_mat(cid,:),col_mat(cid,:)})
%     hold on;
%     scatter(mean(temp_release(:,use_cells_all),2), mean(temp_press(:,use_cells_all),2), 'MarkerEdgeColor', col_mat(cid,:))
%     subplot(1,2,2)
%     errorbarxy(mean(temp_release_RS(:,use_cells_RS),2), mean(temp_press_RS(:,use_cells_RS),2), std(temp_release_RS(:,use_cells_RS),[],2)./sqrt(length(use_cells_RS)), std(temp_press_RS(:,use_cells_RS),[],2)./sqrt(length(use_cells_RS)),{'o', col_mat(cid,:), col_mat(cid,:),col_mat(cid,:)})
%     hold on;
%     scatter(mean(temp_release_RS(:,use_cells_RS),2), mean(temp_press_RS(:,use_cells_RS),2),'MarkerEdgeColor', col_mat(cid,:))
% end
% subplot(1,2,1)
% xlim([0.5 1.5])
% ylim([0.5 1.5])
% hline(1, '--k')
% hold on;
% vline(1, '--k')
% xlabel('Release amplitude')
% ylabel('Press amplitude')
% title('All cells')
% subplot(1,2,2)
% xlim([0 0.6])
% ylim([0 0.6])
% hline(1, '--k')
% hold on;
% vline(1, '--k')
% xlabel('Release amplitude')
% ylabel('Press amplitude')
% title('Responsive cells')
% supertitle('Avg event waveform normalized to spont')
% saveas(fig, [out_base 'Summary_avgpr_scatter_norm2spont.fig']);
% print([out_base 'Summary_avgpr_scatter_norm2spont.eps'], '-depsc');
% print([out_base 'Summary_avgpr_scatter_norm2spont.pdf'], '-dpdf');
% 
% 
% %press and release alone
% fig=figure;
% lever = strvcat('release', 'press');
% % cid = 1;
% % for id = 1:size(mouseID,2) 
% %     if id > 1 && ~strcmp(mouseID{id},mouseID{id-1})
% %             
% %             cid = cid + 1;
% %              
% %     end
% %     subplot(2,2,1)
% %     scatter(ones(size(release_norm_peak{id})), release_norm_peak{id}, 'MarkerEdgeColor',col_mat(cid,:))
% %     hold on;
% %     scatter(ones(size(press_norm_peak{id})).*2, press_norm_peak{id}, 'MarkerEdgeColor',col_mat(cid,:))
% %     hold on;
% %     subplot(2,2,2)
% %     scatter(ones(size(release_norm_peak_RS{id})), release_norm_peak_RS{id}, 'MarkerEdgeColor',col_mat(cid,:))
% %     hold on;
% %     scatter(ones(size(press_norm_peak_RS{id})).*2, press_norm_peak_RS{id}, 'MarkerEdgeColor',col_mat(cid,:))
% %     hold on;
% % end
% % subplot(2,2,1)
% % xlim([0.5 2.5])
% % ylim([0 2.5])
% % hline(1, '--k')
% % set(gca, 'XTick', 1:2, 'XTickLabels', lever);
% % [h_rall p_rall] = ttest(r_all,1);
% % [h_pall p_pall] = ttest(p_all,1);
% % ylabel('Normalized amplitude')
% % title(['Release p = ' num2str(chop(p_rall,2)) '; Press p = ' num2str(chop(p_pall,2))])
% % subplot(2,2,2)
% % xlim([0.5 2.5])
% % ylim([0 2.5])
% % hline(1, '--k')
% % set(gca, 'XTick', 1:2,'XTickLabels', lever);
% % [h_rRS p_rRS] = ttest(r_RS,1);
% % [h_pRS p_pRS] = ttest(p_RS,1);
% % ylabel('Normalized amplitude')
% % title(['Release p = ' num2str(chop(p_rRS,2)) '; Press p = ' num2str(chop(p_pRS,2))])
% cid = 1;
% %average across mice
% for id = 1:size(mouseID,2) 
%     if id > 1 && ~strcmp(mouseID{id},mouseID{id-1})
%             
%             cid = cid + 1;
%              
%     end
%     subplot(1,2,1)
%     temp_release = release_norm_peak{id};
%     temp_press = press_norm_peak{id};
%     temp_release_RS = release_norm_peak_RS{id};
%     temp_press_RS = press_norm_peak_RS{id};
%     errorbar(1, nanmean(temp_release,2),nanstd(temp_release,[],2)./sqrt(size(temp_release,2)), 'o', 'color', col_mat(cid,:))
%     hold on;
%     errorbar(2, nanmean(temp_press,2), nanstd(temp_press,[],2)./sqrt(size(temp_press,2)), 'o', 'color', col_mat(cid,:))
%     hold on;
%     subplot(1,2,2)
%     errorbar(1, nanmean(temp_release_RS,2), nanstd(temp_release_RS,[],2)./sqrt(size(temp_release_RS,2)), 'o', 'color', col_mat(cid,:))
%     hold on;
%     errorbar(2, nanmean(temp_press_RS,2), nanstd(temp_press_RS,[],2)./sqrt(size(temp_press_RS,2)), 'o', 'color',  col_mat(cid,:))
%     hold on;
% end
% subplot(1,2,1)
% xlim([0.5 2.5])
% ylim([0 1.5])
% hline(1, '--k')
% set(gca, 'XTick', 1:2, 'XTickLabel', lever);
% ylabel('Normalized amplitude')
% title(['All cells- Release n = ' num2str(total_cells) ' cells; Press n = ' num2str(sum(~isnan(p_all),2))])
% subplot(1,2,2)
% xlim([0.5 2.5])
% ylim([0 1.5])
% hline(1, '--k')
% set(gca, 'XTick', 1:2, 'XTickLabel', lever);
% ylabel('Normalized amplitude')
% title(['Responsive cells- Release n = ' num2str(total_resp) ' cells; Press n = ' num2str(sum(~isnan(p_RS),2))])
% supertitle('Avg event waveform normalized to spont')
% saveas(fig, [out_base 'Summary_avg_event_norm2spont.fig']);
% print([out_base 'Summary_avg_event_norm2spont.eps'], '-depsc');
% print([out_base 'Summary_avg_event_norm2spont.pdf'], '-dpdf');
% 
% % %% need to correct for acquisition rate
% % %spont waveform comparison
% % press_norm_all = cellfun(@transpose, press_norm_all, 'UniformOutput', 0);
% % press_norm_RS = cellfun(@transpose, press_norm_RS, 'UniformOutput', 0);
% % success_norm_all = cellfun(@transpose, success_norm_all, 'UniformOutput', 0);
% % success_norm_RS = cellfun(@transpose, success_norm_RS, 'UniformOutput', 0);
% % spont_norm_all = cellfun(@transpose, spont_norm_all, 'UniformOutput', 0);
% % spont_norm_RS = cellfun(@transpose, spont_norm_RS, 'UniformOutput', 0);
% % 
% % frame_size = cell2mat(cellfun(@size, success_norm_all, 'UniformOutput', 0));
% % max_frame = max(frame_size(2:2:end));
% % 
% % all_press_norm = interp_frame(press_norm_all, max_frame);
% % RS_press_norm    = interp_frame(press_norm_RS, max_frame);
% % all_release_norm   = interp_frame(success_norm_all, max_frame);
% % RS_success_norm = interp_frame(success_norm_RS, max_frame);
% % all_spont_norm    = interp_frame(spont_norm_all, max_frame);
% % RS_spont_norm   = interp_frame(spont_norm_RS, max_frame);
% % 
% % fig=figure;
% % all_spont_renorm = bsxfun(@rdivide, all_spont_norm, max(all_spont_norm,[],2));
% % all_success_renorm = bsxfun(@rdivide, all_release_norm, max(all_spont_norm,[],2));
% % all_press_renorm = bsxfun(@rdivide, all_press_norm, max(all_spont_norm,[],2));
% % RS_spont_renorm = bsxfun(@rdivide, RS_spont_norm, max(RS_spont_norm,[],2));
% % RS_success_renorm = bsxfun(@rdivide, RS_success_norm, max(RS_spont_norm,[],2));
% % RS_press_renorm = bsxfun(@rdivide, RS_press_norm, max(RS_spont_norm,[],2));
% % cell_use = find(~isnan(mean([all_spont_renorm all_success_renorm all_press_renorm],2)));
% % tth = [1-(size(all_spont_renorm,2)/2):(size(all_spont_renorm,2)/2)].*double(min((TC_ifi)));
% % subplot(1,2,1)
% % errorbar(tth, nanmean(all_spont_renorm(cell_use,:),1), std(all_spont_renorm(cell_use,:),[],1)./sqrt(length(cell_use)), '-r')
% % hold on;
% % errorbar(tth, nanmean(all_success_renorm(cell_use,:),1), std(all_success_renorm(cell_use,:),[],1)./sqrt(length(cell_use)), '-k')
% % hold on;
% % errorbar(tth, nanmean(all_press_renorm(cell_use,:),1), std(all_press_renorm(cell_use,:),[],1)./sqrt(length(cell_use)), '-c')
% % hold on;
% % xlabel('Time (ms)')
% % ylabel('dF/F')
% % title(['Average event- all cells- ' num2str(total_cells) ' cells'])
% % 
% % subplot(1,2,2)
% % cell_use = find(~isnan(mean([RS_spont_renorm RS_success_renorm RS_press_renorm],2)));
% % errorbar(tth, nanmean(RS_spont_renorm(cell_use,:),1), std(RS_spont_renorm(cell_use,:),[],1)./sqrt(length(cell_use)), '-r')
% % hold on;
% % errorbar(tth, nanmean(RS_success_renorm(cell_use,:),1), std(RS_success_renorm(cell_use,:),[],1)./sqrt(length(cell_use)), '-k')
% % hold on;
% % errorbar(tth, nanmean(RS_press_renorm(cell_use,:),1), std(RS_press_renorm(cell_use,:),[],1)./sqrt(length(cell_use)), '-c')
% % hold on;
% % xlabel('Time (ms)')
% % ylabel('dF/F')
% % title(['Average event- responsive cells- ' num2str(total_resp) ' cells'])
% % supertitle('Average event - Red: spont; Black: release; Cyan: press')
% % saveas(fig, [out_base 'Summary_event_TC_norm2spont.fig']);
% % print([out_base 'Summary_event_TC_norm2spont.eps'], '-depsc');
% % print([out_base 'Summary_event_TC_norm2spont.pdf'], '-dpdf');
% 
% 
%spont waveform comparison
normalR_norm1_all_D1 = cellfun(@transpose, normalR_norm1_all_D1, 'UniformOutput', 0);
normalR_norm1_RS_D1 = cellfun(@transpose, normalR_norm1_RS_D1, 'UniformOutput', 0);
normalR_norm1_all_P1 = cellfun(@transpose, normalR_norm1_all_P1, 'UniformOutput', 0);
normalR_norm1_RS_P1 = cellfun(@transpose, normalR_norm1_RS_P1, 'UniformOutput', 0);
omitR_norm1_all_D1 = cellfun(@transpose, omitR_norm1_all_D1, 'UniformOutput', 0);
omitR_norm1_RS_D1 = cellfun(@transpose, omitR_norm1_RS_D1, 'UniformOutput', 0);
normalCue_norm1_all_P1 = cellfun(@transpose, normalCue_norm1_all_P1, 'UniformOutput', 0);
normalCue_norm1_RS_P1 = cellfun(@transpose, normalCue_norm1_RS_P1, 'UniformOutput', 0);
omitCue_norm1_all_P1 = cellfun(@transpose, omitCue_norm1_all_P1, 'UniformOutput', 0);
omitCue_norm1_RS_P1 = cellfun(@transpose, omitCue_norm1_RS_P1, 'UniformOutput', 0);
spont_norm1_all_D1 = cellfun(@transpose, spont_norm1_all_D1, 'UniformOutput', 0);
spont_norm1_RS_D1 = cellfun(@transpose, spont_norm1_RS_D1, 'UniformOutput', 0);
spont_norm1_all_P1 = cellfun(@transpose, spont_norm1_all_P1, 'UniformOutput', 0);
spont_norm1_RS_P1 = cellfun(@transpose, spont_norm1_RS_P1, 'UniformOutput', 0);

frame_size = cell2mat(cellfun(@size, spont_norm1_all_D1, 'UniformOutput', 0));
max_frame = max(frame_size(2:2:end));

all_normalR_norm1_D1 = interp_frame(normalR_norm1_all_D1, max_frame);
RS_normalR_norm1_D1    = interp_frame(normalR_norm1_RS_D1, max_frame);
all_normalR_norm1_P1 = interp_frame(normalR_norm1_all_P1, max_frame);
RS_normalR_norm1_P1    = interp_frame(normalR_norm1_RS_P1, max_frame);
all_omitR_norm1_D1   = interp_frame(omitR_norm1_all_D1, max_frame);
RS_omitR_norm1_D1 = interp_frame(omitR_norm1_RS_D1, max_frame);
all_normalCue_norm1_P1   = interp_frame(normalCue_norm1_all_P1, max_frame);
RS_normalCue_norm1_P1 = interp_frame(normalCue_norm1_RS_P1, max_frame);
all_omitCue_norm1_P1   = interp_frame(omitCue_norm1_all_P1, max_frame);
RS_omitCue_norm1_P1 = interp_frame(omitCue_norm1_RS_P1, max_frame);
all_spont_norm1_D1    = interp_frame(spont_norm1_all_D1, max_frame);
RS_spont_norm1_D1   = interp_frame(spont_norm1_RS_D1, max_frame);
all_spont_norm1_P1    = interp_frame(spont_norm1_all_P1, max_frame);
RS_spont_norm1_P1   = interp_frame(spont_norm1_RS_P1, max_frame);

fig=figure;
normScale = max(all_spont_norm1_D1,[],2);
normScale = 1;
all_spont_renorm1_D1 = bsxfun(@rdivide, all_spont_norm1_D1, normScale);
all_spont_renorm1_P1 = bsxfun(@rdivide, all_spont_norm1_P1, normScale);
all_normalR_renorm1_D1 = bsxfun(@rdivide, all_normalR_norm1_D1, normScale);
all_normalR_renorm1_P1 = bsxfun(@rdivide, all_normalR_norm1_P1, normScale);
all_omitR_renorm1_D1 = bsxfun(@rdivide, all_omitR_norm1_D1, normScale);
all_normalCue_renorm1_P1 = bsxfun(@rdivide, all_normalCue_norm1_P1, normScale);
all_omitCue_renorm1_P1 = bsxfun(@rdivide, all_omitCue_norm1_P1, normScale);
RS_spont_renorm1_D1 = bsxfun(@rdivide, RS_spont_norm1_D1, normScale);
RS_spont_renorm1_P1 = bsxfun(@rdivide, RS_spont_norm1_P1, normScale);
RS_normalR_renorm1_D1 = bsxfun(@rdivide, RS_normalR_norm1_D1, normScale);
RS_normalR_renorm1_P1 = bsxfun(@rdivide, RS_normalR_norm1_P1, normScale);
RS_omitR_renorm1_D1 = bsxfun(@rdivide, RS_omitR_norm1_D1, normScale);
RS_normalCue_renorm1_P1 = bsxfun(@rdivide, RS_normalCue_norm1_P1, normScale);
RS_omitCue_renorm1_P1 = bsxfun(@rdivide, RS_omitCue_norm1_P1, normScale);

tth = [1-(size(all_spont_renorm1_D1,2)/2):(size(all_spont_renorm1_D1,2)/2)].*double(min((ifi)));
subplot(1,2,1)
errorbar(tth, nanmean(all_spont_renorm1_D1,1), nanstd(all_spont_renorm1_D1,[],1)./sqrt(sum(~isnan(all_normalR_norm1_D1(:,1)))), '-m')
hold on;
errorbar(tth, nanmean(all_normalCue_renorm1_P1,1), nanstd(all_normalCue_renorm1_P1,[],1)./sqrt(sum(~isnan(all_normalCue_renorm1_P1(:,1)))), '-g')
hold on;
errorbar(tth, nanmean(all_normalR_renorm1_D1,1), nanstd(all_normalR_renorm1_D1,[],1)./sqrt(sum(~isnan(all_normalR_renorm1_D1(:,1)))), '-k')
% hold on;
% errorbar(tth, nanmean(all_normalR_renorm_P1,1), std(all_normalR_renorm_P1,[],1)./sqrt(length(all_normalR_renorm_P1)), '-r')

xlim([ -500 500]);
xlabel('Time (ms)')
ylabel('dF/F')
axis square
title('All cells');
% title(['Average event- all cells- ' num2str(total_cells) ' cells'])

subplot(1,2,2)
% cell_use = find(~isnan(mean([RS_spont_renorm RS_success_renorm RS_fail_renorm RS_press_renorm],2)));
errorbar(tth, nanmean(RS_spont_renorm1_D1,1), nanstd(RS_spont_renorm1_D1,[],1)./sqrt(sum(~isnan(RS_spont_renorm1_D1(:,1)))), '-m')
hold on;
% errorbar(tth, nanmean(RS_spont_renorm1_P1,1), nanstd(RS_spont_renorm1_D1,[],1)./sqrt(sum(~isnan(RS_spont_renorm1_P1(:,1)))), '-r')
errorbar(tth, nanmean(RS_normalCue_renorm1_P1,1), nanstd(RS_normalCue_renorm1_P1,[],1)./sqrt(sum(~isnan(RS_normalCue_renorm1_P1(:,1)))), '-g')
hold on;
errorbar(tth, nanmean(RS_normalR_renorm1_D1,1), nanstd(RS_normalR_renorm1_D1,[],1)./sqrt(sum(~isnan(RS_normalR_renorm1_D1(:,1)))), '-k')

xlim([ -500 500]);
xlabel('Time (ms)')
ylabel('dF/F')
axis square
title('Responsive cells');
% title(['Average event- responsive cells- ' num2str(total_resp) ' cells'])
suptitle('Average event- Green: Day5 Cue; Black: Day1 Reward; Magenta: spont')
saveas(fig, [out_base 'Summary_event_TC_exptavg.fig']);
print([out_base 'Summary_event_TC_exptavg.eps'], '-depsc');
print([out_base 'Summary_event_TC_exptavg.pdf'], '-dpdf');

% %% plot events amplitude timecourse
% all_success_event = cellfun(@transpose, success_event_TC, 'UniformOutput', 0);
% all_fail_event = cellfun(@transpose, fail_event_TC, 'UniformOutput', 0);
% all_press_event = cellfun(@transpose, press_event_TC, 'UniformOutput', 0);
% all_spont_event = cellfun(@transpose, spont_event_TC, 'UniformOutput', 0);
% RS_success_event = cellfun(@transpose, success_event_TC_RS, 'UniformOutput', 0);
% RS_fail_event = cellfun(@transpose, fail_event_TC_RS, 'UniformOutput', 0);
% RS_press_event = cellfun(@transpose, press_event_TC_RS, 'UniformOutput', 0);
% RS_spont_event = cellfun(@transpose, spont_event_TC_RS, 'UniformOutput', 0);
% 
% frame_size = cell2mat(cellfun(@size, success_event_TC, 'UniformOutput', 0));
% max_frame = max(frame_size(2:2:end));
% 
% all_success_event = interp_frame(all_success_event, max_frame);
% all_fail_event    = interp_frame(all_fail_event, max_frame);
% all_spont_event    = interp_frame(all_spont_event, max_frame);
% all_press_event    = interp_frame(all_press_event, max_frame);
% RS_success_event   = interp_frame(RS_success_event, max_frame);
% RS_fail_event = interp_frame(RS_fail_event, max_frame);
% RS_spont_event = interp_frame(RS_spont_event, max_frame);
% RS_press_event = interp_frame(RS_press_event, max_frame);
% 
% cell_use = find(~isnan(mean([all_success_event all_fail_event all_spont_event all_press_event],1)));
% 
% % tth = [1-(sz1/2):(sz1/2)].*double(min((TC_ifi)));
% tth = [1-(size(all_success_event,2)/2):(size(all_success_event,2)/2)].*double(min((TC_ifi)));
% fig=figure;
% subplot(2,1,1)
% errorbar(tth, nanmean(all_success_event(cell_use,:),1), std(all_success_event(cell_use,:),[],1)./sqrt(length(cell_use)), '-k')
% hold on;
% errorbar(tth, nanmean(all_fail_event(cell_use,:),1), std(all_fail_event(cell_use,:),[],1)./sqrt(length(cell_use)), '-r')
% hold on;
% errorbar(tth, nanmean(all_press_event(cell_use,:),1), std(all_press_event(cell_use,:),[],1)./sqrt(length(cell_use)), '-c')
% hold on;
% errorbar(tth, nanmean(all_spont_event(cell_use,:),1), std(all_spont_event(cell_use,:),[],1)./sqrt(length(cell_use)), '-m')
% xlabel('Time (ms)')
% ylabel('dF/F')
% title(['Average event- all cells- ' num2str(length(cell_use)) ' cells'])
% 
% subplot(1,2,2)
% cell_use = find(~isnan(mean([RS_success_event RS_fail_event RS_press_event],2)));
% errorbar(tth, nanmean(RS_success_event(cell_use,:),1), std(RS_success_event(cell_use,:),[],1)./sqrt(length(cell_use)), '-k')
% hold on;
% errorbar(tth, nanmean(RS_fail_event(cell_use,:),1), std(RS_fail_event(cell_use,:),[],1)./sqrt(length(cell_use)), '-r')
% hold on;
% errorbar(tth, nanmean(RS_press_event(cell_use,:),1), std(RS_press_event(cell_use,:),[],1)./sqrt(length(cell_use)), '-c')
% hold on;
% errorbar(tth, nanmean(RS_spont_event(cell_use,:),1), std(RS_spont_event(cell_use,:),[],1)./sqrt(length(cell_use)), '-c')
% xlabel('Time (ms)')
% ylabel('dF/F')
% title(['Average event- responsive cells- ' num2str(length(cell_use)) ' cells'])
% suptitle('Average event- Black: success; Red: early; Cyan: press; Megneta: spont')
% print([out_base 'Summary_event_TC_exptavg.eps'], '-depsc');
% print([out_base 'Summary_event_TC_exptavg.pdf'], '-dpdf');