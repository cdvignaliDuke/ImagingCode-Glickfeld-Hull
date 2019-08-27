%script for plotting population summaries for the CRP data. 

clear
file_info_CRP_all;
out_base = fullfile('Z:\Analysis\Cue_reward_pairing_analysis\2P\PCA_ICA_recalibration\Summary_folder\');
crp_dir = 'Z:\Analysis\Cue_reward_pairing_analysis\2P\PCA_ICA_recalibration\';
ITI_dir = 'Z:\Analysis\Cue_reward_pairing_analysis\2P\PCA_ICA_recalibration\ITI_licking\';
lick_psth_dir = ['Z:\Analysis\Cue_reward_pairing_analysis\2P\PCA_ICA_recalibration\licking_PSTHs\'];
%load([out_base 'cell_count.mat']);
normalR_nevents_D1 = 0; spont_nevents_D1 = 0;
sess_subset =  [1:15];  % [1:14]%1:6%  [1,2,3,4,5,6,7,8,9,10,11,13,15]; %1:size(days_1,2) %
cell_cats = 'NR_Rew_resp_cells_pos'; %ITI_bout_resp_cells  ITI_lick_resp_cells  allresp_cells_pos
%      allresp_cells_neg  rew_cells_neg  cue_cells_neg rew_cells_pos  cue_cells_pos   OR_Rew_resp_cells_pos
ifi = 33;

%licking ind [-30:76]   cue onset @ 31
%Ca event ind [-31:61]   cue onset @ 32    all_NR_hist dim1=frames  dim2=cells
cue_onset_lick_ind = 31;
%load and group variables for DAY 1
pre_buffer_all = [];
post_buffer_all = [];
for id = sess_subset
    if isempty(days_1{id})
        continue
    end
    %determine session info     %set pathnames
    [dest_sub, dest_sub_spikes, iti_lick_dir,~,~,~,~,~,~] = get_sess_info(days_1{id}, runID, crp_dir, ITI_dir);
    
    %load variables
    all_cell_cats = load([dest_sub '_cell_categories.mat']);
    %load([dest_sub_spikes, '_cell_categories_spike.mat']);
    load([dest_sub_spikes '_evoked_events.mat']);
    load([dest_sub '_cue_movies_lick.mat']);
    load([dest_sub, '_exclude_trials.mat']);
    load([dest_sub, '_pre_cue_ex_cell.mat']);
    %load([iti_lick_dir, 'lick_resp_cells.mat']);
    
    pre_buffer_all = [pre_buffer_all, pre_buffer_cue];
    post_buffer_all = [post_buffer_all, post_buffer_cue];
    
    %Define lick responsive cells
    if exist('ITI_lick_resp_cells', 'var')
        output_cells = cell_cat_selector(cell_cats, all_cell_cats, ITI_lick_resp_cells, ITI_bout_resp_cells);
    else
        output_cells = cell_cat_selector(cell_cats, all_cell_cats, [], []);
    end
    RS_cells_D1{id} = [output_cells & ~pre_cue_ex_cell];
    nCells = length(output_cells);
    
    %average across trials for each cell
    RS_cell_num = 0;
    for cell_num=1:nCells
        %determine which trials to use for this cell
        this_NR_ex = logical(ones(1, size(normalCue(1).hist,2)));
        this_NR_ex([NR_ex_trials{cell_num}']) = 0;
        this_OR_ex = logical(ones(1, size(omitCue(1).hist,2)));
        this_OR_ex([OR_ex_trials{cell_num}']) = 0;
        NR_no_lick = NR_lick_info.no_lick_cue_to_500;
        OR_no_lick = OR_lick_info.no_lick_cue_to_500;
        
        %Group event hists for all NR/OR trials. all/RS. Aligned to cue. 
        normalR_hist_D1{id}(:,cell_num) = mean(normalCue(cell_num).hist(:,[this_NR_ex]),2).*(1000/double(ifi));
        omitR_hist_D1{id}(:,cell_num) = mean(omitCue(cell_num).hist(:,[this_OR_ex]),2).*(1000/double(ifi));
        normalR_hist_nolick_D1{id}(:,cell_num) = mean(normalCue(cell_num).hist(:,[this_NR_ex & NR_no_lick]),2).*(1000/double(ifi));
        omitR_hist_nolick_D1{id}(:,cell_num) = mean(omitCue(cell_num).hist(:,[this_OR_ex & OR_no_lick]),2).*(1000/double(ifi));
        
        %Group event hists for NO LICK NR/OR trials. all/RS
        if ismember(cell_num, find(RS_cells_D1{id}))
            RS_cell_num = RS_cell_num+1;
            normalR_hist_nolick_RS_D1{id}(:,RS_cell_num) = mean(normalCue(cell_num).hist(:,[this_NR_ex & NR_no_lick]),2).*(1000/double(ifi));
            omitR_hist_nolick_RS_D1{id}(:,RS_cell_num) = mean(omitCue(cell_num).hist(:,[this_OR_ex & OR_no_lick]),2).*(1000/double(ifi));
            normalR_hist_RS_D1{id}(:,RS_cell_num) = mean(normalCue(cell_num).hist(:,[this_NR_ex]),2).*(1000/double(ifi));
            omitR_hist_RS_D1{id}(:,RS_cell_num) = mean(omitCue(cell_num).hist(:,[this_OR_ex]),2).*(1000/double(ifi));
        end
    end
   
    %remove pre_cue_ex_cells from all cells condition since they are not useable
    normalR_hist_D1{id} = normalR_hist_D1{id}(:,[~pre_cue_ex_cell]);
    omitR_hist_D1{id} = omitR_hist_D1{id}(:,[~pre_cue_ex_cell]);
    normalR_hist_nolick_D1{id} = normalR_hist_nolick_D1{id}(:,[~pre_cue_ex_cell]);
    omitR_hist_nolick_D1{id} = omitR_hist_nolick_D1{id}(:,[~pre_cue_ex_cell]);
    
    %collect licking data
    normalR_licks_D1{id} = lick_trace_NR;
    omitR_licks_D1{id} = lick_trace_OR;
    normalR_licks_nolick_D1{id} = lick_trace_NR([find(NR_lick_info.no_lick_cue_to_500)],:);
    omitR_licks_nolick_D1{id} = lick_trace_OR([find(OR_lick_info.no_lick_cue_to_500)],:);
end
assert(min(pre_buffer_all)==max(pre_buffer_all));

%Day N
normalCue_nevents_P1 = 0;  spont_nevents_P1 =0;
for id = sess_subset %  [1:14]%  [1,2,3,4,5,6,7,8,9,10,11,13,15]; % 1:size(days_post,2) %1:6 %
    %determine session info     %set pathnames
    [dest_sub, dest_sub_spikes, iti_lick_dir,~,~,~,~,~,~] = get_sess_info(days_post{id}, runID, crp_dir, ITI_dir);
    
    %load variables
    all_cell_cats = load([dest_sub '_cell_categories.mat']);
    %load([dest_sub_spikes, '_cell_categories_spike.mat']);
    load([dest_sub_spikes '_evoked_events.mat']);
    load([dest_sub '_cue_movies_lick.mat']);
    load([dest_sub, '_exclude_trials.mat']);
    load([dest_sub, '_pre_cue_ex_cell.mat']);
    %load([iti_lick_dir, 'lick_resp_cells.mat']);
    load([dest_sub_spikes '_rew_event_trials.mat']);
    load([dest_sub_spikes '_cue_event_trials.mat']);
    
    %Define responsive cells
    if exist('ITI_lick_resp_cells', 'var')
        output_cells = cell_cat_selector(cell_cats, all_cell_cats, ITI_lick_resp_cells, ITI_bout_resp_cells);
    else
        output_cells = cell_cat_selector(cell_cats, all_cell_cats, [], []);
    end
    RS_cells_P1{id} = [output_cells & ~pre_cue_ex_cell];
    nCells = length(output_cells);
    num_OR= size(omitCue(1).hist,2);
    num_NR = size(normalCue(1).hist,2);
    
    %average across trials for each cell
    RS_cell_num = 0;
    for cell_num=1:nCells
        %determine which trials to use for this cell
        this_NR_ex = logical(ones(1, size(normalCue(1).hist,2)));
        this_NR_ex([NR_ex_trials{cell_num}']) = 0;
        this_OR_ex = logical(ones(1, size(omitCue(1).hist,2)));
        this_OR_ex([OR_ex_trials{cell_num}']) = 0;
        NR_no_lick = NR_lick_info.no_lick_cue_to_500;
        OR_no_lick = OR_lick_info.no_lick_cue_to_500;
        %only include trials with an event in the rew window
%          this_NR_ex = this_NR_ex & logical(NR_rew_event_trials(cell_num,:)); %NR_cue_event_trials
%          this_OR_ex = this_OR_ex & logical(OR_rew_event_trials(cell_num,:));  %OR_cue_event_trials
        
        %Group event hists for all NR/OR trials. all/RS. Aligned to cue. 
        normalR_hist_P1{id}(:,cell_num) = mean(normalCue(cell_num).hist(:,[this_NR_ex]),2).*(1000/double(ifi));
        omitR_hist_P1{id}(:,cell_num) = mean(omitCue(cell_num).hist(:,[this_OR_ex]),2).*(1000/double(ifi));
        normalR_hist_nolick_P1{id}(:,cell_num) = mean(normalCue(cell_num).hist(:,[this_NR_ex & NR_no_lick]),2).*(1000/double(ifi));
        omitR_hist_nolick_P1{id}(:,cell_num) = mean(omitCue(cell_num).hist(:,[this_OR_ex & OR_no_lick]),2).*(1000/double(ifi));
        
        %Group event hists for NO LICK NR/OR trials. all/RS
        if ismember(cell_num, find(RS_cells_P1{id}))
            RS_cell_num = RS_cell_num+1;
            normalR_hist_nolick_RS_P1{id}(:,RS_cell_num) = mean(normalCue(cell_num).hist(:,[this_NR_ex & NR_no_lick]),2).*(1000/double(ifi));
            omitR_hist_nolick_RS_P1{id}(:,RS_cell_num) = mean(omitCue(cell_num).hist(:,[this_OR_ex & OR_no_lick]),2).*(1000/double(ifi));
            normalR_hist_RS_P1{id}(:,RS_cell_num) = mean(normalCue(cell_num).hist(:,[this_NR_ex]),2).*(1000/double(ifi));
            omitR_hist_RS_P1{id}(:,RS_cell_num) = mean(omitCue(cell_num).hist(:,[this_OR_ex]),2).*(1000/double(ifi));
        end
    end
    
    %remove pre_cue_ex_cells from all cells condition since they are not useable
    normalR_hist_P1{id} = normalR_hist_P1{id}(:,[~pre_cue_ex_cell]);
    omitR_hist_P1{id} = omitR_hist_P1{id}(:,[~pre_cue_ex_cell]);
    normalR_hist_nolick_P1{id} = normalR_hist_nolick_P1{id}(:,[~pre_cue_ex_cell]);
    omitR_hist_nolick_P1{id} = omitR_hist_nolick_P1{id}(:,[~pre_cue_ex_cell]);
    
    %collect licking data
    normalR_licks_P1{id} = lick_trace_NR;
    omitR_licks_P1{id} = lick_trace_OR;
    normalR_licks_nolick_P1{id} = lick_trace_NR([find(NR_lick_info.no_lick_cue_to_500)],:);
    omitR_licks_nolick_P1{id} = lick_trace_OR([find(OR_lick_info.no_lick_cue_to_500)],:);
    
%     %collect most/least active OR/NR trials in rew/cue windows
%     if sum(RS_cells_P1{id})>1
%         [~,NR_rew_event_trials_sort_ind] = sort(sum(logical(NR_rew_event_trials([RS_cells_P1{id}],:))));
%         NR_rew_low_quart{id} = NR_rew_event_trials_sort_ind(1:round(num_NR/4));
%         NR_rew_high_quart{id} = NR_rew_event_trials_sort_ind(end-round(num_NR/4):end);
%         [~,OR_rew_event_trials_sort_ind] = sort(sum(logical(OR_rew_event_trials([RS_cells_P1{id}],:))));
%         OR_rew_low_quart{id} = OR_rew_event_trials_sort_ind(1:round(num_OR/4));
%         OR_rew_high_quart{id} = OR_rew_event_trials_sort_ind(end-round(num_OR/4):end);
%         [~,NR_cue_event_trials_sort_ind] = sort(sum(logical(NR_cue_event_trials([RS_cells_P1{id}],:))));
%         NR_cue_low_quart{id} = NR_cue_event_trials_sort_ind(1:round(num_NR/4));
%         NR_cue_high_quart{id} = NR_cue_event_trials_sort_ind(end-round(num_NR/4):end);
%         [~,OR_cue_event_trials_sort_ind] = sort(sum(logical(OR_cue_event_trials([RS_cells_P1{id}],:))));
%         OR_cue_low_quart{id} = OR_cue_event_trials_sort_ind(1:round(num_OR/4));
%         OR_cue_high_quart{id} = OR_cue_event_trials_sort_ind(end-round(num_OR/4):end);
%     end
end

%Unexpected Reward Condition
for id = 1:size(days_UR,2)
    %determine session info     %set pathnames
    [dest_sub, dest_sub_spikes, iti_lick_dir,~,~,~,~,~,~] = get_sess_info(days_UR{id}, runID, crp_dir, ITI_dir);
    
    %load data for this session
    all_cell_cats = load([dest_sub '_cell_categories.mat']);
    %load([dest_sub_spikes, '_cell_categories_spike.mat']);
    load([dest_sub_spikes '_evoked_events.mat']);
    load([dest_sub '_cue_movies_lick.mat']);
    load([dest_sub, '_exclude_trials.mat']);
    load([dest_sub, '_pre_cue_ex_cell.mat']);
    
    %fill cells with data
    output_cells = cell_cat_selector(cell_cats, all_cell_cats, NaN, NaN);
    RS_cells_PN1{id} = [output_cells & ~pre_cue_ex_cell];
    nCells = length(output_cells);
    
    %average across trials for each cell
    RS_cell_num = 0;
    for cell_num=1:nCells
        %determine which trials to use for this cell
        this_NR_ex = logical(ones(1, size(normalCue(1).hist,2)));
        this_NR_ex([NR_ex_trials{cell_num}']) = 0;
        this_UR_ex = logical(ones(1, size(unexpCue(1).hist,2)));
        this_UR_ex([UR_ex_trials{cell_num}']) = 0;
        NR_no_lick = NR_lick_info.no_lick_cue_to_500;
        UR_no_lick = UR_lick_info.no_lick_cue_to_500;
        
        %Group event hists for all NR/OR trials. all/RS. Aligned to cue. 
        normalR_hist_PN1{id}(:,cell_num) = mean(normalCue(cell_num).hist(:,[this_NR_ex]),2).*(1000/double(ifi));
        unexpR_hist_PN1{id}(:,cell_num) = mean(unexpCue(cell_num).hist(:,[this_UR_ex]),2).*(1000/double(ifi));
        normalR_hist_nolick_PN1{id}(:,cell_num) = mean(normalCue(cell_num).hist(:,[this_NR_ex & NR_no_lick]),2).*(1000/double(ifi));
        unexpR_hist_nolick_PN1{id}(:,cell_num) = mean(unexpCue(cell_num).hist(:,[this_UR_ex & UR_no_lick]),2).*(1000/double(ifi));
        
        %Group event hists for NO LICK NR/OR trials. all/RS
        if ismember(cell_num, find(RS_cells_PN1{id}))
            RS_cell_num = RS_cell_num+1;
            normalR_hist_nolick_RS_PN1{id}(:,RS_cell_num) = mean(normalCue(cell_num).hist(:,[this_NR_ex & NR_no_lick]),2).*(1000/double(ifi));
            unexpR_hist_nolick_RS_PN1{id}(:,RS_cell_num) = mean(unexpCue(cell_num).hist(:,[this_UR_ex & UR_no_lick]),2).*(1000/double(ifi));
            normalR_hist_RS_PN1{id}(:,RS_cell_num) = mean(normalCue(cell_num).hist(:,[this_NR_ex]),2).*(1000/double(ifi));
            unexpR_hist_RS_PN1{id}(:,RS_cell_num) = mean(unexpCue(cell_num).hist(:,[this_UR_ex]),2).*(1000/double(ifi));
        end
    end

    %remove pre_cue_ex_cells from all cells condition since they are not useable
    normalR_hist_PN1{id} = normalR_hist_PN1{id}(:,[~pre_cue_ex_cell]);
    unexpR_hist_PN1{id} = unexpR_hist_PN1{id}(:,[~pre_cue_ex_cell]);
    normalR_hist_nolick_PN1{id} = normalR_hist_nolick_PN1{id}(:,[~pre_cue_ex_cell]);
    unexpR_hist_nolick_PN1{id} = unexpR_hist_nolick_PN1{id}(:,[~pre_cue_ex_cell]);
    
    %collect licking data
    normalR_licks_PN1{id} = lick_trace_NR;
    unexpR_licks_PN1{id} = lick_trace_UR;
    normalR_licks_nolick_PN1{id} = lick_trace_NR([find(NR_lick_info.no_lick_cue_to_500)],:);
    unexpR_licks_nolick_PN1{id} = lick_trace_UR([find(UR_lick_info.no_lick_cue_to_500)],:);
end

%1000ms delay
for id = 1:size(days_1000,2)
    if isempty(days_1000{id})
        continue
    end
    %determine session info     %set pathnames
    [dest_sub, dest_sub_spikes, iti_lick_dir,~,~,~,~,~,~] = get_sess_info(days_1000{id}, runID, crp_dir, ITI_dir);
    
    %load data for this session
    %load([dest_sub_spikes, '_cell_categories_spike.mat']);
    load([dest_sub_spikes '_event_summary.mat']);
    all_cell_cats = load([dest_sub '_cell_categories.mat']);
    load([dest_sub '_cue_movies_lick.mat']);
    load([dest_sub, '_exclude_trials.mat']);
    load([dest_sub, '_pre_cue_ex_cell.mat']);
    
    %fill cells with data
    output_cells = cell_cat_selector(cell_cats, all_cell_cats, NaN, NaN);
    RS_cells_PN2{id} = [output_cells & ~pre_cue_ex_cell];
    nCells = length(output_cells);
    
    %average across trials for each cell
    RS_cell_num = 0;
    for cell_num=1:nCells
        %determine which trials to use for this cell
        this_NR_ex = logical(ones(1, size(normalCue(1).hist,2)));
        this_NR_ex([NR_ex_trials{cell_num}']) = 0;
        this_OR_ex = logical(ones(1, size(omitCue(1).hist,2)));
        this_OR_ex([OR_ex_trials{cell_num}']) = 0;
        NR_no_lick = NR_lick_info.no_lick_cue_to_500;
        OR_no_lick = OR_lick_info.no_lick_cue_to_500;
        
        %Group event hists for all NR/OR trials. all/RS. Aligned to cue. 
        normalR_hist_PN2{id}(:,cell_num) = mean(normalCue(cell_num).hist(:,[this_NR_ex]),2).*(1000/double(ifi));
        omitR_hist_PN2{id}(:,cell_num) = mean(omitCue(cell_num).hist(:,[this_OR_ex]),2).*(1000/double(ifi));
        normalR_hist_nolick_PN2{id}(:,cell_num) = mean(normalCue(cell_num).hist(:,[this_NR_ex & NR_no_lick]),2).*(1000/double(ifi));
        omitR_hist_nolick_PN2{id}(:,cell_num) = mean(omitCue(cell_num).hist(:,[this_OR_ex & OR_no_lick]),2).*(1000/double(ifi));
        
        %Group event hists for NO LICK NR/OR trials. all/RS
        if ismember(cell_num, find(RS_cells_P1{id}))
            RS_cell_num = RS_cell_num+1;
            normalR_hist_nolick_RS_PN2{id}(:,RS_cell_num) = mean(normalCue(cell_num).hist(:,[this_NR_ex & NR_no_lick]),2).*(1000/double(ifi));
            omitR_hist_nolick_RS_PN2{id}(:,RS_cell_num) = mean(omitCue(cell_num).hist(:,[this_OR_ex & OR_no_lick]),2).*(1000/double(ifi));
            normalR_hist_RS_PN2{id}(:,RS_cell_num) = mean(normalCue(cell_num).hist(:,[this_NR_ex]),2).*(1000/double(ifi));
            omitR_hist_RS_PN2{id}(:,RS_cell_num) = mean(omitCue(cell_num).hist(:,[this_OR_ex]),2).*(1000/double(ifi));
        end
    end
    
    %remove pre_cue_ex_cells from all cells condition since they are not useable
    normalR_hist_PN2{id} = normalR_hist_PN2{id}(:,[~pre_cue_ex_cell]);
    omitR_hist_PN2{id} = omitR_hist_PN2{id}(:,[~pre_cue_ex_cell]);
    normalR_hist_nolick_PN2{id} = normalR_hist_nolick_PN2{id}(:,[~pre_cue_ex_cell]);
    omitR_hist_nolick_PN2{id} = omitR_hist_nolick_PN2{id}(:,[~pre_cue_ex_cell]);
    
    %collect licking data
    normalR_licks_PN2{id} = lick_trace_NR;
    omitR_licks_PN2{id} = lick_trace_OR;
    normalR_licks_nolick_PN2{id} = lick_trace_NR([find(NR_lick_info.no_lick_cue_to_500)],:);
    omitR_licks_nolick_PN2{id} = lick_trace_OR([find(OR_lick_info.no_lick_cue_to_500)],:);
end

tth = [-pre_buffer_cue(1):post_buffer_cue].*double(min(ifi));
tth_lick = [-30:76].*double(min(ifi));

%% Plotting Firing Rate PSTHS

%all cell, all trials
output_fig = PSTH_summary_plotter(normalR_hist_D1, omitR_hist_D1, ...
    normalR_hist_P1, omitR_hist_P1, ...
    normalR_hist_PN1, unexpR_hist_PN1, ...
    normalR_hist_PN2, omitR_hist_PN2, ...
    tth, 'all', 'all', []); %out_base

%all cells, no-lick condition
output_fig = PSTH_summary_plotter(normalR_hist_nolick_D1, omitR_hist_nolick_D1, ...
    normalR_hist_nolick_P1, omitR_hist_nolick_P1, ...
    normalR_hist_nolick_PN1, unexpR_hist_nolick_PN1, ...
    normalR_hist_nolick_PN2, omitR_hist_nolick_PN2, ...
    tth, 'no-lick', 'all', []);

%resp cells, all trials  
output_fig = PSTH_summary_plotter(normalR_hist_RS_D1, omitR_hist_RS_D1, ...
    normalR_hist_RS_P1, omitR_hist_RS_P1, ...
    normalR_hist_RS_PN1, unexpR_hist_RS_PN1, ...
    normalR_hist_RS_PN2, omitR_hist_RS_PN2, ...
    tth, 'all', cell_cats, []);

%resp cells, no-lick condition
output_fig = PSTH_summary_plotter(normalR_hist_nolick_RS_D1, omitR_hist_nolick_RS_D1, ...
    normalR_hist_nolick_RS_P1, omitR_hist_nolick_RS_P1, ...
    normalR_hist_nolick_RS_PN1, unexpR_hist_nolick_RS_PN1, ...
    normalR_hist_nolick_RS_PN2, omitR_hist_nolick_RS_PN2, ...
    tth, 'no-lick', cell_cats, []);

%% for each session plot the mean PSTH and the PSTH for each cell across trials 
% pdf_dir = '\\crash.dhe.duke.edu\data\home\jake\Analysis\Cue_reward_pairing_analysis\2P\PCA_ICA_recalibration\Summary PDF\session_PSTH';
% plot_ops.x_axis = tth; plot_ops.cue = 0;   plot_ops.reward = 600;
% %Day 1 
% plot_ops.day = 'Day 1';
% for session_num = sess_subset
%     plot_ops.mouse = days_1{session_num}(strfind(days_1{session_num},'i'):end);
%     plot_session_PSTH(normalR_hist_D1{session_num}, omitR_hist_D1{session_num}, plot_ops, pdf_dir);
% end
% 
% %Day N+1 
% plot_ops.day = 'Post Learning';
% for session_num = sess_subset
%     plot_ops.mouse = days_post{session_num}(strfind(days_post{session_num},'i'):end);
%     plot_session_PSTH(normalR_hist_P1{session_num}, omitR_hist_P1{session_num}, plot_ops,  pdf_dir);
% end
% 
% %Unexpected Reward days_UR
% plot_ops.day = 'Unexp Reward';
% for session_num = 1:length(days_UR) 
%     plot_ops.mouse = days_UR{session_num}(strfind(days_UR{session_num},'i'):end);
%     plot_session_PSTH(normalR_hist_PN1{session_num}, unexpR_hist_PN1{session_num}, plot_ops, pdf_dir);
% end
% 
% %1000ms delay 
% plot_ops.cue = 0;   plot_ops.reward = 1100;  plot_ops.day = '1000ms';
% for session_num = 1:length(days_1000)
%     plot_ops.mouse = days_1000{session_num}(strfind(days_1000{session_num},'i'):end);
%     plot_session_PSTH(normalR_hist_PN2{session_num}, omitR_hist_PN2{session_num}, plot_ops, pdf_dir);
% end

%% Plotting Licking PSTHs

%All trials
output_fig = licking_summary_plotter(normalR_licks_D1, omitR_licks_D1, ...
    normalR_licks_P1, omitR_licks_P1, ...
    normalR_licks_PN1, unexpR_licks_PN1, ...
    normalR_licks_PN2, omitR_licks_PN2, ...
    tth_lick, 'all', []); %lick_psth_dir

%no-lick condition
output_fig = licking_summary_plotter(normalR_licks_nolick_D1, omitR_licks_nolick_D1, ...
    normalR_licks_nolick_P1, omitR_licks_nolick_P1, ...
    normalR_licks_nolick_PN1, unexpR_licks_nolick_PN1, ...
    normalR_licks_nolick_PN2, omitR_licks_nolick_PN2, ...
    tth_lick, 'no-lick', []);

%% plot licking PSTHs for only the most or least active cells in rew/cue windows for OR resp cells
% normalR_licks_P12 = normalR_licks_P1;
% omitR_licks_P12 = omitR_licks_P1;
% % NR_rew_low_quart{sess_num}  NR_rew_high_quart{sess_num}  NR_cue_high_quart{sess_num}   NR_cue_low_quart{sess_num}
% % OR_rew_low_quart{sess_num}  OR_rew_high_quart{sess_num}  OR_cue_high_quart{sess_num}   OR_cue_low_quart{sess_num} 
% for sess_num = 1:length(normalR_licks_P1)
%     normalR_licks_P12{sess_num} = normalR_licks_P1{sess_num}([NR_cue_low_quart{sess_num}],:);
%     omitR_licks_P12{sess_num} = omitR_licks_P1{sess_num}([OR_cue_low_quart{sess_num}],:);
% end
% 
% output_fig = licking_summary_plotter(normalR_licks_D1, omitR_licks_D1, ...
%     normalR_licks_P12, omitR_licks_P12, ...
%     normalR_licks_PN1, unexpR_licks_PN1, ...
%     normalR_licks_PN2, omitR_licks_PN2, ...
%     tth_lick, 'all', []); %lick_psth_dir

%% Find probability of a Ca event in REW or CUE windows

figure;
%use the same response windows as 
cell_rew_prob = [];
cell_cue_prob = [];
rew_window = [ pre_buffer_cue+1+round(1080/double(ifi))+1 : pre_buffer_cue+1+round(1080/double(ifi))+1+floor(380/ifi)]; 
cue_window = [pre_buffer_cue+1+round(200/ifi): pre_buffer_cue+1+floor(580/ifi)];
for sess_num = 1:length(normalR_hist_RS_P1)   %normalR_hist_P1 normalR_hist_D1
    %check for empty session
    if isempty(normalR_hist_RS_P1{sess_num})
        continue
    end
    this_normalR_hist_PN2 = normalR_hist_RS_P1{sess_num}.*(double(ifi)/1000);
    
    %sum the probability of a Ca event for each frame in the window
    this_rew_prob = sum(this_normalR_hist_PN2([rew_window],:));
    this_cue_prob = sum(this_normalR_hist_PN2([cue_window],:));
    cell_rew_prob = [cell_rew_prob, this_rew_prob];
    cell_cue_prob = [cell_cue_prob, this_cue_prob];
end

subplot(1,2,2);
scatter(cell_rew_prob, cell_cue_prob); hold on;
plot([0 1], [0 1], 'k');
errorbarxy(mean(cell_rew_prob), mean(cell_cue_prob), std(cell_rew_prob)/sqrt(length(cell_rew_prob)), std(cell_cue_prob)/sqrt(length(cell_cue_prob)))
xlabel('Reward Window');
ylabel('Cue Window');
title('Post-learning: NRrew+ Responsive cells: probability of a calcium event');
xlim([0 1]);
ylim([0 1]);

suptitle('post learning');


%% plot mean event prob across neurons for pre vs post learning in different conditions

figure;
D1_this_psth_means = []; 
PL_this_psth_means = [];
for sess_num = 1:15  
    %check for empty session
    if isempty(normalR_hist_RS_P1{sess_num}) | isempty(normalR_hist_RS_D1{sess_num})
        continue
    end
    D1_this_psth = normalR_hist_RS_D1{sess_num}.*(double(ifi)/1000);
    PL_this_psth = normalR_hist_RS_P1{sess_num}.*(double(ifi)/1000);
    
    D1_this_psth = sum(D1_this_psth([cue_window],:));
    PL_this_psth = sum(PL_this_psth([cue_window],:));
    
    D1_this_psth_means(end+1) = mean(D1_this_psth);
    PL_this_psth_means(end+1) = mean(PL_this_psth);
    
end

scatter(D1_this_psth_means, PL_this_psth_means); hold on;
plot([0, 1], [0,1], 'k');
xlim([0 1]);
ylim([0 1]);






