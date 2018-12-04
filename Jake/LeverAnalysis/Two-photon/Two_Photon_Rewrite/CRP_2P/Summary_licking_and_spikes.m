%script for measuring onset of Ca events vs licking relative to cue/reward in the CRP data. 

clear
file_info_CRP_all;
out_base = fullfile('Z:\Analysis\Cue_reward_pairing_analysis\2P\PCA_ICA_recalibration\Summary_folder\first_last_quart_licks');
crp_dir = 'Z:\Analysis\Cue_reward_pairing_analysis\2P\PCA_ICA_recalibration\';
%load([out_base 'cell_count.mat']);
normalR_nevents_D1 = 0; spont_nevents_D1 = 0;
sess_subset = [1:15];  % [1:14]%1:6%  [1,2,3,4,5,6,7,8,9,10,11,13,15]; %1:size(days_1,2) %
ifi = 33;
pre_frames=30;
post_frames=76;
%licking ind [-30:76]   cue onset @ 31
%Ca event ind [-31:61]   cue onset @ 32    all_NR_hist dim1=frames  dim2=cells
x_axis = [-pre_frames:post_frames]*ifi;
cue_fr_ca_hist = 32;
cue_fr_lick_hist = 31;
peak_FR_win_D1 = cue_fr_ca_hist+1:cue_fr_ca_hist+round(1500/ifi);  %0:1500ms following cue onset
peak_FR_win_PL = cue_fr_ca_hist+1:cue_fr_ca_hist+round(900/ifi);  %0:1500ms following cue onset

%load and group variables for DAY 1
bout_aligned_ca_hist_D1 = [];
bout_aligned_ca_hist_sem_D1 = [];
for id = sess_subset
    %determine session info  %set pathnames
    if isempty(days_1{id})
        continue
    end
    [dest_sub, dest_sub_spikes, ~, ~, session, session_date, mouse_num, mouse_ID, rID] ...
        = get_sess_info(days_1{id}, runID, crp_dir, []);
    
    %load variables
    load([dest_sub '_cue_movies_lick.mat']);
    load([dest_sub '\spike_outputs\_event_hist.mat']);
    load([dest_sub_spikes '_evoked_events.mat']);
    load([dest_sub '_cell_categories.mat'], 'allresp_cells');
    assert(pre_buffer_cue == 31); assert(pre_frames==30);
    RS_ind = find(allresp_cells);
    
    %LICKING Group event hists for all NR/OR trials. 
    NR_lick_hist_D1_all{id} = lick_trace_NR;
    %omitR_hist_D1_licks{id} = lick_trace_OR;
    
    %LICKING Group event hists for NO LICK NR/OR trials. all/RS
    NR_lick_hist_D1_nolick{id} = lick_trace_NR(logical(NR_lick_info.no_lick_cue_to_500),:);
    %omitR_hist_nolick_D1{id} = lick_trace_OR(logical(OR_lick_info.no_lick_cue_to_500),:);
    
    %Find the mean lick bout onset time on a trial by trial basis 
    bout_onset_D1_mean(id) = nanmean(NR_lick_info.bout_onset) - cue_fr_lick_hist;
    bout_onset_D1_sem(id) = (nanstd(NR_lick_info.bout_onset) - cue_fr_lick_hist)/sqrt(sum(~isnan(NR_lick_info.bout_onset)));
    
    %Find the time of the peak of the mean firing rate 
    all_NR_hist_avg = mean(all_NR_hist(:,[RS_ind]),2)';
    NR_t_peak_FR_D1_all(id) =  find( all_NR_hist_avg(peak_FR_win_D1) == max(all_NR_hist_avg(peak_FR_win_D1)) );  %find index relative to cue onset frame. idx=1 is the frame after the cue frame (33ms)
    
    %Ca EVENTS group hists for all NR trials. All/RS
    NR_sp_hist_D1_all{id} = all_NR_hist(:,[RS_ind]);  %all_NR_hist has a single mean histogram for each cell   2D
    
    %Find mean time of first Ca Event on a trial by trial basis
    
    %adjust indeces of calcium events to get histograms aligned to lick bout onset.
    all_trial_hist = [];
    for trial_num = 1:size(NR_lick_info.bout_onset,2)
        if NR_lick_info.bout_onset(trial_num) == NaN
            continue
        else
%         this_lick_ind = find(lick_trace_NR(trial_num,cue_fr_lick_hist+1:end),1,'first'); %index of the first lick after the cue (relative to the cue)
%         if isempty(this_lick_ind)
%             continue
%         else
            this_bout_ind = NR_lick_info.bout_onset(trial_num) - cue_fr_lick_hist;  %ind of zero = cue onset
            this_trial_hist = zeros(1, length([-pre_buffer_cue:post_buffer_cue]));  %histogram of Ca spiking across cells for this trial
            for cell_num = [RS_ind]
                this_trial_ind = normalCue(cell_num).ind{trial_num}'-this_bout_ind;   %construct the histogram so that bout onset occurs when cue onset would have occurred.  cue_fr_ca_hist = 32;
                %this_trial_ind = normalCue(cell_num).ind{trial_num}'-this_lick_ind;
                this_trial_ind = this_trial_ind(find(this_trial_ind>0 & this_trial_ind< length([-pre_buffer_cue:post_buffer_cue]) ));
                this_trial_hist_temp = zeros(1, length([-pre_buffer_cue:post_buffer_cue]));
                this_trial_hist_temp(this_trial_ind) = 1;
                this_trial_hist = this_trial_hist + this_trial_hist_temp;
            end
            all_trial_hist = cat(1, all_trial_hist, this_trial_hist);  %for each trial get the histogram of Ca events across all cells
        end
    end
    bout_aligned_ca_hist_D1 = [bout_aligned_ca_hist_D1; mean(all_trial_hist,1)];
    bout_aligned_ca_hist_sem_D1 = [bout_aligned_ca_hist_sem_D1; std(all_trial_hist,1)/sqrt(size(all_trial_hist,1))];
end

%load and group variables for Post Learning
bout_aligned_ca_hist_PL = [];
bout_aligned_ca_hist_sem_PL = [];
cue_aligned_ca_hist_PL = [];
cue_aligned_ca_hist_sem_PL = [];
f_quart_bout_ca_hist_PL = [];
l_quart_bout_ca_hist_PL = [];
f_quart_bout_ca_hist_PL_OR = [];
l_quart_bout_ca_hist_PL_OR = [];
f_quart_FL_NR_PL_Ca_hist = [];
l_quart_FL_NR_PL_Ca_hist = [];
cue_FL_ca_hist_PL = [];
rew_FL_ca_hist_PL = [];
cue_FL_lick_hists_NR_PL = nan(length(sess_subset), 93);;
rew_FL_lick_hists_NR_PL = nan(length(sess_subset), 93);;
first_quart_lick_hists_NR_PL = zeros(length(sess_subset), 93);
last_quart_lick_hists_NR_PL = zeros(length(sess_subset), 93);
first_quart_lick_hists_OR_PL = zeros(length(sess_subset), 93);
last_quart_lick_hists_OR_PL = zeros(length(sess_subset), 93);
NR_FL_f_quart_hist = zeros(length(sess_subset), 93);
NR_FL_l_quart_hist = zeros(length(sess_subset), 93);
NR_FL_PL_late_or_no_lick_hist = zeros(length(sess_subset), 93);

%POST LEARNING
for id = sess_subset
    %determine session info  set pathnames
    [dest_sub, dest_sub_spikes, ~, ~, session, session_date, mouse_num, mouse_ID, rID] ...
        = get_sess_info(days_post{id}, runID, crp_dir, []);
    
    %load variables
    load([dest_sub '_cue_movies_lick.mat']);
    load([dest_sub_spikes '_event_hist.mat']);
    load([dest_sub_spikes '_evoked_events.mat']);
    load([dest_sub, '_pre_cue_ex_cell.mat']);
    %load([dest_sub '_cell_categories.mat'], 'allresp_cells');
    all_cell_cats = load([dest_sub '_cell_categories.mat']);
    assert(pre_buffer_cue == 31); assert(pre_frames==30);
    %RS_ind = find(allresp_cells);
    cell_cats = 'OR_Rew_resp_cells_pos';
    output_cells = cell_cat_selector(cell_cats, all_cell_cats, [], []);  %need to replace RS_ind with output_cells===============================================
    RS_ind = find([output_cells & ~pre_cue_ex_cell]);
    if isempty(RS_ind)
        continue
    end
    
    %LICKING Group event hists for all NR/OR trials. 
    NR_lick_hist_PL_all{id} = lick_trace_NR;
    %omitR_hist_D1_licks{id} = lick_trace_OR;
    
    %LICKING Group event hists for NO LICK NR/OR trials. all/RS
    NR_lick_hist_PL_nolick{id} = lick_trace_NR(logical(NR_lick_info.no_lick_cue_to_500),:);
    %omitR_hist_nolick_PL{id} = lick_trace_OR(logical(OR_lick_info.no_lick_cue_to_500),:);
    
    %Find the mean lick bout onset time on a trial by trial basis 
    bout_onset_PL_mean(id) = nanmean(NR_lick_info.bout_onset) - cue_fr_lick_hist;
    bout_onset_PL_sem(id) = (nanstd(NR_lick_info.bout_onset) - cue_fr_lick_hist)/sqrt(sum(~isnan(NR_lick_info.bout_onset)));
    
    %Find the time of the peak of the mean firing rate 
    all_NR_hist_avg = mean(all_NR_hist(:,[RS_ind]),2)';  %
    NR_t_peak_FR_PL_all(id) =  find( all_NR_hist_avg(peak_FR_win_PL) == max(all_NR_hist_avg(peak_FR_win_PL)), 1, 'first');  %find index relative to cue onset frame. idx=1 is the frame after the cue frame (33ms)
    
    %Ca EVENTS group hists for all NR trials. All/RS
    NR_sp_hist_PL_all{id} = all_NR_hist(:,[RS_ind]);  %all_NR_hist has a single mean histogram for each cell 2D 
    
    %index of the first lick after the cue (relative to the cue)  cue=0
    first_lick_ind = NaN(1,size(lick_trace_NR,1));
    no_lick_trial_num = [];
    for trial_num = 1:size(lick_trace_NR,1)
            this_lick_ind = find(lick_trace_NR(trial_num,cue_fr_lick_hist+1+round(200/ifi):end),1,'first');  %start looking for licks 150ms after the first frame after the cue
            if isempty(this_lick_ind)
                no_lick_trial_num = [no_lick_trial_num, trial_num];
            else
                first_lick_ind(trial_num) = this_lick_ind;
            end
    end
    first_lick_ind = first_lick_ind + round(200/ifi);
    
    %+======================================================================================================
    %find trials with no licks for 1200ms after cue
    no_lick_int = round(1000/ifi);
    late_licks_trials = find(first_lick_ind>=no_lick_int);
    NR_FL_PL_late_and_no_lick_trials = unique(sort([late_licks_trials, no_lick_trial_num]));
    %get licking histograms of late/no lick trials
    NR_FL_PL_late_or_no_lick_hist(id,[2:end]) = mean(lick_trace_NR([NR_FL_PL_late_and_no_lick_trials],[1:92])) *33/10; 
    %get Firing rrate histograms of late/no lick trials 
    late_no_ca_hist = [];
    for trial_num = NR_FL_PL_late_and_no_lick_trials;
        this_trial_hist = zeros(1, length([-pre_buffer_cue:post_buffer_cue]));
        for cell_num = [RS_ind]
            this_trial_ind = normalCue(cell_num).ind{trial_num}';
            this_trial_ind = this_trial_ind(find(this_trial_ind>0 & this_trial_ind< length([-pre_buffer_cue:post_buffer_cue]) ));
            this_trial_hist(this_trial_ind) =  this_trial_hist(this_trial_ind) +1;
        end
        late_no_ca_hist = cat(1, late_no_ca_hist, this_trial_hist);
    end
    NR_FL_PL_late_no_lick_ca_hist{id}= late_no_ca_hist;
    %+======================================================================================================
    
    %% get calcium hists of specific trials based on lick time
    
    % SINGLE LICKS - find earliest 1/4 and latest 1/4 of single licks
    num_first_licks = sum(~isnan(first_lick_ind));
    [~,ind_of_NR_FL_inds] = sort(first_lick_ind);
    f_quart_FL_NR_ind = ind_of_NR_FL_inds(1:round(num_first_licks/4));
    l_quart_FL_NR_ind =  ind_of_NR_FL_inds(end-round(num_first_licks/4):end);
    %store time of first lick
    NR_FL_times_f_quart{id} = first_lick_ind( f_quart_FL_NR_ind )*ifi; 
    NR_FL_times_l_quart{id} = first_lick_ind( l_quart_FL_NR_ind )*ifi;
    %get lick histograms
    NR_FL_f_quart_hist(id,[2:end]) = mean(lick_trace_NR([f_quart_FL_NR_ind],[1:92])) *33/10; 
    NR_FL_l_quart_hist(id,[2:end]) = mean(lick_trace_NR([l_quart_FL_NR_ind],[1:92])) *33/10; 
    
    %get Ca Histogram of trials with earliest 1/4 of SINGLE LICKS
    this_quart_hist = [];
    for cell_num = [RS_ind] %create a matrix of mean hists for each cell dim1=cells dim2=frames
        this_quart_hist = [this_quart_hist; mean(normalCue(cell_num).hist(:,[f_quart_FL_NR_ind]),2)'*30]; %gets the average number of events in a given frame num relative to cue for a subset of trials
    end
    f_quart_FL_NR_PL_Ca_hist = [f_quart_FL_NR_PL_Ca_hist; this_quart_hist];  % need f_quart_FL_NR_PL_Ca_hist to include all cells 
    
    %get Ca Histogram of trials with earliest 1/4 of SINGLE LICKS
    this_quart_hist = [];
    for cell_num = [RS_ind]   
         this_quart_hist = [this_quart_hist; mean(normalCue(cell_num).hist(:,[l_quart_FL_NR_ind]),2)'*30]; 
    end
    l_quart_FL_NR_PL_Ca_hist = [l_quart_FL_NR_PL_Ca_hist; this_quart_hist];
    
    % LICK BOUTS - find earliest 1/4 and latest 1/4 of lick bouts
    NR_bouts_ind = NR_lick_info.bout_onset;
    num_NR_bouts = sum(~isnan(NR_bouts_ind));
    [NR_bouts_ind_sort, sorted_ind] = sort(NR_bouts_ind);
    sorted_ind = sorted_ind(1:num_NR_bouts); %get rid of the NaNs
    first_quart_bout_ind = sorted_ind(1:round(num_NR_bouts/4));
    last_quart_bout_ind = sorted_ind(end-round(num_NR_bouts/4):end);
    %do the same for OR trials but for 1/3 of bouts
    OR_bouts_ind = OR_lick_info.bout_onset;
    num_OR_bouts = sum(~isnan(OR_bouts_ind));
    [OR_bouts_ind_sort, sorted_ind] = sort(OR_bouts_ind);
    sorted_ind = sorted_ind(1:num_OR_bouts); %get rid of the NaNs
    first_quart_bout_ind_OR = sorted_ind(1:round(num_OR_bouts/3));
    last_quart_bout_ind_OR = sorted_ind(end-round(num_OR_bouts/3):end);
    %store times of bout onsets
    first_quart_bout_times_NR_PL{id} = (NR_bouts_ind(first_quart_bout_ind) - (pre_buffer_cue +1) )*ifi;
    last_quart_bout_times_NR_PL{id}  = (NR_bouts_ind(last_quart_bout_ind)  - (pre_buffer_cue +1) )*ifi;
    first_quart_bout_times_OR_PL{id} = (OR_bouts_ind(first_quart_bout_ind_OR) - (pre_buffer_cue +1) )*ifi;
    last_quart_bout_times_OR_PL{id}  = (OR_bouts_ind(last_quart_bout_ind_OR)  - (pre_buffer_cue +1) )*ifi;
    
    %get lick histograms
    first_quart_lick_hists_NR_PL(id,[2:end]) = mean(lick_trace_NR([first_quart_bout_ind],[1:92])) *33/10;  %adjust to Hz and then divide by ten so it fits on the plot
    last_quart_lick_hists_NR_PL(id,[2:end]) = mean(lick_trace_NR([last_quart_bout_ind],[1:92])) *33/10;
    %get lick histograms
    first_quart_lick_hists_OR_PL(id,[2:end]) = mean(lick_trace_OR([first_quart_bout_ind_OR],[1:92])) *33/10;  %adjust to Hz and then divide by ten so it fits on the plot
    last_quart_lick_hists_OR_PL(id,[2:end]) = mean(lick_trace_OR([last_quart_bout_ind_OR],[1:92])) *33/10;
    
    % find Ca histogram of FIRST quarter of lick bouts
    this_quart_hist = [];
    this_quart_hist_OR = [];
    for cell_num = [RS_ind]   
        this_quart_hist = [this_quart_hist; mean(normalCue(cell_num).hist(:,[first_quart_bout_ind]),2)'*30]; 
        this_quart_hist_OR = [this_quart_hist_OR; mean(omitCue(cell_num).hist(:,[first_quart_bout_ind_OR]),2)'*30];
    end
    f_quart_bout_ca_hist_PL = [f_quart_bout_ca_hist_PL; this_quart_hist];
    f_quart_bout_ca_hist_PL_OR = [f_quart_bout_ca_hist_PL_OR; this_quart_hist_OR];
    
    % find Ca histogram of LAST quarter of lick bouts
    this_quart_hist = [];
    this_quart_hist_OR = [];
    for cell_num = [RS_ind]  
        this_quart_hist = [this_quart_hist; mean(normalCue(cell_num).hist(:,[last_quart_bout_ind]),2)'*30];
        this_quart_hist_OR = [this_quart_hist_OR; mean(omitCue(cell_num).hist(:,[last_quart_bout_ind_OR]),2)'*30];
    end
    l_quart_bout_ca_hist_PL = [l_quart_bout_ca_hist_PL; this_quart_hist];   
    l_quart_bout_ca_hist_PL_OR = [l_quart_bout_ca_hist_PL_OR; this_quart_hist_OR];  
    
        
    %CUE vs REW FLs   find trials first lick within a given cue window and reward window
    FL_cue_trial_num = first_lick_ind([find(first_lick_ind>=37 & first_lick_ind<=47)]);
    FL_rew_trial_num = first_lick_ind([find(first_lick_ind>=51 & first_lick_ind<=61)]);
    %get lick histograms
    if ~isempty(FL_cue_trial_num)
        cue_FL_lick_hists_NR_PL(id,[2:end]) = mean(lick_trace_NR([FL_cue_trial_num],[1:92])) *33/10;  %adjust to Hz and then divide by ten so it fits on the plot
        % find Ca histogram of CUE FLs
        this_win_hist = [];
        for cell_num = [RS_ind]
            this_win_hist = [this_win_hist; mean(normalCue(cell_num).hist(:,[FL_cue_trial_num]),2)'*30];
        end
        cue_FL_ca_hist_PL = [cue_FL_ca_hist_PL; this_win_hist];
    end
    if ~isempty(FL_rew_trial_num)
        rew_FL_lick_hists_NR_PL(id,[2:end]) = mean(lick_trace_NR([FL_rew_trial_num],[1:92])) *33/10;
        % find Ca histogram of REW FLs
        this_win_hist = [];
        for cell_num = [RS_ind]
            this_win_hist = [this_win_hist; mean(normalCue(cell_num).hist(:,[FL_rew_trial_num]),2)'*30];
        end
        rew_FL_ca_hist_PL = [rew_FL_ca_hist_PL; this_win_hist];
    end
    
    %% adjust indeces of calcium events to get histograms aligned to lick bout onset. 
    all_trial_hist = [];
    all_trial_hist_cue = [];
    for trial_num = 1:size(NR_lick_info.bout_onset,2)
        if isnan(NR_lick_info.bout_onset(trial_num)) 
            continue
        else
%         this_lick_ind = find(lick_trace_NR(trial_num,cue_fr_lick_hist+1:end),1,'first'); %index of the first lick after the cue (relative to the cue)
%         if isempty(this_lick_ind)
%             continue
%         else

            %bout aligned
            this_bout_ind = NR_lick_info.bout_onset(trial_num) - cue_fr_lick_hist;  %ind of zero = cue onset
            this_trial_hist = zeros(1, length([-pre_buffer_cue:post_buffer_cue]));  %histogram of Ca spiking across cells for this trial
            for cell_num = [RS_ind]
                this_trial_ind = normalCue(cell_num).ind{trial_num}'-this_bout_ind;   %construct the histogram so that bout onset occurs when cue onset would have occurred.  cue_fr_ca_hist = 32;
                %this_trial_ind = normalCue(cell_num).ind{trial_num}'-this_lick_ind;
                this_trial_ind = this_trial_ind(find(this_trial_ind>0 & this_trial_ind< length([-pre_buffer_cue:post_buffer_cue]) ));
                this_trial_hist_temp = zeros(1, length([-pre_buffer_cue:post_buffer_cue]));
                this_trial_hist_temp(this_trial_ind) = 1;
                this_trial_hist = this_trial_hist + this_trial_hist_temp;
            end
            all_trial_hist = cat(1, all_trial_hist, this_trial_hist);  %for each trial get the histogram of Ca events across all cells
            
            %cue aligned
            this_trial_hist_cue = zeros(1, length([-pre_buffer_cue:post_buffer_cue]));
            for cell_num = [RS_ind]
                this_trial_ind = normalCue(cell_num).ind{trial_num}';
                this_trial_ind = this_trial_ind(find(this_trial_ind>0 & this_trial_ind< length([-pre_buffer_cue:post_buffer_cue]) ));
                this_trial_hist_cue_temp = zeros(1, length([-pre_buffer_cue:post_buffer_cue]));
                this_trial_hist_cue_temp(this_trial_ind) = 1;
                this_trial_hist_cue = this_trial_hist_cue + this_trial_hist_cue_temp;
            end
            all_trial_hist_cue = cat(1, all_trial_hist_cue, this_trial_hist_cue);
        end
    end
    bout_aligned_ca_hist_PL = [bout_aligned_ca_hist_PL; mean(all_trial_hist,1)];
    bout_aligned_ca_hist_sem_PL = [bout_aligned_ca_hist_sem_PL; std(all_trial_hist,1)/sqrt(size(all_trial_hist,1))];
    cue_aligned_ca_hist_PL = [cue_aligned_ca_hist_PL; mean(all_trial_hist_cue,1)];
    cue_aligned_ca_hist_sem_PL = [cue_aligned_ca_hist_sem_PL; std(all_trial_hist_cue,1)/sqrt(size(all_trial_hist_cue,1))];
end

%calculate sem for each variable
f_quart_FL_NR_PL_Ca_hist_sem = [];
l_quart_FL_NR_PL_Ca_hist_sem = [];
f_quart_bout_ca_hist_PL_sem = [];
l_quart_bout_ca_hist_PL_sem = [];
f_quart_bout_ca_hist_PL_sem_OR = [];
l_quart_bout_ca_hist_PL_sem_OR = [];
for fr_num=1:size(f_quart_FL_NR_PL_Ca_hist,2)
    f_quart_FL_NR_PL_Ca_hist_sem(fr_num) = std(f_quart_FL_NR_PL_Ca_hist([isfinite(f_quart_FL_NR_PL_Ca_hist(:,fr_num))],fr_num)) / sqrt(size(f_quart_FL_NR_PL_Ca_hist,1));
    l_quart_FL_NR_PL_Ca_hist_sem(fr_num) = std(l_quart_FL_NR_PL_Ca_hist([isfinite(l_quart_FL_NR_PL_Ca_hist(:,fr_num))],fr_num)) / sqrt(size(l_quart_FL_NR_PL_Ca_hist,1));
    f_quart_bout_ca_hist_PL_sem(fr_num) =  std(f_quart_bout_ca_hist_PL([isfinite(f_quart_bout_ca_hist_PL(:,fr_num))],fr_num)) / sqrt(size(f_quart_bout_ca_hist_PL,1));
    l_quart_bout_ca_hist_PL_sem(fr_num) = std(l_quart_bout_ca_hist_PL([isfinite(l_quart_bout_ca_hist_PL(:,fr_num))],fr_num)) / sqrt(size(l_quart_bout_ca_hist_PL,1));
    f_quart_bout_ca_hist_PL_sem_OR(fr_num) =  std(f_quart_bout_ca_hist_PL_OR([isfinite(f_quart_bout_ca_hist_PL_OR(:,fr_num))],fr_num)) / sqrt(size(f_quart_bout_ca_hist_PL_OR,1));
    l_quart_bout_ca_hist_PL_sem_OR(fr_num) = std(l_quart_bout_ca_hist_PL_OR([isfinite(l_quart_bout_ca_hist_PL_OR(:,fr_num))],fr_num)) / sqrt(size(l_quart_bout_ca_hist_PL_OR,1));
end  

%% PLOTTING 

% %mean peak FR and mean bout onset time day1 vs post learning 
% figure; scatter([(NR_t_peak_FR_D1_all*ifi)], [bout_onset_D1_mean*ifi]); hold on;
% scatter(mean(NR_t_peak_FR_D1_all*ifi), mean(bout_onset_D1_mean*ifi), 'filled', 'b');
% scatter([(NR_t_peak_FR_PL_all*ifi)], [bout_onset_PL_mean*ifi], 'r');
% scatter(mean(NR_t_peak_FR_PL_all*ifi), mean(bout_onset_PL_mean*ifi), 'filled', 'r');
% xlabel('time of peak FR relative to cue onset');
% ylabel('time of lick bout onset relative to cue onset');
% title('day 1: blue;   post-learning: red;  RS cells');
% plot([1:100:1500], [1:100:1500]);
% xlim([0 2000]); ylim([0 2000]);
% vline(580); hline(600);
% 
% %plot boutonset aligned histograms
% figure;
% x_axis = [-31:61]*ifi;
% shadedErrorBar( x_axis, mean(bout_aligned_ca_hist_D1), std(bout_aligned_ca_hist_D1,1)/sqrt(size(bout_aligned_ca_hist_D1,1)), 'k'); hold on;
% shadedErrorBar( x_axis, mean(bout_aligned_ca_hist_PL), std(bout_aligned_ca_hist_PL,1)/sqrt(size(bout_aligned_ca_hist_PL,1)), 'b');
% shadedErrorBar( x_axis, mean(cue_aligned_ca_hist_PL), std(cue_aligned_ca_hist_PL,1)/sqrt(size(cue_aligned_ca_hist_PL,1)), 'm');
% vline(0); xlim([-1000 1000]);
% xlabel('time (ms) from lick bout onset');
% ylabel('Ca event firing rate');
% title('Normal Rewarded trials with lick bouts; RS cells; black=day1,  blue=post-learning');

%----------------------------------------------
%LICK BOUTS NR plot earliest 1/4 of bout onsets vs slowest 1/4 of bout onsets
x_axis = [-31:61]*ifi;
figure;
subplot(2,1,1);
shadedErrorBar( x_axis, nanmean(f_quart_bout_ca_hist_PL), f_quart_bout_ca_hist_PL_sem, 'b'); hold on;
shadedErrorBar( x_axis, nanmean(l_quart_bout_ca_hist_PL), l_quart_bout_ca_hist_PL_sem, 'k');
xlabel('time (ms) from cue onset');
ylabel('Ca event firing rate');
title('Normal Rewarded trials with lick bouts; RS cells; All post-learning; blue= 1st 1/4,  black= last 1/4');
xlim([-1000 1500]); ylim([0 2]);

subplot(2,1,2);
plot( x_axis, nanmean(f_quart_bout_ca_hist_PL),  'b'); hold on;
plot( x_axis, nanmean(l_quart_bout_ca_hist_PL),  'k');
xlim([-1000 1500]); ylim([0 2]);
xlabel('time (ms) from cue onset');
ylabel('Ca event firing rate');
title('Normal Rewarded trials with lick bouts; RS cells; All post-learning; blue= 1st 1/4,  black= last 1/4');
bar(x_axis, mean(first_quart_lick_hists_NR_PL), 'b');
bar(x_axis, mean(last_quart_lick_hists_NR_PL), 'k');
errorbar(x_axis, mean(first_quart_lick_hists_NR_PL), std(first_quart_lick_hists_NR_PL)/sqrt(size(first_quart_lick_hists_NR_PL,1)), 'b', 'LineStyle', 'none');
errorbar(x_axis, mean(last_quart_lick_hists_NR_PL), std(last_quart_lick_hists_NR_PL)/sqrt(size(last_quart_lick_hists_NR_PL,1)), 'k', 'LineStyle', 'none');

first_quart_bout_times_NR_PL_mat = [];
last_quart_bout_times_NR_PL_mat = [];
num_trials_per_quart = [];
for id=sess_subset
    first_quart_bout_times_NR_PL_mat = [first_quart_bout_times_NR_PL_mat, mean(first_quart_bout_times_NR_PL{id})];
    last_quart_bout_times_NR_PL_mat = [last_quart_bout_times_NR_PL_mat, mean(last_quart_bout_times_NR_PL{id})];
    num_trials_per_quart = [num_trials_per_quart,  size(last_quart_bout_times_NR_PL{id},2)];
end
[aa, bb] =get_mean_and_sem(first_quart_bout_times_NR_PL_mat)
[aa, bb] =get_mean_and_sem(last_quart_bout_times_NR_PL_mat)
[aa, bb] =get_mean_and_sem(num_trials_per_quart)

%LICK BOUTS OR 
figure;
subplot(2,1,1);
shadedErrorBar( x_axis, nanmean(f_quart_bout_ca_hist_PL_OR), f_quart_bout_ca_hist_PL_sem_OR, 'b'); hold on; %f_quart_bout_ca_hist_PL_sem_OR
shadedErrorBar( x_axis, nanmean(l_quart_bout_ca_hist_PL_OR), l_quart_bout_ca_hist_PL_sem_OR, 'k');
xlabel('time (ms) from cue onset');
ylabel('Ca event firing rate');
title('Omitted Reward trials with lick bouts; RS cells; All post-learning; blue= 1st 1/4,  black= last 1/4');
xlim([-1000 1500]); ylim([0 3]);

subplot(2,1,2);
plot( x_axis, nanmean(f_quart_bout_ca_hist_PL_OR),  'b'); hold on;
plot( x_axis, nanmean(l_quart_bout_ca_hist_PL_OR),  'k');
xlim([-1000 1500]); ylim([0 3]);
xlabel('time (ms) from cue onset');
ylabel('Ca event firing rate');
title('Omissed Reward trials with lick bouts; RS cells; All post-learning; blue= 1st 1/4,  black= last 1/4');
bar(x_axis, nanmean(first_quart_lick_hists_OR_PL), 'b');
bar(x_axis, mean(last_quart_lick_hists_OR_PL), 'k');
errorbar(x_axis, nanmean(first_quart_lick_hists_OR_PL), std(first_quart_lick_hists_OR_PL)/sqrt(size(first_quart_lick_hists_OR_PL,1)), 'b', 'LineStyle', 'none');
errorbar(x_axis, mean(last_quart_lick_hists_OR_PL), std(last_quart_lick_hists_OR_PL)/sqrt(size(last_quart_lick_hists_OR_PL,1)), 'k', 'LineStyle', 'none');

first_quart_bout_times_OR_PL_mat = [];
last_quart_bout_times_OR_PL_mat = [];
num_trials_per_third = [];
for id=sess_subset
    first_quart_bout_times_OR_PL_mat = [first_quart_bout_times_OR_PL_mat, mean(first_quart_bout_times_OR_PL{id})];
    last_quart_bout_times_OR_PL_mat = [last_quart_bout_times_OR_PL_mat, mean(last_quart_bout_times_OR_PL{id})];
    num_trials_per_third = [num_trials_per_third,  size(last_quart_bout_times_OR_PL{id},2)];
end
[aa, bb] =get_mean_and_sem(first_quart_bout_times_OR_PL_mat)
[aa, bb] =get_mean_and_sem(last_quart_bout_times_OR_PL_mat)
[aa, bb] =get_mean_and_sem(num_trials_per_third)
%-----------------------------------

%SINGLE LICKS plot earliest 1/4 of bout onsets vs slowest 1/4 of SINGLE LICK onsets
figure;
subplot(2,1,1);
shadedErrorBar( x_axis, nanmean(f_quart_FL_NR_PL_Ca_hist), f_quart_FL_NR_PL_Ca_hist_sem, 'b'); hold on;
shadedErrorBar( x_axis, nanmean(l_quart_FL_NR_PL_Ca_hist), l_quart_FL_NR_PL_Ca_hist_sem, 'k');
xlabel('time (ms) from cue onset');
ylabel('Ca event firing rate');
title('Normal Rewarded trials with licks; RS cells; All post-learning; blue= 1st 1/4,  black= last 1/4');
xlim([-1000 1500]); ylim([0 2]);

subplot(2,1,2);
plot( x_axis, nanmean(f_quart_FL_NR_PL_Ca_hist),  'b'); hold on;
plot( x_axis, nanmean(l_quart_FL_NR_PL_Ca_hist),  'k');
xlim([-1000 1500]); ylim([0 2]);
xlabel('time (ms) from cue onset');
ylabel('Ca event firing rate');
title('Normal Rewarded trials with licks; RS cells; All post-learning; blue= 1st 1/4,  black= last 1/4');
bar(x_axis, mean(NR_FL_f_quart_hist), 'b');  
bar(x_axis, mean(NR_FL_l_quart_hist), 'k');
errorbar(x_axis, mean(NR_FL_f_quart_hist), std(NR_FL_f_quart_hist)/sqrt(size(NR_FL_f_quart_hist,1)), 'b', 'LineStyle', 'none');
errorbar(x_axis, mean(NR_FL_l_quart_hist), std(NR_FL_l_quart_hist)/sqrt(size(NR_FL_l_quart_hist,1)), 'k', 'LineStyle', 'none');

%mean values
first_quart_FL_times_NR_PL_mat = [];
last_quart_FL_times_NR_PL_mat = [];
num_trials_per_quart = [];
for id=sess_subset
    first_quart_FL_times_NR_PL_mat = [first_quart_FL_times_NR_PL_mat, nanmean(NR_FL_times_f_quart{id})];
    last_quart_FL_times_NR_PL_mat = [last_quart_FL_times_NR_PL_mat, nanmean(NR_FL_times_l_quart{id})];
    num_trials_per_quart = [num_trials_per_quart,  size(NR_FL_times_l_quart{id},2)];
end
[aa, bb] =get_mean_and_sem(first_quart_FL_times_NR_PL_mat)
[aa, bb] =get_mean_and_sem(last_quart_FL_times_NR_PL_mat)
[aa, bb] =get_mean_and_sem(num_trials_per_quart)

%no licking or late licking
NR_FL_PL_Ca_hist_mat = [];
NR_FL_PL_late_or_no_lick_hist2 = [];
for session_num = 1:length(NR_FL_PL_late_no_lick_ca_hist)
    if size(NR_FL_PL_late_no_lick_ca_hist{session_num},1) > 6
        NR_FL_PL_Ca_hist_mat = cat(1, NR_FL_PL_Ca_hist_mat, mean(NR_FL_PL_late_no_lick_ca_hist{session_num},1));
        NR_FL_PL_late_or_no_lick_hist2 = [NR_FL_PL_late_or_no_lick_hist; NR_FL_PL_late_or_no_lick_hist(session_num,:)];
    end
end
[NR_FL_PL_Ca_hist_mat_mean, NR_FL_PL_Ca_hist_mat_sem] = get_mean_and_sem(NR_FL_PL_Ca_hist_mat);
[late_lick_mean, late_lick_sem] = get_mean_and_sem(NR_FL_PL_late_or_no_lick_hist2);
figure;
shadedErrorBar(x_axis, NR_FL_PL_Ca_hist_mat_mean, NR_FL_PL_Ca_hist_mat_sem); hold on;
bar(x_axis, late_lick_mean);
errorbar(x_axis, late_lick_mean, late_lick_sem, 'LineStyle', 'none');
title('trials with no licks or no licks between 200 and 1000ms after cue onset');
xlabel('time (ms) from first lick');
ylabel('Ca event firing rate');
xlim([-1000 2000]);

% 
% figure; 
% for id = 1:size(bout_aligned_ca_hist_D1,1)
%     figure; 
%     plot(x_axis, bout_aligned_ca_hist_D1(id,:), 'k'); hold on;
%     plot(x_axis, bout_aligned_ca_hist_PL(id,:), 'b'); 
%     
% end


