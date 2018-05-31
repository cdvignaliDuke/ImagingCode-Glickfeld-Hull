%ACROSS CELLS

clear

file_info;
dateID = date;
window_width = 200;
figure;
for session_num = 1:30
    session_num
    %check to make sure file exists
    for rID = 1:2
        tc_dir  = fullfile('\\crash.dhe.duke.edu\data\home\ziye\2P_Analysis\2P_Analysis',[dateID{session_num}, '_', runID{rID}, '_', mouseID{session_num}],'\');
        behav_dir = '\\crash.dhe.duke.edu\data\home\andrew\Behavior\Data\';
        dest = tc_dir;
        if exist(dest)
            break
        end
    end
    
    %load variables
    load([dest 'parse_behavior.mat'])
    load([dest '_cell_categories2.mat']);
    load([dest '_event_hist.mat']);
    load([dest '_event_summary.mat']);
    dataName   = dir([behav_dir, '*', behavID{session_num}, '-', date{session_num}(1:6), '*']);
    load([behav_dir, dataName(end).name]);
    b_data = input; clear input;
    
    %define responsive cells
    RS_cells = unique([release_resp_cells success_resp_cells fail_resp_cells press_resp_cells tooFast_resp_cells]);
    success_RS = success(RS_cells);
    fail_RS = fail(RS_cells);
    
    %make an empty matrix to fill in with the mean responses for each cell after adjusting individual trials for RT
    resp_cue_hist = NaN(length(RS_cells), size(success_RS(1).hist,1));  %cells by frames   
    resp_cue_hist(:, [(end+1):(end+(1000/ifi))]) = NaN; %add on 1s of time to the end of the plot
    lever_frame(session_num) =  ceil(size(success_RS(1).hist,1)/2);
    std_window = [lever_frame(session_num)-ceil(window_width/ifi):lever_frame(session_num)+ceil(window_width/ifi)];
    lever_frame_this_window = find(std_window==lever_frame(session_num));
    
    %define some indeces
    lever_frame(session_num) =  ceil(size(success_RS(1).hist,1)/2);
    std_window = [lever_frame(session_num)-ceil(window_width/ifi):lever_frame(session_num)+ceil(window_width/ifi)];
    lever_frame_this_window = find(std_window==lever_frame(session_num));
    
    %remove unimaged trials?
    num_of_trials = size(success_RS(1).hist,2);
    corr_trials_bdata_inx = find(trial_outcome.corr_inx);
    if  ismember(session_num, [25,27,29]) & length(trial_outcome.succ_rm_event_idx) ==1 & trial_outcome.succ_rm_event_idx == length(corr_trials_bdata_inx); %sub =[25,27,29];
        RT_in_frames = cell2mat( b_data.cLeverUp(corr_trials_bdata_inx)) - cell2mat(b_data.cTargetOn(corr_trials_bdata_inx) );
    else
        corr_trials_bdata_inx([trial_outcome.succ_rm_event_idx]) = [];
        RT_in_frames = cell2mat( b_data.cLeverUp(corr_trials_bdata_inx)) - cell2mat(b_data.cTargetOn(corr_trials_bdata_inx) );
    end
    assert(num_of_trials == length(RT_in_frames));
    edges = [200 275 350 600];
    RTs = celleqel2mat_padded(b_data.reactTimesMs(corr_trials_bdata_inx));
    [n_trials bins] = histc(RTs,edges);
    nbin = length(n_trials)-1;
    RT_avg = nan(1,nbin);
    RT_sem = nan(1,nbin);
    cue_spike_lat = nan(nbin,1);
    lever_spike_lat = nan(nbin,1);
    this_cell_avg_cue_hists = nan(nbin,length(RS_cells),size(resp_cue_hist,2));
    this_cell_avg_lever_hists = nan(nbin,length(RS_cells),size(success_RS(1).hist,1));
    all_cell_avg_cue_hists = nan(nbin,size(resp_cue_hist,2));
    all_cell_avg_lever_hists = nan(nbin,size(success_RS(1).hist,1));
    for ibin = 1:nbin
        if n_trials(ibin)>20
            ind = find(bins == ibin);
            RT_avg(ibin) = mean(RTs(ind),2);
            RT_sem(ibin) = std(RTs(ind),[],2)./sqrt(n_trials(ibin));
            all_cell_cue_hists = NaN(length(ind), size(resp_cue_hist,2), length(RS_cells)); %trial#  by frames by cells
            all_cell_lever_hists = NaN(length(ind),size(success_RS(1).hist,1), length(RS_cells));
            for cell_num = 1:length(RS_cells)
%                 %correct trials
%                 this_cell_cue_hists = NaN(length(ind), size(resp_cue_hist,2)); %trial#  by frames

                for trial_num = 1:length(ind)
                    it = ind(trial_num);
                    %adjust the psths for cue aligned
                    this_trial_RT = RT_in_frames(it);
                    this_trial_hist_length = length( success_RS(cell_num).hist(:,it)' );
                    %this_cell_cue_hists(it, [this_trial_RT+1:(this_trial_RT+this_trial_hist_length)]) = success_RS(cell_num).hist(:,it)';  %adjust the indeces of this trial for this cell
                    all_cell_cue_hists(it, [this_trial_RT+1:(this_trial_RT+this_trial_hist_length)], cell_num) = success_RS(cell_num).hist(:,it)';  %adjust the indeces of this trial for this cell
                    
                end
                all_cell_lever_hists(:,:,cell_num) = success_RS(cell_num).hist(:,ind)';
            end
        
%                 this_cell_avg_cue_hists(ibin,cell_num,:) = smooth(nanmean(this_cell_cue_hists,1));
%                 this_cell_avg_lever_hists(ibin,cell_num,:) = smooth(mean(success_RS(cell_num).hist(:,ind),2));
                all_cell_avg_cue_hists(ibin,:) = nanmean(nanmean(all_cell_cue_hists,1),3);
                all_cell_avg_lever_hists(ibin,:) = mean(mean(all_cell_lever_hists,1),3);
                [all_cell_cue_peak_val, all_cell_cue_peak_ind] = max(all_cell_avg_cue_hists(ibin,ceil(1000./double(ifi))+1:(ceil(1000./double(ifi)) + ceil(RT_avg(ibin)./double(ifi)) + ceil(300./double(ifi)))),[],2);
                [all_cell_lever_peak_val, all_cell_lever_peak_ind] = max(all_cell_avg_lever_hists(ibin,ceil(1000./double(ifi))-ceil(RT_avg(ibin)./double(ifi))+1:ceil(1000./double(ifi)) + ceil(300./double(ifi))),[],2);
                
                cue_spike_lat(ibin,:) = squeeze(all_cell_cue_peak_ind .* double(ifi))';
                lever_spike_lat(ibin,:) = squeeze((all_cell_lever_peak_ind -ceil(RT_avg(ibin)./double(ifi))).* double(ifi))'; 
%                 subplot(5,6,session_num)
%                 plot(all_cell_avg_lever_hists')
                subplot(1,2,1)
                scatter(RT_avg(ibin),cue_spike_lat(ibin,:),'o')
                hold on
                subplot(1,2,2)
                scatter(RT_avg(ibin),lever_spike_lat(ibin,:),'o')
                hold on
        end
    end
end
% subplot(1,2,1)
% title('Cue')
% ylim([-200 700])
% subplot(1,2,2)
% title('Lever')
% ylim([-200 700])

%                 [this_cell_cue_peak_val, this_cell_cue_peak_ind] = max(this_cell_avg_cue_hists(ibin,cell_num,ceil(1000./double(ifi)):(ceil(1000./double(ifi)) + ceil(RT_avg(ibin)./double(ifi)) + ceil(300./double(ifi)))),[],3);
%                 if length(find(this_cell_avg_cue_hists(ibin,cell_num,ceil(1000./double(ifi)):ceil(1000./double(ifi)) + ceil(RT_avg(ibin)./double(ifi)) + ceil(300./double(ifi))) == this_cell_cue_peak_val)) < 3
%                     if ~diff(find(this_cell_avg_cue_hists(ibin,cell_num,ceil(1000./double(ifi)):ceil(1000./double(ifi)) + ceil(RT_avg(ibin)./double(ifi)) + ceil(300./double(ifi))) == this_cell_cue_peak_val)) == 1
%                         cue_spike_lat_cell(ibin,cell_num) = nan;
%                     else
%                         cue_spike_lat_cell(ibin,cell_num) = squeeze(this_cell_cue_peak_ind .* double(ifi))'; 
%                     end
%                 else
%                     cue_spike_lat_cell(ibin,cell_num) = squeeze(this_cell_cue_peak_ind .* double(ifi))';
%                 end
%                 [this_cell_lever_peak_val, this_cell_lever_peak_ind] = max(this_cell_avg_lever_hists(ibin,cell_num,ceil(1000./double(ifi))-ceil(RT_avg(ibin)./double(ifi)):ceil(1000./double(ifi)) + ceil(300./double(ifi))),[],3);
%                 if length(find(this_cell_avg_lever_hists(ibin,cell_num,ceil(1000./double(ifi))-ceil(RT_avg(ibin)./double(ifi)):ceil(1000./double(ifi))+ ceil(300./double(ifi))) == this_cell_lever_peak_val)) == 2
%                     if find(this_cell_avg_lever_hists(ibin,cell_num,ceil(1000./double(ifi)) -ceil(RT_avg(ibin)./double(ifi)):ceil(1000./double(ifi)) + ceil(300./double(ifi))) == this_cell_lever_peak_val)) == 1
%                         lever_spike_lat_cell(ibin,cell_num) = nan;
%                     else
%                         lever_spike_lat_cell(ibin,cell_num) = squeeze((this_cell_lever_peak_ind .* double(ifi))- RT_avg(ibin))'; 
%                     end
%                 else
%                     lever_spike_lat_cell(ibin,cell_num) = squeeze((this_cell_lever_peak_ind .* double(ifi))- RT_avg(ibin))'; 
% %                 end
%             end
%         end
%     end
% end

                
%                 %find the index of all spikes in a 100ms window around release
%                 [ind_frame, ind_trial] = find(success_RS(cell_num).hist([std_window],:));  %row indeces = frame num     column indeces = trial_num     
%                 %convert ind_frame to ms. Log the mean latency and std
%                 lever_ind_frame = (ind_frame - lever_frame_this_window)*double(ifi);
%                 lever_spike_std_cell(cell_num) = std(lever_ind_frame);
%                 lever_spike_lat_cell(cell_num) = mean(lever_ind_frame); 
% 
%                 %early trials: find the index of spikes in window
%                 [ind_frame, ind_trial] = find(fail_RS(cell_num).hist([std_window],:));  %row indeces = frame num     column indeces = trial_num     
%                 %convert ind_frame to ms. Log the mean latency and std
%                 lever_ind_frame_fail = (ind_frame - lever_frame_this_window)*double(ifi);
%                 lever_spike_std_cell_fail(cell_num) = std(lever_ind_frame_fail);
%                 lever_spike_lat_cell_fail(cell_num) = mean(lever_ind_frame_fail); 
% 
%                 %ACROSS CELLS
%                 for trial_num = 1:size(success_RS(1).hist,2)
%                     success_RS(:).hist(std_window,trial_num)
%                 end
%                 for trial_num = 1:size(fail_RS(1).hist,2)
% 
%                 end
%     end
%     
%     %store the average values of spike std for each session
%     lever_spike_std_cell_mean(session_num) =  nanmean(lever_spike_std_cell);
%     lever_spike_std_cell_sem(session_num) =  nanstd(lever_spike_std_cell)/sqrt(length(lever_spike_std_cell));
%     lever_spike_std_cell_fail_mean(session_num) = nanmean(lever_spike_std_cell_fail);
%     lever_spike_std_cell_fail_sem(session_num) = nanstd(lever_spike_lat_cell_fail)/sqrt(length(lever_spike_lat_cell_fail));
%     
%     %store the average values of spike latency for each session
%     lever_spike_lat_cell_mean(session_num) =  nanmean(lever_spike_lat_cell);
%     lever_spike_lat_cell_sem(session_num) =  nanstd(lever_spike_lat_cell)/sqrt(length(lever_spike_lat_cell));
%         
%     %store 
%     TC_ifi(session_num) = ifi;
%     
%     clear lever_spike_std_cell cue_spike_std_cell  lever_spike_lat_cell cue_spike_lat_cell lever_spike_std_cell_fail lever_spike_lat_cell_fail
% end


