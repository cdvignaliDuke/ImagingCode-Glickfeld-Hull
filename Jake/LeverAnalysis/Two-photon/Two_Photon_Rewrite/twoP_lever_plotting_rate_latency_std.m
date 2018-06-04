% New plotting script to caclulate and plot spike rate, latency, and std
%ACROSS EVENTS
clear

file_info;
dateID = date;
window_width = 200;
%min_hold_dur = 1000;
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
 
    %log mean correct trial reaction time for each mouse   
    assert(length(cell2mat(b_data.reactTimesMs)) == length(b_data.reactTimesMs));
    RTs_ms_this_mouse = cell2mat(b_data.reactTimesMs(corr_trials_bdata_inx));
    cLeverUp = cell2mat_padded(b_data.cLeverUp);
    cTargetOn = celleqel2mat_padded(b_data.cTargetOn);
    cTimes = frame_info.times;
    counterRTs = cTimes(cLeverUp(corr_trials_bdata_inx))-cTimes(cTargetOn(corr_trials_bdata_inx));
    RTs_ms_mean(session_num) = mean((counterRTs));
    RTs_ms_sem(session_num)= std(double((counterRTs)))/sqrt(length(counterRTs));
%      RTs_ms_mean(session_num) = mean(RTs_ms_this_mouse);
%      RTs_ms_sem(session_num)= std([double(RTs_ms_this_mouse)])/sqrt(length(RTs_ms_this_mouse));

    for cell_num = 1:length(RS_cells)
        %correct trials
        this_cell_cue_hists = NaN(size(RT_in_frames,1), size(resp_cue_hist,2)); %trial#  by frames
        for trial_num = 1:length(RT_in_frames)
            %adjust the psths for cue aligned
            this_trial_RT = RT_in_frames(trial_num);
            this_trial_hist_length = length( success_RS(cell_num).hist(:,trial_num)' );
            this_cell_cue_hists(trial_num, [this_trial_RT+1:(this_trial_RT+this_trial_hist_length)]) = success_RS(cell_num).hist(:,trial_num)';  %adjust the indeces of this trial for this cell
        end
        
        %average across trials for this cell and store the mean
        resp_cue_hist(cell_num, :) = nanmean(this_cell_cue_hists,1);  
        
        %find the index of all spikes in a X ms window around release
        [ind_frame, ind_trial] = find(success_RS(cell_num).hist([std_window],:));  %row indeces = frame num     column indeces = trial_num     
        %convert ind_frame to ms. Log the mean latency and std
        lever_ind_frame = (ind_frame - lever_frame_this_window)*double(ifi);
        lever_spike_std_cell(cell_num) = std(lever_ind_frame);
        lever_spike_lat_cell(cell_num) = mean(lever_ind_frame); 
        
        %adjust index of spike by RT and store std and latency
        %cue_adj_vec = repmat(lever_frame_this_window, [1,length(ind_frame)]);
        cue_adj_vec = double(RT_in_frames([ind_trial]));
        cue_ind_frame = ((ind_frame' - lever_frame_this_window) + cue_adj_vec)*double(ifi);
        cue_spike_std_cell(cell_num) = std(cue_ind_frame);
        cue_spike_lat_cell(cell_num) = mean(cue_ind_frame); 
        
        %early trials: find the index of spikes in window
        [ind_frame, ind_trial] = find(fail_RS(cell_num).hist([std_window],:));  %row indeces = frame num     column indeces = trial_num     
        %convert ind_frame to ms. Log the mean latency and std
        lever_ind_frame_fail = (ind_frame - lever_frame_this_window)*double(ifi);
        lever_spike_std_cell_fail(cell_num) = std(lever_ind_frame_fail);
        lever_spike_lat_cell_fail(cell_num) = mean(lever_ind_frame_fail); 
    end
    
    %store the average values of spike std for each session
    lever_spike_std_cell_mean(session_num) =  nanmean(lever_spike_std_cell);
    lever_spike_std_cell_sem(session_num) =  nanstd(lever_spike_std_cell)/sqrt(length(lever_spike_std_cell));
    cue_spike_std_cell_mean(session_num) =  nanmean(cue_spike_std_cell);
    cue_spike_std_cell_sem(session_num) =  nanstd(cue_spike_std_cell)/sqrt(length(cue_spike_std_cell));
    lever_spike_std_cell_fail_mean(session_num) = nanmean(lever_spike_std_cell_fail);
    lever_spike_std_cell_fail_sem(session_num) = nanstd(lever_spike_std_cell_fail)/sqrt(length(lever_spike_std_cell_fail));
    
    %store the average values of spike latency for each session
    lever_spike_lat_cell_mean(session_num) =  nanmean(lever_spike_lat_cell);
    lever_spike_lat_cell_sem(session_num) =  nanstd(lever_spike_lat_cell)/sqrt(length(lever_spike_lat_cell));
    cue_spike_lat_cell_mean(session_num) =  nanmean(cue_spike_lat_cell);
    cue_spike_lat_cell_sem(session_num) =  nanstd(cue_spike_lat_cell)/sqrt(length(cue_spike_lat_cell));
    lever_spike_lat_cell_fail_mean(session_num) = nanmean(lever_spike_lat_cell_fail);
    lever_spike_lat_cell_fail_sem(session_num) = nanstd(lever_spike_lat_cell_fail)/sqrt(length(lever_spike_lat_cell_fail));
        
    %store 
    TC_ifi(session_num) = ifi;
    
    clear lever_spike_std_cell cue_spike_std_cell  lever_spike_lat_cell cue_spike_lat_cell lever_spike_std_cell_fail lever_spike_lat_cell_fail
    
%     %cell arrays with the average response for each cell in Hz
%     success_hist_RS{session_num} = resp_success_hist.*(1000/double(ifi));
%     fail_hist_RS{session_num} = resp_fail_hist.*(1000/double(ifi));
%     cue_hist_RS{session_num} = resp_cue_hist.*(1000/double(ifi));

end



%% spike std scatterplot
figure; subplot(1,3,1);
errorbarxy(lever_spike_std_cell_mean, cue_spike_std_cell_mean, lever_spike_std_cell_sem, cue_spike_std_cell_sem, {'o', 'k', 'k'});
xx= [1:300]; hold on;
plot(xx,xx);
ylim([0 300]); xlim([0 300]);
axis square
xlabel('lever aligned trials');
ylabel('cue aligned trials');
title(['std of spike latencies for spike -', num2str(window_width), ':', num2str(window_width), ' around release']);
[hh,pp, ~, stats_spikes] = ttest(lever_spike_std_cell_mean, cue_spike_std_cell_mean);

% spike latency scatterplot
% subplot(1,2,2); 
% errorbarxy(lever_spike_lat_cell_mean, cue_spike_lat_cell_mean, lever_spike_lat_cell_sem, cue_spike_lat_cell_sem, {'o', 'k', 'k'});
% xx= [-50:500]; hold on;
% plot(xx,xx);
% axis square
% %ylim([0 215]); xlim([0 215]);
% xlabel('lever aligned trials');
% ylabel('cue aligned trials');
% title('spike latencies for spike -100:100 around release');

%RT vs spike latency relative to the cue
subplot(1,3,2);
errorbarxy(RTs_ms_mean, cue_spike_lat_cell_mean, RTs_ms_sem, cue_spike_lat_cell_sem, {'o', 'k', 'k'});
xx= [0:500]; hold on;
plot(xx,xx);
axis square
%ylim([0 215]); xlim([0 215]);
xlabel('Reaction times (ms)');
ylabel('cue aligned trials (ms)');
title(['spike latencies for spikes -',num2str(window_width), ':',num2str(window_width), ' around release']);

%RT vs spike latency relative to the lever
subplot(1,3,3);
errorbarxy(RTs_ms_mean, lever_spike_lat_cell_mean, RTs_ms_sem, lever_spike_lat_cell_sem, {'o', 'k', 'k'});
%xx= [-100:300]; hold on;
%plot(xx,xx);
axis square
ylim([-100 400]); xlim([0 500]);
xlabel('Reaction times (ms)');
ylabel('lever aligned trials (ms)');
title(['spike latencies for spikes -',num2str(window_width), ':',num2str(window_width), ' around release']);

%%  std of latency for corrects vs earlies across cells vs across events   (figure 6 d
%ACROSS TRIALS
figure; ;
errorbarxy(lever_spike_std_cell_mean, lever_spike_std_cell_fail_mean, lever_spike_std_cell_sem, lever_spike_std_cell_fail_sem, {'o', 'k', 'k'});
xx= [-100:250]; hold on;
plot(xx,xx);
ylim([0 215]); xlim([0 215]);
axis square
xlabel('correct trials');
ylabel('early trials');
title(['Across trials: std of spike times -',num2str(window_width), ':',num2str(window_width), ' around release']);
[hh,pp, ~, stats_spikes] = ttest(lever_spike_std_cell_mean, lever_spike_std_cell_fail_mean);

%spike latency across trials correct vs early
figure; ;
errorbarxy(lever_spike_lat_cell_mean, lever_spike_lat_cell_fail_mean, lever_spike_lat_cell_sem, lever_spike_lat_cell_fail_sem, {'o', 'k', 'k'});
xx= [-100:250]; hold on;
plot(xx,xx);
%ylim([0 215]); xlim([0 215]);
axis square
xlabel('correct trials');
ylabel('early trials');
title(['Across trials: latency of spike times -',num2str(window_width), ':',num2str(window_width), ' around release']);
%[hh,pp, ~, stats_spikes] = ttest(lever_spike_std_cell_mean, lever_spike_std_cell_fail_mean);


%% PSTH
% %transpose cell arrays. Cue hist is already transposed.
% success_hist_RS = cellfun(@transpose, success_hist_RS, 'UniformOutput', 0);
% fail_hist_RS = cellfun(@transpose, fail_hist_RS, 'UniformOutput', 0);
% 
% frame_size = cell2mat(cellfun(@size, success_hist_RS, 'UniformOutput', 0));
% max_frame = max(frame_size(2:2:end));
% frame_size_cue = cell2mat(cellfun(@size, cue_hist_RS, 'UniformOutput', 0));
% max_frame_cue = max(frame_size_cue(2:2:end));
% 
% %interpolate the psths so they can be plotted on the same axes
% success_hist_RS_all   = interp_frame(success_hist_RS, max_frame);
% fail_hist_RS_all = interp_frame(fail_hist_RS, max_frame);
% cue_hist_RS_all   = interp_frame(cue_hist_RS, max_frame_cue);
% 
% %get means and sems
% succ_mean_psth = nanmean(success_hist_RS_all,1);
% succ_sem_psth = nanstd(success_hist_RS_all,[],1)./sqrt(size(success_hist_RS_all,1));
% fail_mean_psth = nanmean(fail_hist_RS_all,1);
% fail_sem_psth = nanstd(fail_hist_RS_all,[],1)./sqrt(size(fail_hist_RS_all,1));
% cue_mean_psth = nanmean(cue_hist_RS_all,1);
% cue_sem_psth = nanstd(cue_hist_RS_all,[],1)./sqrt(size(cue_hist_RS_all,1));
% 
% % PLOTTING
% %tth = [1-(size(fail_hist_RS_all,2)/2):(size(fail_hist_RS_all,2)/2)].*double(min((TC_ifi)));
% tth = [-ceil(1000/(double(min(TC_ifi)))):ceil(1000/(double(min(TC_ifi))))]*double(min(TC_ifi));
% ttc = [-ceil(1000/(double(min(TC_ifi)))):ceil(2000/(double(min(TC_ifi))))]*double(min(TC_ifi));
% fig=figure;
% %corrects vs earlies
% subplot(2,1,1)
% shadedErrorBar(tth, succ_mean_psth, succ_sem_psth, 'k');
% hold on;
% shadedErrorBar(tth, fail_mean_psth, fail_sem_psth, 'r');
% xlim([-1000 1000])
% ylim([0 3])
% xlabel('Time (ms)')
% ylabel('Firing rate (Hz)')
% title(['PSTH- responsive cells- n = ' num2str(size(success_hist_RS_all,1))]);
% %lever vs cue aligned
% subplot(2,1,2)
% tth = [1-(size(cue_hist_RS_all,2)/2):(size(cue_hist_RS_all,2)/2)].*double(min((TC_ifi)));
% %shadedErrorBar(tth, nanmean(success_hist_RS_all,1), std(success_hist_RS_all,[],1)./sqrt(size(success_hist_RS_all,1)), 'k');
% %hold on;
% shadedErrorBar(ttc, cue_mean_psth, cue_sem_psth, 'b');
% xlim([-1000 1000])
% ylim([0 3])
% xlabel('Time (ms)')
% ylabel('Firing rate (Hz)')
% title(['PSTH- responsive cells- n = ' num2str(size(success_hist_RS_all,1))])
% supertitle(['PSTH- responsive cells- Black: success; Red: Failure; Blue: cue aligned'])
% 
% %lever vs cue aligned corrects overlayed
% subplot(2,1,1);
% x_axis = [1:34]*30-525; hold on;
% shadedErrorBar(x_axis, cue_mean_psth([(40-16):(40+17)]), cue_sem_psth([(40-16):(40+17)]), 'b');
% vline(-195 ,'k')

%% 
% for session_num = 1:30
%     
% succ_latency_cell_mean(session_num) = nanmean(success_syn_cell_RS_spike_ind_all{session_num});
% succ_latency_cell_sem(session_num) = nanstd(success_syn_cell_RS_spike_ind_all{session_num})/sqrt(length(success_syn_cell_RS_spike_ind_all{session_num}));
% fail_latency_cell_mean(session_num) = nanmean(fail_syn_cell_RS_spike_ind_all{session_num});
% fail_latency_cell_sem(session_num) = nanstd(fail_syn_cell_RS_spike_ind_all{session_num})/sqrt(length(fail_syn_cell_RS_spike_ind_all{session_num}));
% cue_latency_cell_mean(session_num) = nanmean(cue_syn_cell_RS_spike_ind_all{session_num});
% cue_latency_cell_sem(session_num) = nanstd(cue_syn_cell_RS_spike_ind_all{session_num})/sqrt(length(cue_syn_cell_RS_spike_ind_all{session_num}));
% end
% figure;
% errorbarxy(succ_latency_cell_mean, fail_latency_cell_mean, succ_latency_cell_sem, fail_latency_cell_sem); hold on;
% xx = [-1:0.1:1];
% plot(xx,xx);
% xlabel('correct trials mean latency in frame nums');
% ylabel('fail trials mean latency in frame nums');
% title('mean latency to the nearest spike for each cell across trials');
% 
% figure;
% errorbarxy(succ_latency_cell_mean, cue_latency_cell_mean, succ_latency_cell_sem, cue_latency_cell_sem); hold on;
% xx = [-1:0.1:1];
% plot(xx,xx);
% xlabel('lever aligned trials mean latency in frame nums');
% ylabel('cue aligned trials mean latency in frame nums');
% title('mean latency to the nearest spike for each cell across trials');
% 
% % figure;   
% % for session_num = 1:30
% %     scatter(nanmean(success_syn_all_RS{session_num}),   nanmean(fail_syn_all_RS{session_num}));%,    nanstd(success_syn_all_RS{session_num})/sum(~isnan(success_syn_all_RS{session_num})),    nanstd(fail_syn_all_RS{session_num})/sum(~isnan(fail_syn_all_RS{session_num})));
% %     hold on;
% % end
% % 
% % 
% % %plot latencies
% % success_syn_cell_RS_spike_ind_all_mean = 






