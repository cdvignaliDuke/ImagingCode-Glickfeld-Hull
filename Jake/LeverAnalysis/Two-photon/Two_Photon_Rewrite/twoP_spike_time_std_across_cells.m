%ACROSS CELLS spike latencies and st.dev
clear

file_info;
dateID = date;
window_width = 200;
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
    
    %define some indeces
    lever_frame(session_num) =  ceil(size(success_RS(1).hist,1)/2);
    std_window = [lever_frame(session_num)-ceil(window_width/ifi):lever_frame(session_num)+ceil(window_width/ifi)];
    lever_frame_this_window = find(std_window==lever_frame(session_num));
    
    %correct trials
    for trial_num = 1:size(success_RS(1).hist,2)
        this_trial_hist = NaN( length(std_window), length(success_RS)); %  frame#(window) by   cell# 
        for cell_num = 1:length(success_RS)
            this_trial_hist(:, cell_num) = success_RS(cell_num).hist([std_window],trial_num);
        end
        [frame_ind, cell_ind] = find(this_trial_hist); 
        lever_ind_frame = (frame_ind - lever_frame_this_window)*double(ifi); %center the frame indeces around lever frame and then convert to ms 
        lever_spike_std(trial_num) = std(lever_ind_frame);
        lever_spike_lat(trial_num) = mean(lever_ind_frame); 
    end
    
    %early trials
    for trial_num = 1:size(fail_RS(1).hist,2)
        this_trial_hist = NaN( length(std_window), length(fail_RS)); %  frame#(window) by   cell# 
        for cell_num = 1:length(fail_RS)
            this_trial_hist(:, cell_num) = fail_RS(cell_num).hist([std_window],trial_num);
        end
        [frame_ind, cell_ind] = find(this_trial_hist); 
        lever_ind_frame = (frame_ind - lever_frame_this_window)*double(ifi); %center the frame indeces around lever frame and then convert to ms 
        lever_spike_std_fail(trial_num) = std(lever_ind_frame);
        lever_spike_lat_fail(trial_num) = mean(lever_ind_frame); 
    end
    
    %store the average values of spike std for each session
    lever_spike_std_mean(session_num) =  nanmean(lever_spike_std);
    lever_spike_std_sem(session_num) =  nanstd(lever_spike_std)/sqrt(length(lever_spike_std));
    lever_spike_std_fail_mean(session_num) = nanmean(lever_spike_std_fail);
    lever_spike_std_fail_sem(session_num) = nanstd(lever_spike_std_fail)/sqrt(length(lever_spike_std_fail));
    [corr_fail_ttest_h(session_num), corr_fail_ttest_p(session_num)] = ttest2(lever_spike_std, lever_spike_std_fail);
    
    %store the average values of spike latency for each session
    lever_spike_lat_mean(session_num) =  nanmean(lever_spike_lat);
    lever_spike_lat_sem(session_num) =  nanstd(lever_spike_lat)/sqrt(length(lever_spike_lat));
    lever_spike_lat_fail_mean(session_num) =  nanmean(lever_spike_lat_fail);
    lever_spike_lat_fail_sem(session_num) =  nanstd(lever_spike_lat_fail)/sqrt(length(lever_spike_lat_fail));
        
    %store 
    TC_ifi(session_num) = ifi;
    
    clear lever_spike_std      lever_spike_lat     lever_spike_std_fail      lever_spike_lat_fail
    
%     %cell arrays with the average response for each cell in Hz
%     success_hist_RS{session_num} = resp_success_hist.*(1000/double(ifi));
%     fail_hist_RS{session_num} = resp_fail_hist.*(1000/double(ifi));
%     cue_hist_RS{session_num} = resp_cue_hist.*(1000/double(ifi));

end


%%  std of latency for corrects vs earlies across cells    (figure 6 e
%ACROSS CELLS
figure; 
errorbarxy(lever_spike_std_mean, lever_spike_std_fail_mean, lever_spike_std_sem, lever_spike_std_fail_sem, {'o', 'k', 'k'});
sig_x = [lever_spike_std_mean(find(corr_fail_ttest_h))];
sig_y = [lever_spike_std_fail_mean(find(corr_fail_ttest_h))]
hold on; scatter(sig_x, sig_y, 'MarkerFaceColor', 'r' );
xx= [-100:250]; hold on;
plot(xx,xx);
ylim([0 160]); xlim([0 160]);
axis square
xlabel('correct trials');
ylabel('early trials');
title(['Across cells: std of spike times -',num2str(window_width), ':',num2str(window_width), ' around release']);
[hh,pp, ~, stats_spikes] = ttest(lever_spike_std_mean, lever_spike_std_fail_mean);


%spike latency across cells correct vs early
% figure; ;
% errorbarxy(lever_spike_lat_mean, lever_spike_lat_fail_mean, lever_spike_lat_sem, lever_spike_lat_fail_sem, {'o', 'k', 'k'});
% xx= [-100:250]; hold on;
% plot(xx,xx);
% %ylim([0 215]); xlim([0 215]);
% axis square
% xlabel('correct trials');
% ylabel('early trials');
% title(['Across cells: latency of spike times -',num2str(window_width), ':',num2str(window_width), ' around release']);

