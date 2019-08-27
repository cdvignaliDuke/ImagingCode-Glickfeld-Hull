%licking summary script 

clear
file_info_CRP_all;
out_base = fullfile('Z:\Analysis\Cue_reward_pairing_analysis\2P\PCA_ICA_recalibration\Summary_folder\');
crp_dir = 'Z:\Analysis\Cue_reward_pairing_analysis\2P\PCA_ICA_recalibration\';
%load([out_base 'cell_count.mat']);
normalR_nevents_D1 = 0; spont_nevents_D1 = 0;
sess_subset = [1:15];%[1:14];  % [1:14]%1:6%  [1,2,3,4,5,6,7,8,9,10,11,13,15]; %1:size(days_1,2) %
ifi = 33;
pre_frames=30;
post_frames=76;
x_axis = [-pre_frames:post_frames]*ifi;

%load and group variables for DAY 1
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
    dest_sub  = fullfile(crp_dir, [session_date, '_', runID{rID}, '_', mouse_ID],'\');
    dest_sub_spikes  = fullfile(dest_sub, 'spike_outputs', '\');
    
    %load variables
    load([dest_sub '_cue_movies_lick.mat']);

    %Group event hists for all NR/OR trials. all/RS
    normalR_hist_D1_licks{id} = lick_trace_NR;
    omitR_hist_D1_licks{id} = lick_trace_OR;
    
    %Group event hists for NO LICK NR/OR trials. all/RS
    normalR_hist_nolick_D1{id} = lick_trace_NR(logical(NR_lick_info.no_lick_cue_to_500),:);
    omitR_hist_nolick_D1{id} = lick_trace_OR(logical(OR_lick_info.no_lick_cue_to_500),:);
end

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
    dest_sub  = fullfile(crp_dir, [session_date, '_', runID{rID}, '_', mouse_ID],'\');
    dest_sub_spikes  = fullfile(dest_sub, 'spike_outputs', '\');
    
    %load variables
    load([dest_sub '_cue_movies_lick.mat']);
    
    %Group event hists for all NR/OR trials. all/RS
    normalR_hist_PL_licks{id} = lick_trace_NR;
    omitR_hist_PL_licks{id} = lick_trace_OR;
    
    %Group event hists for NO LICK NR/OR trials. all/RS
    normalR_hist_nolick_PL{id} = lick_trace_NR(logical(NR_lick_info.no_lick_cue_to_500),:);
    omitR_hist_nolick_PL{id} = lick_trace_OR(logical(OR_lick_info.no_lick_cue_to_500),:);
end

%calculate mean lick trace for each session
normalR_avg_D1_licks = NaN(length(normalR_hist_nolick_PL), size(normalR_hist_nolick_PL{1},2)); %dim1=animal  dim2=frame/time
omitR_avg_D1_licks = NaN(length(normalR_hist_nolick_PL), size(normalR_hist_nolick_PL{1},2));
normalR_avg_D1_nolick = NaN(length(normalR_hist_nolick_PL), size(normalR_hist_nolick_PL{1},2));
omitR_avg_D1_nolick = NaN(length(normalR_hist_nolick_PL), size(normalR_hist_nolick_PL{1},2));
normalR_avg_PL_licks = NaN(length(normalR_hist_nolick_PL), size(normalR_hist_nolick_PL{1},2));
omitR_avg_PL_licks = NaN(length(normalR_hist_nolick_PL), size(normalR_hist_nolick_PL{1},2));
normalR_avg_PL_nolick = NaN(length(normalR_hist_nolick_PL), size(normalR_hist_nolick_PL{1},2));
omitR_avg_PL_nolick = NaN(length(normalR_hist_nolick_PL), size(normalR_hist_nolick_PL{1},2));

for id = 1:length(normalR_hist_nolick_PL)
    normalR_avg_D1_licks(id,:) = mean(normalR_hist_D1_licks{id},1);
    omitR_avg_D1_licks(id,:) = mean(omitR_hist_D1_licks{id},1);
    normalR_avg_D1_nolick(id,:) = mean(normalR_hist_nolick_D1{id},1);
    omitR_avg_D1_nolick(id,:) = mean(omitR_hist_nolick_D1{id},1);
    
    normalR_avg_PL_licks(id,:) = mean(normalR_hist_PL_licks{id},1);
    omitR_avg_PL_licks(id,:) = mean(omitR_hist_PL_licks{id},1);
    normalR_avg_PL_nolick(id,:) = mean(normalR_hist_nolick_PL{id},1);
    omitR_avg_PL_nolick(id,:) = mean(omitR_hist_nolick_PL{id},1);
end

%plot the mean lick trace for each session
for id = 1:length(normalR_hist_D1_licks)
    figure;
    subplot(2,3,1); hold on;
    errorbar(x_axis, mean(normalR_hist_D1_licks{id},1), std(normalR_hist_D1_licks{id},1)/sqrt(size(normalR_hist_D1_licks{id},1)), 'k'); 
    errorbar(x_axis, mean(normalR_hist_PL_licks{id},1), std(normalR_hist_PL_licks{id},1)/sqrt(size(normalR_hist_PL_licks{id},1)), 'b');
    title('Normal reward: All trials');
    
    subplot(2,3,4); hold on;
    errorbar(x_axis, mean(omitR_hist_D1_licks{id},1), std(omitR_hist_D1_licks{id},1)/sqrt(size(omitR_hist_D1_licks{id},1)), 'k');
    errorbar(x_axis, mean(omitR_hist_PL_licks{id},1), std(omitR_hist_PL_licks{id},1)/sqrt(size(omitR_hist_PL_licks{id},1)), 'b');
    title('Omitted reward: All trials');

    suptitle([days_1{id}(end-2:end), ' lick traces. Black: day1;  Blue: post learning']);
end

%plot cumulative histogram for each session
for id = 1:length(normalR_hist_D1_licks)
      figure('rend', 'painters', 'pos', [100 50 1250 500]); 
    subplot(2,2,1); hold on;
    plot(x_axis, cumsum(mean(normalR_hist_D1_licks{id},1)),  'k'); 
    plot(x_axis, cumsum(mean(normalR_hist_PL_licks{id},1)),  'b');
    title('Normal reward: All trials: Black: D1;   Blue: post-learning');
    vline(600, 'b'); vline(0, 'k');
    xlim([min(x_axis) max(x_axis)]); ylim([0 12]);
    
    subplot(2,2,3); hold on;
    plot(x_axis, cumsum(mean(omitR_hist_D1_licks{id},1)),  'k'); 
    plot(x_axis, cumsum(mean(omitR_hist_PL_licks{id},1)),  'b');
    title('Omitted reward: All trials.  Black: D1;   Blue: post-learning');
    vline(600, '--'); vline(0, 'k');
    xlim([min(x_axis) max(x_axis)]); ylim([0 6]);
    
    subplot(2,2,2); hold on;
    num_trials = floor(size(normalR_hist_D1_licks{id},1)*(1/2));
    num_trials2 = floor(size(normalR_hist_D1_licks{id},1)*(2/3));
    num_trials3 = floor(size(normalR_hist_D1_licks{id},1)*(3/4));
    plot(x_axis, cumsum(mean(normalR_hist_D1_licks{id}([1:num_trials],:),1)),  'k'); 
    plot(x_axis, cumsum(mean(normalR_hist_D1_licks{id}([1:num_trials2],:),1)),  'b'); 
    plot(x_axis, cumsum(mean(normalR_hist_D1_licks{id}([1:num_trials3],:),1)),  'g'); 
    title('Normal reward D1: black: 1/2 trials;   blue: 2/3 trials;   green: 3/4 trials');
    vline(600, 'b'); vline(0, 'k');
    xlim([min(x_axis) max(x_axis)]); ylim([0 12]);
    
    subplot(2,2,4); hold on;
    num_trials = floor(size(omitR_hist_D1_licks{id},1)*(1/2));
    num_trials2 = floor(size(omitR_hist_D1_licks{id},1)*(2/3));
    num_trials3 = floor(size(omitR_hist_D1_licks{id},1)*(3/4));
    plot(x_axis, cumsum(mean(omitR_hist_D1_licks{id}([1:num_trials],:),1)),  'k'); 
    plot(x_axis, cumsum(mean(omitR_hist_D1_licks{id}([1:num_trials2],:),1)),  'b'); 
    plot(x_axis, cumsum(mean(omitR_hist_D1_licks{id}([1:num_trials3],:),1)),  'g'); 
    title('Omitted reward D1: black: 1/2 trials;   blue: 2/3 trials;   green: 3/4 trials');
    vline(600, 'b'); vline(0, 'k');
    xlim([min(x_axis) max(x_axis)]); ylim([0 6]);

    suptitle([days_1{id}(end-5:end)]);
end

%plot mean cumulative histograms of Day1 vs post learning 500ms   All trials
figure('rend', 'painters', 'pos', [100 50 1250 500]);
subplot(2,1,1); hold on;
plot(x_axis, cumsum(mean(normalR_avg_D1_licks,1)),  'k');
plot(x_axis, cumsum(mean(normalR_hist_PL_licks{id},1)),  'b');
title('Normal reward: All trials: Black: D1;   Blue: post-learning');
vline(600, 'b'); vline(0, 'k');
xlim([min(x_axis) max(x_axis)]); ylim([0 12]);

subplot(2,1,2); hold on;
plot(x_axis, cumsum(mean(omitR_avg_D1_licks,1)),  'k');
plot(x_axis, cumsum(mean(omitR_avg_PL_licks,1)),  'b');
title('Omitted reward: All trials.  Black: D1;   Blue: post-learning');
vline(600, '--'); vline(0, 'k');
xlim([min(x_axis) max(x_axis)]); ylim([0 6]);
suptitle('across animals mean cumulative histogram of licking');

%plot mean cumulative histograms of Day1 vs post learning 500ms  nolick trials
figure('rend', 'painters', 'pos', [100 50 1250 500]);
subplot(2,1,1); hold on;
plot(x_axis, cumsum(mean(normalR_avg_D1_nolick,1)),  'k');
plot(x_axis, cumsum(mean(normalR_avg_PL_nolick,1)),  'b');
title('Normal reward: All trials: Black: D1;   Blue: post-learning');
vline(600, 'b'); vline(0, 'k');
xlim([min(x_axis) max(x_axis)]); ylim([0 12]);

subplot(2,1,2); hold on;
plot(x_axis, cumsum(mean(omitR_avg_D1_nolick,1)),  'k');
plot(x_axis, cumsum(mean(omitR_avg_PL_nolick,1)),  'b');
title('Omitted reward: All trials.  Black: D1;   Blue: post-learning');
vline(600, '--'); vline(0, 'k');
xlim([min(x_axis) max(x_axis)]); ylim([0 6]);
suptitle('Across animals mean cumulative histogram of licking; no licking cue:cue+500ms');








