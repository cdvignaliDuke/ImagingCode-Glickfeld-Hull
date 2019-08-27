%% plots the peak magnitude for the mean df/f trace. Find the mean and sem across animals 
%in this script day N actually refers to Day N+1 the post learning imaging day
% DayN+1 actually refers to day N+2 the unexpected reward condition

clear; 
WF_CRP_list_of_days;
mean_TC_dir = 'Y:\home\jake\Analysis\Cue_reward_pairing_analysis\WF\onset latencies and peak mags\mean interpolated TCs\';
color_pallette = {'b', 'r', 'k', 'm', 'g', 'c', 'b'};
peak_window = [800:1501];
baseline_window = [400:-1:200];
cue_time_interp = 501;
peak_diff_window= 500;
no_peak_diff_window = 1000;
yy = [1:2000];
pre_rew_win = [cue_time_interp+600-351 : cue_time_interp+600-1];
post_rew_win = [cue_time_interp + 600 +100 : cue_time_interp+600+100+350];

day1_rew_peak = [];
day1_om_peak = [];
dayN_rew_peak = [];
dayN_om_peak = [];
dayUR_exp_peak = [];
dayUR_unexp_peak = [];

day1_rew_baseline = [];
day1_om_baseline = [];
dayN_rew_baseline = [];
dayN_om_baseline = [];
dayUR_exp_baseline = [];
dayUR_unexp_baseline = [];

day1_pre_rew_mags = [];
dayN_pre_rew_mags = [];
day1_post_rew_mags = [];
dayN_post_rew_mags = [];
day1_pre_om_mags = [];
dayN_pre_om_mags = [];
day1_post_om_mags = [];
dayN_post_om_mags = [];

dayU_pre_rew_mags = []; 
dayU_post_rew_mags = []; 
dayU_pre_ur_mags = []; 
dayU_post_ur_mags = []; 

%% Day 1 rewarded and omission trials
figure; 
for session_num = 1:length(days_1);
    if exist([mean_TC_dir, 'day_1_vars_', days_1{session_num}, '_all.mat']);
        
        %load TC data and check number of trials per condition
        load([mean_TC_dir, 'day_1_vars_', days_1{session_num}, '_all']);
        if size(rew_om_roi_interp,2)<7 || size(rew_roi_interp,2) < 7 
            continue
        end
        
        %log the mean and sem values for pre-reward and post-reward windows
        day1_pre_rew_mags{session_num} = mean(rew_roi_interp([pre_rew_win],:)); 
        day1_post_rew_mags{session_num} = mean(rew_roi_interp([post_rew_win+150],:)); 
%         subplot(2,length(days_1),session_num); 
%         plot(rew_roi_interp, 'k'); hold on; vline(501, 'k'); vline(1101, 'k');
%         subplot(2,length(days_1),session_num+length(days_1)); 
%         plot(rew_om_roi_interp, 'r'); hold on; vline(501, 'k'); vline(1101, 'k');
        

        %find the peak times in the mean TC
        rew_roi_interp = mean(rew_roi_interp, 2)';
        this_peak = findpeaks(rew_roi_interp(peak_window));
        if size(this_peak) == 1 & this_peak < max(rew_roi_interp(1:cue_time_interp)) %if the peak magnitude was less than the max value before the cue then do not count the peak
            this_peak = [];
        end
        if isempty(this_peak) %if there was no peak in peak_window then simply take the max in the peak window
            no_peak = 1;
            this_peak = max(rew_roi_interp(peak_window));
        else
            no_peak = 0;
        end
        this_peak = max(this_peak);  %if there are two peaks then take the later one
        day1_rew_peak = [day1_rew_peak, this_peak];
        
        %find a baseline window 
        day1_rew_baseline = [day1_rew_baseline, mean(rew_roi_interp(1:cue_time_interp-100))];
        
        %plot each animal's mean trace
        subplot(2,2,3); hold on;
        plot(rew_roi_interp, color_pallette{session_num}); 
        
        %Day 1 omission trials -------------------------------------
        
        %log the mean and sem values for pre-omission and post-omission windows
        day1_pre_om_mags{session_num} = mean(rew_om_roi_interp([pre_rew_win],:)); 
        day1_post_om_mags{session_num} = mean(rew_om_roi_interp([post_rew_win+150],:)); 
        
        %find the peak df/f value
        rew_om_roi_interp = mean(rew_om_roi_interp, 2)';
        this_peak = findpeaks(rew_om_roi_interp(peak_window));
        if size(this_peak) == 1 & this_peak < max(rew_om_roi_interp(1:cue_time_interp))
            this_peak = [];
        end
        if isempty(this_peak)
            this_peak = [this_peak, max(rew_om_roi_interp(peak_window));];
            no_peak = 1;
        else
            no_peak = 0;
        end
        this_peak = max(this_peak);
        day1_om_peak = [day1_om_peak, this_peak];
        
        %find a baseline window 
        day1_om_baseline = [day1_om_baseline, mean(rew_om_roi_interp(1:cue_time_interp-100))];
        
        %plot each animal's mean trace
        subplot(2,2,4); hold on;
        plot(rew_om_roi_interp, color_pallette{session_num}); 
    end
end
%calculate means and sems
% [day1_rew_peak_mean, day1_rew_peak_sem] = get_mean_and_sem(day1_rew_peak);
% [day1_om_peak_mean, day1_om_peak_sem] = get_mean_and_sem(day1_om_peak);
% [day1_rew_pb_diff_mean, day1_rew_pb_diff_sem] = get_mean_and_sem(day1_rew_peak-day1_rew_baseline);
% [day1_om_pb_diff_mean, day1_om_pb_diff_sem] = get_mean_and_sem(day1_om_peak-day1_om_baseline);
% 
% subplot(2,2,1);
% scatter(day1_rew_peak, day1_om_peak); hold on;
% plot(day1_rew_peak_mean, day1_om_peak_mean, 'or'); 
% plot([day1_rew_peak_mean-day1_rew_peak_sem, day1_rew_peak_mean+day1_rew_peak_sem], [day1_om_peak_mean,day1_om_peak_mean], 'r');
% plot([day1_rew_peak_mean, day1_rew_peak_mean], [day1_om_peak_mean-day1_om_peak_sem,day1_om_peak_mean+day1_om_peak_sem], 'r');
% xlabel('rewarded trials'); ylabel('omission trials');
% title('day 1 peak vals');
% plot([-0.1, 0.3], [-0.1, 0.3], 'k');
% xlim([-0.1 0.3]); vline(0, 'k'); hline(0,'k');
% 
% subplot(2,2,2);
% scatter(day1_rew_peak-day1_rew_baseline, day1_om_peak-day1_om_baseline); hold on;
% plot(day1_rew_pb_diff_mean, day1_om_pb_diff_mean, 'or');
% plot([day1_rew_pb_diff_mean-day1_rew_pb_diff_sem, day1_rew_pb_diff_mean+day1_rew_pb_diff_sem], [day1_om_pb_diff_mean,day1_om_pb_diff_mean], 'r');
% plot([day1_rew_pb_diff_mean, day1_rew_pb_diff_mean], [day1_om_pb_diff_mean-day1_om_pb_diff_sem,day1_om_pb_diff_mean+day1_om_pb_diff_sem], 'r');
% xlabel('rewarded trials'); ylabel('omission trials');
% title('day 1 peak to baseline diff');
% plot([-0.1, 0.3], [-0.1, 0.3], 'k');
% xlim([-0.1 0.3]); vline(0, 'k'); hline(0,'k');
% 
% subplot(2,2,3); 
% title('Day 1 rewarded mean trace for each animal'); ylabel('df/f');
% ylim([-0.1 0.3]);
% vline(501, 'k'); vline(1101, 'k');
% subplot(2,2,4);
% title('Day 1 omission mean trace for each animal'); ylabel('df/f');
% ylim([-0.1 0.3]);
% vline(501, 'k'); vline(1101, 'k');

% %make sure all sessions have same number of trials
% day1_post_rew_mags = min_trial_finder(day1_post_rew_mags);
% day1_pre_rew_mags = min_trial_finder(day1_pre_rew_mags);
% day1_post_om_mags = min_trial_finder(day1_post_om_mags);
% day1_pre_om_mags = min_trial_finder(day1_pre_om_mags);

%calculate mean and sem in the pre and post reward windows
day1_post_rew_mags_mean = cellfun(@mean, day1_post_rew_mags);
day1_pre_rew_mags_mean = cellfun(@mean, day1_pre_rew_mags);
day1_post_om_mags_mean = cellfun(@mean, day1_post_om_mags);
day1_pre_om_mags_mean = cellfun(@mean, day1_pre_om_mags);

%remove nans 
day1_post_rew_mags_mean = day1_post_rew_mags_mean(~isnan(day1_post_rew_mags_mean));
day1_pre_rew_mags_mean = day1_pre_rew_mags_mean(~isnan(day1_pre_rew_mags_mean));
day1_post_om_mags_mean = day1_post_om_mags_mean(~isnan(day1_post_om_mags_mean));
day1_pre_om_mags_mean = day1_pre_om_mags_mean(~isnan(day1_pre_om_mags_mean));

%plot pre- post-reward magnitudes
figure; subplot(1,2,1);
scatter(day1_pre_rew_mags_mean, day1_post_rew_mags_mean); hold on;
plot(mean(day1_pre_rew_mags_mean), mean(day1_post_rew_mags_mean), 'or');
errorbarxy(mean(day1_pre_rew_mags_mean), mean(day1_post_rew_mags_mean), std(day1_pre_rew_mags_mean)/sqrt(length(day1_pre_rew_mags_mean)), std(day1_post_rew_mags_mean)/sqrt(length(day1_post_rew_mags_mean)));
ylim([-0.1 0.1]); xlim([-0.1 0.1]); hold on;
plot([-0.1, 0.1], [-0.1, 0.1], 'k');
hline(0, 'k'); vline(0,'k');
xlabel('df/f magnitude in pre-reward window');
ylabel('df/f magnitude in post-reward window');
title('day 1: rewarded trials');

subplot(1,2,2);
scatter(day1_pre_om_mags_mean, day1_post_om_mags_mean); hold on;
plot(mean(day1_pre_om_mags_mean), mean(day1_post_om_mags_mean), 'or');
errorbarxy(mean(day1_pre_om_mags_mean), mean(day1_post_om_mags_mean), std(day1_pre_om_mags_mean)/sqrt(length(day1_pre_om_mags_mean)), std(day1_post_om_mags_mean)/sqrt(length(day1_post_om_mags_mean)));
ylim([-0.1 0.1]); xlim([-0.1 0.1]); hold on;
plot([-0.1, 0.1], [-0.1, 0.1], 'k');
hline(0, 'k'); vline(0,'k');
xlabel('df/f magnitude in pre-reward window');
ylabel('df/f magnitude in post-reward window');
title('day 1: omission trials');

%calculate ttest
% [hh, pp] = ttest(day1_pre_rew_mags_mean, day1_post_rew_mags_mean);
% disp(['mean of day1 pre reward window magnitude = ', num2str(mean(day1_pre_rew_mags_mean)), ' +/- ', num2str(std(day1_pre_rew_mags_mean)/sqrt(length(day1_pre_rew_mags_mean)))]);
% disp(['mean of day1 post reward window magnitude = ', num2str(mean(day1_post_rew_mags_mean)), ' +/- ', num2str(std(day1_post_rew_mags_mean)/sqrt(length(day1_post_rew_mags_mean)))]);
% disp(['ttest for diff between pre and post reward windows on day 1 rewarded trials: h=', num2str(hh), ' p=', num2str(pp)]);
% 
% [hh, pp] = ttest(day1_pre_om_mags_mean, day1_post_om_mags_mean);
% disp(['mean of day1 pre reward window magnitude = ', num2str(mean(day1_pre_om_mags_mean)), ' +/- ', num2str(std(day1_pre_om_mags_mean)/sqrt(length(day1_pre_om_mags_mean)))]);
% disp(['mean of day1 post reward window magnitude = ', num2str(mean(day1_post_om_mags_mean)), ' +/- ', num2str(std(day1_post_om_mags_mean)/sqrt(length(day1_post_om_mags_mean)))]);
% disp(['ttest for diff between pre and post reward windows on day 1 omission trials: h=', num2str(hh), ' p=', num2str(pp)]);

%try calculating mean from all trials
[day1_post_rew_mags_ttmean, day1_post_rew_mags_ttsem] = get_mean_and_sem(cell2mat(day1_post_rew_mags));
[day1_pre_rew_mags_ttmean, day1_pre_rew_mags_ttsem] = get_mean_and_sem(cell2mat(day1_pre_rew_mags));
[day1_post_om_mags_ttmean, day1_post_om_mags_ttsem] = get_mean_and_sem(cell2mat(day1_post_om_mags));
[day1_pre_om_mags_ttmean, day1_pre_om_mags_ttsem] = get_mean_and_sem(cell2mat(day1_pre_om_mags));
[rew_h, rew_p] = ttest(cell2mat(day1_post_rew_mags), cell2mat(day1_pre_rew_mags));
[om_h, om_p] = ttest(cell2mat(day1_post_om_mags), cell2mat(day1_pre_om_mags));
disp(['mean of day1 pre reward window magnitude = ', num2str(day1_pre_rew_mags_ttmean), ' +/- ', num2str(day1_pre_rew_mags_ttsem)]);
disp(['mean of day1 post reward window magnitude = ', num2str(day1_post_rew_mags_ttmean), ' +/- ', num2str(day1_post_rew_mags_ttsem)]);
disp(['ttest for diff between pre and post reward windows on day 1 rewarded trials: h=', num2str(rew_h), ' p=', num2str(rew_p)]);
disp(['mean of day1 pre omission window magnitude = ', num2str(day1_pre_om_mags_ttmean), ' +/- ', num2str(day1_pre_om_mags_ttsem)]);
disp(['mean of day1 post omission window magnitude = ', num2str(day1_post_om_mags_ttmean), ' +/- ', num2str(day1_post_om_mags_ttsem)]);
disp(['ttest for diff between pre and post omission windows on day 1 omission trials: h=', num2str(om_h), ' p=', num2str(om_p)]);

%% Day N: Rewarded vs Omission
figure;
for session_num = 1:length(days_post);
    if  exist([mean_TC_dir, 'day_N_vars_', days_post{session_num}, '_all.mat'])
        
        load([mean_TC_dir, 'day_N_vars_', days_post{session_num}, '_all']);
        if size(rew_roi_interp,2) < 7 || size(rew_om_roi_interp,2) < 7  % must have at least two trials per condition
            continue
        end
        
        %log the mean and sem values for pre-reward and post-reward windows
        dayN_pre_rew_mags{session_num} = mean(rew_roi_interp([pre_rew_win],:)); 
        dayN_post_rew_mags{session_num} = mean(rew_roi_interp([post_rew_win],:)); 
        
        %get the mean and the derivative of the TC
        rew_roi_interp = mean(rew_roi_interp, 2)';
        
        %find the peaks of the mean TC
        this_peak = findpeaks(rew_roi_interp(peak_window));
        if size(this_peak) == 1 & this_peak < max(rew_roi_interp(1:cue_time_interp)) %if the peak mag is less than a value in the pre-cue period then there is no peak
            this_peak = [];
        end
        if isempty(this_peak)
            this_peak = [this_peak, max(rew_roi_interp(peak_window));];
            no_peak = 1;
        else
            no_peak = 0;
        end
        this_peak = max(this_peak);
        dayN_rew_peak = [dayN_rew_peak, this_peak];
        
        %find a baseline window
        dayN_rew_baseline = [dayN_rew_baseline, mean(rew_roi_interp(1:cue_time_interp-100))];
        
        %plot each animal's mean trace
        subplot(2,2,3); hold on;
        plot(rew_roi_interp, color_pallette{session_num});
        
        %Reward omission Day N+1 onset latencies ----------------------
        
        %log the mean and sem values for pre-omission and post-omission windows
        dayN_pre_om_mags{session_num} = mean(rew_om_roi_interp([pre_rew_win],:)); 
        dayN_post_om_mags{session_num} = mean(rew_om_roi_interp([post_rew_win],:));
        
        rew_om_roi_interp = mean(rew_om_roi_interp, 2)';
        
        %find the peaks of the mean TC
        this_peak = findpeaks(rew_om_roi_interp(peak_window));
        if size(this_peak) == 1 & this_peak < max(rew_om_roi_interp(1:cue_time_interp))
            this_peak = [];
        end
        if isempty(this_peak);
            this_peak = [this_peak, max(rew_om_roi_interp(peak_window));];
            no_peak = 1;
        else
            no_peak = 0;
        end
        this_peak = max(this_peak);
        dayN_om_peak = [dayN_om_peak, this_peak];
        
        %find a baseline window
        dayN_om_baseline = [dayN_om_baseline, mean(rew_om_roi_interp(1:cue_time_interp-100))];
        
        %plot each animal's mean trace
        subplot(2,2,4); hold on;
        plot(rew_om_roi_interp, color_pallette{session_num});
        
    end
end
%calculate means and sems
% [dayN_rew_peak_mean, dayN_rew_peak_sem] = get_mean_and_sem(dayN_rew_peak);
% [dayN_om_peak_mean, dayN_om_peak_sem] = get_mean_and_sem(dayN_om_peak);
% [dayN_rew_pb_diff_mean, dayN_rew_pb_diff_sem] = get_mean_and_sem(dayN_rew_peak-dayN_rew_baseline);
% [dayN_om_pb_diff_mean, dayN_om_pb_diff_sem] = get_mean_and_sem(dayN_om_peak-dayN_om_baseline);
% 
% subplot(2,2,1);
% scatter(dayN_rew_peak, dayN_om_peak); hold on;
% plot(dayN_rew_peak_mean, dayN_om_peak_mean, 'or'); 
% plot([dayN_rew_peak_mean-dayN_rew_peak_sem, dayN_rew_peak_mean+dayN_rew_peak_sem], [dayN_om_peak_mean,dayN_om_peak_mean], 'r');
% plot([dayN_rew_peak_mean, dayN_rew_peak_mean], [dayN_om_peak_mean-dayN_om_peak_sem, dayN_om_peak_mean+dayN_om_peak_sem], 'r');
% xlabel('rewarded trials'); ylabel('omission trials');
% title('day N peak vals');
% plot([-0.1, 0.3], [-0.1, 0.3], 'k');
% xlim([-0.1 0.3]); vline(0, 'k'); hline(0,'k');
% 
% subplot(2,2,2);
% scatter(dayN_rew_peak-dayN_rew_baseline, dayN_om_peak-dayN_om_baseline); hold on;
% plot(dayN_rew_pb_diff_mean, dayN_om_pb_diff_mean, 'or');
% plot([dayN_rew_pb_diff_mean-dayN_rew_pb_diff_sem, dayN_rew_pb_diff_mean+dayN_rew_pb_diff_sem], [dayN_om_pb_diff_mean,dayN_om_pb_diff_mean], 'r');
% plot([dayN_rew_pb_diff_mean, dayN_rew_pb_diff_mean], [dayN_om_pb_diff_mean-dayN_om_pb_diff_sem,dayN_om_pb_diff_mean+dayN_om_pb_diff_sem], 'r');
% xlabel('rewarded trials'); ylabel('omission trials');
% title('day N peak to baseline diff');
% plot([-0.1, 0.3], [-0.1, 0.3], 'k');
% xlim([-0.1 0.3]); vline(0, 'k'); hline(0,'k');
% 
% subplot(2,2,3); 
% title('Day N rewarded mean trace for each animal'); ylabel('df/f');
% ylim([-0.1 0.3]);
% vline(501, 'k'); vline(1101, 'k');
% subplot(2,2,4);
% title('Day N omission mean trace for each animal'); ylabel('df/f');
% ylim([-0.1 0.3]);
% vline(501, 'k'); vline(1101, 'k');
% 
% %make sure all sessions have same number of trials
% dayN_post_rew_mags = min_trial_finder(dayN_post_rew_mags);
% dayN_pre_rew_mags = min_trial_finder(dayN_pre_rew_mags);
% dayN_post_om_mags = min_trial_finder(dayN_post_om_mags);
% dayN_pre_om_mags = min_trial_finder(dayN_pre_om_mags);

%calculate mean and sem in the pre and post reward windows
dayN_post_rew_mags_mean = cellfun(@mean, dayN_post_rew_mags);
dayN_pre_rew_mags_mean = cellfun(@mean, dayN_pre_rew_mags);
dayN_post_om_mags_mean = cellfun(@mean, dayN_post_om_mags);
dayN_pre_om_mags_mean = cellfun(@mean, dayN_pre_om_mags);

%remove nans 
dayN_post_rew_mags_mean = dayN_post_rew_mags_mean(~isnan(dayN_post_rew_mags_mean));
dayN_pre_rew_mags_mean  = dayN_pre_rew_mags_mean(~isnan(dayN_pre_rew_mags_mean));
dayN_post_om_mags_mean  = dayN_post_om_mags_mean(~isnan(dayN_post_om_mags_mean));
dayN_pre_om_mags_mean   = dayN_pre_om_mags_mean(~isnan(dayN_pre_om_mags_mean));

%plot pre- post-reward magnitudes
figure; subplot(1,2,1);
scatter(dayN_pre_rew_mags_mean, dayN_post_rew_mags_mean); hold on;
plot(mean(dayN_pre_rew_mags_mean), mean(dayN_post_rew_mags_mean), 'or');
errorbarxy(mean(dayN_pre_rew_mags_mean), mean(dayN_post_rew_mags_mean), std(dayN_pre_rew_mags_mean)/sqrt(length(dayN_pre_rew_mags_mean)), std(dayN_post_rew_mags_mean)/sqrt(length(dayN_post_rew_mags_mean)));
ylim([-0.1 0.1]); xlim([-0.1 0.1]); hold on;
plot([-0.1, 0.1], [-0.1, 0.1], 'k');
hline(0, 'k'); vline(0,'k');
xlabel('df/f magnitude in pre-reward window');
ylabel('df/f magnitude in post-reward window');
title('day N: rewarded trials');
subplot(1,2,2);
scatter(dayN_pre_om_mags_mean, dayN_post_om_mags_mean); hold on;
plot(mean(dayN_pre_om_mags_mean), mean(dayN_post_om_mags_mean), 'or');
errorbarxy(mean(dayN_pre_om_mags_mean), mean(dayN_post_om_mags_mean), std(dayN_pre_om_mags_mean)/sqrt(length(dayN_pre_om_mags_mean)), std(dayN_post_om_mags_mean)/sqrt(length(dayN_post_om_mags_mean)));
ylim([-0.1 0.1]); xlim([-0.1 0.1]); hold on;
plot([-0.1, 0.1], [-0.1, 0.1], 'k');
hline(0, 'k'); vline(0,'k');
xlabel('df/f magnitude in pre-reward window');
ylabel('df/f magnitude in post-reward window');
title('day N: omission trials');

%calculate ttest
% [hh, pp] = ttest(dayN_pre_rew_mags_mean, dayN_post_rew_mags_mean);
% disp(['mean of PL pre reward window magnitude = ', num2str(mean(dayN_pre_rew_mags_mean)), ' +/- ', num2str(std(dayN_pre_rew_mags_mean)/sqrt(length(dayN_pre_rew_mags_mean)))]);
% disp(['mean of PL post reward window magnitude = ', num2str(mean(dayN_post_rew_mags_mean)), ' +/- ', num2str(std(dayN_post_rew_mags_mean)/sqrt(length(dayN_post_rew_mags_mean)))]);
% disp(['ttest for diff between pre and post reward windows on PL rewarded trials: h=', num2str(hh), ' p=', num2str(pp)]);
% 
% [hh, pp] = ttest(dayN_pre_om_mags_mean, dayN_post_om_mags_mean);
% disp(['mean of PL pre omission window magnitude = ', num2str(mean(dayN_pre_om_mags_mean)), ' +/- ', num2str(std(dayN_pre_om_mags_mean)/sqrt(length(dayN_pre_om_mags_mean)))]);
% disp(['mean of PL post omission window magnitude = ', num2str(mean(dayN_post_om_mags_mean)), ' +/- ', num2str(std(dayN_post_om_mags_mean)/sqrt(length(dayN_post_om_mags_mean)))]);
% disp(['ttest for diff between pre and post reward windows on PL omission trials: h=', num2str(hh), ' p=', num2str(pp)]);

%try calculating mean from all trials
[dayN_post_rew_mags_ttmean, dayN_post_rew_mags_ttsem] = get_mean_and_sem(cell2mat(dayN_post_rew_mags));
[dayN_pre_rew_mags_ttmean, dayN_pre_rew_mags_ttsem] = get_mean_and_sem(cell2mat(dayN_pre_rew_mags));
[dayN_post_om_mags_ttmean, dayN_post_om_mags_ttsem] = get_mean_and_sem(cell2mat(dayN_post_om_mags));
[dayN_pre_om_mags_ttmean, dayN_pre_om_mags_ttsem] = get_mean_and_sem(cell2mat(dayN_pre_om_mags));
[rew_h, rew_p] = ttest(cell2mat(dayN_post_rew_mags), cell2mat(dayN_pre_rew_mags));
[om_h, om_p] = ttest(cell2mat(dayN_post_om_mags), cell2mat(dayN_pre_om_mags));
disp(['mean of dayN pre reward window magnitude = ', num2str(dayN_pre_rew_mags_ttmean), ' +/- ', num2str(dayN_pre_rew_mags_ttsem)]);
disp(['mean of dayN post reward window magnitude = ', num2str(dayN_post_rew_mags_ttmean), ' +/- ', num2str(dayN_post_rew_mags_ttsem)]);
disp(['ttest for diff between pre and post reward windows on day N rewarded trials: h=', num2str(rew_h), ' p=', num2str(rew_p)]);
disp(['mean of dayN pre omission window magnitude = ', num2str(dayN_pre_om_mags_ttmean), ' +/- ', num2str(dayN_pre_om_mags_ttsem)]);
disp(['mean of dayN post omission window magnitude = ', num2str(dayN_post_om_mags_ttmean), ' +/- ', num2str(dayN_post_om_mags_ttsem)]);
disp(['ttest for diff between pre and post omission windows on day N omission trials: h=', num2str(om_h), ' p=', num2str(om_p)]);

%% Day N+1: expected vs unexpecteds
figure;
for session_num = 1:length(days_UR);
    if  exist([mean_TC_dir, 'day_N+1_vars_', days_UR{session_num}, '_all.mat'])
        
        load([mean_TC_dir, 'day_N+1_vars_', days_UR{session_num}, '_all']);
        if size(rew_roi_interp,2) < 7 || size(unexp_rew_roi_interp,2) < 7  %must have a minimum of two trials to count...  %actual minimum trial number is 7 though
            continue
        end
        
        %log the mean and sem values for pre-reward and post-reward windows
        dayU_pre_rew_mags{session_num} = mean(rew_roi_interp([pre_rew_win],:)); 
        dayU_post_rew_mags{session_num} = mean(rew_roi_interp([post_rew_win],:)); 
        
        %get the mean TC and the derivative of that mean
        rew_roi_interp = mean(rew_roi_interp, 2)';
        
        %find the peak mag in the mean TC
        this_peak = findpeaks(rew_roi_interp(peak_window));
        if size(this_peak) == 1 & this_peak < max(rew_roi_interp(1:cue_time_interp)) %if the peak magnitude was less than the max value before the cue then do not count the peak
            this_peak = [];
        end
        if isempty(this_peak) %if there was no peak in peak_window then simply take the max in the peak window
            this_peak = max(rew_roi_interp(peak_window));
            no_peak = 1;
        else
            no_peak = 0;
        end
        this_peak = max(this_peak);  %if there are two peaks then take the later one
        dayUR_exp_peak = [dayUR_exp_peak, this_peak];
        
        %find a baseline window
        dayUR_exp_baseline = [dayUR_exp_baseline, mean(rew_roi_interp(1:cue_time_interp-100))];
        
        %plot each animal's mean trace
        subplot(2,2,3); hold on;
        plot(rew_roi_interp, color_pallette{session_num});
        
        %Unexpecte reward -----------------------------------------------
        
        %log the mean and sem values for pre-reward and post-reward windows
        dayU_pre_ur_mags{session_num} = mean(unexp_rew_roi_interp([pre_rew_win],:)); 
        dayU_post_ur_mags{session_num} = mean(unexp_rew_roi_interp([post_rew_win],:)); 
        
        %get the mean TC and the derivative of that mean
        unexp_rew_roi_interp = mean(unexp_rew_roi_interp, 2)';
        
        %find the peak times in the mean TC
        this_peak = findpeaks(unexp_rew_roi_interp(peak_window));
        if size(this_peak) == 1 & this_peak < max(unexp_rew_roi_interp(1:cue_time_interp)) %if the peak magnitude was less than the max value before the cue then do not count the peak
            this_peak = [];
        end
        if isempty(this_peak) %if there was no peak in peak_window then simply take the max in the peak window
            this_peak = max(unexp_rew_roi_interp(peak_window));
            no_peak = 1;
        else
            no_peak = 0;
        end
        this_peak = max(this_peak);  %if there are two peaks then take the later one
        dayUR_unexp_peak = [dayUR_unexp_peak, this_peak];
        
        %find a baseline window
        dayUR_unexp_baseline = [dayUR_unexp_baseline, mean(unexp_rew_roi_interp(1:cue_time_interp-100))];
        
        %plot each animal's mean trace
        subplot(2,2,4); hold on;
        plot(unexp_rew_roi_interp, color_pallette{session_num});
    end
end
% %calculate means and sems
% [dayUR_rew_peak_mean, dayUR_rew_peak_sem] = get_mean_and_sem(dayUR_exp_peak);
% [dayUR_om_peak_mean, dayUR_om_peak_sem] = get_mean_and_sem(dayUR_unexp_peak);
% [dayUR_rew_pb_diff_mean, dayUR_rew_pb_diff_sem] = get_mean_and_sem(dayUR_exp_peak-dayUR_exp_baseline);
% [dayUR_om_pb_diff_mean, dayUR_om_pb_diff_sem] = get_mean_and_sem(dayUR_unexp_peak-dayUR_unexp_baseline);
% 
% subplot(2,2,1);
% scatter(dayUR_exp_peak, dayUR_unexp_peak); hold on;
% plot(dayUR_rew_peak_mean, dayUR_om_peak_mean, 'or'); 
% plot([dayUR_rew_peak_mean-dayUR_rew_peak_sem, dayUR_rew_peak_mean+dayUR_rew_peak_sem], [dayUR_om_peak_mean,dayUR_om_peak_mean], 'r');
% plot([dayUR_rew_peak_mean, dayUR_rew_peak_mean], [dayUR_om_peak_mean-dayUR_om_peak_sem, dayUR_om_peak_mean+dayUR_om_peak_sem], 'r');
% xlabel('cued trials'); ylabel('unexpected trials');
% title('day UR peak vals');
% plot([-0.1, 0.3], [-0.1, 0.3], 'k');
% xlim([-0.1 0.3]); vline(0, 'k'); hline(0,'k');
% 
% subplot(2,2,2);
% scatter(dayUR_exp_peak-dayUR_exp_baseline, dayUR_unexp_peak-dayUR_unexp_baseline); hold on;
% plot(dayUR_rew_pb_diff_mean, dayUR_om_pb_diff_mean, 'or');
% plot([dayUR_rew_pb_diff_mean-dayUR_rew_pb_diff_sem, dayUR_rew_pb_diff_mean+dayUR_rew_pb_diff_sem], [dayUR_om_pb_diff_mean,dayUR_om_pb_diff_mean], 'r');
% plot([dayUR_rew_pb_diff_mean, dayUR_rew_pb_diff_mean], [dayUR_om_pb_diff_mean-dayUR_om_pb_diff_sem,dayUR_om_pb_diff_mean+dayUR_om_pb_diff_sem], 'r');
% xlabel('cued trials'); ylabel('unexpected trials');
% title('day UR peak to baseline diff');
% plot([-0.1, 0.3], [-0.1, 0.3], 'k');
% xlim([-0.1 0.3]); vline(0, 'k'); hline(0,'k');
% 
% subplot(2,2,3); hold on;
% title('Day UR cued mean trace for each animal'); ylabel('df/f');
% vline(501, 'k'); vline(1101, 'k');
% ylim([-0.1 0.3]);
% subplot(2,2,4); hold on;
% title('Day UR unexpected mean trace for each animal'); ylabel('df/f');
% vline(501, 'k'); vline(1101, 'k');
% ylim([-0.1 0.3]);

%----------------------------------------------------------------------
% %make sure all sessions have same number of trials
% dayU_post_rew_mags = min_trial_finder(dayU_post_rew_mags);
% dayU_pre_rew_mags = min_trial_finder(dayU_pre_rew_mags);
% dayU_post_ur_mags = min_trial_finder(dayU_post_ur_mags);
% dayU_pre_ur_mags = min_trial_finder(dayU_pre_ur_mags);

%calculate mean and sem in the pre and post reward windows
dayUR_post_rew_mags_mean = cellfun(@mean, dayU_post_rew_mags);
dayUR_pre_rew_mags_mean = cellfun(@mean, dayU_pre_rew_mags);
dayUR_post_om_mags_mean = cellfun(@mean, dayU_post_ur_mags);
dayUR_pre_om_mags_mean = cellfun(@mean, dayU_pre_ur_mags);

%remove nans 
dayUR_post_rew_mags_mean = dayUR_post_rew_mags_mean(~isnan(dayUR_post_rew_mags_mean));
dayUR_pre_rew_mags_mean = dayUR_pre_rew_mags_mean(~isnan(dayUR_pre_rew_mags_mean));
dayUR_post_om_mags_mean = dayUR_post_om_mags_mean(~isnan(dayUR_post_om_mags_mean));
dayUR_pre_om_mags_mean = dayUR_pre_om_mags_mean(~isnan(dayUR_pre_om_mags_mean));

%plot pre- post-reward magnitudes
figure; subplot(1,2,1);
scatter(dayUR_pre_rew_mags_mean, dayUR_post_rew_mags_mean); hold on;
plot(mean(dayUR_pre_rew_mags_mean), mean(dayUR_post_rew_mags_mean), 'or');
errorbarxy(mean(dayUR_pre_rew_mags_mean), mean(dayUR_post_rew_mags_mean), std(dayUR_pre_rew_mags_mean)/sqrt(length(dayUR_pre_rew_mags_mean)), std(dayUR_post_rew_mags_mean)/sqrt(length(dayUR_post_rew_mags_mean)));
ylim([-0.1 0.1]); xlim([-0.1 0.1]); hold on;
plot([-0.1, 0.1], [-0.1, 0.1], 'k');
hline(0, 'k'); vline(0,'k');
xlabel('df/f magnitude in pre-reward window');
ylabel('df/f magnitude in post-reward window');
title('day UR: rewarded trials');

subplot(1,2,2);
scatter(dayUR_pre_om_mags_mean, dayUR_post_om_mags_mean); hold on;
plot(mean(dayUR_pre_om_mags_mean), mean(dayUR_post_om_mags_mean), 'or');
errorbarxy(mean(dayUR_pre_om_mags_mean), mean(dayUR_post_om_mags_mean), std(dayUR_pre_om_mags_mean)/sqrt(length(dayUR_pre_om_mags_mean)), std(dayUR_post_om_mags_mean)/sqrt(length(dayUR_post_om_mags_mean)));
ylim([-0.1 0.1]); xlim([-0.1 0.1]); hold on;
plot([-0.1, 0.1], [-0.1, 0.1], 'k');
hline(0, 'k'); vline(0,'k');
xlabel('df/f magnitude in pre-reward window');
ylabel('df/f magnitude in post-reward window');
title('day UR: omission trials');

%calculate ttest
% [hh, pp] = ttest(dayUR_pre_rew_mags_mean, dayUR_post_rew_mags_mean);
% disp(['mean of dayUR pre reward window magnitude = ', num2str(mean(dayUR_pre_rew_mags_mean)), ' +/- ', num2str(std(dayUR_pre_rew_mags_mean)/sqrt(length(dayUR_pre_rew_mags_mean)))]);
% disp(['mean of dayUR post reward window magnitude = ', num2str(mean(dayUR_post_rew_mags_mean)), ' +/- ', num2str(std(dayUR_post_rew_mags_mean)/sqrt(length(dayUR_post_rew_mags_mean)))]);
% disp(['ttest for diff between pre and post reward windows on day UR rewarded trials: h=', num2str(hh), ' p=', num2str(pp)]);
% 
% [hh, pp] = ttest(dayUR_pre_om_mags_mean, dayUR_post_om_mags_mean);
% disp(['mean of dayUR pre reward window magnitude = ', num2str(mean(dayUR_pre_om_mags_mean)), ' +/- ', num2str(std(dayUR_pre_om_mags_mean)/sqrt(length(dayUR_pre_om_mags_mean)))]);
% disp(['mean of dayUR post reward window magnitude = ', num2str(mean(dayUR_post_om_mags_mean)), ' +/- ', num2str(std(dayUR_post_om_mags_mean)/sqrt(length(dayUR_post_om_mags_mean)))]);
% disp(['ttest for diff between pre and post reward windows on day UR omission trials: h=', num2str(hh), ' p=', num2str(pp)]);

%try calculating mean from all trials
[dayU_post_rew_mags_ttmean, dayU_post_rew_mags_ttsem] = get_mean_and_sem(cell2mat(dayU_post_rew_mags));
[dayU_pre_rew_mags_ttmean, dayU_pre_rew_mags_ttsem] = get_mean_and_sem(cell2mat(dayU_pre_rew_mags));
[dayU_post_ur_mags_ttmean, dayU_post_ur_mags_ttsem] = get_mean_and_sem(cell2mat(dayU_post_ur_mags));
[dayU_pre_ur_mags_ttmean, dayU_pre_ur_mags_ttsem] = get_mean_and_sem(cell2mat(dayU_pre_ur_mags));
[rew_h, rew_p] = ttest(cell2mat(dayU_post_rew_mags), cell2mat(dayU_pre_rew_mags));
[om_h, om_p] = ttest(cell2mat(dayU_post_ur_mags), cell2mat(dayU_pre_ur_mags));
disp(['mean of dayU pre reward window magnitude = ', num2str(dayU_pre_rew_mags_ttmean), ' +/- ', num2str(dayU_pre_rew_mags_ttsem)]);
disp(['mean of dayU post reward window magnitude = ', num2str(dayU_post_rew_mags_ttmean), ' +/- ', num2str(dayU_post_rew_mags_ttsem)]);
disp(['ttest for diff between pre and post reward windows on day U rewarded trials: h=', num2str(rew_h), ' p=', num2str(rew_p)]);
disp(['mean of dayU pre unexp window magnitude = ', num2str(dayU_pre_ur_mags_ttmean), ' +/- ', num2str(dayU_pre_ur_mags_ttsem)]);
disp(['mean of dayU post unexp window magnitude = ', num2str(dayU_post_ur_mags_ttmean), ' +/- ', num2str(dayU_post_ur_mags_ttsem)]);
disp(['ttest for diff between pre and post unexp windows on day UR unexp trials: h=', num2str(om_h), ' p=', num2str(om_p)]);




