%% plot onset times of the mean interpolated TCs which were saved in WF_CRP_across_animal_TCs
clear; 
WF_CRP_list_of_days;
mean_TC_dir = 'Y:\home\jake\Analysis\Cue_reward_pairing_analysis\WF\onset latencies and peak mags\mean interpolated TCs\';
days = days_UR;
ROI_cell = days_UR_ROI;
peak_window = [800:1501];
baseline_window = [400:-1:200];
cue_time_interp = 501;
peak_diff_window= 500;
no_peak_diff_window = 1000;
yy = [1:2000];

day_1_onsets = [];
day_N_onsets = [];

day1_rew_peak = [];
day1_om_peak = [];
dayN_rew_peak = [];
dayN_om_peak = [];
dayUR_exp_peak = [];
dayUR_unexp_peak = [];

%% Rewarded trails: Day1 vs DayN
%figure; plot_num=1;
for session_num = 1:length(days_1);
    if exist([mean_TC_dir, 'day_1_vars_', days_1{session_num}, '_all.mat']) && exist([mean_TC_dir, 'day_N_vars_', days_post{session_num}, '_all.mat']);
        %load TC data
        load([mean_TC_dir, 'day_1_vars_', days_1{session_num}, '_all']);
        
        %plot each trial individually
        %figure; subplot(1,2,1); plot(rew_roi_interp); 
        %title(['day 1 first half of trials. lick excemption. mean of ROIs. n=', size(rew_roi_interp,2), ' ', days_1(session_num)]); 
        %xlabel('cue at 500'); ylabel('df/f'); vline(cue_time_interp); vline(cue_time_interp+600, 'b');
        %subplot(2,3,plot_num); plot(diff(mean(rew_roi_interp, 2)')); plot_num=plot_num+1; vline(cue_time_interp); vline(cue_time_interp+600, 'b'); title(['day 1 first half of trials. lick excemption. n=', size(rew_roi_interp,2), ' ', days_1(session_num)]);
        
        %find the mean TC and the derivative of that mean TC
        rew_roi_interp = mean(rew_roi_interp, 2)';
        rew_roi_interp_diff = diff(rew_roi_interp);
        %subplot(1,2,1); hold on; plot(rew_roi_interp); title('day 1 rewarded');
        %vline(cue_time_interp); vline(cue_time_interp+600, 'b'); title(['day 1 mean timecourse', ' ', days_1(session_num)]); %plot the mean TC
        
        %find the peak times in the mean TC
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
        
        %find the peak latency
        this_peak_lat = find(rew_roi_interp(peak_window) == this_peak, 1, 'first') + min(peak_window)-1;
        
%         %find a baseline window (this variable is not used again in this forloop)
%         if this_peak_lat-max(baseline_window) < cue_time_interp   %if the peak latency - (duration+offset) of the baseline window 
%             this_baseline = min(rew_roi_interp([cue_time_interp:this_peak_lat]));
%         else
%             this_baseline = min(rew_roi_interp([this_peak_lat-baseline_window]));
%         end

        %find the onset latency
        if no_peak == 0
            peak_win_use = peak_diff_window;
        else
            peak_win_use = no_peak_diff_window;
        end
        this_onset_lat = find(rew_roi_interp_diff == max(rew_roi_interp_diff([this_peak_lat-peak_win_use:this_peak_lat])));  %find the maximum rate of change in df/f in a 500ms (or 1000ms) window before the peak
        this_onset_lat = this_onset_lat+1; %adjust for the diff
        this_onset_lat = this_onset_lat(this_onset_lat>this_peak_lat-peak_win_use); %just in case there were two copies of the same rate. Must have occured within the window
        this_onset_lat = this_onset_lat(1); %this value will be equal to the first value where the rate reached peak
         day_1_onsets = [day_1_onsets, this_onset_lat]; 
        
        %Find onset latencies for Day N rewarded trials
        load([mean_TC_dir, 'day_N_vars_', days_post{session_num}, '_all']);
        %figure; subplot(1,2,1); plot(rew_roi_interp); title(['day N rewarded trials. lick excemption. n=', size(rew_roi_interp,2), ' ', days_post(session_num)]); xlabel('cue at 500'); ylabel('df/f'); vline(cue_time_interp); vline(cue_time_interp+600, 'b');
            %subplot(2,3,plot_num); plot(diff(mean(rew_roi_interp, 2)')); plot_num=plot_num+1; vline(cue_time_interp); vline(cue_time_interp+600, 'b'); title(['day N rewarded trials. lick excemption. n=', size(rew_roi_interp,2), ' ', days_post(session_num)]);
        rew_roi_interp = mean(rew_roi_interp, 2)';
        rew_roi_interp_diff = diff(rew_roi_interp);
        %subplot(1,2,2); hold on; plot(rew_roi_interp);  title('day N rewarded');
        %vline(cue_time_interp); vline(cue_time_interp+600, 'b'); title('mean timecourse');
        %find the peak df/f value
        this_peak = findpeaks(rew_roi_interp(peak_window));
        if size(this_peak) == 1 & this_peak < max(rew_roi_interp(1:cue_time_interp))
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
        
        %find the peak latency
        this_peak_lat = find(rew_roi_interp(peak_window) == this_peak, 1, 'first') + min(peak_window)-1;
%         if this_peak_lat-max(baseline_window) < cue_time_interp
%             this_baseline = min(rew_roi_interp([cue_time_interp:this_peak_lat]));
%         else
%             this_baseline = min(rew_roi_interp([this_peak_lat-baseline_window]));
%         end

        %find the onset latency
        if no_peak == 0
            peak_win_use = peak_diff_window;
        else
            peak_win_use = no_peak_diff_window;
        end
        this_onset_lat = find(rew_roi_interp_diff == max(rew_roi_interp_diff([this_peak_lat-peak_win_use:this_peak_lat])));  
        this_onset_lat = this_onset_lat+1; %adjust for the diff
        this_onset_lat = this_onset_lat(this_onset_lat>this_peak_lat-peak_win_use); %just in case there were two copies of the same rate. Must have occured within the window
        this_onset_lat = this_onset_lat(1); %this value will be equal to the first value where the rate reached peak
        day_N_onsets = [day_N_onsets, this_onset_lat];
    end
end
figure; hold on; 
scatter(day_1_onsets-cue_time_interp, day_N_onsets-cue_time_interp, 'k', 'filled');
plot(yy, 'k'); xlim([0 1000]); ylim([0 1000]);
xlabel('day 1'); ylabel('day N'); title('rewarded trials onset latency relative to cue onset');
[day_1_onsets_mean, day_1_onsets_sem] = get_mean_and_sem(day_1_onsets-cue_time_interp);
[day_N_onsets_mean, day_N_onsets_sem] = get_mean_and_sem(day_N_onsets-cue_time_interp);
scatter(day_1_onsets_mean, day_N_onsets_mean, 'r', 'filled');
plot([day_1_onsets_mean, day_1_onsets_mean], [day_N_onsets_mean-day_N_onsets_sem, day_N_onsets_mean+day_N_onsets_sem], 'r');
plot([day_1_onsets_mean-day_1_onsets_sem, day_1_onsets_mean+day_1_onsets_sem], [day_N_onsets_mean, day_N_onsets_mean], 'r');
[hh, pp] = ttest(day_1_onsets, day_N_onsets);


%% Day N: Rewarded vs Omission
day_N_rew_onsets = [];
day_N_om_onsets = [];
%figure;
for session_num = 1:length(days_post);
    if  exist([mean_TC_dir, 'day_N_vars_', days_post{session_num}, '_all.mat'])
        load([mean_TC_dir, 'day_N_vars_', days_post{session_num}, '_all']);
        
        if size(rew_roi_interp,2) > 2 && size(rew_om_roi_interp,2) > 2  % must have at least two trials per condition
            %get the mean and the derivative of the TC
            rew_roi_interp = mean(rew_roi_interp, 2)';
            rew_roi_interp_diff = diff(rew_roi_interp);
           % figure; plot(rew_roi_interp); title(['rewarded', days_post(session_num)]); vline(cue_time_interp); vline(cue_time_interp+600);
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
            this_peak_lat = find(rew_roi_interp(peak_window) == this_peak, 1, 'first') + min(peak_window)-1;
            if this_peak_lat-max(baseline_window) < cue_time_interp
                this_baseline = min(rew_roi_interp([cue_time_interp:this_peak_lat]));
            else
                this_baseline = min(rew_roi_interp([this_peak_lat-baseline_window]));
            end
            if no_peak == 0
                peak_win_use = peak_diff_window;
            else
                peak_win_use = no_peak_diff_window;
            end
            this_onset_lat = find(rew_roi_interp_diff == max(rew_roi_interp_diff([this_peak_lat-peak_win_use:this_peak_lat])));
            this_onset_lat = this_onset_lat+1; %adjust for the diff
            this_onset_lat = this_onset_lat(this_onset_lat>this_peak_lat-peak_win_use); %just in case there were two copies of the same rate. Must have occured within the window
            this_onset_lat = this_onset_lat(1); %this value will be equal to the first value where the rate reached peak
            day_N_rew_onsets = [day_N_rew_onsets, this_onset_lat]; 

            %Reward omission Day N+1 onset latencies
            rew_om_roi_interp = mean(rew_om_roi_interp, 2)';
            %figure; 
            %plot(rew_om_roi_interp); hold on; 
            %title(['omission', days_post(session_num)]); vline(cue_time_interp); vline(cue_time_interp+600);
            rew_om_roi_interp_diff = diff(rew_om_roi_interp);
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
            this_peak_lat = find(rew_om_roi_interp(peak_window) == this_peak, 1, 'first') + min(peak_window)-1;
            if this_peak_lat-max(baseline_window) < cue_time_interp
                this_baseline = min(rew_om_roi_interp([cue_time_interp:this_peak_lat]));
            else
                this_baseline = min(rew_om_roi_interp([this_peak_lat-baseline_window]));
            end
            
            if no_peak == 0
                peak_win_use = peak_diff_window;
            else
                peak_win_use = no_peak_diff_window;
            end
            this_onset_lat = find(rew_om_roi_interp_diff == max(rew_om_roi_interp_diff([this_peak_lat-peak_win_use:this_peak_lat])));
            this_onset_lat = this_onset_lat+1; %adjust for the diff
            this_onset_lat = this_onset_lat(this_onset_lat>this_peak_lat-peak_win_use); %just in case there were two copies of the same rate. Must have occured within the window
            this_onset_lat = this_onset_lat(1); %this value will be equal to the first value where the rate reached peak
            day_N_om_onsets = [day_N_om_onsets, this_onset_lat];
        end
    end
end
figure; hold on;
scatter(   day_N_rew_onsets-cue_time_interp,  day_N_om_onsets-cue_time_interp, 'k', 'filled');
plot(yy, 'k'); xlim([0 800]); ylim([0 800]);
xlabel('rewarded'); ylabel('omission'); title('Day N rewarded vs omission onset latencies');

[day_N_rew_onsets_mean, day_N_rew_onsets_sem] = get_mean_and_sem(day_N_rew_onsets-cue_time_interp);
[day_N_om_onsets_mean, day_N_om_onsets_sem] = get_mean_and_sem(day_N_om_onsets-cue_time_interp);
scatter(day_N_rew_onsets_mean, day_N_om_onsets_mean, 'r', 'filled');
plot([day_N_rew_onsets_mean, day_N_rew_onsets_mean], [day_N_om_onsets_mean-day_N_om_onsets_sem, day_N_om_onsets_mean+day_N_om_onsets_sem], 'r');
plot([day_N_rew_onsets_mean-day_N_rew_onsets_sem, day_N_rew_onsets_mean+day_N_rew_onsets_sem], [day_N_om_onsets_mean, day_N_om_onsets_mean], 'r');
[hh, pp] = ttest(day_N_rew_onsets, day_N_om_onsets);


%% Day N+1: expected vs unexpecteds
day_UR_rew_onsets = [];
day_UR_une_onsets = [];
figure; plot_num=1;
for session_num = 1:length(days_UR);
    if  exist([mean_TC_dir, 'day_N+1_vars_', days_UR{session_num}, '_all.mat'])
        load([mean_TC_dir, 'day_N+1_vars_', days_UR{session_num}, '_all']);
        if size(rew_roi_interp,2) > 2 && size(unexp_rew_roi_interp,2) > 2  %must have a minimum of two trials to count...  %actual minimum trial number is 7 though
               %figure; plot(rew_roi_interp); title(['expected reward. lick excemption. n=', size(rew_roi_interp,2), ' ', days_UR(session_num)]); xlabel('cue at 500'); ylabel('df/f'); vline(cue_time_interp); vline(cue_time_interp+600, 'b');
               %figure; plot(mean(rew_roi_interp, 2)');  vline(cue_time_interp); vline(cue_time_interp+600, 'b'); title(['Mean: expected reward. lick excemption. n=', size(rew_roi_interp,2), ' ', days_UR(session_num)]); xlabel('cue at 500');
              %subplot(2,3,plot_num); plot(diff(mean(rew_roi_interp, 2)')); plot_num=plot_num+1; vline(cue_time_interp); vline(cue_time_interp+600, 'b'); title(['expected reward. lick excemption. n=', size(rew_roi_interp,2), ' ', days_UR(session_num)]);
            
            %get the mean TC and the derivative of that mean
            rew_roi_interp = mean(rew_roi_interp, 2)';
            rew_roi_interp_diff = diff(rew_roi_interp);
            subplot(1,2,1); hold on; plot(rew_roi_interp); title('UR day expected rewards');
            
%             %find peak magnitude and latency in the mean TC
%             this_peak = findpeaks(rew_roi_interp(peak_window));
%             this_peak = [this_peak, max(rew_roi_interp(peak_window));]; %peak must be within the peak window
%             this_peak = max(this_peak); %take the peak with the largest magnitude 
            
            %==========================================
            %find the peak times in the mean TC
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
            %-----------------------------------------------
            
            %find the peak latency
            this_peak_lat = find(rew_roi_interp(peak_window) == this_peak, 1, 'first') + min(peak_window)-1;
            
            %find the onset latency before the peak
            if no_peak == 0
                peak_win_use = peak_diff_window;
            else
                peak_win_use = no_peak_diff_window;
            end
            this_onset_lat = find(rew_roi_interp_diff == max(rew_roi_interp_diff([this_peak_lat-peak_win_use:this_peak_lat])));
            this_onset_lat = this_onset_lat+1; %adjust for the diff
            this_onset_lat = this_onset_lat(this_onset_lat>this_peak_lat-peak_win_use); %just in case there were two copies of the same rate. Must have occured within the window
            this_onset_lat = this_onset_lat(1); %this value will be equal to the first value where the rate reached peak
            day_UR_rew_onsets = [day_UR_rew_onsets, this_onset_lat];
            
            
               %figure; subplot(1,2,1); plot(mean(unexp_rew_roi_interp, 2)');  vline(cue_time_interp+600, 'b'); title(['Mean: unexpected reward. lick excemption. n=', size(unexp_rew_roi_interp,2), ' ', days_UR(session_num)]); 
               %subplot(1,2,2); plot(diff(mean(unexp_rew_roi_interp, 2)')); plot_num=plot_num+1; vline(cue_time_interp); vline(cue_time_interp+600, 'b'); title(['unexpected reward. lick excemption. n=', size(rew_roi_interp,2), ' ', days_UR(session_num)]);
            
            %get the mean TC and the derivative of that mean
            unexp_rew_roi_interp = mean(unexp_rew_roi_interp, 2)';
            unexp_rew_roi_interp_diff = diff(unexp_rew_roi_interp);
            subplot(1,2,2); hold on; plot(unexp_rew_roi_interp); title('UR day Unexpected');
            
%             %find peak magnitude and latency in the mean TC
%             this_peak = findpeaks(unexp_rew_roi_interp(peak_window));
%             this_peak = [this_peak, max(unexp_rew_roi_interp(peak_window));];
%             this_peak = max(this_peak);
            
            %====================================
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
            %=======================================
            
            %find the peak latency
            this_peak_lat = find(unexp_rew_roi_interp(peak_window) == this_peak, 1, 'first') + min(peak_window)-1;
            if this_peak_lat < 1102
                this_peak = max(unexp_rew_roi_interp(peak_window));
                this_peak_lat = find(unexp_rew_roi_interp(peak_window) == this_peak, 1, 'first') + min(peak_window)-1;
            end

            %find the onset latency before the peak
            if no_peak == 0
                peak_win_use = peak_diff_window;
            else
                peak_win_use = no_peak_diff_window;
            end
            this_onset_lat = find(unexp_rew_roi_interp_diff == max(unexp_rew_roi_interp_diff([this_peak_lat-peak_win_use:this_peak_lat])));
            this_onset_lat = this_onset_lat+1; %adjust for the diff
            this_onset_lat = this_onset_lat(this_onset_lat>this_peak_lat-peak_win_use); %just in case there were two copies of the same rate. Must have occured within the window
            this_onset_lat = this_onset_lat(1); %this value will be equal to the first value where the rate reached peak
            day_UR_une_onsets = [day_UR_une_onsets, this_onset_lat];
        end
    end
end
figure; hold on;
scatter(   day_UR_rew_onsets-cue_time_interp-600,  day_UR_une_onsets-cue_time_interp-600, 'k', 'filled');
plot([-500:500],[-500:500], 'k'); xlim([-500 500]); ylim([-500 500]);
xlabel('rewarded'); ylabel('unexpected'); title('Day N+1 rewarded vs unexpected onset latencies');

[day_UR_rew_onsets_mean, day_UR_rew_onsets_sem] = get_mean_and_sem(day_UR_rew_onsets-cue_time_interp-600);
[day_UR_une_onsets_mean, day_UR_une_onsets_sem] = get_mean_and_sem(day_UR_une_onsets-cue_time_interp-600);
scatter(day_UR_rew_onsets_mean, day_UR_une_onsets_mean, 'r', 'filled');
plot([day_UR_rew_onsets_mean, day_UR_rew_onsets_mean], [day_UR_une_onsets_mean-day_UR_une_onsets_sem, day_UR_une_onsets_mean+day_UR_une_onsets_sem], 'r');
plot([day_UR_rew_onsets_mean-day_UR_rew_onsets_sem, day_UR_rew_onsets_mean+day_UR_rew_onsets_sem], [day_UR_une_onsets_mean, day_UR_une_onsets_mean], 'r');
hline(0, 'k'); vline(0, 'k');
[hh, pp] = ttest(day_UR_rew_onsets, day_UR_une_onsets);



