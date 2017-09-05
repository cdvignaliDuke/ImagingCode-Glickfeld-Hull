%% plot onset times of the mean interpolated TCs which were saved in WF_CRP_across_animal_TCs
clear; 
WF_CRP_list_of_days;
mean_TC_dir = 'Z:\Analysis\Cue_reward_pairing_analysis\WF\onset latencies and peak mags\mean interpolated TCs\';
days = days_UR;
ROI_cell = days_UR_ROI;
peak_window = [800:1501];
baseline_window = [400:-1:200];
cue_time_interp = 501;
yy = [1:2000];

% day_1_onsets = [];
% day_N_onsets = [];
% %Rewarded trails: Day1 vs DayN
% for session_num = 1:length(days_1);
%     if exist([mean_TC_dir, 'day_1_vars_', days_1{session_num}, '.mat']) && exist([mean_TC_dir, 'day_N_vars_', days_post{session_num}, '.mat']);
%         load([mean_TC_dir, 'day_1_vars_', days_1{session_num}]);
%         rew_roi_interp = mean(rew_roi_interp, 2)';
%         this_peak = findpeaks(rew_roi_interp(peak_window));
%         if size(this_peak) == 1 & this_peak < max(rew_roi_interp(1:cue_time_interp))
%             this_peak = [];
%         end
%         if isempty(this_peak)
%             this_peak = max(rew_roi_interp(peak_window));
%             this_peak_lat = find(rew_roi_interp(peak_window) == this_peak, 1, 'first') + min(peak_window)-1;
%             this_baseline = min(rew_roi_interp([this_peak_lat-baseline_window]));
%         else 
%             this_peak = max(this_peak);
%             this_peak_lat = find(rew_roi_interp(peak_window) == this_peak, 1, 'first') + min(peak_window)-1;
%             if this_peak_lat-max(baseline_window) < cue_time_interp
%                 this_baseline = min(rew_roi_interp([cue_time_interp:this_peak_lat]));
%             else
%                 this_baseline = min(rew_roi_interp([this_peak_lat-baseline_window]));
%             end
%         end
%         days_1(session_num)
%         this_peak_lat
%         this_peak_base_diff = this_peak - this_baseline; assert(this_peak_base_diff > 0);
%         this_ten_val = 0.1*this_peak_base_diff+this_baseline;
%         this_onset_lat = find( rew_roi_interp(1:this_peak_lat)<this_ten_val, 1, 'last' );
%         day_1_onsets = [day_1_onsets, this_onset_lat]; this_onset_lat
%         
%         
%         
%         load([mean_TC_dir, 'day_N_vars_', days_post{session_num}]);
%         load([mean_TC_dir, 'day_N_vars_', days_post{session_num}]);
%         rew_roi_interp = mean(rew_roi_interp, 2)';
%         this_peak = findpeaks(rew_roi_interp(peak_window));
%         if size(this_peak) == 1 & this_peak < max(rew_roi_interp(1:cue_time_interp))
%             this_peak = [];
%         end
%         if isempty(this_peak)
%             this_peak = max(rew_roi_interp(peak_window));
%             this_peak_lat = find(rew_roi_interp(peak_window) == this_peak, 1, 'first') + min(peak_window)-1;
%             this_baseline = min(rew_roi_interp([this_peak_lat-baseline_window]));
%         else
%             this_peak = max(this_peak);
%             this_peak_lat = find(rew_roi_interp(peak_window) == this_peak, 1, 'first') + min(peak_window)-1;
%             if this_peak_lat-max(baseline_window) < cue_time_interp
%                 this_baseline = min(rew_roi_interp([cue_time_interp:this_peak_lat]));
%             else
%                 this_baseline = min(rew_roi_interp([this_peak_lat-baseline_window]));
%             end
%         end
%         this_peak_base_diff = this_peak - this_baseline; assert(this_peak_base_diff > 0);
%         this_ten_val = 0.1*this_peak_base_diff+this_baseline;
%         this_onset_lat = find( rew_roi_interp(1:this_peak_lat)<this_ten_val, 1, 'last' );
%         day_N_onsets = [day_N_onsets, this_onset_lat];
%     end
% end
% figure; hold on; 
% scatter(day_1_onsets-cue_time_interp, day_N_onsets-cue_time_interp, 'k', 'filled');
% plot(yy, 'k'); xlim([0 800]); ylim([0 800]);
% xlabel('day 1'); ylabel('day N'); title('rewarded trials onset latency relative to cue onset');
% [day_1_onsets_mean, day_1_onsets_sem] = get_mean_and_sem(day_1_onsets-cue_time_interp);
% [day_N_onsets_mean, day_N_onsets_sem] = get_mean_and_sem(day_N_onsets-cue_time_interp);
% scatter(day_1_onsets_mean, day_N_onsets_mean, 'r', 'filled');
% plot([day_1_onsets_mean, day_1_onsets_mean], [day_N_onsets_mean-day_N_onsets_sem, day_N_onsets_mean+day_N_onsets_sem], 'r');
% plot([day_1_onsets_mean-day_1_onsets_sem, day_1_onsets_mean+day_1_onsets_sem], [day_N_onsets_mean, day_N_onsets_mean], 'r');
% 


% Day N: Rewarded vs Omission
day_N_rew_onsets = [];
day_N_om_onsets = [];
for session_num = 1:length(days_post);
    if  exist([mean_TC_dir, 'day_N_vars_', days_post{session_num}, '.mat'])
        load([mean_TC_dir, 'day_N_vars_', days_post{session_num}]);
        
        if size(rew_roi_interp,2) > 2 && size(rew_om_roi_interp,2) > 2
            rew_roi_interp = mean(rew_roi_interp, 2)';
            this_peak = findpeaks(rew_roi_interp(peak_window));
            if size(this_peak) == 1 & this_peak < max(rew_roi_interp(1:cue_time_interp))
                this_peak = [];
            end
            if isempty(this_peak)
                this_peak = max(rew_roi_interp(peak_window));
                this_peak_lat = find(rew_roi_interp(peak_window) == this_peak, 1, 'first') + min(peak_window)-1;
                this_baseline = min(rew_roi_interp([this_peak_lat-baseline_window]));
            else
                this_peak = max(this_peak);
                this_peak_lat = find(rew_roi_interp(peak_window) == this_peak, 1, 'first') + min(peak_window)-1;
                if this_peak_lat-max(baseline_window) < cue_time_interp
                    this_baseline = min(rew_roi_interp([cue_time_interp:this_peak_lat]));
                else
                    this_baseline = min(rew_roi_interp([this_peak_lat-baseline_window]));
                end
            end
            this_peak_base_diff = this_peak - this_baseline; assert(this_peak_base_diff > 0);
            this_ten_val = 0.1*this_peak_base_diff+this_baseline;
            this_onset_lat = find( rew_roi_interp(1:this_peak_lat)<this_ten_val, 1, 'last' );
            day_N_rew_onsets = [day_N_rew_onsets, this_onset_lat];
            
            
            rew_om_roi_interp = mean(rew_om_roi_interp, 2)';
            this_peak = findpeaks(rew_om_roi_interp(peak_window));
            if size(this_peak) == 1 & this_peak < max(rew_om_roi_interp(1:cue_time_interp))
                this_peak = [];
            end
            if isempty(this_peak)
                this_peak = max(rew_om_roi_interp(peak_window));
                this_peak_lat = find(rew_om_roi_interp(peak_window) == this_peak, 1, 'first') + min(peak_window)-1;
                this_baseline = min(rew_om_roi_interp([this_peak_lat-baseline_window]));
            else
                this_peak = max(this_peak);
                this_peak_lat = find(rew_om_roi_interp(peak_window) == this_peak, 1, 'first') + min(peak_window)-1;
                if this_peak_lat-max(baseline_window) < cue_time_interp
                    this_baseline = min(rew_om_roi_interp([cue_time_interp:this_peak_lat]));
                else
                    this_baseline = min(rew_om_roi_interp([this_peak_lat-baseline_window]));
                end
            end
            this_peak_base_diff = this_peak - this_baseline; assert(this_peak_base_diff > 0);
            this_ten_val = 0.1*this_peak_base_diff+this_baseline;
            this_onset_lat = find( rew_om_roi_interp(1:this_peak_lat)<this_ten_val, 1, 'last' );
            day_N_om_onsets = [day_N_om_onsets, this_onset_lat];
            days_post(session_num)
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



% Day N+1: expected vs unexpecteds
% day_UR_rew_onsets = [];
% day_UR_une_onsets = [];
% for session_num = 1:length(days_UR);
%     if  exist([mean_TC_dir, 'day_N+1_vars_', days_UR{session_num}, '.mat'])
%         load([mean_TC_dir, 'day_N+1_vars_', days_UR{session_num}]);
%         if size(rew_roi_interp,2) > 2 && size(unexp_rew_roi_interp,2) > 2
%             rew_roi_interp = mean(rew_roi_interp, 2)';
%             this_peak = findpeaks(rew_roi_interp(peak_window));
%             if isempty(this_peak)
%                 this_peak = max(rew_roi_interp(peak_window));
%                 this_peak_lat = find(rew_roi_interp(peak_window) == this_peak, 1, 'first') + min(peak_window)-1;
%                 this_baseline = min(rew_roi_interp([this_peak_lat-baseline_window]));
%             else
%                 this_peak = max(this_peak);
%                 this_peak_lat = find(rew_roi_interp(peak_window) == this_peak, 1, 'first') + min(peak_window)-1;
%                 if this_peak_lat-max(baseline_window) < cue_time_interp
%                     this_baseline = min(rew_roi_interp([cue_time_interp:this_peak_lat]));
%                 else
%                     this_baseline = min(rew_roi_interp([this_peak_lat-baseline_window]));
%                 end
%             end
%             this_peak_base_diff = this_peak - this_baseline; assert(this_peak_base_diff > 0);
%             this_ten_val = 0.1*this_peak_base_diff+this_baseline;
%             this_onset_lat = find( rew_roi_interp(1:this_peak_lat)<this_ten_val, 1, 'last' );
%             day_UR_rew_onsets = [day_UR_rew_onsets, this_onset_lat];
%             
%             
%             
%             unexp_rew_roi_interp = mean(unexp_rew_roi_interp, 2)';
%             this_peak = findpeaks(unexp_rew_roi_interp(peak_window));
%             if isempty(this_peak)
%                 this_peak = max(unexp_rew_roi_interp(peak_window));
%                 this_peak_lat = find(unexp_rew_roi_interp(peak_window) == this_peak, 1, 'first') + min(peak_window)-1;
%                 this_baseline = min(unexp_rew_roi_interp([this_peak_lat-baseline_window]));
%             else
%                 this_peak = max(this_peak);
%                 this_peak_lat = find(unexp_rew_roi_interp(peak_window) == this_peak, 1, 'first') + min(peak_window)-1;
%                 if this_peak_lat-max(baseline_window) < cue_time_interp
%                     this_baseline = min(unexp_rew_roi_interp([cue_time_interp:this_peak_lat]));
%                 else
%                     this_baseline = min(unexp_rew_roi_interp([this_peak_lat-baseline_window]));
%                 end
%             end
%             this_peak_base_diff = this_peak - this_baseline; assert(this_peak_base_diff > 0);
%             this_ten_val = 0.1*this_peak_base_diff+this_baseline;
%             this_onset_lat = find( unexp_rew_roi_interp(1:this_peak_lat)<this_ten_val, 1, 'last' );
%             day_UR_une_onsets = [day_UR_une_onsets, this_onset_lat];
%         end
%     end
% end
% figure; hold on;
% scatter(   day_UR_rew_onsets-cue_time_interp,  day_UR_une_onsets-cue_time_interp, 'k', 'filled');
% plot(yy, 'k'); xlim([0 800]); ylim([0 800]);
% xlabel('rewarded'); ylabel('unexpected'); title('Day N+1 rewarded vs unexpected onset latencies');
% 
% [day_UR_rew_onsets_mean, day_UR_rew_onsets_sem] = get_mean_and_sem(day_UR_rew_onsets-cue_time_interp);
% [day_UR_une_onsets_mean, day_UR_une_onsets_sem] = get_mean_and_sem(day_UR_une_onsets-cue_time_interp);
% scatter(day_UR_rew_onsets_mean, day_UR_une_onsets_mean, 'r', 'filled');
% plot([day_UR_rew_onsets_mean, day_UR_rew_onsets_mean], [day_UR_une_onsets_mean-day_UR_une_onsets_sem, day_UR_une_onsets_mean+day_UR_une_onsets_sem], 'r');
% plot([day_UR_rew_onsets_mean-day_UR_rew_onsets_sem, day_UR_rew_onsets_mean+day_UR_rew_onsets_sem], [day_UR_une_onsets_mean, day_UR_une_onsets_mean], 'r');
