% WF_CRP_analysis_overview

clear; 
WF_CRP_list_of_days;
F_TC_dir    = 'Z:\Analysis\WF Lever Analysis\LeverSummaryFolder\';
lick_TC_dir = 'Z:\Analysis\WF Lever Analysis\BxAndAnalysisOutputs\BxOutputs\'; 
out_dir = 'Z:\Analysis\Cue_reward_pairing_analysis\WF\onset_latencies\';
BEHAVE_DIR = 'Z:\Data\WidefieldImaging\GCaMP\behavior\';
cue_fr = 6;
reward_fr = 12;
pre_fr = 5;
post_fr = 20;
x_axis = [-pre_fr:post_fr].*100;

% days_subset = days_post;
% ROIs_subset = days_post_ROI;
days_subset = days_1;
ROIs_subset = days_1_ROI;
NR_all_trials = [];
OR_all_trials = [];
figure;
for sess_num=find(~cellfun(@isempty,days_subset))
    
    %load TCs and and licking -------------------------------------
    b_data = get_bx_data(BEHAVE_DIR, days_subset{sess_num});  %find the correct behavior file and loads it.
    load([lick_TC_dir, days_subset{sess_num}, '_bx_outputs'], 'licking_data');
    lick_trace_rew = licking_data.lick_trace_rew;
    load([F_TC_dir, days_subset{sess_num}, '_reward_trials']);
    if exist([F_TC_dir, days_subset{sess_num}, '_rew_om_trials.mat'])
        %if strcmp(days_subset{1}, days_1{1}) | strcmp(days_subset{1}, days_post{1}) | strcmp(days_subset{1}, days_1000{1});
        if strcmp(days_subset{1}, days_1{1}) | strcmp(days_subset{1}, days_post{1}) | strcmp(days_subset{1}, days_1000{1});
            load([F_TC_dir, days_subset{sess_num}, '_rew_om_trials']);
            lick_trace_rew_om = licking_data.lick_trace_rew_om;
        end
    end
    if exist([F_TC_dir, days_subset{sess_num}, '_unexp_trials.mat']) && strcmp(days_subset{1}, days_UR{1});
        load([F_TC_dir, days_subset{sess_num}, '_unexp_trials']);
        lick_trace_unexp_rew = licking_data.lick_trace_unexp_rew;
    end
    
    %only include LS ROIs 
    rewarded_roi = rewarded_roi(:,[ROIs_subset{sess_num}],:);
    rew_om_roi = rew_om_roi(:,[ROIs_subset{sess_num}],:);
    plot(x_axis, squeeze(mean(mean(rewarded_roi,2),1)), 'k'); hold on;
    plot(x_axis, squeeze(mean(mean(rew_om_roi,2),1)), 'r');
    %get first half of trials
    rew_1st_half = round(size(rewarded_roi,1)/2);
    om_1st_half = round(size(rew_om_roi,1)/2);
    rewarded_roi = rewarded_roi([1:rew_1st_half],:,:);
    rew_om_roi = rew_om_roi([1:om_1st_half],:,:);
    %average together ROIs within a session  %store all trials
    NR_all_trials = [NR_all_trials; squeeze(mean(rewarded_roi,2))];
    if om_1st_half >2
        OR_all_trials = [OR_all_trials; squeeze(mean(rew_om_roi,2))];
    end
end

%% a grand mean plotting 

figure; 
subplot(1,2,1);
errorbar(x_axis, mean(NR_all_trials), std(NR_all_trials)/sqrt(size(NR_all_trials,1)), 'k'); hold on;
errorbar(x_axis, mean(OR_all_trials), std(OR_all_trials)/sqrt(size(OR_all_trials,1)), 'r');
vline(0, 'k');
vline(600, 'b');
ylabel('df/f');
xlabel('time (ms) relative to cue onset');
title(['NR: ', num2str(size(NR_all_trials,1)), ' trials,   OR: ', num2str(size(OR_all_trials,1)), 'trials,  ', num2str(sum(~cellfun(@isempty,days_subset))), ' sessions']);
hline(0,'k');
xlim([-500 2000]);

subplot(1,2,2);
plot(x_axis, mean(OR_all_trials)-mean(NR_all_trials), 'k');
ylim([-0.03 0.07]);
vline(0, 'k');
vline(600, 'b');
ylabel('df/f');
xlabel('time (ms) relative to cue onset');
title(['Omission - Normal rewarded trials']);
hline(0,'k');


