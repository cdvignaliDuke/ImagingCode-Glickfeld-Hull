%script for measuring the df/f response just before vs after the cue for
%naive CRP animals. Used in response to NN manuscript reviews. 
clear; 
WF_CRP_list_of_days;
F_TC_dir    = 'Z:\Analysis\WF Lever Analysis\LeverSummaryFolder\';
lick_TC_dir = 'Z:\Analysis\WF Lever Analysis\BxAndAnalysisOutputs\BxOutputs\'; 
out_dir = 'Z:\Analysis\Cue_reward_pairing_analysis\WF\onset_latencies\';
days = days_1([1,4,5,6,7]);
pre_cue_window = [2:5];
post_cue_window = [6:9];

for session_num = 1:length(days)
    %load TCs and and licking -------------------------------------
    load([lick_TC_dir, days{session_num}, '_bx_outputs'], 'licking_data');
    lick_trace_rew = licking_data.lick_trace_rew;
    load([F_TC_dir, days{session_num}, '_reward_trials']);
    if exist([F_TC_dir, days{session_num}, '_rew_om_trials.mat'])
        if strcmp(days{1}, days_1{1}) | strcmp(days{1}, days_post{1}) | strcmp(days{1}, days_1000{1});
            load([F_TC_dir, days{session_num}, '_rew_om_trials']);
            lick_trace_rew_om = licking_data.lick_trace_rew_om;
        end
    end
    
    %average together the ROIs ------------------------
    rew_om_roi = squeeze(mean(rew_om_roi,2));
    
    %find the mean df/f amplitude from the 400ms before vs after the cue
    pre_cue_trial_mean = mean(rew_om_roi(:,pre_cue_window),2);
    post_cue_trial_mean = mean(rew_om_roi(:,post_cue_window),2);
    [pre_cue_session_mean(session_num), pre_cue_session_sem(session_num)] = get_mean_and_sem(pre_cue_trial_mean);
    [post_cue_session_mean(session_num), post_cue_session_sem(session_num)] = get_mean_and_sem(post_cue_trial_mean);
end

figure; errorbar([1:5], post_cue_session_mean, post_cue_session_sem);
ylim([-0.1 0.1]);


%% plot an example TC from img95
days = days_1(4);
ROI_num =4;
x_axis = [-500:100:2000];

session_num = 1
    %load TCs and and licking -------------------------------------
    lick_trace_rew = licking_data.lick_trace_rew;
    load([F_TC_dir, days{session_num}, '_rew_om_trials']);
    
    
    %make a session average TC from an ROI ------------------------
    %rew_om_roi = rew_om_roi([1:ceil(size(rew_om_roi,1)/2)],:,:);
    rew_om_roi = squeeze(rew_om_roi(:,ROI_num,:));
    rew_om_roi_mean = mean(rew_om_roi,1);
    rew_om_roi_sem = std(rew_om_roi)/sqrt(size(rew_om_roi,1));
    figure; 
    errorbar(x_axis, rew_om_roi_mean, rew_om_roi_sem);
    ylim([-0.1 0.1]);
    hline(0); hline(0.01,'k');


