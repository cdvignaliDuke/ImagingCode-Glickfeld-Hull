%find the grand mean of the df/f and licking TCs for the widefield animals
%on day1 dayN dayN+1 and 1000ms

%alter the script to get day1 vs DayN TCs. Exclude second half of trials
%for day1 to get naive mice. 

WF_CRP_list_of_days;
TC_dir = 'Z:\Analysis\WF Lever Analysis\LeverSummaryFolder\';
licking_dir = 'Z:\Analysis\WF Lever Analysis\BxAndAnalysisOutputs\BxOutputs\';
days = days_1;
ROI_cell = days_1_ROI;
reward_frame_num = 12;
no_lick_window = [9:11];

rewarded_roi_all = [];
rew_om_roi_all = [];
unexpected_roi_all = [];
lick_trace_rew_all = [];
lick_trace_rew_om_all = [];

%figure; subplot(2,2,3); hold on;  
colors = {'k', 'r', 'b', 'm', 'g', 'c', 'y', 'k', 'r'}; %title('all sessions df/f TCs omission');
for ii = 5 % 1:length(days)
    if isempty(days{ii});
        continue;
    end
    load([TC_dir, days{ii}, '_rew_om_trials']);
    load([TC_dir, days{ii}, '_reward_trials']);
    load([licking_dir, days{ii}, '_bx_outputs'], 'licking_data');
    lick_trace_rew = licking_data.lick_trace_rew;
    lick_trace_rew_om = licking_data.lick_trace_rew_om;
    
    %only use first half of trials on day 1
    if strcmp(days{1}, days_1{1})
        rewarded_roi = rewarded_roi([1:round(size(rewarded_roi,1)/2)], :, :);
        rew_om_roi = rew_om_roi([1:round(size(rew_om_roi,1)/2)], :, :);
        lick_trace_rew = lick_trace_rew([1:round(size(rewarded_roi,1)/2)], :);
        lick_trace_rew_om = lick_trace_rew_om([1:round(size(rew_om_roi,1)/2)],:);
    end
    
    %identify trials where lick bout starts before reward delivery 
    rew_no_lick_inx = remove_trials_with_licks(lick_trace_rew, no_lick_window, days{ii});
    rew_om_no_lick_inx = remove_trials_with_licks(lick_trace_rew_om, no_lick_window, days{ii});
    disp([days{ii}, ' rewarded n=', num2str(length(rew_no_lick_inx)), '   omission n=', num2str(length(rew_om_no_lick_inx))]);
    
    %remove trials with licking between cue and reward
%     lick_trace_rew = lick_trace_rew([rew_no_lick_inx], :);
%     lick_trace_rew_om = lick_trace_rew_om([rew_om_no_lick_inx], :);
%     rewarded_roi = rewarded_roi([rew_no_lick_inx], :, :);
%     rew_om_roi = rew_om_roi([rew_om_no_lick_inx], :, :);
    
    %store avg of remaining trials
    lick_trace_rew = mean(lick_trace_rew,1);
    lick_trace_rew_om = mean(lick_trace_rew_om,1);
    lick_trace_rew_all = [lick_trace_rew_all; lick_trace_rew];
    lick_trace_rew_om_all = [lick_trace_rew_om_all; lick_trace_rew_om];
    
    %====================================================================
    %plot the lick trace and df/f traces for each animal
    x_axis = [1:size(rewarded_roi,3)]*100-600;
    temp_n = size(rewarded_roi,1);
    rewarded_roi_error = std(squeeze(rewarded_roi(:,3,:)),[],1)/sqrt(size(rewarded_roi,1));
    rewarded_roi = squeeze(mean(rewarded_roi,1));
    figure; subplot(2,1,1);
    bar(x_axis, lick_trace_rew/10); hold on;
    for iii = 3% 1:size(rewarded_roi,1)
        plot(x_axis, rewarded_roi(iii, :)', 'Color', colors{iii});
        errorbar(x_axis, rewarded_roi(iii, :)', rewarded_roi_error);
    end
    title([days{ii}, ' day 1 rewarded trials, no lick window, n=', num2str(temp_n)]); xlabel('ms relative to cue onset');
    %--------- reward omission
    subplot(2,1,2); hold on;
    temp_n = size(rew_om_roi,1);
    rew_om_roi_error = std(squeeze(rew_om_roi(:,3,:)),[],1)/sqrt(size(rew_om_roi,1));
    rew_om_roi = squeeze(mean(rew_om_roi,1));
    bar(x_axis, lick_trace_rew_om/10);
    for iii = 3% 1:size(rew_om_roi,1)
        plot(x_axis, rew_om_roi(iii, :)', 'Color', colors{iii});
        errorbar(x_axis, rew_om_roi(iii, :)', rew_om_roi_error);
    end
    title([days{ii}, ' day 1 omission trials, no lick window, n=', num2str(temp_n)]); xlabel('ms relative to cue onset');
    %continue
    %=================================================================
    
    %avg together the ROIs (dim2) and trials (dim1)
    rewarded_roi = squeeze(mean(rewarded_roi(:,ROI_cell{ii},:),2));
    %errorbar( [1:size(rewarded_roi,2)]*100-600,  mean(rewarded_roi),  std(rewarded_roi)/sqrt(length(rewarded_roi)), colors{ii});
    rewarded_roi = mean(rewarded_roi);
    rewarded_roi_all = [rewarded_roi_all; rewarded_roi];
    
    %do the same for omission and unexpected trials if the exist
    if exist('rew_om_roi') && size(rew_om_roi,1) >2
        rew_om_roi = squeeze(mean(rew_om_roi(:,ROI_cell{ii},:),2));
        %errorbar( [1:size(rew_om_roi,2)]*100-600,  mean(rew_om_roi),  std(rew_om_roi)/sqrt(length(rew_om_roi)), colors{ii});
        rew_om_roi= mean(rew_om_roi);
        rew_om_roi_all = [rew_om_roi_all; rew_om_roi];
    end
    if exist('unexpected_roi');
        unexpected_roi = squeeze(mean(unexpected_roi(:,ROI_cell{ii},:),1));
        unexpected_roi = mean(unexpected_roi);
        unexpected_roi_all = [unexpected_roi_all; unexpected_roi];
    end
end

[lick_trace_rew_all_mean, lick_trace_rew_all_sem] = get_mean_and_sem(lick_trace_rew_all);
[lick_trace_rew_om_all_mean, lick_trace_rew_om_all_sem] = get_mean_and_sem(lick_trace_rew_om_all);
[rewarded_roi_all_mean, rewarded_roi_all_sem] = get_mean_and_sem(rewarded_roi_all);
if ~isempty(rew_om_roi_all)
    [rew_om_roi_all_mean, rew_om_roi_all_sem] = get_mean_and_sem(rew_om_roi_all);
end
if ~isempty(unexpected_roi_all)
    [unexpected_roi_all_mean, unexpected_roi_all_sem] = get_mean_and_sem(unexpected_roi_all);
end

%plot the timecourses
x_axis = [1:length(rewarded_roi_all_mean)]*100-600;
figure; 
subplot(2,2,1); hold on;
bar(x_axis, lick_trace_rew_all_mean/10);
errorbar(x_axis, lick_trace_rew_all_mean/10, lick_trace_rew_all_sem/10, '.');
errorbar(x_axis, rewarded_roi_all_mean, rewarded_roi_all_sem, 'g'); 
title('day 1 rewarded trials grand mean');
xlabel('ms relative to cue onset');
ylabel('df/f  and avg lick rate Hz/100');
xlim([min(x_axis)-50 max(x_axis)+50]); ylim([-0.07 0.14]);
vline(600);

% figure; 
subplot(2,2,2); hold on;
bar(x_axis, lick_trace_rew_om_all_mean/10);
errorbar(x_axis, lick_trace_rew_om_all_mean/10, lick_trace_rew_om_all_sem/10, '.');
errorbar(x_axis, rew_om_roi_all_mean, rew_om_roi_all_sem, 'g'); 
title('day 1 omission trials grand mean');
xlabel('ms relative to cue onset');
ylabel('df/f  and avg lick rate Hz/100');
ylim([-0.08 0.1]);
xlim([min(x_axis)-50 max(x_axis)+50]); ylim([-0.07 0.14]);
vline(600);

