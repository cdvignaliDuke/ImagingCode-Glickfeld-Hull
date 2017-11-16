clear;
%trying to compare failed trials with licking to failed trials without
%licking. 

%sessions to be analyzed
%days = {'160904_img55', '160905_img55', '160916_img61', '160918_img61', '160920_img61', '160921_img61', '161030_img62', '160904_img55'};
%days = {'161031_img68','161101_img68', '161030_img69', '161030_img70', '161101_img69', '161101_img70'};

bx_source     = ['Z:\Data\WidefieldImaging\GCaMP\behavior\'];
image_source_base  = ['Z:\Data\WidefieldImaging\GCaMP\']; %location of permanently stored image files for retreiving meta data
image_dest_base    = ['Z:\Analysis\WF Lever Analysis\BxAndAnalysisOutputs\']; %stores the data on crash in the lever analysis folder
bx_outputs_dir = ['Z:\Analysis\WF Lever Analysis\BxAndAnalysisOutputs\BxOutputs\'];
TC_dir = ['Z:\Analysis\WF Lever Analysis\LeverSummaryFolder\'];
output_dir_early = ['Z:\Analysis\WF Lever Analysis\licking_investigation\early_trials_lick_v_no_lick\'];
output_dir_corr = ['Z:\Analysis\WF Lever Analysis\licking_investigation\correct_trials_lick_v_no_lick\'];
colors = [1,0,0; 0,1,0; 0,0,1; 0.5,0.5,0.5; 1,0,1; 1,1,0; 0,1,1]; %sets up the color scheme for plotting multiple ROIs with errorbar     
old_cd = cd; %save old cd so I can restore it later
WF_plotting_lists_of_days;
no_lick_window = [6:9];
lever_frame = 6;
percent_fails_lick_all = [];
percent_corr_lick_all = [];

%% SECTION TWO  EARLY TRIALS
for ii = 1:length(days);
    %load bx and failed trials TCs
    b_data = get_bx_data(bx_source, days{ii});  %find the correct behavior file and loads it.
    bx_out_dir  = [bx_outputs_dir days{ii} '_bx_outputs'];
    load(bx_out_dir);
    load([TC_dir, days{ii}, '_fail']);
    lick_trace_fail = licking_data.lick_trace_fail; 
    
    %segregate trials based on licking output
    lick_fail_trials = find(sum(lick_trace_fail(:, no_lick_window),2))';
    no_lick_fail_trials = [1:size(lick_trace_fail,1)];
    no_lick_fail_trials(lick_fail_trials) = [];
    save([output_dir_early, days{ii}, '_trial_indeces'], 'lick_fail_trials', 'no_lick_fail_trials', 'no_lick_window');
    
    %quantify percent of failed trials w/ licks in the window
    num_failed_trials = size(lick_trace_fail,1);
    if isempty(lick_trace_fail)
        percent_fails_lick = 0;
    else
        percent_fails_lick = length(lick_fail_trials)/num_failed_trials;
    end
    percent_fails_lick_all = [percent_fails_lick_all, percent_fails_lick];
    
    %skip plotting if nothing to plot
    if isempty(no_lick_fail_trials)
        disp([days{ii}, ' no trials WITHOUT lickiing in window']);
        continue;
    elseif isempty(lick_fail_trials)
        disp([days{ii}, ' no trials WITH lickiing in window']);
        continue;
    end
    
    %plot failed trials with vs without licks in window
    offset = 0;
    x_axis = [-5:10];
    figure; hold on;
   
    for ROI_num = 1:size(fail_roi,2)
        no_lick_TC = squeeze(mean(fail_roi(no_lick_fail_trials, ROI_num, :),1));
        lick_TC = squeeze(mean(fail_roi(lick_fail_trials, ROI_num, :),1));
        shift_no_lick= mean(no_lick_TC(2:4));
        shift_lick = mean(lick_TC(2:4));
        no_lick_TC = no_lick_TC - shift_no_lick;
        lick_TC = lick_TC- shift_lick;
        no_lick_sem = squeeze(std(fail_roi(no_lick_fail_trials, ROI_num, :),1))'/sqrt(length(no_lick_fail_trials));
        lick_sem = squeeze(std(fail_roi(lick_fail_trials, ROI_num, :),1))'/sqrt(length(lick_fail_trials));
        if length(no_lick_sem) == length(no_lick_TC)
            errorbar(x_axis, no_lick_TC+offset, no_lick_sem, 'Color','r');
        else
            plot(x_axis, no_lick_TC+offset, 'Color','r');
        end
        if length(lick_sem) == length(lick_TC)
             errorbar(x_axis, lick_TC+offset, lick_sem, 'Color','m');
        else
             plot(x_axis, lick_TC+offset, 'Color','m');
        end
        offset = offset +0.2;
    end
    axis tight;
    vline(0);
    title([['no licking in frames ', num2str(no_lick_window-lever_frame), days(ii)], ['failed trials with licks n=', num2str(length(lick_fail_trials)), ' (magenta).'] ['Without licks n=', num2str(length(no_lick_fail_trials)), '(red)']]);
    xlabel('frame number relative to lever release');
    ylabel('df/f segregated by ROI');
    savefig([output_dir_early, days{ii}]);
   
    %save the time courses for use in other scripts
    early_TCs.no_lick_TCs = fail_roi(no_lick_fail_trials, :, :);
    early_TCs.lick_TCs    = fail_roi(lick_fail_trials, :, :);
    save([output_dir_early, days{ii}], 'early_TCs');
    
end

%% Section THREE CORRECT TRIALS
for ii = 1:length(days);
    %load bx and correct trials TCs
    b_data = get_bx_data(bx_source, days{ii});  %find the correct behavior file and loads it.
    bx_out_dir  = [bx_outputs_dir days{ii} '_bx_outputs'];
    load(bx_out_dir);
    load([TC_dir, days{ii}, '_success']);
    lick_trace_corr = licking_data.lick_trace_succ;
    
    %segregate trials based on licking output
    lick_corr_trials = find(sum(lick_trace_corr(:, no_lick_window),2))';
    no_lick_corr_trials = [1:size(lick_trace_corr,1)];
    no_lick_corr_trials(lick_corr_trials) = [];
    save([output_dir_corr, days{ii}, '_trial_indeces'], 'lick_corr_trials', 'no_lick_corr_trials', 'no_lick_window');
    
    %quantify percent of correct trials w/ licks in the window
    num_corr_trials = size(lick_trace_corr,1);
    if isempty(lick_trace_corr)
        percent_corr_lick = 0;
    else
        percent_corr_lick = length(lick_corr_trials)/num_corr_trials;
    end
    percent_corr_lick_all = [percent_corr_lick_all, percent_corr_lick];
    
    %skip plotting if nothing to plot
    if isempty(no_lick_corr_trials)
        disp([days{ii}, ' no trials WITHOUT lickiing in window']);
        continue;
    elseif isempty(lick_corr_trials)
        disp([days{ii}, ' no trials WITH lickiing in window']);
        continue;
    end
    
    %plot correct trials with vs without licks in window
    offset = 0;
    x_axis = [-5:10];
    figure; hold on;
    
    for ROI_num = 1:size(success_roi,2)
        no_lick_TC = squeeze(mean(success_roi(no_lick_corr_trials, ROI_num, :),1));
        lick_TC = squeeze(mean(success_roi(lick_corr_trials, ROI_num, :),1));
        shift_no_lick= mean(no_lick_TC(2:4));
        shift_lick = mean(lick_TC(2:4));
        no_lick_TC = no_lick_TC - shift_no_lick;
        lick_TC = lick_TC- shift_lick;
        no_lick_sem = squeeze(std(success_roi(no_lick_corr_trials, ROI_num, :),1))'/sqrt(length(no_lick_corr_trials));
        lick_sem = squeeze(std(success_roi(lick_corr_trials, ROI_num, :),1))'/sqrt(length(lick_corr_trials));
        if length(no_lick_sem) == length(no_lick_TC)
            errorbar(x_axis, no_lick_TC+offset, no_lick_sem, 'Color','r');
        else
            plot(x_axis, no_lick_TC+offset, 'Color','r');
        end
        if length(lick_sem) == length(lick_TC)
             errorbar(x_axis, lick_TC+offset, lick_sem, 'Color','m');
        else
             plot(x_axis, lick_TC+offset, 'Color','m');
        end
        offset = offset +0.2;
    end
    axis tight;
    vline(0);
    title([['no licking in frames ', num2str(no_lick_window-lever_frame), days(ii)], ['correct trials with licks n=', num2str(length(lick_corr_trials)), ' (magenta)'], ['Without licks n=', num2str(length(no_lick_corr_trials)), '(red)']]);
    xlabel('frame number relative to lever release');
    ylabel('df/f segregated by ROI');
    savefig([output_dir_corr, days{ii}]);
    
    %save the time courses for use in other scripts
    corr_TCs.no_lick_TCs = success_roi(no_lick_corr_trials, :, :);
    corr_TCs.lick_TCs    = success_roi(lick_corr_trials, :, :);
    save([output_dir_corr, days{ii}], 'corr_TCs');
end

%% section four - plot percent fail vs corr trials with licks
figure; scatter(percent_corr_lick_all, percent_fails_lick_all, 'o'); 
hold on; 
[corr_lick_mean, corr_lick_sem] = get_mean_and_sem(percent_corr_lick_all); 
[fail_lick_mean, fail_lick_sem] = get_mean_and_sem(percent_fails_lick_all); 
errorbarxy(corr_lick_mean, fail_lick_mean, corr_lick_sem, fail_lick_sem);
xlabel('correct trials: percent with lick in window');
ylabel('early trials: percent with lick in window');
title(['WF data % of trials with lick(s) in window ' num2str(no_lick_window(1)-lever_frame), ' to ', num2str(no_lick_window(end)-lever_frame)]);
xlim([0 1]);
ylim([0 1]); hold on; 
plot([0:.1:1], [0:.1:1]);
%savefig(['Z:\Analysis\WF Lever Analysis\licking_investigation\percent_with_licks_scatter']);









