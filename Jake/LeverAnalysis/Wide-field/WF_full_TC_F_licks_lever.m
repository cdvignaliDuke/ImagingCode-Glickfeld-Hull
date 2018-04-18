%WF_full_TC_F_licks_lever
%script which collects the necessary data to plot a single time course for
%each session encompassing the entire session. X axis will be in frames. Y
%plots will include df/f for each ROI and licks per frame. Black vertical
%lines will indicate lever presses. Green vertical lines will indicate rewarded lever
%releases. Red vertical lines will represent unrewarded lever releases. 
clear;

pre_frames = 9;
post_frames = 50;

days = {'161107_img68', '161107_img69', '161030_img68', '161030_img70', '161030_img69', '161031_img68', '161101_img69', '161101_img70'};  %cue-reward pairing delayed reward. 

bx_source     = ['Z:\Data\WidefieldImaging\GCaMP\behavior\'];
image_source_base  = ['Z:\Data\WidefieldImaging\GCaMP\']; %location of permanently stored image files for retreiving meta data
image_dest_base    = ['Z:\Analysis\WF Lever Analysis\BxAndAnalysisOutputs\']; %stores the data on crash in the lever analysis folder
bx_outputs_dir = ['Z:\Analysis\WF Lever Analysis\BxAndAnalysisOutputs\BxOutputs\'];
colors = [1,0,0; 0,1,0; 0,0,1; 0.5,0.5,0.5; 1,0,1; 1,1,0; 0,1,1]; %sets up the color scheme for plotting multiple ROIs with errorbar     
old_cd = cd; %save old cd so I can restore it later
WF_plotting_lists_of_days;

%datasets with good licking traces
days = {[], [], [], [], '151211_img32', '151212_img32', '160129_img35', '160131_img35', '160129_img36', '160131_img36', '160314_img38', '160315_img38', '160319_img41', '160320_img41', '160606_img46', '160722_img53', '160904_img55'}; %'150718_img27', '150719_img27', '150716_img28', '150717_img28', 

summary_bout_TCs = [];
summary_bout_TCs_early = [];
%% SECTION TWO  
for session_num = 5% 1:length(days);
    if isempty(days{session_num})
        lick_triggered_ttest{session_num}= [];
        continue
    end
    b_data = get_bx_data(bx_source, days{session_num});  %find the correct behavior file and loads it.
    bx_out_dir  = [bx_outputs_dir days{session_num} '_bx_outputs'];
    load(bx_out_dir);
    
    imaging_start_MW_T = frame_info.imaging_start_MW_T;
    licksByFrame = licking_data.licksByFrame; 
    leverPress = cell2mat(b_data.holdStartsMs)-imaging_start_MW_T;
    holdDur = cell2mat(b_data.holdTimesMs);
    leverRelease = leverPress + holdDur;
    reqHoldTime = cell2mat(b_data.tTotalReqHoldTimeMs);
    changeOrientation = trial_outcome.change_orientation;
    
    %find and trim the trials which were not imaged. 
    trials_to_trim = [(1:frame_info.f_frame_trial_num), (frame_info.l_frame_trial_num:length(leverPress))];
    
    if ~isempty(trials_to_trim);
        leverPress(trials_to_trim) = [];
        holdDur(trials_to_trim) = [];
        leverRelease(trials_to_trim) = [];
        reqHoldTime(trials_to_trim) = [];
        changeOrientation(trials_to_trim) = [];
    end
    leverPressFrame = [];
    for iii = 1:length(leverPress)
        leverPressFrame = [leverPressFrame, frame_info.counter(leverPress(iii))];
    end
    leverReleaseFrame = [];
    for iii = 1:length(leverRelease)
        leverReleaseFrame = [leverReleaseFrame, frame_info.counter(leverRelease(iii))];
    end
    corrInx = zeros(1, length(holdDur));
    for iii = 1:length(holdDur)
        if holdDur(iii) > reqHoldTime(iii)
            corrInx(iii)=1;
        end
    end
    earlyInx = find(corrInx == 0);
    corrInx  = find(corrInx == 1);
    
    %changeOrientation = changeOrientation(find(corrInx)); %trim changeOr so it only includes corrects
    changeOrientation = round(changeOrientation(find(~isnan(changeOrientation)))); %get rid of the NaNs
    cueByFrame = [];
    for iii = 1:length(changeOrientation) 
        cueByFrame = [cueByFrame, frame_info.counter(changeOrientation(iii))];
    end
    
    no_fidget_earlyInx = zeros(1, length(holdDur));
    for iii = 1:length(holdDur)
        if holdDur(iii) < reqHoldTime(iii) & holdDur(iii) > 200
            no_fidget_earlyInx(iii)=1;
        end
    end
    no_fidget_earlyInx = find(no_fidget_earlyInx == 1);
    
    no_TFC_corrInx = zeros(1, length(holdDur));
    for iii = 1:length(holdDur)
        if holdDur(iii) > reqHoldTime(iii) & holdDur(iii) > reqHoldTime(iii)+150
            no_TFC_corrInx(iii)=1;
        end
    end
    no_TFC_corrInx = find(no_TFC_corrInx == 1);
    
    
    figure; 
    bar(licksByFrame*0.01);
    alpha(.25);
    hold on; 
    for iii = 1%:size(tc_dfoverf,1)
        plot([1:length(licksByFrame)], tc_dfoverf(iii,[1:length(licksByFrame)]), 'Color', colors(iii,:))
    end
    for iii = 1:length(leverPressFrame)
        vline(leverPressFrame(iii),'k');
    end
    for iii = 1:length(earlyInx)
        vline(leverReleaseFrame(earlyInx(iii)),'r');
    end
    for iii = 1:length(corrInx)
        vline(leverReleaseFrame(corrInx(iii)),'g');
    end
    
    %plot individual trials
    frame_int = [1003:1021];    %correct examples:   [459:483] [5097:5121]  [9508:9545]    early examples: [8470:8509]   [7787:7819]  [5011:5030]
    lever_press_plot = intersect(frame_int, leverPressFrame)-frame_int(1)+1;
    corr_release_plot = intersect(frame_int, leverReleaseFrame(corrInx))-frame_int(1)+1;
    early_release_plot = intersect(frame_int, leverReleaseFrame(earlyInx))-frame_int(1)+1;
    cue_plot = intersect(frame_int, cueByFrame)-frame_int(1)+1;
    figure; hold on;
    iii=1;
    bar(licksByFrame([frame_int])*0.01);
    plot([1:length(frame_int)], tc_dfoverf(iii,[frame_int]), 'r');
    ylim([-0.05 0.15])
    vline(lever_press_plot,'k');
    if ~isempty(early_release_plot)
        vline(early_release_plot,'r');
    end
    if ~isempty(corr_release_plot)
        vline(corr_release_plot,'g');
        vline(cue_plot,'m');
    end
    title(['early trial examlpe: frame range ', num2str(frame_int(1)), ' ', num2str(frame_int(end))]);
    xlim([1 40]);
    
    hold off
    
    % make figure for example traces
    %early lever release example
    early_ex_press_inx = [2026, 7793, 8475];
    early_ex_hold_inx = [8 26, 13];
    figure; 
    plot()


%%  corrInx is the trial # of correct trials not the frame number
    %plot y frames before and x frames after the lever release   Write a
    %different one for the cue-reward pairing sessions
%     corr_release_by_frame = leverReleaseFrame(no_TFC_corrInx);
%     corr_release_by_frame =  corr_release_by_frame(corr_release_by_frame > pre_frames & corr_release_by_frame < length(tc_dfoverf)-post_frames); %trimming to make sure we have enough frames before and after
%     corr_TCs = NaN(pre_frames+post_frames+1,length(corr_release_by_frame)-1,size(tc_dfoverf,1));
%     corr_licks = NaN(pre_frames+post_frames+1,length(corr_release_by_frame)-1);
%     for  iii = 1:length(corr_release_by_frame)-1;  %corr_release_by_frame are the frames on which correct lever releases occurred 
%         for iv = 1:size(tc_dfoverf,1)
%             corr_TCs(:,iii,iv) = tc_dfoverf(iv, [corr_release_by_frame(iii)-pre_frames]:[corr_release_by_frame(iii)+post_frames]);
%             corr_licks(:,iii) = licksByFrame([corr_release_by_frame(iii)-pre_frames]:[corr_release_by_frame(iii)+post_frames]);
%         end
%     end
%     
%     assert(isempty(find(isnan(corr_TCs))));
%     assert(isempty(find(isnan(corr_licks))));
%     corr_TCs_mean = squeeze(mean(corr_TCs,2)); %dim 1 will be frames and dim 2 will be ROI
%     corr_TCs_sem = squeeze(std(corr_TCs,[],2)/sqrt(size(corr_TCs,2)));
%     corr_licks_mean = mean(corr_licks,2);
%     
%     figure; hold on;
%     bar([1:size(corr_TCs_mean,1)]-(pre_frames+1)-10, corr_licks_mean/10);
%     alpha(0.25);
%     for iii = 1:size(corr_TCs_mean,2)
%         plot([1:size(corr_TCs_mean,1)]-(pre_frames+1)-10, corr_TCs_mean(:,iii), 'Color', colors(iii,:)); hold on;
%         errorbar([1:size(corr_TCs_mean,1)]-(pre_frames+1)-10, corr_TCs_mean(:,iii), corr_TCs_sem(:,iii), 'Color', colors(iii,:));
%     end
%     title([days{ii}(1:6), ' ', days{ii}(8:end), 'reward aligned']);
%     xlabel('frame # relative to reward delivery');
%     ylabel('df/f');
%     vline(-14, 'c');
%     vline(-10, 'r');
%     vline(0, 'b');
%     savefig(['\\crash\data\public\presentations\SfN\licking_investigation\cue_reward_pairing_5000ms_TCs\', days{ii}, '_reward_aligned']);
%     summary_bout_TCs(ii,:) = corr_licks_mean;
    
    %% plot TCs for incorrect trials
%     fail_release_by_frame = leverReleaseFrame(no_fidget_earlyInx);
%     fail_release_by_frame =  fail_release_by_frame(fail_release_by_frame > pre_frames & fail_release_by_frame < length(tc_dfoverf)-post_frames); %trimming to make sure we have enough frames before and after
%     early_TCs = NaN(pre_frames+post_frames+1,length(fail_release_by_frame)-1,size(tc_dfoverf,1));
%     early_licks = NaN(pre_frames+post_frames+1,length(fail_release_by_frame)-1);
%     for  iii = 1:length(fail_release_by_frame)-1;  %corr_release_by_frame are the frames on which correct lever releases occurred 
%         for iv = 1:size(tc_dfoverf,1)
%             early_TCs(:,iii,iv) = tc_dfoverf(iv, [fail_release_by_frame(iii)-pre_frames]:[fail_release_by_frame(iii)+post_frames]);
%             early_licks(:,iii) = licksByFrame([fail_release_by_frame(iii)-pre_frames]:[fail_release_by_frame(iii)+post_frames]);
%         end
%     end
%     
%     assert(isempty(find(isnan(early_TCs))));
%     assert(isempty(find(isnan(early_licks))));
%     early_TCs_mean = squeeze(mean(early_TCs,2)); %dim 1 will be frames and dim 2 will be ROI
%     early_licks_mean = mean(early_licks,2);
%     summary_bout_TCs_early(ii,:) = early_licks_mean;
    
    %% plot TCs for cue aligned corrects
%     cueByFrame =  cueByFrame(cueByFrame > pre_frames & cueByFrame < length(tc_dfoverf)-post_frames); %trimming to make sure we have enough frames before and after
%     cue_TCs = NaN(pre_frames+post_frames+1,length(cueByFrame)-1,size(tc_dfoverf,1));
%     cue_licks = NaN(pre_frames+post_frames+1,length(cueByFrame)-1);
%     for  iii = 1:length(cueByFrame)-1;  %corr_release_by_frame are the frames on which correct lever releases occurred 
%         for iv = 1:size(tc_dfoverf,1)
%             cue_TCs(:,iii,iv) = tc_dfoverf(iv, [cueByFrame(iii)-pre_frames]:[cueByFrame(iii)+post_frames]);
%             cue_licks(:,iii) = licksByFrame([cueByFrame(iii)-pre_frames]:[cueByFrame(iii)+post_frames]);
%         end
%     end
%     
%     assert(isempty(find(isnan(cue_TCs))));
%     assert(isempty(find(isnan(cue_licks))));
%     cue_TCs_mean = squeeze(mean(cue_TCs,2)); %dim 1 will be frames and dim 2 will be ROI
%     cue_TCs_sem = squeeze(std(cue_TCs,[],2)/sqrt(size(cue_TCs,2)));
%     cue_licks_mean = mean(cue_licks,2);
%     
%     figure; hold on;
%     bar([1:size(cue_TCs_mean,1)]-(pre_frames+1)+4, cue_licks_mean/10);
%     alpha(0.25);
%     for iii = 1:size(cue_TCs_mean,2)
%         plot([1:size(cue_TCs_mean,1)]-(pre_frames+1)+4, cue_TCs_mean(:,iii), 'Color', colors(iii,:)); hold on;
%         errorbar([1:size(cue_TCs_mean,1)]-(pre_frames+1)+4, cue_TCs_mean(:,iii), cue_TCs_sem(:,iii), 'Color', colors(iii,:));
%     end
%     title([days{ii}(1:6), ' ', days{ii}(8:end), ' cue aligned']);
%     xlabel('frame # relative to cue');
%     ylabel('df/f');
%     vline(0, 'c');
%     vline(4, 'r');
%     vline(14, 'b');
%     savefig(['\\crash\data\public\presentations\SfN\licking_investigation\cue_reward_pairing_5000ms_TCs\', days{ii},'_cue_aligned']);
%     %summary_bout_TCs(ii,:) = cue_licks_mean;
%     
end


% figure; 
% bar([1:size(corr_TCs_mean,1)]-(pre_frames+1), mean(summary_bout_TCs)); hold on;
% alpha(0.25);
% bar([1:size(early_TCs_mean,1)]-(pre_frames+1), mean(summary_bout_TCs_early),'r');
% alpha(0.25);
% errorbar([1:size(corr_TCs_mean,1)]-(pre_frames+1), mean(summary_bout_TCs), std(summary_bout_TCs)/sqrt(ii));
% errorbar([1:size(early_TCs_mean,1)]-(pre_frames+1), mean(summary_bout_TCs_early), std(summary_bout_TCs_early)/sqrt(ii),'r');
% title(['average licking traces across animals n=', num2str(ii), ' sessions']);
% ylabel('average number of licks per frame');
% xlabel('frames number relative to lever release');



