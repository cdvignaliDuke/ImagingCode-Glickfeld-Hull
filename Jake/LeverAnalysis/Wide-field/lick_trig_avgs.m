% script for isolating lick triggered averages
clear;

pre_frames = 8;
post_frames = 4;
lever_frame = pre_frames+1;
peak_window = [pre_frames+1:pre_frames+2];

bx_source      = ['Z:\Data\WidefieldImaging\GCaMP\behavior\'];
bx_outputs_dir = ['Z:\Analysis\WF Lever Analysis\BxAndAnalysisOutputs\BxOutputs\'];  
old_cd = cd; %save old cd so I can restore it later
WF_plotting_lists_of_days;

%datasets with good licking traces
days = {'151021_img29', '151022_img29', [], [], '151211_img32', '151212_img32', '160129_img35', '160131_img35', '160129_img36', '160131_img36', '160314_img38', '160315_img38', '160319_img41', '160320_img41', '160606_img46', '160722_img53', '160904_img55'}; %'150718_img27', '150719_img27', '150716_img28', '150717_img28', 
colors_roi = {'b', 'r', 'g', 'k', 'c', 'm'};

lick_triggered_ttest = {};
lick_trig_mags = {}; 
all_ttest = {};

for session_num = 1:length(days);
    if isempty(days{session_num})
        lick_triggered_ttest{session_num}= [];
        lick_trig_mags{session_num} = [];
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
    
    %collect all the lever events and index correct trials
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
    
    %find licking bouts not near lever events
    early_release_frame = leverReleaseFrame(earlyInx);
    corr_release_frame =  leverReleaseFrame(corrInx);
    bout_start_inx = [];
    bout_middle = [];
    for iii = 1:length(early_release_frame)-1
        this_lick_int_inx = [(early_release_frame(iii)+8):(leverPressFrame(earlyInx(iii)+1))-post_frames]; %look at the interval of licking from a given lever release to the subsequent lever press
        for iv = 1:length(this_lick_int_inx) %look at each frame in the selected interval
            if sum(licksByFrame((this_lick_int_inx(iv)-3):(this_lick_int_inx(iv)-1))) == 0     %if the three frames before this frame have no licks.
                if licksByFrame(this_lick_int_inx(iv))>0 & sum(licksByFrame([this_lick_int_inx(iv)+1:this_lick_int_inx(iv)+3]))>1  %And this frame has atleast one lick and there is at least 2 lick in the next 3 frames
                    bout_start_inx = [bout_start_inx, this_lick_int_inx(iv)];
                end
            end
        end
    end
    for iii = 1:length(corr_release_frame)-1
        this_lick_int_inx = [(corr_release_frame(iii)+12):(leverPressFrame(corrInx(iii)+1))-post_frames]; %look at the interval of licking from a given lever release to the subsequent lever press
        for iv = 1:length(this_lick_int_inx) %look at each frame in the selected interval
            if sum(licksByFrame((this_lick_int_inx(iv)-3):(this_lick_int_inx(iv)-1))) == 0     %if the three frames before this frame have no licks.
                if licksByFrame(this_lick_int_inx(iv))>0 & sum(licksByFrame([this_lick_int_inx(iv)+1:this_lick_int_inx(iv)+3]))>1 %And this frame has atleast one lick and there is at least 2 lick in the next 3 frames
                    bout_start_inx = [bout_start_inx, this_lick_int_inx(iv)];
                end
            end
        end
    end
    
    %if there are valid licks 
    if ~isempty(bout_start_inx)
        %extract traces from each bout
        lick_start_trig = [];
        for bout_num = 1:length(bout_start_inx)
            for iv = 1:length(LS_ROIs{session_num}) % 1:size(tc_dfoverf,1)
                lick_start_trig(bout_num,iv,:) = tc_dfoverf(LS_ROIs{session_num}(iv), [(bout_start_inx(bout_num)-pre_frames):(bout_start_inx(bout_num)+post_frames)]);   %dim1 bout#  dim2 ROI #   dim3 f values
            end
        end
        
        %shift each TC so 0 is the mean of...
        %shift = mean(lick_start_trig(:,:,[5:6]),3);
        shift = mean(lick_start_trig(:,:,[1:6]),3);
        shift = repmat(shift, 1, 1, post_frames+pre_frames+1);
        lick_start_trig = lick_start_trig-shift;
        
        %optional - plot individual trials
        figure; 
        for ROI_num = 1:size(lick_start_trig,2)
            subplot(1, size(lick_start_trig,2), ROI_num);
            this_baseline = squeeze(lick_start_trig(:,ROI_num,2));
            ttest_results = [];
            for bout_num = 1:size(lick_start_trig,1)
                plot([-pre_frames:post_frames], squeeze(lick_start_trig(bout_num,ROI_num,:))); hold on;
                ttest_results = [ttest_results, ttest( squeeze(mean(lick_start_trig(:,ROI_num,[1:5]),3)), squeeze(mean(lick_start_trig(bout_num,ROI_num,[lever_frame-1:lever_frame+1]))), 'Tail', 'left'   )];
            end
            title(['ROI#', num2str(ROI_num), '. Fraction of trials with a significant increase ', num2str(sum(ttest_results)/length(ttest_results))]);
            ylabel('df/f');
            xlabel('frames relative to lick');
            all_ttest{session_num}(ROI_num) = sum(ttest_results)/length(ttest_results);
        end
        suptitle([days{session_num}, ' 3 lick bout traces. n=', num2str(size(lick_start_trig,1))]);
        savefig(['Z:\Analysis\WF Lever Analysis\licking_investigation\lick_triggered_averages\individual lick bout TCs\', days{session_num}, '_bout']);
        continue
        
        %get mean and sem across trials for each ROI
        lick_start_sem = squeeze(std(lick_start_trig,1)/sqrt(size(lick_start_trig,1)));
        lick_start_trig_avg = squeeze(mean(lick_start_trig));
        if length(LS_ROIs{session_num}) ==1
            lick_start_trig_avg = lick_start_trig_avg';
            lick_start_sem = lick_start_sem';
        end
        
        %find and store the peak values
        if size(lick_start_trig,1) >= 4
            lick_trig_mags{session_num} = max(lick_start_trig_avg(:,peak_window),[],2)';
        else
            lick_trig_mags{session_num} = [];
        end
        
        %plot the mean TCs for each ROI
        figure; hold on;
        for iii = 1:length(LS_ROIs{session_num})
            errorbar([-pre_frames:post_frames], lick_start_trig_avg(iii,:), lick_start_sem(iii,:));
            %plot([-4:4], tc_dfoverf([(bout_start_inx(iii)-4):(bout_start_inx(iii)+4)]))
            title([days(session_num), ['df/f response to licking bout start n=' num2str(length(bout_start_inx))], 'minimum 3 lick per bout', 'no licks in 3 frames before zero']);
            xlabel('frame number relative to licking bout start');
            ylabel('df/f');
        end
        savefig(['Z:\Analysis\WF Lever Analysis\licking_investigation\lick_triggered_averages\', days{session_num},'_bout']);
        
        %ttest for significant licking response
        session_ttest = NaN(2,length(LS_ROIs{session_num}));
        for ROI_num = 1:length(LS_ROIs{session_num})
            [session_ttest(1,ROI_num), session_ttest(2,ROI_num)] = ttest(squeeze(mean(lick_start_trig(:,ROI_num,peak_window),3)) ,  squeeze(mean(lick_start_trig(:,ROI_num,[1:2]),3)) ,  'tail', 'right' );
        end
        lick_triggered_ttest{session_num} = session_ttest;
    end
    
end
all_ttest_mean = [];
for session_num = [1:length(days)]
    all_ttest_mean = [all_ttest_mean, mean(all_ttest{session_num})];
end
figure;
plot(ones(1, length(find(~isnan(all_ttest_mean)))), all_ttest_mean(find(~isnan(all_ttest_mean))), 'o');
ylim([0 1]); ylabel('percent of lick bouts with a significant positive df/f response');
title('ttest on individual ITI lick bouts averaged across trials and ROIs');
%savefig(['Z:\Analysis\WF Lever Analysis\licking_investigation\lick_triggered_averages\individual lick bout TCs\ttest_results_fig']);
%save(['Z:\Analysis\WF Lever Analysis\licking_investigation\lick_triggered_averages\individual lick bout TCs\all_ttest'], 'all_ttest', 'all_ttest_mean');
%save(['Z:\Analysis\WF Lever Analysis\licking_investigation\lick_triggered_averages\lick_trig_mags'], 'lick_trig_mags');
%save(['Z:\Analysis\WF Lever Analysis\licking_investigation\lick_triggered_averages\lick_triggered_ttest_bout'], 'lick_triggered_ttest');


