% %align licking and F to lever press. 
% use_ev_press = round(lever.press);
% [press_roi, use_times_press, lick_trace_press] = trigger_movie_by_event_licks(tc_dfoverf, frame_info, licksByFrame, use_ev_press, pre_frames, post_frames+35);
% ts_press = (-500:100:4500);
% 
% figure; 
% avg_press_roi = squeeze(func(press_roi,1));
% std_success = squeeze(std(press_roi,1));
% sm_press = std_success./sqrt(size(press_roi,1));
% for i = 1:cluster.num_cluster  %baseline each curve so it passes through zero
%     shift = (-1)*avg_press_roi(i,3);
%     avg_press_roi(i,:) = avg_press_roi(i,:)+shift;
% end
% bar(ts_press, mean(lick_trace_press)/10); hold on
% for i = 1:size(ts,1);
%     errorbar(ts_press, avg_press_roi(i,:), sm_press(i,:), 'Color', colors(i,:)); hold on;
% end
% xlabel('Time from press (ms)');
% ylabel('dF/F');
% title([days ' presses n=', num2str(length(use_ev_press))]);
% axis tight;
% 

clear
WRITE_VEDIO = 0;
BIN_SIZE = 1;   % if the ROI is large use large bin size (10)

DATA_DIR =  'C:\Users\jake\TempData\';
FRAME_TIME_DIR = 'C:\Users\jake\TempData\';
BEHAVE_DIR = 'Z:\Data\WidefieldImaging\GCaMP\behavior\';
ANALYSIS_DIR ='Z:\Analysis\LeverAnalysis\';
mask = 0;      %set to 1 in order to take all ROIs as one combined mask.  Set to 0 to take TCs of each ROI separately 
% -----------
%days = {'160208_img35', '160209_img36'}; %NO-LEVER DATASETS
%days = {'160208_img36', '160207_img35'}; '160208_img36',  %NO Vis Stim DATASETS
%days = {'160207_img36', '160205_img35'}; %NO Aud Stim DATASETS

%days = {'150514_img24', '150518_img24', '150519_img24', '150521_img24' '150706_img24'} ;
%days = {'150517_img25', '150518_img25', '150514_img25', '150515_img25'};
%days = {'150716_img27', '150718_img27', '150719_img27'};  %'150717_img27',
%days = {'150716_img28', '150717_img28', '150719_img28'}; 
%days = {'151021_img29', '151022_img29'}; %'151015_img29'
%days = {'151009_img30', '151011_img30'};
%days = {'160129_img35', '160131_img35', '160129_img36','160131_img36'};

%days = {'150718_img27', '150719_img27', '150716_img28', '150717_img28', '151021_img29', '151022_img29', '151009_img30', '151011_img30', '151211_img32', '151212_img32', '160129_img35', '160131_img35', '160129_img36','160131_img36', '150518_img24', '150519_img24', '150518_img25', '150517_img25'};
days = {'160131_img36'};

for kk=1:length(days)
    session = '';
    ROI_name  =  days{kk};
    
    frame_info_dest  = [FRAME_TIME_DIR days{kk} '\' ROI_name '_ROI_frame_times.mat'];
    image_dest  = [DATA_DIR days{kk} '\' ROI_name '_ROIshrink.tif'];
    %frame_info_dest = [image_dest(1:end-4) '_frame_times.mat'];
    
    %behave_dest = [BEHAVE_DIR behave_name];
    bfile = dir([BEHAVE_DIR 'data-*i9' days{kk}(end-1:end) '-' days{kk}(1:6) '*' ]);
    behave_dest = [BEHAVE_DIR bfile.name];
    assert(length(bfile)) =1;
    b_data = load(behave_dest);
    ftimes = load(frame_info_dest);
    ifi = (ftimes.frame_times(end)-ftimes.frame_times(1))/length(ftimes.frame_times);
    Sampeling_rate = 1000/ifi;
    info = imfinfo(image_dest);
    
    if(~exist('first_frame', 'var'))
        f_frame =1;
    else
        f_frame = first_frame;
    end
    
    if(~exist('last_frame', 'var'))
        l_frame =min(length(ftimes.frame_times), length(info));
    else
        l_frame = last_frame;
    end
    
    %ftimes.frame_times = []; %load(frame_info_dest);
    
    % [mouse_run, airpuff] = parse_behavior_for_running(b_data.input); % get running and airpuff times
    holdT_min = 400000;
    [lever, frame_info, trial_outcome, lickTimes] = parse_behavior_for_HAD(b_data.input, ...
        f_frame, l_frame, ftimes.frame_times, holdT_min);
    [img, sz]  = get_movie_by_ROI(image_dest, info,[], [], BIN_SIZE, f_frame, l_frame);
    

    % calculate some stats for the raw movie
    avg_img = mean(img,2);
    std_img = std(img,[], 2);
    all_sd = std(img(:));
    
    clear cluster;    %clear existing clusters 
    cluster_file = [image_dest(1:end-4) 'cluster.mat'];
    load(cluster_file); %load cluster
    
    roi_sz = sum(cluster.roi_mask,2); %log the #of pixels in each ROI (binned pixles) 
    
    % ----- create mask of all ROIs combined
    roi_mask_combined = sum(cluster.roi_mask,1);
    roi_combined_sz = sum(roi_mask_combined);
    
    % ----- cluster data by ROI
    if mask == 0; 
        disp('extracting ROIs separately');
        data_tc = [];
        for i = 1:cluster.num_cluster;
            data_tc(i,:) = stackGetTimeCourses(reshape(img,sz(1),sz(2),size(img,2)), reshape(cluster.roi_mask(i,:),[sz(1) sz(2)]));
        end
    else    %----combines all ROIs into one single maske to take a single TC from that mask
        disp('extracting ROIs as a single mask');
        data_tc = [];
        data_tc = stackGetTimeCourses(reshape(img,sz(1),sz(2),size(img,2)), reshape(roi_mask_combined,[sz(1) sz(2)]));
        data_tc= data_tc';
    end
    
    %----------------------------------------------------------
    %Take df/f of TC(s)
    first_baseline = find(~isnan(lever.baseline_timesMs(1,:)),1, 'first');    %find the first trial / baseline_timesMs window that is not NaN
    %first_baseline_hold = find(~isnan(lever.baseline_times_holdMs(1,:)),1, 'first'); 
    for i = 1:length(b_data.input.counterTimesUs);   %finds the MWtime of the first counter
        if find(~isempty(cell2mat(b_data.input.counterTimesUs(i))))==1;
            if length(cell2mat(b_data.input.counterTimesUs(i)))>1;
                StartT = b_data.input.counterTimesUs{i}(1)/1000;
                break
            end
        end
    end
    tc_dfoverf = nan(size(data_tc));
    F_range = [];    
    for iT=frame_info.f_frame_trial_num+1: frame_info.l_frame_trial_num-1;    %only looks at the first and last fully imaged trials
        if ~isnan(lever.baseline_timesMs(1,iT));   %if there is a valid baseline interval then make that the new F
            F_range = frame_info.counter(lever.baseline_timesMs(1,iT)):frame_info.counter(lever.baseline_timesMs(2,iT));
        elseif isempty(F_range)   %if these are trials before there was ever a valid F_range then use the first valid F_range as this trials F_range also. 
            F_range = frame_info.counter(lever.baseline_timesMs(1,first_baseline)):frame_info.counter(lever.baseline_timesMs(2,first_baseline));
        end    %if there was no valid F_range but there was previously a valid F_range then F_range will remain unchanged and the most recent one will be use. 
        if iT == frame_info.f_frame_trial_num+1; %if this is the first fully imaged trial then t_range includes all frames up to this point
            t_range = 1:(frame_info.counter(cell2mat(b_data.input.tThisTrialStartTimeMs(iT+1))-StartT)-1); 
        else   %t_range will be the time interval encompasing the whole trial. Used to find the frames to which we will apply this df/f. 
            t_range = frame_info.counter(cell2mat(b_data.input.tThisTrialStartTimeMs(iT))-StartT):(frame_info.counter(cell2mat(b_data.input.tThisTrialStartTimeMs(iT+1))-StartT)-1);
        end
        if iT == frame_info.l_frame_trial_num-1;
            t_range = (t_range(1)+4):length(data_tc);
        elseif iT == frame_info.f_frame_trial_num+1;
            t_range = 1:(t_range(end)+4); 
        else t_range = t_range + 4;    %added this shift because we have a 1s anaylsis window post release but trial ends 600ms after release.
        end
        for i = 1:size(data_tc,1); %find the avf_f for F_range and apply it to all the frames in that trial
            F_avg= mean(data_tc(i,F_range)); 
            t_df = bsxfun(@minus, data_tc(i, t_range), F_avg);   %use bsxfun because sometimes F_avg is a vector
            t_dfoverf = bsxfun(@rdivide, t_df, F_avg);
            tc_dfoverf(i,t_range) = t_dfoverf;
        end
    end
    %%
    % ---- do simple movie analysis.
    func = @median;
    %func = @mean;
    %func = @std;
    pre_frames = 5;
    post_frames = 10;
    
    
    %  --- first plot  clustering
    colors = [1,0,0;0,1,0;0,0,1;0.5,0.5,0.5; 1,0,1;1,1,0;0,1,1]; %sets up the color scheme for plotting multiple ROIs with errorbar     
    f2 = figure;
    subplot(1,2,1); imagesc(reshape(avg_img, sz(1), sz(2)));
    if b_data.input.doLever == 0;
        title(['CONTROL: no lever ' days{kk}]);
    else
        title(days{kk});
    end
    shading flat; hold on;
    for i=1:cluster.num_cluster
        line(  cluster.roi_position{i}(:,1),   cluster.roi_position{i}(:,2) ,'color', 'w', 'linewidth', 2)
        text(mean(cluster.roi_position{i}(:,1)),mean(cluster.roi_position{i}(:,2)), ...
            [num2str(i)], 'color', 'k', 'FontSize', 30);
    end  
    xlabel(['req=' num2str(b_data.input.fixedReqHoldTimeMs) ' rand=' num2str(b_data.input.randReqHoldMaxMs)]);
    ts = (-pre_frames:post_frames)*1000/round(Sampeling_rate);
    tot_frame = pre_frames + post_frames+1;
   
    %bin licking times according to frame number
    frameTimes = frame_info.times; %start by adjusting frameTimes so it accounts for missing frames
    frameMisses = find(diff(frameTimes)>150);
    frameShift=0;
    for i = frameMisses;
        missedFrame = ((frameTimes(i+frameShift+1)-frameTimes(i+frameShift))/2) + frameTimes(i+frameShift); %approximate the time of the missed frame 
        frameTimes = [frameTimes(1:i+frameShift), missedFrame, frameTimes(i+frameShift+1:end)];
        frameShift=frameShift+1; %necessary in order to account for the fact that frames are being inserted into the TC as the forloop progresses
    end
    licksByFrame = [];
    for i = 1:length(frameTimes);
        if i == length(frameTimes);
            licksThisFrame = sum(lickTimes>=frameTimes(i) & lickTimes<frameTimes(i)+diff(fliplr(frameTimes(i-1:i))));
            licksByFrame = [licksByFrame licksThisFrame];
            break
        end
        licksThisFrame = sum(lickTimes>=frameTimes(i) & lickTimes<frameTimes(i+1));
        licksByFrame = [licksByFrame licksThisFrame];
    end
   
    %STANDARDIZE YLIMS, SAVE VARIABLES, REPORT Ns---------------------------------
%     subplot(2,3,2); ylim([min(YL(:,1)) max(YL(:,2))]);
%     subplot(2,3,3); ylim([min(YL(:,1)) max(YL(:,2))]);
%     subplot(2,3,4); ylim([min(YL(:,1)) max(YL(:,2))]); 
%     subplot(2,3,5); ylim([min(YL(:,1)) max(YL(:,2))]); 
    
%      success_roi = squeeze(success_roi);
%      fail_roi = squeeze(fail_roi);
%      destySucc = strcat(ANALYSIS_DIR, 'LeverSummaryFolder\', days{kk}, '_success');
%      destyFail = strcat(ANALYSIS_DIR, 'LeverSummaryFolder\', days{kk}, '_fail');
%      destyFig = strcat(ANALYSIS_DIR, 'LeverFigureFolder\', days{kk}, '_fig');
%      save([destySucc], 'success_roi');
%      save([destyFail], 'fail_roi');
%      savefig([destyFig]);
%     
%     disp(['day/animal: ' num2str(ROI_name)])
%     disp(['# of successful trials = ' num2str(size(success_roi,1))])
%     disp(['# of failed trials = ' num2str(size(fail_roi,1))])
%     disp(['# of tooFast_successes = ' num2str(size(tooFast_roi,1))])
    
    %PLOT ZEROED AT LEVER PRESS
    use_ev_press = round(lever.press);
    [press_roi, use_times_press, lick_trace_press] = trigger_movie_by_event_licks(tc_dfoverf, frame_info, licksByFrame, use_ev_press, pre_frames, post_frames+35);
    ts_press = (-500:100:4500);
  
    avg_press_roi = squeeze(func(press_roi,1));
    std_success = squeeze(std(press_roi,1));
    sm_press = std_success./sqrt(size(press_roi,1));
    for i = 1:cluster.num_cluster  %baseline each curve so it passes through zero
        shift = (-1)*avg_press_roi(i,3);
        avg_press_roi(i,:) = avg_press_roi(i,:)+shift;
    end
    subplot(1,2,2); bar(ts_press, mean(lick_trace_press)/10); hold on
    for i = 1:cluster.num_cluster;
        subplot(1,2,2); errorbar(ts_press, avg_press_roi(i,:), sm_press(i,:), 'Color', colors(i,:)); hold on;
    end
    xlabel('Time from press(ms)');
    ylabel('dF/F');
    title([days{kk} ' presses n=', num2str(length(use_ev_press))]);
    legend(['avg licks/frame'; cellstr(num2str([1:cluster.num_cluster]'))]);
    axis tight;
    
    %PLOT TC of licking and lever activity
    %find frame numbers for lever presses and releases
    pressTimes=lever.press;
    pressFrameNum=[];
    for i = 1:length(pressTimes);
        thisPressFrame = find(frameTimes<=pressTimes(i),1,'Last');
        pressFrameNum = [pressFrameNum thisPressFrame];
    end
    
    releaseTimes=lever.release;
    releaseFrameNum=[];
    for i = 1:length(releaseTimes);
        thisReleaseFrame = find(frameTimes<=releaseTimes(i),1,'Last');
        releaseFrameNum = [releaseFrameNum thisReleaseFrame];
    end
    figure; bar(licksByFrame, 'k')
    hold on
    scatter(pressFrameNum,ones(1,length(pressFrameNum)),[],'r')
    scatter(releaseFrameNum,ones(1,length(releaseFrameNum)),[],'b','d')
    legend('number of licks', 'lever press', 'lever release')
    xlabel('frame number')
    ylabel('licks per frame')
    title(days{kk})
end