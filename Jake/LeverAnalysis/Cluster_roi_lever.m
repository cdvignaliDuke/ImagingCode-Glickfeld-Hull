%CLUSTER ROI
clear
WRITE_VEDIO = 0;
BIN_SIZE = 1;   % if the ROI is large use large bin size (10)

DATA_DIR =  'C:\Users\jake\TempData\';
FRAME_TIME_DIR = 'C:\Users\jake\TempData\';
BEHAVE_DIR = 'C:\Users\jake\TempData\behavior\';
mask = 1;      %set to 1 in order to take all ROIs as one combined mask.  Set to 0 to take TCs of each ROI separately 
% -----------
%days = {'150725_img27'};
%days = {'150716_img27', '150718_img27', '150719_img27'};  %'150717_img27',
%days = {'150514_img24', '150518_img24', '150519_img24', '150521_img24' '150706_img24'} ;
%days = {'150716_img28', '150717_img28', '150719_img28'}; 
%days = {'150517_img25', '150518_img25', '150514_img25', '150515_img25'};
days = {'150514_img24', '150518_img24', '150519_img24', '150521_img24' '150706_img24', '150517_img25', '150518_img25', '150514_img25', '150515_img25', '150716_img27', '150718_img27', '150719_img27', '150716_img28', '150717_img28', '150719_img28'}; 

success_trials = [];
fail_trials = [];

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
    %
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
    [lever, frame_info, trial_outcome] = parse_behavior_for_HAD(b_data.input, ...
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
    for i = 1:length(b_data.input.counterTimesUs);   %finds the MWtime of the first counter
        if find(~isempty(cell2mat(b_data.input.counterTimesUs(i))))==1;
            StartT = b_data.input.counterTimesUs{i}(1)/1000;
            break
        end
    end
    tc_dfoverf = nan(size(data_tc));
    F_range = [];  
    for iT=frame_info.f_frame_trial_num+1: frame_info.l_frame_trial_num-1;    %only looks at the first and last fully imaged trials
        %F_range is the # of each frame which we will use to generate f
        if ~isnan(lever.baseline_timesMs(1,iT));   %if there is a valid baseline interval then make that the new F
            F_range = frame_info.counter(lever.baseline_timesMs(1,iT)):frame_info.counter(lever.baseline_timesMs(2,iT));
        elseif isempty(F_range)    
            F_range = frame_info.counter(lever.baseline_timesMs(1,first_baseline)):frame_info.counter(lever.baseline_timesMs(2,first_baseline));
        end
        %assign the frame numbers which correspond to this trial
        if iT == frame_info.f_frame_trial_num+1; %if this is the first fully imaged trial then t_range includes all frames up to this point
            t_range = 1:(frame_info.counter(cell2mat(b_data.input.tThisTrialStartTimeMs(iT+1))-StartT)-1);
        else
            t_range = frame_info.counter(cell2mat(b_data.input.tThisTrialStartTimeMs(iT))-StartT):(frame_info.counter(cell2mat(b_data.input.tThisTrialStartTimeMs(iT+1))-StartT)-1);
        end
        if iT == frame_info.l_frame_trial_num-1;
            t_range = (t_range(1)+4):length(data_tc);
        elseif iT == frame_info.f_frame_trial_num+1;
            t_range = 1:(t_range(end)+4);
        else t_range = t_range + 4;    %added this shift because we have a 1s anaylsis window post release but trial ends 600ms after release.
        end
        for i = 1:size(data_tc,1);
            F_avg= mean(data_tc(i,F_range));
            t_df = bsxfun(@minus, data_tc(i, t_range), F_avg);   %use bsxfun because sometimes F_avg is a vector
            t_dfoverf = bsxfun(@rdivide, t_df, F_avg);
            tc_dfoverf(i,t_range) = t_dfoverf;
        end
    end  %---------------- %remember to cut nondf/f frames (beginning and end) from movie
    disp([mat2str(sum(isnan(tc_dfoverf))) ' NaNs remaining in tc_dfoverf']);
    %----------------------------------------------------------------------------------------------------------------------
    
    
    %plot the corr coef
    if mask == 0   
        cnames = cell([1,size(tc_dfoverf,1)]);
        rnames = cell([1,size(tc_dfoverf,1)]);
    coefMat = corrcoef(tc_dfoverf');
    for i = 1:size(tc_dfoverf,1);
        cnames{1,i} = ['ROI' mat2str(i)];
        rnames{1,i} = ['ROI' mat2str(i)];
    end
    f = figure;
    title(days{kk});
    t = uitable(f, 'Data', coefMat,...
        'ColumnName', cnames,...
        'RowName', rnames);
    end     
    %
    % ---- do simple movie analysis.
    func = @median;
    %func = @mean;
    %func = @std;
    pre_frames = 5;
    post_frames = 20;
    




    %  --- first plot  clustering
    f = figure;      
    subplot(2,2,1); imagesc(reshape(avg_img, sz(1), sz(2)));
    title(days{kk});
    shading flat; hold on;
    for i=1:cluster.num_cluster
        line(  cluster.roi_position{i}(:,1),   cluster.roi_position{i}(:,2) ,'color', 'w', 'linewidth', 2)
        text(mean(cluster.roi_position{i}(:,1)),mean(cluster.roi_position{i}(:,2)), ...
            [num2str(i)], 'color', 'k', 'FontSize', 30);
    end  
    
    ts = (-pre_frames:post_frames)*1000/round(Sampeling_rate);
    tot_frame = pre_frames + post_frames+1;
    use_ev_success = trial_outcome.success_time;
    
    
    time_before = 500; % in ms, time before event w/o release
    time_after =1000; % in ms, time after event w/o press
    
    %----- uncomment to use only events w/o lever press after release
    %     use_ev_success = remove_events_by_lever_state(use_ev_success,  ...
    %         lever.state, 10,time_after, 0);
    %     %------ uncomment to use event only w/o lever press before release time
    %     use_ev_success = remove_events_by_lever_state(use_ev_success,  ...
    %         lever.state, -time_before,0, 1);
    success_roi = trigger_movie_by_event(tc_dfoverf, frame_info, ...
        use_ev_success, pre_frames, post_frames);
    if mask == 0;   %plotting ROI condition  
        ts = repmat(ts,[cluster.num_cluster 1]);
        avg_success_roi = squeeze(func(success_roi,1));
        std_success = squeeze(std(success_roi,1));
        sm_success = std_success./sqrt(size(success_roi,1));
        for i = 1:size(ts,1);
        subplot(2,2,2); errorbar(ts(i,:), avg_success_roi(i,:), sm_success(i,:)); hold on;
        end
        xlabel('Time from release (ms)');
        ylabel('dF/F');
        title('success');
        %axis tight;
        %legend(cellstr(num2str([1:cluster.num_cluster]')));
    else    %plotting mask condition
        avg_success_roi = squeeze(func(success_roi,1));
        std_success= std(squeeze(success_roi));
        sm_success = std_success./sqrt(size(success_roi,1));
        subplot(2,2,2); errorbar(ts,avg_success_roi,sm_success);
        xlabel('Time from release (ms)');
        ylabel('dF/F');
        title('success');
        axis tight;
        %legend(cellstr(num2str([1:cluster.num_cluster]')));
    end
    hold off   
    use_ev_fail = trial_outcome.early_time;
    %----- uncomment to use only events w/o lever press after release
    %   use_ev_fail = remove_events_by_lever_state(use_ev_fail,  ...
    %     lever.state, 10,time_after, 0);
    %------ uncomment to use event only w/o lever press before release time
    %   use_ev_fail = remove_events_by_lever_state(use_ev_fail,  ...
    %      lever.state, -time_before,0, 1);
    fail_roi = trigger_movie_by_event(tc_dfoverf, frame_info, ...
        use_ev_fail, pre_frames, post_frames);
    if mask ==0;        
        avg_fail_roi = squeeze(func(fail_roi,1));
        std_fail = squeeze(std(fail_roi,1));
        sm_fail = std_fail./sqrt(size(fail_roi,1));
        for i = 1:size(ts,1);
        hold on; subplot(2,2,3); errorbar(ts(i,:), avg_fail_roi(i,:), sm_fail(i,:));
        end
        xlabel('Time from release (ms)');
        ylabel('dF/F');
        title('fail');
        axis tight;
        legend(cellstr(num2str([1:cluster.num_cluster]')));
    else
        avg_fail_roi = squeeze(func(fail_roi,1));
        std_fail= std(squeeze(fail_roi));
        sm_fail = std_fail./sqrt(size(fail_roi,1));
        subplot(2,2,3); errorbar(ts, avg_fail_roi, sm_fail);
        xlabel('Time from release (ms)');
        ylabel('dF/F');
        title('fail');
        axis tight;
    end
    hold off        
    sub_sm = sqrt(sm_fail.^2+sm_success.^2); 
    if mask == 0
        for i = 1:size(ts,1);
            hold on; subplot(2,2,4); errorbar(ts(i,:), avg_success_roi(i,:) - avg_fail_roi(i,:), sm_fail(i,:));
        end
    else
        subplot(2,2,4); errorbar(ts, (avg_success_roi - avg_fail_roi), sub_sm);
        axis tight;
    end
    xlabel('Time from release (ms)');
    ylabel('dF/F');
    title('success - fail');
    hold off
    if mask ==1;
        success_trials(kk,:) = avg_success_roi;
        fail_trials(kk,:) = avg_fail_roi;
    end
    %set y scale to be equal for all 3 subplots
    YL = []
    for i =2:4
        subplot(2,2,i); YL(i-1,:) = ylim;
    end
    subplot(2,2,2); ylim([min(YL(:,1)) max(YL(:,2))]);
    subplot(2,2,3); ylim([min(YL(:,1)) max(YL(:,2))]);
    subplot(2,2,4); ylim([min(YL(:,1)) max(YL(:,2))]); 
    success_roi = squeeze(success_roi);
    fail_roi = squeeze(fail_roi);
    destySucc = strcat(DATA_DIR, 'summaryFolder\', days{kk}, '_success')
    destyFail = strcat(DATA_DIR, 'summaryFolder\', days{kk}, '_fail')
    save([destySucc], 'success_roi');
    save([destyFail], 'fail_roi');
end
if mask == 1;
    %  ----- plot summary
    figure;
    ts = (-pre_frames:post_frames)*1000/round(Sampeling_rate);
    
    subplot(2,2,1);
    m =  squeeze(mean(squeeze(success_trials),1));
    sd = squeeze(std(squeeze(success_trials),[], 1));
    sm  = sd./sqrt(size(success_trials,1));
    errorbar(ts, m, sm); axis tight;
    xlabel('Time from lever release (ms)');
    ylabel('df/f');
    %legend(cellstr(num2str([1:cluster.num_cluster]')));
    m2 =  squeeze(mean(fail_trials,1));
    sd2 = squeeze(std(fail_trials,[], 1));
    sm2  = sd2./sqrt(size(fail_trials,1));
    title('success summary');
    subplot(2,2,2);
    errorbar(ts,m2,sm2);  axis tight;
    xlabel('Time from lever release (ms)');
    ylabel('df/f');
    title('fail summary');
    m2 =  squeeze(mean(fail_trials,1));
    
    subplot(2,2,3);
    errorbar(ts,m-m2,sqrt(sm.^2+sm2.^2));  axis tight;
    xlabel('Time from lever release (ms)');
    ylabel('df/f');
    title('success - fail summary');
end