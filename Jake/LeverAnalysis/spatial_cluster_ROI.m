%clear;
clear
WRITE_VEDIO = 0;
BIN_SIZE = 1;   % if the ROI is large use large bin size (10)

DATA_DIR =  'C:\Users\jake\TempData\';
FRAME_TIME_DIR = 'C:\Users\jake\TempData\';
BEHAVE_DIR = 'C:\Users\jake\TempData\behavior\';
mask = 1;      %set to 1 in order to take all ROIs as one combined mask.  Set to 0 to take TCs of each ROI separately 
% -----------
days = {'150718_img28'};
%days = {'150206_img16','150207_img16', '150209_img16','150210_img16','150211_img16','150212_img16','150213_img16', '150214_img16'} ;


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
    
    %------------Obtain a df/f movie using Lindsey's baseline_times
    first_baseline = find(~isnan(lever.baseline_timesMs(1,:)),1, 'first');    %find the first trial / baseline_timesMs window that is not NaN
    for i = 1:length(b_data.input.counterTimesUs);   %finds the MWtime of the first counter
        if find(~isempty(cell2mat(b_data.input.counterTimesUs(i))))==1;
            StartT = b_data.input.counterTimesUs{i}(1)/1000;
            break
        end
    end
    img_dfoverf = zeros(size(img));
    F_range = [];
    for iT=frame_info.f_frame_trial_num+1: frame_info.l_frame_trial_num-1;    %only looks at the first and last fully imaged trials
        if ~isnan(lever.baseline_timesMs(1,iT));
            F_range = frame_info.counter(lever.baseline_timesMs(1,iT)):frame_info.counter(lever.baseline_timesMs(2,iT));
        elseif isempty(F_range)
            F_range = frame_info.counter(lever.baseline_timesMs(1,first_baseline)):frame_info.counter(lever.baseline_timesMs(2,first_baseline));
        end
        F_avg= mean(img(:,F_range),2); %there are zeroes in F_avg which result in a 0/0 situation later on
        %F_avg(F_avg==0) = 0.1;   %if F_avg has values that =0 then this will force the code to divide by zero when calculating df/f. Replaced zeroes with a very low value (0.1 on a scale of 1:256)
        %need to write a smart t_range which includes all fully imaged trials
        %before the first trial with a valid baseline
        t_range = frame_info.counter(cell2mat(b_data.input.tThisTrialStartTimeMs(iT))-StartT):(frame_info.counter(cell2mat(b_data.input.tThisTrialStartTimeMs(iT+1))-StartT)-1);
        if iT == frame_info.l_frame_trial_num-1;
            t_range = t_range(1:(end-4)) + 4;
        else t_range = t_range + 4;    %added this shift because we have a 1s anaylsis window post release but trial ends 600ms after release.
        end
        %problematic bc it looks at times before the first counter then
        %subtracts the time of the first counter
        t_df = bsxfun(@minus, double(img(:,t_range)), F_avg);
        t_dfoverf = bsxfun(@rdivide, t_df, F_avg);   %this is where the NaNs happen. CAUSE 0/0 = NaN     There are no NaNs in img_dfover when I run HAD_comp
        img_dfoverf(:,t_range) = t_dfoverf;
    end  %---------------- %remember to cut nondf/f frames (beginning and end) from movie
    img_mean = mean(img,2);
    img = img_dfoverf;



    % remove average
    avg_img = mean(img,2);
    std_img = std(img,[], 2);
    all_sd = std(img(:));
    
    clear cluster;
    cluster_file = [image_dest(1:end-4) 'cluster.mat'];
    load(cluster_file); % load cluster
    
    % -----  plot results
    roi_sz = sum(cluster.roi_mask,2);
    
    % ----- create mask of all ROIs combined
    roi_mask_combined = sum(cluster.roi_mask,1);
    roi_combined_sz = sum(roi_mask_combined);
    if mask == 0;
        % ----- cluster data by ROI
        signal = [];
        disp('plotting ROIs separately');
        for i=1:size(img,2)
            for j=1:size(cluster.roi_mask,1)
                signal(i,j) = nansum(img(:,i).* cluster.roi_mask(j,:)')/roi_sz(j);
            end
        end
    else    %----combines all ROIs into one single maske to take a single TC from that mask
        signal = [];
        disp('plotting ROIs as a single mask');
        for i=1:size(img,2)
            signal(i) = nansum(img(:,i).* roi_mask_combined')/roi_combined_sz;
        end
        signal = signal'; %for some reason signal comes out inverted in this condition so you have to fix it. 
    end
    
    % ---- do simple movie analysis.
    func = @median;
    %func = @mean;
    %func = @std;
    pre_frames = 5;
    post_frames = 20;
    
    %  --- first plot  clustering
    figure;
    subplot(2,2,1);
    title(days{kk});
    figure; imagesc(reshape(img_mean, sz(1), sz(2)));
    im = imagesc(reshape(avg_img', sz)); axis ij; colormap jet;   % it looks weird bc I am viewing the df/f
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
    success_roi = trigger_movie_by_event(signal', frame_info, ...
        use_ev_success, pre_frames, post_frames);
    if mask == 0;   %plotting ROI condition
        %[df_success, df_f_success] = adapt_baselines(success_roi, pre_frames);
        avg_success_roi = squeeze(func(df_f_success,1));
        subplot(2,2,2); plot(ts, avg_success_roi);
        xlabel('Time from release (ms)');
        ylabel('dF/F');
        title('success');
        axis tight;
        legend(cellstr(num2str([1:cluster.num_cluster]')));
    else    %plotting mask condition
        avg_success_roi = squeeze(func(success_roi,1));
        std_success= std(squeeze(success_roi));
        sm_success = std_success./sqrt(size(success_roi,1))
        subplot(2,2,2); errorbar(ts,avg_success_roi,sm_success);
        xlabel('Time from release (ms)');
        ylabel('dF/F');
        title('success');
        axis tight;
        %legend(cellstr(num2str([1:cluster.num_cluster]')));
    end
    
    use_ev_fail = trial_outcome.early_time;
    %----- uncomment to use only events w/o lever press after release
    %   use_ev_fail = remove_events_by_lever_state(use_ev_fail,  ...
    %     lever.state, 10,time_after, 0);
    %------ uncomment to use event only w/o lever press before release time
    %   use_ev_fail = remove_events_by_lever_state(use_ev_fail,  ...
    %      lever.state, -time_before,0, 1);
    fail_roi = trigger_movie_by_event(signal', frame_info, ...
        use_ev_fail, pre_frames, post_frames);
    if mask ==0;
        %[df_fail, df_f_fail] = adapt_baselines(fail_roi, pre_frames);
        avg_fail_roi = squeeze(func(df_f_fail,1));
        subplot(2,2,3); plot(ts, avg_fail_roi)
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
    
    sub_sm = sqrt(sm_fail.^2+sm_success.^2);
    subplot(2,2,4); errorbar(ts, (avg_success_roi - avg_fail_roi), sub_sm);
    axis tight;
    xlabel('Time from release (ms)');
    ylabel('dF/F');
    title('success - fail');
    success_trials(kk,:,:) = avg_success_roi;
    fail_trials(kk,:,:) = avg_fail_roi;
end

%  ----- plot summary
figure;
ts = (-pre_frames:post_frames)*1000/round(Sampeling_rate);

subplot(2,2,1);
m =  squeeze(mean(squeeze(success_trials),1));
sd = squeeze(std(squeeze(success_trials),[], 1));
sm  = sd./sqrt(size(success_trials,1));
errorbar(ts, m', sm'); axis tight;
xlabel('Time from lever release (ms)');
ylabel('df/f');
legend(cellstr(num2str([1:cluster.num_cluster]')));
m2 =  squeeze(mean(fail_trials,1));
sd2 = squeeze(std(fail_trials,[], 1));
sm2  = sd2./sqrt(size(fail_trials,1));
title('success summary');
subplot(2,2,2);
plot(ts,m2);  axis tight;
xlabel('Time from lever release (ms)');
ylabel('df/f');
title('fail summary');
m2 =  squeeze(mean(fail_trials,1));

subplot(2,2,3);
plot(ts,m-m2);  axis tight;
xlabel('Time from lever release (ms)');
ylabel('df/f');
title('success - fail summary');

