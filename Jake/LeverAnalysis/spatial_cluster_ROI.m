%clear;
clear
WRITE_VEDIO = 0;
BIN_SIZE = 1;   % if the ROI is large use large bin size (10)

DATA_DIR =  'C:\Users\jake\TempData\';
FRAME_TIME_DIR = 'C:\Users\jake\TempData\';
BEHAVE_DIR = 'C:\Users\jake\TempData\behavior\';
% -----------
days = {'150514_img25'};
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
    [lever, frame_info, trial_outcome] = parse_behavior_for_HAD(b_data.input, ...
        f_frame, l_frame, ftimes.frame_times);
    
    
    [img, sz]  = get_movie_by_ROI(image_dest, info,[], [], BIN_SIZE, f_frame, l_frame);
    
    % remove avergae
    avg_img = mean(img,2);
    std_img = std(img,[], 2);
    all_sd = std(img(:));
    
    clear cluster;
    cluster_file = [image_dest(1:end-4) 'cluster.mat'];
    load(cluster_file); % load cluster
    
    % -----  plot results
    roi_sz = sum(cluster.roi_mask,2);
    % ----- cluster data by ROI
    signal = [];
    for i=1:size(img,2)
        for j=1:size(cluster.roi_mask,1)
            signal(i,j) = sum(img(:,i).* cluster.roi_mask(j,:)')/roi_sz(j);
        end
    end
    
    % ---- do simple movie analysis.
    func = @median;
    %func = @mean;
    %func = @std;
    pre_frames = 5;
    post_frames = 30;
    
    %  --- first plot  clustering
    figure;
    subplot(2,2,1);
    title(days{kk});
    im = imagesc(reshape(std_img', sz)); axis ij; colormap jet;
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
    
    [df_success, df_f_success] = adapt_baselines(success_roi, pre_frames);
    avg_success_roi = squeeze(func(df_f_success,1));
    
    subplot(2,2,2); plot(ts, avg_success_roi);
    xlabel('Time from release (ms)');
    ylabel('dF/F');
    title('success');
    axis tight;
    legend(cellstr(num2str([1:cluster.num_cluster]')));
    
    use_ev_fail = trial_outcome.early_time;
    %----- uncomment to use only events w/o lever press after release
    %   use_ev_fail = remove_events_by_lever_state(use_ev_fail,  ...
    %     lever.state, 10,time_after, 0);
    %------ uncomment to use event only w/o lever press before release time
    %   use_ev_fail = remove_events_by_lever_state(use_ev_fail,  ...
    %      lever.state, -time_before,0, 1);
    fail_roi = trigger_movie_by_event(signal', frame_info, ...
        use_ev_fail, pre_frames, post_frames);
    
    [df_fail, df_f_fail] = adapt_baselines(fail_roi, pre_frames);
    avg_fail_roi = squeeze(func(df_f_fail,1));
    
    subplot(2,2,3); plot(ts, avg_fail_roi)
    xlabel('Time from release (ms)');
    ylabel('dF/F');
    title('fail');
    axis tight;
    legend(cellstr(num2str([1:cluster.num_cluster]')));
    subplot(2,2,4); plot(ts, avg_success_roi - avg_fail_roi);
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
m =  squeeze(mean(success_trials,1));
sd = squeeze(std(success_trials,[], 1));
sm  = sd./sqrt(size(success_trials,1));

plot(ts,m );  axis tight;
%errorbar(repmat(ts, size( success_trials,2),1)',m',sm' );  axis tight;
xlabel('Time from lever release (ms)');
ylabel('df/f');
legend(cellstr(num2str([1:cluster.num_cluster]')));
m2 =  squeeze(mean(fail_trials,1));
sd2 = squeeze(std(fail_trials,[], 1));
sm2  = sd2./sqrt(size(fail_trials,1));
title('success');
subplot(2,2,2);
%errorbar(repmat(ts, size( fail_trials,2),1)',m2',sm2' );  axis tight;
plot(ts,m2);  axis tight;
xlabel('Time from lever release (ms)');
ylabel('df/f');
title('fail');
m2 =  squeeze(mean(fail_trials,1));

subplot(2,2,3);
plot(ts,m-m2);  axis tight;
xlabel('Time from lever release (ms)');
ylabel('df/f');
title('success - fail');

