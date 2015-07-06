%clear;
clear
WRITE_VEDIO = 0;
BIN_SIZE = 1;   % if the ROI is large use large bin size (10)


% DATA_DIR =  '/Volumes/Promise RAID/mati/imaging_data/';
% BEHAVE_DIR = '/Volumes/Promise RAID/mati/behave_data/';
DATA_DIR =  'C:\Users\jake\TempData\';
BEHAVE_DIR = 'C:\Users\jake\TempData\behavior\';

% ----------
%days = {'150320_img20'};
days = {'150518_img24'};

for kk=1:length(days)
    session = '';
    ROI_name  =  days{kk};
    
    image_dest  = [DATA_DIR days{kk} '\' session '\' ROI_name '_ROI.tif'];
    frame_info_dest = [image_dest(1:end-4) '_frame_times.mat'];
    % ---- load behavior file
    bfile = dir([BEHAVE_DIR 'data-*i9' days{kk}(end-1:end) '-' days{kk}(1:6) '*' ]);
    behave_dest = [BEHAVE_DIR bfile.name];
    assert(length(bfile)) =1; % only one behavior file per day
    b_data = load(behave_dest);
    % --- load frame times
    ftimes = load(frame_info_dest);
    ifi = (ftimes.frame_times(end)-ftimes.frame_times(1))/length(ftimes.frame_times);
    Sampeling_rate = 1000/ifi;
    
    if(~exist('first_frame', 'var'))
        f_frame =1;
    else
        f_frame = first_frame;
    end
    
    if(~exist('last_frame', 'var'))
        l_frame =length(ftimes.frame_times);
    else
        l_frame = last_frame;
    end
    % ---- parse behavior
    [lever, frame_info, trial_outcome] = parse_behavior_for_HAD(b_data.input, ...
        f_frame, l_frame, ftimes.frame_times);
    
    info = imfinfo(image_dest);
    [img, sz]  = get_movie_by_ROI(image_dest, info,[], [], BIN_SIZE, f_frame, l_frame);
   
    
    %Obtain a df/f movie using Lindsey's baseline_times
    img_dfoverf = zeros(size(img));           %this could be problematic due to the frame skipping issue
    for iT=1:length(lever.baseline_timesMs);    %this could be problematic due to unremoved NaNs
        if ~isnan(lever.baseline_timesMs(1,iT));
            F_range = frame_info.counter(lever.baseline_timesMs(1,iT)):frame_info.counter(lever.baseline_timesMs(2,iT));
            t_range = frame_info.counter(cell2mat(b_data.input.tThisTrialStartTimeMs(iT))):frame_info.counter(cell2mat(b_data.input.tThisTrialStartTimeMs(iT+1)));
        end
        F_avg= mean(img(:,:,F_range),3);
        %t_df = img(:,:,t_range)-F_avg;
        t_df = bsxfun(@minus, img(:,:,t_range), F_avg);
        t_dfoverf = bsxfun(@rdivide, t_df, F_avg);
        img_dfoverf(:,:,t_range) = t_dfoverf;
    end
       

        
        
        
        
        
    % remove avergae
    raw_img = img;
    avg_img = mean(img,2);
    std_img = std(img,[], 2);
    all_sd = std(img(:));
    for i=1:size(img,2)  % subtract average from movie
        % df
        img(:,i)  =  (img(:,i) - avg_img); %this is actually dF
        % df/SD
        %img(:,i)  =  (img(:,i) - avg_img)./all_sd; % I have no idea what this is.
        % df/std
        %img(:,i)  =  (img(:,i) - avg_img)./std_img;
        % df/f
        %img(:,i)  =  (img(:,i) - avg_img)./avg_img;
    end
    
    
    % ---- do simple movie analysis
    func = @median;
    % func = @mean;
    %func = @std;
    pre_frames = 5;
    post_frames = 10;
    post_press_frames = 10;   
    rm_baseline_plot = 0; % 1 for removing baseline when plotiing
    
    ts = (-pre_frames:post_frames)*1000/round(Sampeling_rate);
    tot_frame = pre_frames + post_frames+1;
    
    use_ev_success = trial_outcome.success_time;
    %----- uncomment to use only events w/o lever press after release
    %     use_ev_success = remove_events_by_lever_state(use_ev_success,  ...
    %         lever.state, 10,ceil(post_frames*1000/Sampeling_rate), 0);
    %------ uncomment to use event only w/o lever press before release time
    %     use_ev_success = remove_events_by_lever_state(use_ev_success,  ...
    %         lever.state, -ceil(pre_frames*1000/Sampeling_rate),0, 1);
    %
    success_movie = trigger_movie_by_event(img, frame_info, ...
        use_ev_success, pre_frames, post_frames);
    avg_success_move = squeeze(func(success_movie,1));
    fig1 = figure; plot_movie(avg_success_move,sz,rm_baseline_plot, ts);
    title('Hits');
    
    use_ev_fail = trial_outcome.early_time;
    %----- uncomment to use only events w/o lever press after release
    %     use_ev_fail = remove_events_by_lever_state(use_ev_fail,  ...
    %         lever.state, 10,ceil(post_frames*1000/Sampeling_rate), 0);
    %------ uncomment to use event only w/o lever press before release time
    %     use_ev_fail = remove_events_by_lever_state(use_ev_fail,  ...
    %         lever.state, -ceil(pre_frames*1000/Sampeling_rate),0, 1);
    %
    
    % -----trigger movie by early release
    fail_movie = trigger_movie_by_event(img, frame_info, ...
        use_ev_fail, pre_frames, post_frames);
    avg_fail_move = squeeze(func(fail_movie,1));
    fig2 = figure; plot_movie(avg_fail_move,sz,rm_baseline_plot, ts);
    title('False Alarms');
    % --- make color scale the same
    min_val = min([avg_fail_move(:);avg_success_move(:)]);
    max_val = max([avg_fail_move(:);avg_success_move(:)]);
    figure(fig1); caxis([min_val max_val]);
    figure(fig2); caxis([min_val max_val]);
    
    % ----- success  - early!!
    diff_succ_fail = avg_success_move - avg_fail_move;
    fig3 = figure; plot_movie(diff_succ_fail,sz,rm_baseline_plot, ts);    %JH ALTERED   Added "fig3"
    title('Hits  - False Alarms');
    figure(fig3); caxis([min_val max_val]);                             %JH ADDED this entire line to make axes constant
    
 
    % ---- Trigger movie off lever press IF subsequent hold time is greater
    % than 500ms
    post_frames = post_press_frames;
    if lever.release(1) < lever.press(1)
        holdDuration = NaN(1, size((lever.release),2)-1);
        for i = 1:(length(holdDuration))
            holdDuration(i) = lever.release(i+1)-lever.press(i);
        end
    else
        holdDuration = NaN(1, size((lever.release),2));
        for i = 1:(length(holdDuration))
            holdDuration(i) = lever.release(i)-lever.press(i);
        end
    end
    
    longHoldInd = zeros(size(holdDuration));
    for i = 1:size(holdDuration, 2);
        if holdDuration(i) > 500
            longHoldInd(i) = 1;
        end
    end
    
    use_ev_press = NaN(1, sum(longHoldInd));
    aa = 1;
    for i = 1:length(longHoldInd)
        if longHoldInd(i) == 1;
            use_ev_press(aa) = lever.press(i);
            aa=aa+1;
        end
    end
    use_ev_press = round(use_ev_press);
    %----- uncomment to use only events w/o lever press after release
    %     use_ev_success = remove_events_by_lever_state(use_ev_success,  ...
    %         lever.state, 10,ceil(post_frames*1000/Sampeling_rate), 0);
    %------ uncomment to use event only w/o lever press before release time
    %     use_ev_success = remove_events_by_lever_state(use_ev_success,  ...
    %         lever.state, -ceil(pre_frames*1000/Sampeling_rate),0, 1);
    %
    press_movie = trigger_movie_by_event(img, frame_info, ...
        use_ev_press, pre_frames, post_frames);
    avg_press_move = squeeze(func(press_movie,1));
    fig4 = figure; plot_movie(avg_press_move,sz,rm_baseline_plot, ts);
    title('lever press');
    figure(fig4);  caxis([min_val max_val]);     %JH ADDED this line to make axes constant
    
    %PLOT TRIAL BY TRIAL VARIABILITY
    %Hits
    subset1 = round(linspace(1,size(success_movie,1),10));
    success_subset = [];
    success_subset = success_movie(subset1,:,:);
    fig5 = figure; plot_movie(success_subset,sz,rm_baseline_plot, ts);
    title('Hits');
    clim([-5 12])   % NEED a smart way to do scaling
    
    %False Alarms   
    subset2 = round(linspace(1,size(fail_movie,1),10));
    fail_subset = [];
    fail_subset = success_movie(30:39,:,:);
    fig6 = figure; plot_movie(fail_subset,sz,rm_baseline_plot, ts);
    title('False Alarms');
    clim([-12 30])
    
    %Lever Presses
    %Hits-FAsNot sure this will be informative. Also will require a
    %lot of tweeking
%     subset3 = round(linspace(1,size(_movie,1),10));
%     _subset = [];
%     _subset = success_movie(subset3,:,:);
%     fig5 = figure; plot_movie(_subset,sz,rm_baseline_plot, ts);
%     title('Hits - False Alarms');
%     clim([-5 12])


end


