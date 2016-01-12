%generates frame by frame movies of an ROI. Heatmaps. 
%redefines df/f using lever.baslineTimes
clear
WRITE_VEDIO = 0;
BIN_SIZE = 1;   % if the ROI is large use large bin size (10)


% DATA_DIR =  '/Volumes/Promise RAID/mati/imaging_data/';
% BEHAVE_DIR = '/Volumes/Promise RAID/mati/behave_data/';
DATA_DIR =  'C:\Users\jake\TempData\';
BEHAVE_DIR = 'C:\Users\jake\TempData\behavior\';

% ----------
%days = {'150320_img20'};
days = {'150706_img24'};
holdT_min = 500000;

for kk=1:length(days)
    session = '';
    ROI_name  =  days{kk};
    
    image_dest  = [DATA_DIR days{kk} '\' session ROI_name '_ROI.tif'];
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
        f_frame, l_frame, ftimes.frame_times, holdT_min);
    
    info = imfinfo(image_dest);
    [img, sz]  = get_movie_by_ROI(image_dest, info,[], [], BIN_SIZE, f_frame, l_frame);
   
    
  %Obtain a df/f movie using Lindsey's baseline_times
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
    F_avg= mean(img(:,F_range),2);
    
    %need to write a smart t_range which includes all fully imaged trials
    %before the first trial with a valid baseline
    t_range = frame_info.counter(cell2mat(b_data.input.tThisTrialStartTimeMs(iT))-StartT):(frame_info.counter(cell2mat(b_data.input.tThisTrialStartTimeMs(iT+1))-StartT)-1);
    if iT == frame_info.l_frame_trial_num-1;
        t_range = t_range(1:(end-4)) + 4;
    else t_range = t_range + 4;    %added this shift because we have a 1s anaylsis window post release but trial ends 600ms after release.
    end
    %problematic bcimg it looks at times before the first counter then
    %subtracts the time of the first counter
    t_df = bsxfun(@minus, double(img(:,t_range)), F_avg);
    t_dfoverf = bsxfun(@rdivide, t_df, F_avg);
    img_dfoverf(:,t_range) = t_dfoverf;
end 
%remember to cut nondf/f frames (beginning and end) from movie

% PLOTTING 
    % ---- do simple movie analysis
    func = @median;
    % func = @mean;
    %func = @std;
    pre_frames = 5;
    post_frames = 10;  
    rm_baseline_plot = 0; % 1 for removing baseline when plotiing
    
    ts = (-pre_frames:post_frames)*1000/round(Sampeling_rate);
    tot_frame = pre_frames + post_frames+1;
    
    use_ev_success = trial_outcome.success_time;
    %get rid of first and last trials
    if strcmp(b_data.input.trialOutcomeCell{1}, 'success')
        use_ev_success(1) = [];
    elseif strcmp(b_data.input.trialOutcomeCell{end}, 'success')
        use_ev_success(end) = [];
    end
    %----- uncomment to use only events w/o lever press after release
    %     use_ev_success = remove_events_by_lever_state(use_ev_success,  ...
    %         lever.state, 10,ceil(post_frames*1000/Sampeling_rate), 0);
    %------ uncomment to use event only w/o lever press before release time
    %     use_ev_success = remove_events_by_lever_state(use_ev_success,  ...
    %         lever.state, -ceil(pre_frames*1000/Sampeling_rate),0, 1);
    %
    success_movie = trigger_movie_by_event(img_dfoverf, frame_info, ...
        use_ev_success, pre_frames, post_frames);
    avg_success_move = squeeze(func(success_movie,1));
    fig1 = figure; plot_movie(avg_success_move,sz,rm_baseline_plot, ts);
    title('Hits');
    clim([-0.5 0.5]);
    disp(['hits: n = ' num2str(size(success_movie,1))]);
    
    use_ev_fail = trial_outcome.early_time;
    %get rid of first and last trials
    if strcmp(b_data.input.trialOutcomeCell{1}, 'failure')
        use_ev_fail(1) = [];
    elseif strcmp(b_data.input.trialOutcomeCell{end}, 'failure')
        use_ev_fail(end) = [];
    end
    %----- uncomment to use only events w/o lever press after release
    %     use_ev_fail = remove_events_by_lever_state(use_ev_fail,  ...
    %         lever.state, 10,ceil(post_frames*1000/Sampeling_rate), 0);
    %------ uncomment to use event only w/o lever press before release time
    %     use_ev_fail = remove_events_by_lever_state(use_ev_fail,  ...
    %         lever.state, -ceil(pre_frames*1000/Sampeling_rate),0, 1);
    %
    
    % -----trigger movie by early release
    fail_movie = trigger_movie_by_event(img_dfoverf, frame_info, ...
        use_ev_fail, pre_frames, post_frames);
    avg_fail_move = squeeze(func(fail_movie,1));
    fig2 = figure; plot_movie(avg_fail_move,sz,rm_baseline_plot, ts);
    title('False Alarms');
    disp(['false alarms: n = ' num2str(size(fail_movie,1))]);
    clim([-0.5 0.5]);
    
    % ----- success  - early!!
    diff_succ_fail = avg_success_move - avg_fail_move;
    fig3 = figure; plot_movie(diff_succ_fail,sz,rm_baseline_plot, ts);    %JH ALTERED   Added "fig3"
    title('Hits  - False Alarms');
    figure(fig3);                             
    clim([-0.5 0.5]);
 
    % -- Trigger movie off lever press IF subsequent hold time is >1000ms
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
        if holdDuration(i) > 1000
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
    press_movie = trigger_movie_by_event(img_dfoverf, frame_info, ...
        use_ev_press, pre_frames, post_frames);
    avg_press_move = squeeze(func(press_movie,1));
    fig4 = figure; plot_movie(avg_press_move,sz,rm_baseline_plot, ts);
    title('lever press');
    clim([-0.5 0.5]);
    disp(['press: n = ' num2str(size(press_movie,1))]);
    
    %PLOT TRIAL BY TRIAL VARIABILITY
    %Hits
    subset1 = round(linspace(2,size(success_movie,1),10));
    success_subset = [];
    success_subset = success_movie(30:39,:,:);
    fig5 = figure; plot_movie(success_subset,sz,rm_baseline_plot, ts);
    title('Hits');
    clim([-0.5 0.5])   % NEED a smart way to do scaling
    
    %False Alarms   
    subset2 = round(linspace(2,size(fail_movie,1),10));
    fail_subset = [];
    fail_subset = fail_movie(10:19,:,:);
    fig6 = figure; plot_movie(fail_subset,sz,rm_baseline_plot, ts);
    title('False Alarms');
    clim([-0.5 0.5])
    
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


