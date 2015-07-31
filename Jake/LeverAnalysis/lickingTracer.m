%CLUSTER ROI
clear
WRITE_VEDIO = 0;
BIN_SIZE = 1;   % if the ROI is large use large bin size (10)

DATA_DIR =  'C:\Users\jake\TempData\';
FRAME_TIME_DIR = 'C:\Users\jake\TempData\';
BEHAVE_DIR = 'C:\Users\jake\TempData\behavior\';
mask = 1;      %set to 1 in order to take all ROIs as one combined mask.  Set to 0 to take TCs of each ROI separately 
% -----------
%days = {'150519_img24'};
%days = {'150514_img24', '150518_img24', '150519_img24', '150521_img24' '150706_img24'} ;
days = {'150716_img28', '150717_img28', '150719_img28'}; 


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
    
    pre_frames = 5;
    post_frames = 20;
    ts = (-pre_frames:post_frames)*1000/round(Sampeling_rate);
    tot_frame = pre_frames + post_frames+1;
    
    lickVec = ;
    
    use_ev_success = trial_outcome.success_time;
    time_before = 500; % in ms, time before event w/o release
    time_after =1000; % in ms, time after event w/o press
    success_roi = trigger_movie_by_event(lickVec, frame_info, ...
        use_ev_success, pre_frames, post_frames);
    
    %PLOT LICKING TRACES FOR SUCCESSFUL TRIALS HERE. 
    
    
    
    use_ev_fail = trial_outcome.early_time;
    fail_roi = trigger_movie_by_event(lickVec, frame_info, ...
        use_ev_fail, pre_frames, post_frames);
    %PLOT FAIL LICK TRACES HERE
    
    
    
end