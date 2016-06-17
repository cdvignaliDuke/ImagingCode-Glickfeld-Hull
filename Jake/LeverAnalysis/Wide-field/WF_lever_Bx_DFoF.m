%CLUSTER ROI
%this script replaced spatial_cluster_ROI    
%Must run ShrinkMovie first
%requires _ROIshrink.tif file and _ROI_frame_times.mat file
%can plot TCs of individual ROIs or group them as a single mask
%plots correlation coefficient between the ROIs
%plots TCs of the df/f for success, earlies, and press. Plot window defined
%by time_before and time_after. df/f redfined periodically anytime there is
%a period of no lever activity within the iti of a given length. 
clear
WRITE_VEDIO = 0;
BIN_SIZE = 1;   % if the ROI is large use large bin size (10)

DATA_DIR =  'C:\Users\jake\TempData\';
FRAME_TIME_DIR = 'C:\Users\jake\TempData\';
BEHAVE_DIR = 'Z:\Data\WidefieldImaging\GCaMP\behavior\';
ANALYSIS_DIR ='Z:\Analysis\LeverAnalysis\';
% -----------
%days = {'160208_img35', '160209_img36', '151222_img32', '151019_img30'}; %NO-LEVER DATASETS %   150725_img27   150707_img25  150707_img24
%days = {'160208_img36', '160207_img35'}; %151221_img32,  %NO Vis Stim DATASETS 
%days = {'160207_img36', '160205_img35'}; %NO Aud Stim DATASETS

%days = {'150514_img24', '150518_img24', '150519_img24', '150521_img24' '150706_img24'} ;
%days = {'150517_img25', '150518_img25', '150514_img25', '150515_img25'};
%days = {'150716_img27', '150718_img27', '150719_img27'};  %'150717_img27',
%days = {'150716_img28', '150717_img28', '150719_img28'}; 
%days = {'151021_img29', '151022_img29'}; %'151015_img29'
%days = {'151009_img30', '151011_img30'};
%days = {'160129_img35', '160131_img35', '160129_img36','160131_img36'};

%days = {'160131_img36', '160131_img35', '151212_img32', '160315_img38', '160320_img41'};  %rand = 1000
%days = {'160129_img36', '160129_img35', '151009_img30', '151011_img30', '151211_img32', '160314_img38', '160319_img41'};  %rand=4500

%days = {'150718_img27', '150719_img27', '150716_img28', '150717_img28', '151021_img29', '151022_img29', '151009_img30', '151011_img30', '151211_img32', '151212_img32', '160129_img35', '160131_img35', '160129_img36','160131_img36', '160319_img41', '160320_img41', '150518_img24', '150519_img24', '150518_img25', '150517_img25'};
%days = {'150518_img25', '150517_img25', '150716_img27', '150718_img27', '150716_img28', '150717_img28', '151021_img29', '151022_img29', '151009_img30', '151011_img30', '151211_img32', '151212_img32', '160129_img35', '160131_img35', '160129_img36','160131_img36', '160314_img38', '160315_img38', '160319_img41', '160320_img41'}; %'150718_img27', '150719_img27',
%'150518_img24', '150519_img24',
days = {'160319_img41'};

%days = {'160129_img36', '160129_img35', '160131_img35', '151211_img32', '150717_img28', '150716_img28', '150718_img27', '151022_img29', '150719_img27', '160319_img41', '160320_img41'}; %days with usable licking data

success_trials = [];
fail_trials = [];
tooFast_trials = [];

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
    sampling_rate = 1000/ifi;
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
    disp('extracting ROIs separately');
    data_tc = [];
    for i = 1:cluster.num_cluster;
        data_tc(i,:) = stackGetTimeCourses(reshape(img,sz(1),sz(2),size(img,2)), reshape(cluster.roi_mask(i,:),[sz(1) sz(2)]));
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
    %tc_dfoverf_hold = nan(size(data_tc));
    F_range = [];  
    %F_range_hold = [];  
    for iT=frame_info.f_frame_trial_num+1: frame_info.l_frame_trial_num-1;    %only looks at the first and last fully imaged trials
        %F_range is the # of each frame which we will use to generate f
        %iti based F
        if ~isnan(lever.baseline_timesMs(1,iT));   %if there is a valid baseline interval then make that the new F
            F_range = frame_info.counter(lever.baseline_timesMs(1,iT)):frame_info.counter(lever.baseline_timesMs(2,iT));
        elseif isempty(F_range)   %if these are trials before there was ever a valid F_range then use the first valid F_range as this trials F_range also. 
            F_range = frame_info.counter(lever.baseline_timesMs(1,first_baseline)):frame_info.counter(lever.baseline_timesMs(2,first_baseline));
        end    %if there was no valid F_range but there was previously a valid F_range then F_range will remain unchanged and the most recent one will be use. 
%         %hold based F
%         if ~isnan(lever.baseline_times_holdMs(1,iT));   %same concept but for the hold based F
%             F_range_hold = frame_info.counter(lever.baseline_times_holdMs(1,iT)):frame_info.counter(lever.baseline_times_holdMs(2,iT));
%         elseif isempty(F_range_hold)    
%             F_range_hold = frame_info.counter(lever.baseline_times_holdMs(1,first_baseline_hold)):frame_info.counter(lever.baseline_times_holdMs(2,first_baseline_hold));
%         end
        %assign the frame numbers which correspond to this trial
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
%         for i = 1:size(data_tc,1); %find the avf_f for F_range and apply it to all the frames in that trial
%             F_avg_hold= mean(data_tc(i,F_range_hold));
%             t_df = bsxfun(@minus, data_tc(i, t_range), F_avg_hold);   %use bsxfun because sometimes F_avg is a vector
%             t_dfoverf = bsxfun(@rdivide, t_df, F_avg_hold);
%             tc_dfoverf_hold(i,t_range) = t_dfoverf;
%         end
    end  
    %disp([mat2str(sum(isnan(tc_dfoverf))) ' NaNs remaining in tc_dfoverf']);
    
    destBxOutput = strcat(ANALYSIS_DIR, 'BxAndAnalysisOutputs\BxOutputs\', days{kk}, '_bx_outputs');
    save([destBxOutput], 'trial_outcome', 'lever', 'frame_info', 'avg_img', 'sz', 'cluster', 'sampling_rate', 'lickTimes', 'tc_dfoverf');
end