%% SECTION ONE - assign pathnames and datasets to be analyzed/written. 
clear;
%NEED TO UPDATE THIS SO IT ACCESSES SPREADSHEET INSTEAD OF JUST WRITING IN THE NAMES
sessions = {'180414_img1008_1'}; 
days = '1008-180414_1';
bx_source = 'Z:\Data\Behv_MovingDots\behavior_raw';
image_source_base  = 'Z:\Data\WF imaging\'; %location of permanently stored image files for retreiving meta data
image_dest_base    = 'Z:\Analysis\WF_MovingDots_Analysis\BxAndAnalysisOutputs\'; %stores the data on crash in the movingDots analysis folder
behav_dest = ['Z:\Analysis\WF_MovingDots_Analysis\behavioral_analysis\' days];% behavior analysis results 

%% SECTION TWO - Uses a gui to allow user to draw ROIs 
for ii = 1:length(sessions)
    image_source = [image_source_base, sessions{ii}];
    image_dest = [image_dest_base sessions{ii} '\' sessions{ii}];
    if exist([image_dest_base sessions{ii}], 'file') ~= 7;
        mkdir(image_dest_base, sessions{ii});
    end
    WF_draw_ROIs_for_movingDots(sessions{ii}, image_source, image_dest); %automatically saves the ROIs
end

%% SECTION THREE - find/save frame times and meta data collected from the tiff files 
for ii = 1:length(sessions)
    image_source = [image_source_base, sessions{ii}];
    [all_files, meta_data, meta_data2] = obtain_tif_meta_data(image_source);
    frame_times = get_frame_time_by_movie_info(meta_data);
    dest =  [image_dest_base sessions{ii} '\' sessions{ii}];
    if exist([image_dest_base sessions{ii}], 'file') ~= 7;
        mkdir(image_dest_base, sessions{ii});
    end
    save([dest '_frame_times'],  'frame_times');
    save([dest '_meta_data'],  'meta_data', 'meta_data2', 'all_files');
end

%% SECTION FOUR - calculate TC of raw F for each ROI
for ii = 1:length(sessions)
    image_dest = [image_dest_base sessions{ii} '\' sessions{ii}];
    meta_data_dir = [image_dest_base sessions{ii} '\' sessions{ii} '_meta_data'];
    cluster_dir = [image_dest_base sessions{ii} '\' sessions{ii} '_cluster'];
    [avg_img, roi_sz,data_tc] = calculate_data_tc_jin(meta_data_dir, cluster_dir);
    rawF_dir = [image_dest, '_raw_F.mat'];
    save(rawF_dir, 'data_tc', 'roi_sz', 'avg_img');%automatically saves data_tc to bxOutputs
end

%% SECTION FIVE - calculate df/f
% if F is the mean of the stationary states
for ii = 1: length(sessions)
    image_dest = [image_dest_base sessions{ii} '\' sessions{ii}];
    raw_F = load([image_dest, '_raw_F.mat']);
    data_tc = raw_F.data_tc;
    %data_tc is the fluorescence data for this session, and it will be a
    %n*num of frames double, n=number of ROIs.
    behav_output = load([behav_dest '\' days '_behavAnalysis.mat']);
    stay = behav_output.frames_stay_cell;
    for n = 1:size(data_tc,1);
        data_tc_stay = data_tc(n,cell2mat(stay));
        F_staybase = mean(data_tc_stay);
        dfOvF_staybase(n,:) = (data_tc(n,:) - F_staybase)/F_staybase;
    end
   save([image_dest,'_dfOvF_staybase'],'dfOvF_staybase');
end

% if F is the mean of the whole session
%for ii = 1: length(sessions)
 %   image_dest = [image_dest_base sessions{ii} '\' sessions{ii}];
  %  raw_F = load([image_dest, '_raw_F.mat']);
   % data_tc = raw_F.data_tc;
    %data_tc is the fluorescence data for this session, and it will be a
    %n*num of frames double, n=number of ROIs.
    %for n = 1:size(data_tc,1)
     %   F_allbase = mean(data_tc(n,:));
      %  dfOvF_allbase = (data_tc - F_allbase)/F_allbase;
    %end
  % save([image_dest,'_dfOvF_allbase'],'dfOvF_allbase');
%end
