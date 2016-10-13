%WF_lever_Bx_DFoF -- This script replaced spatial_cluster_ROI. Must run ShrinkMovie first
%Requires _ROIshrink.tif file and _ROI_frame_times.mat file

clear
BIN_SIZE = 1;   % if the ROI is large use large bin size (10)
holdT_min = 400000;

image_source = ['Z:\Data\WidefieldImaging\GCaMP\', day]; %location of permanently stored image files for retreiving meta data
image_source_temp = ['C:\Users\jake\TempData\', day]; %looks at temporary storage to find the raw image data
frame_and_cluster_dest   = ['Z:\tempHoldingArea\full_size_analysis\' day]; %stores the data on crash in the lever analysis folder

DATA_DIR = 'Z:\Analysis\WF Lever Analysis\BxAndAnalysisOutputs\';
BEHAVE_DIR = 'Z:\Data\WidefieldImaging\GCaMP\behavior\';
ANALYSIS_DIR = 'Z:\Analysis\LeverAnalysis\';

for kk=1:length(days)
    frame_info_dest  = [frame_and_cluster_dest '_ROI_frame_times.mat'];
    load([frame_and_cluster_dest '_meta_data.mat']);
    image_dest  = [DATA_DIR days{kk} '\' days{kk} '_ROIshrink.tif'];
    
    bfile = dir([BEHAVE_DIR 'data-*i9' days{kk}(end-1:end) '-' days{kk}(1:6) '*' ]);
    behave_dest = [BEHAVE_DIR bfile.name];
    assert(length(bfile)) = 1;
    b_data = load(behave_dest); 
    b_data = b_data.input;
    ftimes = load(frame_info_dest);
    info = imfinfo(image_dest);
    
    if(~exist('first_frame', 'var'))  %need to redo how I select these variables 
        f_frame =1;
    else
        f_frame = first_frame;
    end
    if(~exist('last_frame', 'var'))
        l_frame = min(length(ftimes.frame_times), length(info));
    else
        l_frame = last_frame;
    end
    
    %obtain bx and frame info from MWorks.  Trim movie by ROI (shurnken and rectangular one)
    [lever, frame_info, trial_outcome, lickTimes] = parse_behavior_for_HAD(b_data, ...
        f_frame, l_frame, ftimes.frame_times, holdT_min);             % ftimes is not used again after parse_bx
    %load tiff file
    img = [];
    for ii = 1:length(meta_data2)
        for iii = 1:length(meta_data2{ii});
            this_frame = imread(meta_data2{ii}.Filename, iii);
            img = cat(3, img, this_frame);
        end
    end
    
    clear cluster;    %clear existing clusters 
    cluster_file = [image_dest(1:end-4) 'cluster.mat'];
    load(cluster_file); %load cluster
    roi_sz = sum(cluster.roi_mask,2); %log the #of pixels in each ROI (binned pixles)
    
    % ----- cluster data by ROI
    disp('extracting ROIs separately');
    data_tc = [];
    for i = 1:cluster.num_cluster;
        data_tc(i,:) = stackGetTimeCourses(reshape(img,sz(1),sz(2),size(img,2)), reshape(cluster.roi_mask(i,:),[sz(1) sz(2)]));
    end
    
    %Obtain a dFoverF TC on a trial by trial basis
    tc_dfoverf = tbyt_dfoverf(lever, b_data, data_tc, frame_info); 
    
    %calculate sampling rate
    frame_times = frame_info.times;
    sampling_rate = 1000/mode(diff(frame_times));
    
    destBxOutput = strcat(DATA_DIR, 'BxOutputs\', days{kk}, '_bx_outputs');
    save([destBxOutput], 'trial_outcome', 'lever', 'frame_info', 'sampling_rate','sz', 'avg_img', 'cluster', 'lickTimes', 'tc_dfoverf'); %  
end