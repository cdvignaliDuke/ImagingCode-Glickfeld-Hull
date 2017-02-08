%% SECTION ONE - assign pathnames and datasets to be analyzed/written. 
clear;
%these are the days which I am including in my summary statistics so far. Need to update with newer datasets. 
%all datasets which were included in the scatterplot as of 10/13/16
days = {'151021_img29', '151022_img29', '151009_img30', '151011_img30', '151211_img32', '151212_img32', '160129_img35', '160131_img35', '160129_img36', '160131_img36', '160314_img38', '160315_img38', '160319_img41', '160320_img41', '160606_img46', '160722_img53', '160904_img55'}; %'150718_img27', '150719_img27', '150716_img28', '150717_img28', 

%naive mice 

%sessions to be analyzed
%days = {'161216_img74', '161216_img73', '161216_img65', '161215_img74', '161215_img73',};  

bx_source     = ['Z:\Data\WidefieldImaging\GCaMP\behavior\'];
image_source_base  = ['Z:\Data\WidefieldImaging\GCaMP\']; %location of permanently stored image files for retreiving meta data
image_dest_base    = ['Z:\Analysis\WF Lever Analysis\BxAndAnalysisOutputs\']; %stores the data on crash in the lever analysis folder
bx_outputs_dir = ['Z:\Analysis\WF Lever Analysis\BxAndAnalysisOutputs\BxOutputs\'];
PCA_output_dir_base = ['Z:\Analysis\WF Lever Analysis\PCA_output_dir\'];
old_cd = cd; %save old cd so I can restore it later

%% SECTION TWO - Uses a gui to allow user to draw ROIs 
for ii = 1:length(days)
    days{ii}
    image_source = [image_source_base, days{ii}];
    image_dest = [image_dest_base days{ii} '\' days{ii}];
    if exist([image_dest_base days{ii}], 'file') ~= 7;
        mkdir(image_dest_base, days{ii});
    end
    WF_draw_ROIs_for_lever(days{ii}, image_source, image_dest); %automatically saves the ROIs
end

%% SECTION TWO B - Use PCA to generate functional ROIs 
for ii = 16 %:length(days)  %
    days{ii}
    image_source = [image_source_base, days{ii} '\' days{ii} '_MMStack.ome.tif'];
    meta_data_dir = [image_dest_base days{ii} '\' days{ii} '_meta_data'];
    PCA_output_dir = PCA_output_dir_base;
    load(meta_data_dir, 'meta_data2');
    movie_downsamp = [];
    movie_avg = [];
    for iii = 1:length(meta_data2)-1
        for iv = 1:5:4000;
            img_subset = imread(meta_data2{iii}(1).Filename,iv);
            movie_downsamp = cat(3,movie_downsamp,img_subset);
        end
        iii
    end
    movie_avg = mean(movie_downsamp,3);
    writetiff(movie_downsamp, [PCA_output_dir, days{ii}, '_movies']);
    num_frames = size(movie_downsamp,3);
    clear movie_downsamp;
    
    %crop tif in imageJ mannually
    
    movie_downsamp = readtiff([PCA_output_dir, days{ii}, '_movies.tif']);
    movie_avg = mean(movie_downsamp,3); 
    
    %manually generate a mask in imageJ
    mask_tif = readtiff([PCA_output_dir, days{ii}, '_movies_mask.tif']);
    mask_tif = reshape(mask_tif, 1, size(mask_tif,1)*size(mask_tif,2));
    mask_zeros = find(mask_tif==0);
    movie_downsamp = reshape(movie_downsamp, size(movie_downsamp,1)*size(movie_downsamp,2), num_frames);
    movie_downsamp(mask_zeros,:) = 0;
    movie_downsamp = reshape(movie_downsamp, size(movie_avg,1), size(movie_avg,2), num_frames);
    movie_avg = mean(movie_downsamp,3); 
    
    global rt
    rt = movie_downsamp;
    [mixedsig, mixedfilters, CovEvals, covtrace, movm, movtm] = CellsortPCA2([PCA_output_dir, days{ii}, '_movies.tif'], [], 100,1, [PCA_output_dir], []);
    [PCuse] = CellsortChoosePCs(movie_avg, mixedfilters);
end
%% SECTION THREE - find/save frame times and meta data collected from the tiff files 
for ii = 1:length(days)
    image_source = [image_source_base, days{ii}];
    [all_files, meta_data, meta_data2] = obtain_tif_meta_data(image_source);
    frame_times = get_frame_time_by_movie_info(meta_data);
    dest =  [image_dest_base days{ii} '\' days{ii}];
    if exist([image_dest_base days{ii}], 'file') ~= 7;
        mkdir(image_dest_base, days{ii});
    end
    save([dest '_frame_times'],  'frame_times');
    save([dest '_meta_data'],  'meta_data', 'meta_data2', 'all_files');
end

%% SECTION FOUR - calculate TC of raw F for each ROI
for ii = 1:length(days)
    image_dest = [image_dest_base days{ii} '\' days{ii}];
    meta_data_dir = [image_dest_base days{ii} '\' days{ii} '_meta_data'];
    cluster_dir = [image_dest_base days{ii} '\' days{ii} '_cluster'];
    bx_out_dir  = [bx_outputs_dir days{ii} '_bx_outputs'];
    calculate_data_tc(meta_data_dir, cluster_dir, bx_out_dir); %automatically saves data_tc to bxOutputs
end

%% SECTION FIVE - obtain bx and frame info from MWorks.
for ii = 1:length(days)
    days(ii)
    image_dest = [image_dest_base days{ii} '\' days{ii}];
    meta_data_dir = [image_dest '_meta_data'];
    frame_times_dir = [image_dest  '_frame_times'];
    bx_out_dir  = [bx_outputs_dir days{ii} '_bx_outputs'];
    b_data = get_bx_data(bx_source, days{ii});  %find the correct behavior file and loads it.
    [lever, frame_info, trial_outcome, lickTimes] = parse_behavior_for_HAD(b_data, meta_data_dir, frame_times_dir);
    save(bx_out_dir, 'lever', 'frame_info', 'trial_outcome', 'lickTimes', '-append');
end

%% SECTION SIX - calculate df/f TC. baseline f redefined on a trial by trial basis
for ii = 1:length(days)
    bx_out_dir  = [bx_outputs_dir days{ii} '_bx_outputs'];
    b_data = get_bx_data(bx_source, days{ii});  %find the correct behavior file and loads it.
    tc_dfoverf = tbyt_dfoverf(b_data, bx_out_dir);
    save(bx_out_dir, 'tc_dfoverf', '-append');
end

%% PLOT TCs save TCs for summary plots
 WF_lever_plotting_TCs;
% struct_maker;
% plot_variables_from_structs;

%% SECTION FIVE - population analyses to create scatter plots
% WF_summary_scatters





