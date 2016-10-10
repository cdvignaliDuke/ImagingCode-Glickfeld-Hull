%% SECTION ONE  
% 1)assign pathnames and datasets to be analyzed/written. 
clear;
day = '160129_img35';
%day = '160606_img46';
image_source = ['Z:\Data\WidefieldImaging\GCaMP\', day]; %location of permanently stored image files for retreiving meta data
image_source_temp = ['C:\Users\jake\TempData\', day]; %looks at temporary storage to find the raw image data
image_dest   = ['Z:\Analysis\WF Lever Analysis\BxAndAnalysisOutputs\' day]; %stores the data on crash in the lever analysis folder
bx_source    =  'Z:\Data\WidefieldImaging\GCaMP\behavior\';
if exist(image_dest, 'file') ~= 7;   %check to make sure that image_dest exists as a folder in the correct location. 
    mkdir(image_dest);               %if it does not exist then create that folder. Should make this into a function and insert into the 2P code. 
end
old_cd = cd; %save old cd so I can restore it later

%% SECTION TWO
% 1) Collects the frame times from the tif meta data. 
% 2) Shrinks the movie, reducing resolution so it uses less memory.
%find/save frame times and meta data collected from the tiff files 
[all_files, meta_data, meta_data2] = obtain_tif_meta_data(image_source);
frame_times = get_frame_time_by_movie_info(meta_data);
dest =  [image_dest '\' day '_ROI'];
save([dest '_frame_times'],  'frame_times');  %ALTER PATH so it saves in WF Lever Analysis
save([dest '_meta_data'],  'meta_data', 'meta_data2', 'all_files');

%shrink tiff files and save to image_dest
shrink_movie(day, image_source_temp, image_dest, meta_data2); %automatically saves the shrunken movie

%% SECTION THREE - Uses a gui to allow user to draw ROIs 
WF_draw_ROIs_for_lever(day, image_dest)  %automatically saves the ROIs

%% SECTION FOUR - analyze bx data, find baseline times, make df/f, align frame counter/times, plot TCs
WF_lever_Bx_DFoF
% - load paths, frame_info, bx file
% - calculate info, ifi, sampling_rate, f_frame, l_frame
% - parse_behavior_for_HAD
    % allocates memory for lever, frame, and trial_outcome
    % goes through each trial. concatenates lever state, press/release times, calculates frame counter times and counterbytimes  (needs counter fixer)
    % stores timestamps needed to trim non-imaged baseline times
    % concatenates licktimes. Aligns licktimes to first frame.
    % aligns frame.times to beginning of exp.  Generates counter_by_frame_times
    % collects hold/trial info from bx data. calculates trial outcome categories
    % attemps to correct for spurious counters before camera starts and on camera activation
    % removes all events not between first and last frames
    % determines baseline times. Trims non-imaged baseline times
    % stores trial nums of first and last imaged trials 
    % aligns baseline times to frame.counter
    % removes nonimaged trials from lever hold info
    % stores in three structures
% - get_movie_by_ROI
% - uses baseline times to calculate DFoF TC
% - save outputs

WF_lever_plotting_TCs
% - pathnames. days. 
% - loads bx file and bx_outputs(event times)
% - plots ROIs on a heatmap
% - Bins licking times according to frame # (should increase the resolution for plotting)
% - uses trigger_movie_by_event_licks to get the individual calcium traces for each trial outcome condition
% - plots average TCs for each trial outcome condition
% - standardizes the axes
% - reports # of trials in each category 
% - save figures 
% - plot cue triggered TCs 
% - save indiv TC data
% - commented out - plot lever press aligned TC
% - plot corr coef of licking and df/f
% - saved lapse trials 



%% SECTION FIVE - population analyses to create scatter plots
WF_summary_scatters







