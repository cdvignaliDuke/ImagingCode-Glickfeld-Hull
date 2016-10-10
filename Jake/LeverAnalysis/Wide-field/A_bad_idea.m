%%
% clear
% tic
% day = '151009_img30';
% %day = '160606_img46';
% image_source = ['Z:\Data\WidefieldImaging\GCaMP\', day]; %location of permanently stored image files for retreiving meta data
% image_source_temp = ['C:\Users\jake\TempData\', day]; %looks at temporary storage to find the raw image data
% image_dest   = ['Z:\Analysis\WF Lever Analysis\BxAndAnalysisOutputs\' day]; %stores the data on crash in the lever analysis folder
% bx_source    =  'Z:\Data\WidefieldImaging\GCaMP\behavior\';
% if exist(image_dest, 'file') ~= 7;   %check to make sure that image_dest exists as a folder in the correct location. 
%     mkdir(image_dest);               %if it does not exist then create that folder. Should make this into a function and insert into the 2P code. 
% end
% old_cd = cd; %save old cd so I can restore it later
% [all_files, meta_data, meta_data2] = obtain_tif_meta_data(image_source);
% frame_times = get_frame_time_by_movie_info(meta_data);
% dest =  [image_dest '\' day '_ROI'];
% save([dest '_frame_times'],  'frame_times');  %ALTER PATH so it saves in WF Lever Analysis
% save([dest '_meta_data'],  'meta_data', 'meta_data2', 'all_files');
% shrink_movie(day, image_source_temp, image_dest, meta_data2); %automatically saves the shrunken movie
% toc
% clear
% %%
% clear
% tic
% day = '151011_img30';
% %day = '160606_img46';
% image_source = ['Z:\Data\WidefieldImaging\GCaMP\', day]; %location of permanently stored image files for retreiving meta data
% image_source_temp = ['C:\Users\jake\TempData\', day]; %looks at temporary storage to find the raw image data
% image_dest   = ['Z:\Analysis\WF Lever Analysis\BxAndAnalysisOutputs\' day]; %stores the data on crash in the lever analysis folder
% bx_source    =  'Z:\Data\WidefieldImaging\GCaMP\behavior\';
% if exist(image_dest, 'file') ~= 7;   %check to make sure that image_dest exists as a folder in the correct location. 
%     mkdir(image_dest);               %if it does not exist then create that folder. Should make this into a function and insert into the 2P code. 
% end
% old_cd = cd; %save old cd so I can restore it later
% [all_files, meta_data, meta_data2] = obtain_tif_meta_data(image_source);
% frame_times = get_frame_time_by_movie_info(meta_data);
% dest =  [image_dest '\' day '_ROI'];
% save([dest '_frame_times'],  'frame_times');  %ALTER PATH so it saves in WF Lever Analysis
% save([dest '_meta_data'],  'meta_data', 'meta_data2', 'all_files');
% shrink_movie(day, image_source_temp, image_dest, meta_data2); %automatically saves the shrunken movie
% toc
% clear
% %%
% clear
% tic
% day = '151021_img29';
% %day = '160606_img46';
% image_source = ['Z:\Data\WidefieldImaging\GCaMP\', day]; %location of permanently stored image files for retreiving meta data
% image_source_temp = ['C:\Users\jake\TempData\', day]; %looks at temporary storage to find the raw image data
% image_dest   = ['Z:\Analysis\WF Lever Analysis\BxAndAnalysisOutputs\' day]; %stores the data on crash in the lever analysis folder
% bx_source    =  'Z:\Data\WidefieldImaging\GCaMP\behavior\';
% if exist(image_dest, 'file') ~= 7;   %check to make sure that image_dest exists as a folder in the correct location. 
%     mkdir(image_dest);               %if it does not exist then create that folder. Should make this into a function and insert into the 2P code. 
% end
% old_cd = cd; %save old cd so I can restore it later
% [all_files, meta_data, meta_data2] = obtain_tif_meta_data(image_source);
% frame_times = get_frame_time_by_movie_info(meta_data);
% dest =  [image_dest '\' day '_ROI'];
% save([dest '_frame_times'],  'frame_times');  %ALTER PATH so it saves in WF Lever Analysis
% save([dest '_meta_data'],  'meta_data', 'meta_data2', 'all_files');
% shrink_movie(day, image_source_temp, image_dest, meta_data2); %automatically saves the shrunken movie
% toc
% clear
%%
clear
tic
day = '151022_img29';
%day = '160606_img46';
image_source = ['Z:\Data\WidefieldImaging\GCaMP\', day]; %location of permanently stored image files for retreiving meta data
image_source_temp = ['C:\Users\jake\TempData\', day]; %looks at temporary storage to find the raw image data
image_dest   = ['Z:\Analysis\WF Lever Analysis\BxAndAnalysisOutputs\' day]; %stores the data on crash in the lever analysis folder
bx_source    =  'Z:\Data\WidefieldImaging\GCaMP\behavior\';
if exist(image_dest, 'file') ~= 7;   %check to make sure that image_dest exists as a folder in the correct location. 
    mkdir(image_dest);               %if it does not exist then create that folder. Should make this into a function and insert into the 2P code. 
end
old_cd = cd; %save old cd so I can restore it later
[all_files, meta_data, meta_data2] = obtain_tif_meta_data(image_source);
frame_times = get_frame_time_by_movie_info(meta_data);
dest =  [image_dest '\' day '_ROI'];
save([dest '_frame_times'],  'frame_times');  %ALTER PATH so it saves in WF Lever Analysis
save([dest '_meta_data'],  'meta_data', 'meta_data2', 'all_files');
shrink_movie(day, image_source_temp, image_dest, meta_data2); %automatically saves the shrunken movie
toc
clear
%%
clear
tic
day = '151211_img32';
%day = '160606_img46';
image_source = ['Z:\Data\WidefieldImaging\GCaMP\', day]; %location of permanently stored image files for retreiving meta data
image_source_temp = ['C:\Users\jake\TempData\', day]; %looks at temporary storage to find the raw image data
image_dest   = ['Z:\Analysis\WF Lever Analysis\BxAndAnalysisOutputs\' day]; %stores the data on crash in the lever analysis folder
bx_source    =  'Z:\Data\WidefieldImaging\GCaMP\behavior\';
if exist(image_dest, 'file') ~= 7;   %check to make sure that image_dest exists as a folder in the correct location. 
    mkdir(image_dest);               %if it does not exist then create that folder. Should make this into a function and insert into the 2P code. 
end
old_cd = cd; %save old cd so I can restore it later
[all_files, meta_data, meta_data2] = obtain_tif_meta_data(image_source);
frame_times = get_frame_time_by_movie_info(meta_data);
dest =  [image_dest '\' day '_ROI'];
save([dest '_frame_times'],  'frame_times');  %ALTER PATH so it saves in WF Lever Analysis
save([dest '_meta_data'],  'meta_data', 'meta_data2', 'all_files');
shrink_movie(day, image_source_temp, image_dest, meta_data2); %automatically saves the shrunken movie
toc
clear
%%
clear
tic
day = '151212_img32';
%day = '160606_img46';
image_source = ['Z:\Data\WidefieldImaging\GCaMP\', day]; %location of permanently stored image files for retreiving meta data
image_source_temp = ['C:\Users\jake\TempData\', day]; %looks at temporary storage to find the raw image data
image_dest   = ['Z:\Analysis\WF Lever Analysis\BxAndAnalysisOutputs\' day]; %stores the data on crash in the lever analysis folder
bx_source    =  'Z:\Data\WidefieldImaging\GCaMP\behavior\';
if exist(image_dest, 'file') ~= 7;   %check to make sure that image_dest exists as a folder in the correct location. 
    mkdir(image_dest);               %if it does not exist then create that folder. Should make this into a function and insert into the 2P code. 
end
old_cd = cd; %save old cd so I can restore it later
[all_files, meta_data, meta_data2] = obtain_tif_meta_data(image_source);
frame_times = get_frame_time_by_movie_info(meta_data);
dest =  [image_dest '\' day '_ROI'];
save([dest '_frame_times'],  'frame_times');  %ALTER PATH so it saves in WF Lever Analysis
save([dest '_meta_data'],  'meta_data', 'meta_data2', 'all_files');
shrink_movie(day, image_source_temp, image_dest, meta_data2); %automatically saves the shrunken movie
toc
clear
%%
clear
tic
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
[all_files, meta_data, meta_data2] = obtain_tif_meta_data(image_source);
frame_times = get_frame_time_by_movie_info(meta_data);
dest =  [image_dest '\' day '_ROI'];
save([dest '_frame_times'],  'frame_times');  %ALTER PATH so it saves in WF Lever Analysis
save([dest '_meta_data'],  'meta_data', 'meta_data2', 'all_files');
shrink_movie(day, image_source_temp, image_dest, meta_data2); %automatically saves the shrunken movie
toc
clear
