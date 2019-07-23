% save facial data as tiff

%% SECTION ONE - assign pathnames and datasets to be analyzed/written. 
clear;
%NEED TO UPDATE THIS SO IT ACCESSES SPREADSHEET INSTEAD OF JUST WRITING IN THE NAMES
sessions = '190520_img1021'; 
%ID = '1016';
image_source_base  = 'Z:\Data\2photon\'; %location of permanently stored image files for retreiving meta data
%image_analysis_base    = 'Z:\Analysis\2photon_test\'; %stores the data on crash in the movingDots analysis folder
image_analysis_base    = 'Z:\Analysis\2P_MovingDots_Analysis\imaging_analysis\'; %stores the data on crash in the movingDots analysis folder
%image_source = [image_source_base, sessions,'\',ID,'\'];
image_source = [image_source_base, sessions,'\'];
image_analysis_dest = [image_analysis_base, sessions, '\' 'face\'];

%% 
filename = dir([image_source '*' '_ball.mat']);
facedata = load([image_source filename.name]);
facedata = squeeze(facedata.data);
f1 = 1;
f2 = 3000;
saveastiff(facedata(:,:,f1:f2),[image_analysis_dest sessions '_grooming' num2str(f1) '_' num2str(f2)]);

