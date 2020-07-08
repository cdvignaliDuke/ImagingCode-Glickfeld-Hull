% make the average running onset movies for both 2P and wide field.
% the output movie is the average across trials. 

%% Wide Field
clear;
sessions = '190617_img1021_1'; 
days = '1021-190617_1';
image_dest_base    = 'Z:\Analysis\WF_MovingDots_Analysis\BxAndAnalysisOutputs\'; %stores the data on crash in the movingDots analysis folder
dest =  [image_dest_base sessions '\' sessions];
behav_dest = ['Z:\Analysis\WF_MovingDots_Analysis\behavioral_analysis\' days '\'];

% read all of the tiff files (check if there's any variable saved after loading) % do the average % write into tiff

meta_data_dir = [image_dest_base sessions '\' sessions '_meta_data'];
load(meta_data_dir, 'meta_data2');

%load tiff files
img = [];
for ii = 1:length(meta_data2)
    img_subset = readtiff(meta_data2{ii}(1).Filename);
    img = cat(3,img,img_subset);
    disp(['file #', num2str(ii), ' download complete']);  
end
save([dest '_img'],  'img', '-v7.3');
%load([dest '_img'],'img');

behav_output = load([behav_dest days '_behavAnalysis.mat']);
frames_run_cell = behav_output.frames_run_cell;
speed = double(behav_output.speed);
period = 3; % # of frames before and after running
befoRunStay = 10; %1s, # of frames that speed =0 before running
totalT = 20; %2s, # of frmaes before running onset + # of frames following running onset/ # of frames before running offset + # of frames after running offset
aftRunOff = 10; %1s,# of frames that speed = 0 after running
[~,~,frms_runTrig_mat,~,~]= findFrames_runWindows_2P(speed,frames_run_cell,period,befoRunStay,totalT,aftRunOff);
% [~,~,frms_runTrig_NC_mat,frms_runoff_NC_mat,~]= findFrames_run_no_criteria_2P(speed,...
%     frames_run_cell,period,befoRunStay,totalT,aftRunOff);

% img = 1002x1004x17900 (widthxlengthxframes)
% frms_runTrig_NC_mat = 20x40 (framesxtrials)
% want: width x length x frames x trials

img_run = img(:,:,frms_runTrig_mat(:));
img_run = reshape(img_run, size(img_run,1), size(img_run,2), size(frms_runTrig_mat,1), size(frms_runTrig_mat,2));

% img_run = img(:,:,frms_runTrig_NC_mat(:));
% img_run = reshape(img_run, size(img_run,1), size(img_run,2), size(frms_runTrig_NC_mat,1), size(frms_runTrig_NC_mat,2));

size(img_run)
ave_img_run = squeeze(mean(img_run,4));
writetiff(ave_img_run,[dest '_runtrigave_stayzero']);

%% 2 photon: identify the frames you wanna average, load sbxfile

clear;
%NEED TO UPDATE THIS SO IT ACCESSES SPREADSHEET INSTEAD OF JUST WRITING IN THE NAMES
sessions = '200316_img1064'; 
image_analysis_base    = 'Z:\Analysis\2P_MovingDots_Analysis\imaging_analysis\'; %stores the data on crash in the movingDots analysis folder
%image_source = [image_source_base, sessions,'\',ID,'\'];
image_analysis_dest = [image_analysis_base, sessions, '\' 'getTC\'];
days = '089-190602_1';
behav_dest = ['Z:\Analysis\2P_MovingDots_Analysis\behavioral_analysis\' days '\'];
file = [image_analysis_dest, sessions '_000_img_rgs_ref56.mat' ];
img_rgs = load(file,'img_rgs');
img_rgs = img_rgs.img_rgs;

behav_output = load([behav_dest days '_behavAnalysis.mat']);
frames_run_cell = behav_output.frames_run_cell;
speed = double(behav_output.speed);
period = 9; % # of frames before and after running
befoRunStay = 30; %1s, # of frames that speed =0 before running
totalT = 60; %2s, # of frmaes before running onset + # of frames following running onset/ # of frames before running offset + # of frames after running offset
aftRunOff = 30; %1s,# of frames that speed = 0 after running
[~,~,frms_runTrig_NC_mat,frms_runoff_NC_mat,~]= findFrames_run_no_criteria_2P(speed,...
    frames_run_cell,period,befoRunStay,totalT,aftRunOff);

% img = 1002x1004x17900 (widthxlengthxframes)
% frms_runTrig_NC_mat = 20x40 (framesxtrials)
% want: width x length x frames x trials

img_run_2P = img_rgs(:,:,frms_runTrig_NC_mat(:));
img_run_2P = reshape(img_run_2P, size(img_run_2P,1), size(img_run_2P,2), size(frms_runTrig_NC_mat,1), size(frms_runTrig_NC_mat,2));

size(img_run_2P)
ave_img_run_2P = squeeze(mean(img_run_2P,4));
writetiff(ave_img_run_2P,[image_analysis_base, '\' sessions '\' sessions '_runtrigave']);

