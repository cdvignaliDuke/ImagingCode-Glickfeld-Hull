% make the average running onset movies for 2P during motorized experiment.


%% 2 photon: identify the frames you wanna average, load sbxfile

clear;
%NEED TO UPDATE THIS SO IT ACCESSES SPREADSHEET INSTEAD OF JUST WRITING IN THE NAMES
sessions = '200319_img1064'; 
image_analysis_base = 'Z:\Analysis\motorizedWheel_Analysis\running\imaging_analysis\'; 
%image_source = [image_source_base, sessions,'\',ID,'\'];
image_analysis_dest = [image_analysis_base, sessions, '\' 'getTC\'];
days = '1064-200319_1';
behav_dest = ['Z:\Analysis\motorizedWheel_Analysis\running\behavioral_analysis\' days '\'];
file = [image_analysis_dest, sessions '_000_img_rgs_ref45.mat' ];
img_rgs = load(file,'img_rgs');
img_rgs = img_rgs.img_rgs;

behav_output = load([behav_dest days '_behavAnalysis.mat']);
runtrig_fast_mat = behav_output.runtrig_fast_mat; % trial * frame
runtrig_slow_mat = behav_output.runtrig_slow_mat;
runoff_slow_mat = behav_output.runoff_slow_mat;
runoff_fast_mat = behav_output.runoff_fast_mat;

% img = 264*796*29999 (widthxlengthxframes)
% want: width x length x trials x frames
img_runonset_fast = img_rgs(:,:,runtrig_fast_mat(:));
img_runonset_fast = reshape(img_runonset_fast,size(img_runonset_fast,1),size(img_runonset_fast,2),size(runtrig_fast_mat,1),size(runtrig_fast_mat,2));
ave_img_fastrunOnset = squeeze(mean(img_runonset_fast,3));
size(ave_img_fastrunOnset)
writetiff(ave_img_fastrunOnset,[image_analysis_base, '\' sessions '\' sessions '_aveRunOnsetFast']);

img_runoffset_fast = img_rgs(:,:,runoff_fast_mat(:));
img_runoffset_fast = reshape(img_runoffset_fast,size(img_runoffset_fast,1),size(img_runoffset_fast,2),size(runoff_fast_mat,1),size(runoff_fast_mat,2));
ave_img_fastrunOffset = squeeze(mean(img_runoffset_fast,3));
size(ave_img_fastrunOffset)
writetiff(ave_img_fastrunOffset,[image_analysis_base, '\' sessions '\' sessions '_aveRunOffsetFast']);

img_runonset_slow = img_rgs(:,:,runtrig_slow_mat(:));
img_runonset_slow = reshape(img_runonset_slow,size(img_runonset_slow,1),size(img_runonset_slow,2),size(runtrig_slow_mat,1),size(runtrig_slow_mat,2));
ave_img_slowrunOnset = squeeze(mean(img_runonset_slow,3));
size(ave_img_slowrunOnset)
writetiff(ave_img_slowrunOnset,[image_analysis_base, '\' sessions '\' sessions '_aveRunOnsetSlow']);

img_runoffset_slow = img_rgs(:,:,runoff_slow_mat(:));
img_runoffset_slow = reshape(img_runoffset_slow,size(img_runoffset_slow,1),size(img_runoffset_slow,2),size(runoff_slow_mat,1),size(runoff_slow_mat,2));
ave_img_slowrunOffset = squeeze(mean(img_runoffset_slow,3));
size(ave_img_slowrunOffset)
writetiff(ave_img_slowrunOffset,[image_analysis_base, '\' sessions '\' sessions '_aveRunOffsetSlow']);

