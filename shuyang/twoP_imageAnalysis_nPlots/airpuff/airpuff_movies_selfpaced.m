% make the average running onset movies for 2P during motorized experiment.
%{
clear;
sessions = {'191114_img1040','191115_img1039','191115_img1041','191115_img1042','200316_img1064_airpuff_2'};
image_analysis_base  = 'Z:\Analysis\Airpuff_analysis\imaging_analysis\';%stores the data on crash in the movingDots analysis folder
% behavior analysis results
days = {'1040-191114_1','1039-191115_1','1041-191115_1','1042-191115_1','1064-200316_2'};
%}
%% 2 photon: identify the frames you wanna average, load sbxfile

clear;
%NEED TO UPDATE THIS SO IT ACCESSES SPREADSHEET INSTEAD OF JUST WRITING IN THE NAMES
sessions = '191115_img1042'; 
image_analysis_base  = 'Z:\Analysis\Airpuff_analysis\imaging_analysis\';
%image_source = [image_source_base, sessions,'\',ID,'\'];
image_analysis_dest = [image_analysis_base, sessions, '\' 'getTC\'];
days = '1042-191115_1';
behav_dest = ['Z:\Analysis\Airpuff_analysis\behavioral_analysis\' days];
file = [image_analysis_dest, sessions '_000_img_rgs_ref18.mat' ];
img_rgs = load(file,'img_rgs');
img_rgs = img_rgs.img_rgs;

behav_output = load([behav_dest '\' days '_behavAnalysis.mat']);
frms_airstim_stay = behav_output.frms_airstim_stay; %trials*frames

% img = 264*796*29999 (widthxlengthxframes)
% want: width x length x trials x frames
img_airstim_stay = img_rgs(:,:,frms_airstim_stay(:));
img_airstim_stay = reshape(img_airstim_stay,size(img_airstim_stay,1),size(img_airstim_stay,2),size(frms_airstim_stay,1),size(frms_airstim_stay,2));
ave_img_airstimstay = squeeze(mean(img_airstim_stay,3));
size(ave_img_airstimstay)
writetiff(ave_img_airstimstay,[image_analysis_base, '\' sessions '\' sessions '_aveAirstimStay']);

