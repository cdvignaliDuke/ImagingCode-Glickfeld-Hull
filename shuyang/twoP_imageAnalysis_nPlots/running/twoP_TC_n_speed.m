clear;
sessions = '190507_img1024'; 
image_analysis_base    = 'Z:\Analysis\2P_MovingDots_Analysis\imaging_analysis\'; %stores the data on crash in the movingDots analysis folder
image_analysis_dest = [image_analysis_base, sessions, '\'];

% behavior analysis results 
days = '1024-190507_1';
behav_dest = ['Z:\Analysis\2P_MovingDots_Analysis\behavioral_analysis\' days '\'];
color_code = {'b','r','k','c'};

%%

image_analysis_dest = [image_analysis_base, sessions, '\'];
% behavior analysis results
%behav_dest = ['Z:\Analysis\Airpuff_analysis\behavioral_analysis\' days{i} '\'];
behav_dest = ['Z:\Analysis\2P_MovingDots_Analysis\behavioral_analysis\' days '\'];
%load data
filename = dir([image_analysis_dest 'getTC\' '*' '_TCave.mat']);
TCave = load([image_analysis_dest 'getTC\' filename.name]);
tc_avg = TCave.tc_avg;
behav_output = load([behav_dest days '_behavAnalysis.mat']);
speed = double(behav_output.speed);

plotTC_Jin(tc_avg(25000:28992,:),1,1,[84,82,81,72,69,66,57,51,47,36:38,34,33],1);
figure;plot(speed(25000:end)*2*3.1415926*7.5/128);
