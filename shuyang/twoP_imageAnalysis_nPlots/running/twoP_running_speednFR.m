clear;
sessions = '190507_img1024'; 
days = '1024-190507_1';
image_analysis_base    = 'Z:\Analysis\2P_MovingDots_Analysis\imaging_analysis\'; %stores the data on crash in the movingDots analysis folder
color_code = {'b','r','k','c'};

%%
%behav_dest = ['Z:\Analysis\Airpuff_analysis\behavioral_analysis\' days{ii}];
behav_dest = ['Z:\Analysis\2P_MovingDots_Analysis\behavioral_analysis\' days '\'];
behav_struct = load([behav_dest '\' days '_behavAnalysis.mat']);
speed = behav_struct.speed; speed = double(speed);
image_analysis_dest = [image_analysis_base, sessions, '\'];
threshold = -4;
spk_deconv_output = load([image_analysis_dest sessions,'_spk_deconvolve_threshold' num2str(threshold) '.mat']);
spk_logic = spk_deconv_output.spk_logic_cl;
spk_logic_ave = mean(spk_logic,2);
FR_ave = spk_logic_ave * 30; %firing rate = spike probablity*imaging rate
% generate matrix for spikes

frame = (0:(length(speed)-1));
% convert frames to seconds
x = frame/30;
figure; n1 = 25000; n2 = 28992;
%convert speed to cm/s: each unit = 2pai*r/128 (there're 128 pulses in total for each circle)
speedcm = speed*2*3.1415926*7.5/128;

[hAx,hline1,hline2] = plotyy(x(n1:n2),speedcm(n1:n2),x(n1:n2),FR_ave(n1:n2)');
%[hAx,hline1,hline2] = plotyy(frame,speed,frame,(avedfOvF(1:length(speed))'));

%set(hline1,'color', 'b'); 
%set(hline2,'color',color_code{n});
legend('speed','firing rate','Location','northwest');
ylim(hAx(1),[-15 25]);
ylim(hAx(2),[-1,15]);
xlabel('time (s)');
title(['speed and average FR' days]);
savefig([image_analysis_dest sessions '_aveFR_wspeed']);




