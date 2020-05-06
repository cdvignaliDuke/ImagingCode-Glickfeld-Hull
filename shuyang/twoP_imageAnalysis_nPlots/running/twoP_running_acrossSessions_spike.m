% run triggered average for across all sessions using only trials that mice
% are totally still before/after running

%% assign document paths and experimental sessions
clear;
sessions = {'190429_img1021','190430_img1023','190507_img1024','190603_img1025'};
image_analysis_base    = 'Z:\Analysis\2P_MovingDots_Analysis\imaging_analysis\'; %stores the data on crash in the movingDots analysis folder
% behavior analysis results 
days = {'1021-190429_1','1023-190430_1','1024-190507_1','1025-190603_1'};
color_code = {'b','r','k','c'};

%% 
ave_spike_runtrig_allsession = [];% frames*total number of cells in all sessions
ave_spike_runoff_allsession = [];
aveSpd_runtrig_100ms_all_session = [];% frames*total number of sessions that have runtrig
aveSpd_runoff_100ms_all_session = [];
nPCs_total = 0;
ntrials_runtrig_total = 0;
ntrials_runoff_total = 0;
for ii = 1:length(sessions) 
    %load things
    image_analysis_dest = [image_analysis_base, sessions{ii}, '\'];
    threshold = -4;
    spike_output = load([image_analysis_dest sessions{ii},'_spk_deconvolve_threshold' num2str(threshold) '.mat']);
    spike_runoff_mat = spike_output.spk_runoff_mat; %frames*trials*cells
    spike_runtrig_mat = spike_output.spk_runtrig_mat; %frames*trials*cells
    
    behav_dest = ['Z:\Analysis\2P_MovingDots_Analysis\behavioral_analysis\' days{ii} '\'];
    behav_output = load([behav_dest days{ii} '_behavAnalysis.mat']);
    mean_runtrig_every100ms = behav_output.mean_runtrig_every100ms;
    mean_runoff_every100ms = behav_output.mean_runoff_every100ms;
    
%     aveSpd_runtrig = spike_output.speed_runtrig; %frames*1
%     aveSpd_runoff = spike_output.speed_runoff; %frames*1
    
    nPCs_total = nPCs_total + size(spike_runoff_mat,3);
    ntrials_runtrig_total = ntrials_runtrig_total + size(spike_runoff_mat,2);
    ntrials_runoff_total = ntrials_runoff_total + size(spike_runtrig_mat,2);
    
    %for spike data, average first. average across cells and trials
    mean_spike_runtrig_cells = squeeze(mean(spike_runtrig_mat,2));
    mean_spike_runtrig = squeeze(mean(mean_spike_runtrig_cells,2));
    mean_spike_runoff_cells = squeeze(mean(spike_runoff_mat,2));
    mean_spike_runoff = squeeze(mean(mean_spike_runoff_cells,2));
    
    ave_spike_runtrig_allsession = cat(1,ave_spike_runtrig_allsession,mean_spike_runtrig');%when cat matrix, the matrix with NaN will be cat. And you have matrix with NaN because for some sessions there's no running trials fullfill the runtrig/runoff criteria
    ave_spike_runoff_allsession = cat(1,ave_spike_runoff_allsession,mean_spike_runoff');
    aveSpd_runtrig_100ms_all_session = cat(1,aveSpd_runtrig_100ms_all_session,mean_runtrig_every100ms);%when cat vectors, the empty vectors won't be cat
    aveSpd_runoff_100ms_all_session = cat(1,aveSpd_runoff_100ms_all_session,mean_runoff_every100ms);
    
%     aveSpd_runtrig_all_session = cat(2,aveSpd_runtrig_all_session,aveSpd_runtrig);%when cat vectors, the empty vectors won't be cat
%     aveSpd_runoff_all_session = cat(2,aveSpd_runoff_all_session,aveSpd_runoff);
end

ave_spike_runtrig_allsession = ave_spike_runtrig_allsession(all(~isnan(ave_spike_runtrig_allsession),2),:);% remove rows with NaN
ave_spike_runoff_allsession = ave_spike_runoff_allsession(all(~isnan(ave_spike_runoff_allsession),2),:);

ave_spike_runtrig_acrossessions = mean(ave_spike_runtrig_allsession)*30;
ste_spike_runtrig_acrossessions = std(ave_spike_runtrig_allsession)*30/sqrt(size(ave_spike_runtrig_allsession,1));

ave_spike_runoff_acrossessions = mean(ave_spike_runoff_allsession)*30;
ste_spike_runoff_acrossessions = std(ave_spike_runoff_allsession)*30/sqrt(size(ave_spike_runoff_allsession,1));

aveSpd_runtrig_100ms_acrossessions = mean(aveSpd_runtrig_100ms_all_session,1);
steSpd_runtrig_100ms_acrossessions = std(aveSpd_runtrig_100ms_all_session,0,1)/sqrt(size(aveSpd_runtrig_100ms_all_session,1));
aveSpd_runoff_100ms_acrossessions = mean(aveSpd_runoff_100ms_all_session,1);
steSpd_runoff_100ms_acrossessions = std(aveSpd_runoff_100ms_all_session,0,1)/sqrt(size(aveSpd_runoff_100ms_all_session,1));

% aveSpd_runtrig_acrossessions = mean(aveSpd_runtrig_all_session,2);
% steSpd_runtrig_acrossessions = std(aveSpd_runtrig_all_session,0,2)/sqrt(size(aveSpd_runtrig_all_session,2));
% aveSpd_runoff_acrossessions = mean(aveSpd_runoff_all_session,2);
% steSpd_runoff_acrossessions = std(aveSpd_runoff_all_session,0,2)/sqrt(size(aveSpd_runoff_all_session,2));

%%
%plot
x = (1:length(ave_spike_runtrig_acrossessions))/30;
figure; subplot(2,1,1);
shadedErrorBar(x,ave_spike_runtrig_acrossessions,ste_spike_runtrig_acrossessions,{'color',[0.1922 0.6392 0.3294]});hold on;
ylabel('firing rate');
ylim([0 2.5]);
vline(1,'r');
title(' Firing rate and speed aligned with running onset across all sessions');
subplot(2,1,2); 
%convert speed to cm/s: each unit = 2pai*r/128 (there're 128 pulses in total for each circle)
shadedErrorBar(x,aveSpd_runtrig_100ms_acrossessions*2*3.1415926*7.5/128, ...
    steSpd_runtrig_100ms_acrossessions*2*3.1415926*7.5/128,{'color',[0.0196 0.4392 0.6902]});hold on;
ylabel('speed(cm/s)');xlabel('time(s)');
%text(0.1,min(aveSpd_runtrig*2*3.1415926*7.5/128)+2,['n = ' num2str(size(frms_runTrig_mat,2))]);
ylim([-5 10]);
vline(1,'r');

x = (1:length(ave_spike_runtrig_acrossessions))/30;
figure; subplot(2,1,1);
shadedErrorBar(x,ave_spike_runoff_acrossessions,ste_spike_runoff_acrossessions,{'color',[0.1922 0.6392 0.3294]});hold on;
ylabel('firing rate');
ylim([0 2.5]);
vline(1,'r'); 
title('Firing rate and speed aligned with running offset across all sessions');
subplot(2,1,2); 
%convert speed to cm/s: each unit = 2pai*r/128 (there're 128 pulses in total for each circle)
shadedErrorBar(x,aveSpd_runoff_100ms_acrossessions*2*3.1415926*7.5/128, ...
    steSpd_runoff_100ms_acrossessions*2*3.1415926*7.5/128,{'color',[0.0196 0.4392 0.6902]});hold on;
ylabel('speed(cm/s)');xlabel('time(s)');
ylim([-5 10]);vline(1,'r');

