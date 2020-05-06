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
ave_dfOvF_runtrig_allsession = [];% frames*total number of cells in all sessions
ave_dfOvF_runoff_allsession = [];
aveSpd_runtrig_all_session = [];% frames*total number of sessions that have runtrig
aveSpd_runoff_all_session = [];
nPCs_total = 0;
ntrials_runtrig_total = 0;
ntrials_runoff_total = 0;

for ii = 1:length(sessions) 
    %load things
    image_analysis_dest = [image_analysis_base, sessions{ii}, '\'];
    dfOvF_output = load([image_analysis_dest sessions{ii} '_dfOvF.mat']);
    ave_dfOvF_runtrig_cells = dfOvF_output.ave_dfOvF_runtrig_cells; %frames*cells
    ave_dfOvF_runoff_cells = dfOvF_output.ave_dfOvF_runoff_cells; %frames*cells
    behav_dest = ['Z:\Analysis\2P_MovingDots_Analysis\behavioral_analysis\' days{ii} '\'];
    behav_output = load([behav_dest days{ii} '_behavAnalysis.mat']);
    aveSpd_runtrig = behav_output.aveSpd_runtrig; %frames*1
    aveSpd_runoff = behav_output.aveSpd_runoff; %frames*1
    frms_runTrig_mat = behav_output.frms_runTrig_mat;
    frms_runoff_mat = behav_output.frms_runoff_mat;
    %use avedfOvF_runoff_cells and avedfOvF_runtrig_cells for dfOvF averaging and ste with cell by cell variation n = #of mice
    %use aveSpd_runoff and aveSpd_runtrig for speed averaging, ste will be wrong for now, need to rerun twoP_runing_dfOvF_graphs to save spd in all trials to calculate ste
    
    nPCs_total = nPCs_total + size(ave_dfOvF_runtrig_cells,2);
    ntrials_runtrig_total = ntrials_runtrig_total + size(frms_runTrig_mat,2);
    ntrials_runoff_total = ntrials_runoff_total + size(frms_runoff_mat,2);
    
    ave_dfOvF_runtrig_allsession = cat(2,ave_dfOvF_runtrig_allsession,ave_dfOvF_runtrig_cells);%when cat matrix, the matrix with NaN will be cat. And you have matrix with NaN because for some sessions there's no running trials fullfill the runtrig/runoff criteria
    ave_dfOvF_runoff_allsession = cat(2,ave_dfOvF_runoff_allsession,ave_dfOvF_runoff_cells);
    aveSpd_runtrig_all_session = cat(2,aveSpd_runtrig_all_session,aveSpd_runtrig);%when cat vectors, the empty vectors won't be cat
    aveSpd_runoff_all_session = cat(2,aveSpd_runoff_all_session,aveSpd_runoff);
    
end
ave_dfOvF_runtrig_allsession = ave_dfOvF_runtrig_allsession(:,all(~isnan(ave_dfOvF_runtrig_allsession)));% remove columns with NaN
ave_dfOvF_runoff_allsession = ave_dfOvF_runoff_allsession(:,all(~isnan(ave_dfOvF_runoff_allsession)));

ave_dfOvF_runtrig_acrossessions = mean(ave_dfOvF_runtrig_allsession,2);
ste_dfOvF_runtrig_acrossessions = std(ave_dfOvF_runtrig_allsession,0,2)/sqrt(size(ave_dfOvF_runtrig_allsession,2));

ave_dfOvF_runoff_acrossessions = mean(ave_dfOvF_runoff_allsession,2);
ste_dfOvF_runoff_acrossessions = std(ave_dfOvF_runoff_allsession,0,2)/sqrt(size(ave_dfOvF_runoff_allsession,2));

aveSpd_runtrig_acrossessions = mean(aveSpd_runtrig_all_session,2);
steSpd_runtrig_acrossessions = std(aveSpd_runtrig_all_session,0,2)/sqrt(size(aveSpd_runtrig_all_session,2));
aveSpd_runoff_acrossessions = mean(aveSpd_runoff_all_session,2);
steSpd_runoff_acrossessions = std(aveSpd_runoff_all_session,0,2)/sqrt(size(aveSpd_runoff_all_session,2));

bin = length(aveSpd_runtrig_acrossessions)/3;
mean_runtrig_every100ms_across = zeros(1,length(aveSpd_runtrig_acrossessions));
ste_runtrig_every100ms_across = zeros(1,length(aveSpd_runtrig_acrossessions));

for x1 = 1:bin
    y = 3*x1 -2;
    mean_runtrig_every100ms_across(y:y+2) = mean(aveSpd_runtrig_acrossessions(y:y+2));
    ste_runtrig_every100ms_across(y:y+2) = mean(steSpd_runtrig_acrossessions(y:y+2));
    if abs(mean_runtrig_every100ms_across(y))<2 %on average moves less than 2 units per 0.1s, consider as stationary
     mean_runtrig_every100ms_across(y:y+2) = 0;
     ste_runtrig_every100ms_across(y:y+2) = 0;
    end
end


bin = length(aveSpd_runoff_acrossessions)/3;
mean_runoff_every100ms_across = zeros(1,length(aveSpd_runoff_acrossessions));
ste_runoff_every100ms_across = zeros(1,length(aveSpd_runoff_acrossessions));
a = 1;
for x1 = 1:bin
    y = 3*x1 -2;
    mean_runoff_every100ms_across(y:y+2) = mean(aveSpd_runoff_acrossessions(y:y+2));
    ste_runoff_every100ms_across(y:y+2) = mean(steSpd_runoff_acrossessions(y:y+2));
    if abs(mean_runoff_every100ms_across(y))<2 %on average moves less than 2 units per 0.1s, consider as stationary
     mean_runoff_every100ms_across(y:y+2) = 0;
     ste_runoff_every100ms_across(y:y+2) = 0;
    end
end

% save variable: df/f and speed. speed can be used when plotting spikes.

%%
%plot
x1 = (1:length(ave_dfOvF_runtrig_acrossessions))/30;
figure; subplot(2,1,1);
shadedErrorBar(x1,ave_dfOvF_runtrig_acrossessions,ste_dfOvF_runtrig_acrossessions,{'color',[0.1922 0.6392 0.3294]});hold on;
ylabel('df/f');
ylim([0 0.5]);vline(1,'r');
title(' running triggered average across all sessions');
subplot(2,1,2); 
%convert speed to cm/s: each unit = 2pai*r/128 (there're 128 pulses in total for each circle)
shadedErrorBar(x1,mean_runtrig_every100ms_across*2*3.1415926*7.5/128, ...
    ste_runtrig_every100ms_across*2*3.1415926*7.5/128,{'color',[0.0196 0.4392 0.6902]});hold on;
ylabel('speed(cm/s)');xlabel('time(s)');
%text(0.1,min(aveSpd_runtrig*2*3.1415926*7.5/128)+2,['n = ' num2str(size(frms_runTrig_mat,2))]);
ylim([0 12]);vline(1,'r');

x2 = (1:length(ave_dfOvF_runoff_acrossessions))/30;
figure; subplot(2,1,1);
shadedErrorBar(x2,ave_dfOvF_runoff_acrossessions,ste_dfOvF_runoff_acrossessions,{'color',[0.1922 0.6392 0.3294]});hold on;
ylabel('df/f');
ylim([0 0.5]);vline(1,'r'); 
title('df/f and speed aligned with running offset across all sessions');
subplot(2,1,2); 
%convert speed to cm/s: each unit = 2pai*r/128 (there're 128 pulses in total for each circle)
shadedErrorBar(x2,mean_runoff_every100ms_across*2*3.1415926*7.5/128, ...
    ste_runoff_every100ms_across*2*3.1415926*7.5/128,{'color',[0.0196 0.4392 0.6902]});hold on;
ylabel('speed(cm/s)');xlabel('time(s)');
ylim([0 10]);vline(1,'r');

















