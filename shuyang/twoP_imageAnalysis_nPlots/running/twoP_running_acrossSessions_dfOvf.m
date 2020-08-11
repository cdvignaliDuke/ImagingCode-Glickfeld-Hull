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
speed_runtrig_100ms_all_session = [];% frames*total number of sessions that have runtrig
speed_runoff_100ms_all_session = [];
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
    speed_runtrig_100ms = behav_output.speed_runtrig_100ms; %frames*trials
    speed_runoff_100ms = behav_output.speed_runoff_100ms; %frames*trials
    
    %use avedfOvF_runoff_cells and avedfOvF_runtrig_cells for dfOvF averaging and ste with cell by cell variation n = #of mice
    %use aveSpd_runoff and aveSpd_runtrig for speed averaging
    
    nPCs_total = nPCs_total + size(ave_dfOvF_runtrig_cells,2);
    ntrials_runtrig_total = ntrials_runtrig_total + size(speed_runtrig_100ms,2);
    ntrials_runoff_total = ntrials_runoff_total + size(speed_runoff_100ms,2);
    
    ave_dfOvF_runtrig_allsession = cat(2,ave_dfOvF_runtrig_allsession,ave_dfOvF_runtrig_cells);%when cat matrix, the matrix with NaN will be cat. And you have matrix with NaN because for some sessions there's no running trials fullfill the runtrig/runoff criteria
    ave_dfOvF_runoff_allsession = cat(2,ave_dfOvF_runoff_allsession,ave_dfOvF_runoff_cells);
    
    speed_runtrig_100ms_all_session = cat(2,speed_runtrig_100ms_all_session,speed_runtrig_100ms);
    speed_runoff_100ms_all_session = cat(2,speed_runoff_100ms_all_session,speed_runoff_100ms);
    
end
ave_dfOvF_runtrig_allsession = ave_dfOvF_runtrig_allsession(:,all(~isnan(ave_dfOvF_runtrig_allsession)));% remove columns with NaN
ave_dfOvF_runoff_allsession = ave_dfOvF_runoff_allsession(:,all(~isnan(ave_dfOvF_runoff_allsession)));

ave_dfOvF_runtrig_acrossessions = mean(ave_dfOvF_runtrig_allsession,2); %average across all cells
ste_dfOvF_runtrig_acrossessions = std(ave_dfOvF_runtrig_allsession,0,2)/sqrt(size(ave_dfOvF_runtrig_allsession,2));

ave_dfOvF_runoff_acrossessions = mean(ave_dfOvF_runoff_allsession,2);
ste_dfOvF_runoff_acrossessions = std(ave_dfOvF_runoff_allsession,0,2)/sqrt(size(ave_dfOvF_runoff_allsession,2));

aveSpd100ms_runtrig_acrosssession = mean(speed_runtrig_100ms_all_session,2); %average across all running trials from all sessions
steSpd100ms_runtrig_acrosssession = std(speed_runtrig_100ms_all_session,0,2)/sqrt(size(speed_runtrig_100ms_all_session,2));

aveSpd100ms_runoff_acrosssession = mean(speed_runoff_100ms_all_session,2); 
steSpd100ms_runoff_acrosssession = std(speed_runoff_100ms_all_session,0,2)/sqrt(size(speed_runoff_100ms_all_session,2));



% save variable: df/f and speed. speed can be used when plotting spikes.

%%
%plot
x = (1:length(ave_dfOvF_runtrig_acrossessions))/30-1;
dfOvF_Onset_fig = figure; 
subplot(2,1,1);
shadedErrorBar(x,ave_dfOvF_runtrig_acrossessions,ste_dfOvF_runtrig_acrossessions,{'color',[0.8431 0.0980 0.1098]},{'Linewidth',1});hold on;
xlim([-1.1 1.6]);ylim([0 0.5]);
vline(0,'k'); ylabel('df/f');
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontSize',7);
%title([sessions ' running triggered average']);
subplot(2,1,2); %plot(speed_runtrig);
%convert speed to cm/s: each unit = 2pai*r/128 (there're 128 pulses in total for each circle)
shadedErrorBar(x,aveSpd100ms_runtrig_acrosssession*2*3.1415926*7.5/128, steSpd100ms_runtrig_acrosssession*2*3.1415926*7.5/128,{'color',[0 0 0]},{'Linewidth',1});hold on;
xlim([-1.1 1.6]); ylim([0 12]);
vline(0,'k');ylabel('speed(cm/s)');xlabel('time from running onset(s)');
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontSize',7);
%text(0.1,min(aveSpd_runtrig*2*3.1415926*7.5/128)+2,['n = ' num2str(size(frms_runTrig_mat,2))]);
dfOvF_Onset_fig.Units = 'centimeters';
dfOvF_Onset_fig.Position = [1 1 7 6];
fig_name = 'across_session_runonset';
path = 'Z:\Analysis\figures\figure2_2Prun\';
orient(dfOvF_Onset_fig,'landscape');
print(dfOvF_Onset_fig,[path,fig_name],'-r600','-dpdf');

dfOvF_Offset_fig = figure; 
subplot(2,1,1);
shadedErrorBar(x,ave_dfOvF_runoff_acrossessions,ste_dfOvF_runoff_acrossessions,{'color',[0.8431 0.0980 0.1098]},{'Linewidth',1});hold on;
xlim([-1.1 1.6]);ylim([0 0.5]);
vline(0,'k'); ylabel('df/f');
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontSize',7);
%title([sessions ' running triggered average']);
subplot(2,1,2); %plot(speed_runtrig);
%convert speed to cm/s: each unit = 2pai*r/128 (there're 128 pulses in total for each circle)
shadedErrorBar(x,aveSpd100ms_runoff_acrosssession*2*3.1415926*7.5/128, steSpd100ms_runoff_acrosssession*2*3.1415926*7.5/128,{'color',[0 0 0]},{'Linewidth',1});hold on;
xlim([-1.1 1.6]); ylim([0 12]);
vline(0,'k');ylabel('speed(cm/s)');xlabel('time from running offset(s)');
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontSize',7);
%text(0.1,min(aveSpd_runtrig*2*3.1415926*7.5/128)+2,['n = ' num2str(size(frms_runTrig_mat,2))]);
dfOvF_Offset_fig.Units = 'centimeters';
dfOvF_Offset_fig.Position = [1 1 7 6];
fig_name = 'across_session_runoffset';
path = 'Z:\Analysis\figures\figure2_2Prun\';
orient(dfOvF_Offset_fig,'landscape');
print(dfOvF_Offset_fig,[path,fig_name],'-r600','-dpdf');

















