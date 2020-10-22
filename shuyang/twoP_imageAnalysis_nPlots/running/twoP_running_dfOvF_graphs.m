%% 2 photon: running - dfOvF
%plot average dfOvF for all cells in one trial, aligned with running onset and offset
%ave_dfOvF_runtrig_cells: frame*cells, averaged across trials
%ave_dfOvF_runtrig:     frames*1, averaged across cells
%ste_dfOvF_runtrig:     frames*1, cell by cell variation
%aveSpd_runtrig:        1*frames, averaged across trials
%steSpd_runtrig:        1*frames if more than 1 trial, trial by trial variation.
%                       a string that says no variation if only 1 trial


%% assign document paths and experimental sessions

%sessions = {'190429_img1021','190430_img1023','190507_img1024','190603_img1025'};
%days = {'1021-190429_1','1023-190430_1','1024-190507_1','1025-190603_1'};

clear;
sessions = '190507_img1024'; 
image_analysis_base    = 'Z:\Analysis\2P_MovingDots_Analysis\imaging_analysis\'; %stores the data on crash in the movingDots analysis folder
image_analysis_dest = [image_analysis_base, sessions, '\'];

% behavior analysis results 
days = '1024-190507_1';
behav_dest = ['Z:\Analysis\2P_MovingDots_Analysis\behavioral_analysis\' days '\'];
color_code = {'b','r','k','c'};

%% SECTION I df/f and speed
% load df/f bottom (F using bottom 10 percent of fluorescence values)
dfOvF_output = load([image_analysis_dest sessions '_dfOvF.mat']);
dfOvF_btm = dfOvF_output.dfOvF_btm_cl;
avedfOvF_btm_cl = mean(dfOvF_btm,2);
stedfOvF_btm_cl = std(dfOvF_btm,0,2)/sqrt(size(dfOvF_btm,2));
behav_output = load([behav_dest days '_behavAnalysis.mat']);
speed = behav_output.speed;
frames = 20250:21150;%range of frames you want to plot
speed_plot = speed(frames);
%smooth speed into every 100ms
bin = length(speed_plot)/3;
mean_spd_every100ms = zeros(1,length(speed_plot));
for x1 = 1:bin
    y = 3*x1 -2;
    mean_spd_every100ms(y:y+2) = mean(speed_plot(y:y+2));
    if abs(mean_spd_every100ms(y))<2 %on average moves less than 2 units per 0.1s, consider as stationary
     mean_spd_every100ms(y:y+2) = 0;
    end
end
% convert frames to seconds
x_plot = frames/30 - frames(1)/30;

dfOvF_speed_fig = figure;
yyaxis left; hold on;
plot(x_plot,mean_spd_every100ms*2*3.1415926*7.5/128,'linewidth', 1,'color','k');
trial1_start = 172/30; trial1_length = (330-172)/30;
trial2_start = 625/30; trial2_length = (774-625)/30;
r1 = rectangle('Position',[trial1_start -10.3333*2*3.1415926*7.5/128 trial1_length 35+10.3333*2*3.1415926*7.5/128] , ...
    'EdgeColor',[0.7825    0.7825    0.7825], 'FaceColor', [0.7825    0.7825    0.7825]);
uistack(r1,'bottom');
r2 = rectangle('Position',[trial2_start -10.3333*2*3.1415926*7.5/128 trial2_length 35+10.3333*2*3.1415926*7.5/128] , ...
    'EdgeColor',[0.7825    0.7825    0.7825], 'FaceColor', [0.7825    0.7825    0.7825]);
uistack(r2,'bottom');
ylim([-5 35]);xlim([min(x_plot) max(x_plot)]);
ylabel('speed(cm/s)')
set(gca,'YColor','k');
yyaxis right;
shadedErrorBar(x_plot,avedfOvF_btm_cl(frames),stedfOvF_btm_cl(frames),{'color',[0.8431    0.0980    0.1098]},{'Linewidth',1}); 
%p = plot(x_plot,avedfOvF_btm_cl(x1),'Color',);
ylabel('df/F'); set(gca,'YColor','r');
ylim([-0.1 1.4]);xlim([min(x_plot) max(x_plot)]);

xlabel('time(s)');
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontSize',7);
dfOvF_speed_fig.Units = 'centimeters';
dfOvF_speed_fig.Position = [1 3 7 4];
fig_name = 'eg_speed_and_dfOvF_190603_img1025_20300-21100';
path = 'Z:\Analysis\figures\figure2_2Prun\';
orient(dfOvF_speed_fig,'landscape')
print(dfOvF_speed_fig,[path,fig_name],'-r600','-depsc');

%save([image_analysis_dest sessions '_dfOvF.mat'],'avedfOvF_btm_cl','-append');
% title(['speed and average df/f' days]);
% savefig([image_analysis_dest sessions '_avedfOvF_wspeed']);

%% SECTION II run triggered ave and run offset ave --- df/f
% the mice has to be stationary for all frames before or after running
dfOvF_output = load([image_analysis_dest sessions '_dfOvF.mat']);
dfOvF_btm = dfOvF_output.dfOvF_btm_cl;
behav_output = load([behav_dest days '_behavAnalysis.mat']);
frames_run_cell = behav_output.frames_run_cell;
speed = double(behav_output.speed);

period = 9; % # of frames right before and after running--- probably not useful
befoRunStay = 27; %1s, # of frames that speed =0 before running
befoRun = 30;%# of frames to include before running starts
totalruntrig = 75; %2s # total T = of frames before running onset + # of frames following running onset/ # of frames before running offset + # of frames after running offset
aftRunStay = 42; % # of frames that speed = 0 after running
aftRun = 45;%# of frames to include after running ends
totalrunoff = 75; 
[~,~,frms_runTrig_mat,frms_runoff_mat,frames_runoff_include,~]= findFrames_runWindows_2P(speed,...
    frames_run_cell,period,befoRunStay,totalruntrig,totalrunoff,aftRunStay,befoRun,aftRun);
% save([behav_dest days '_behavAnalysis.mat'],'befoRunStay','befoRun','totalruntrig',...
%     'aftRunStay','aftRun','totalrunoff','frms_runTrig_mat',...
%     'frms_runoff_mat','frames_runoff_include','-append');

% aligned with running onset ----------------------------------------------------------------------------
dfOvF_runtrig_mat = zeros(size(frms_runTrig_mat,1),size(frms_runTrig_mat,2),size(dfOvF_btm,2));
for i = 1: size(frms_runTrig_mat,2)                                    % for each running trig window
    dfOvF_runtrig_mat(:,i,:) = dfOvF_btm(frms_runTrig_mat(:,i),:);         % df/f for all cells during that window, frame*cells 
end
ave_dfOvF_runtrig_cells = squeeze(mean(dfOvF_runtrig_mat,2)); % average across trials. % if there's no trials fullfill the criteria, then ave_dfOvF will be [] because dfOvF_runoff_mat is [] and mean [] is NaN.
ave_dfOvF_runtrig = squeeze(mean(ave_dfOvF_runtrig_cells,2)); % average across cells
ste_dfOvF_runTrig = std(ave_dfOvF_runtrig_cells,0,2)/sqrt(size(ave_dfOvF_runtrig_cells,2)); % ste is reflecting cell by cell variation
% speed
speed_runtrig = speed(frms_runTrig_mat);
if size(frms_runTrig_mat,2) > 1
    aveSpd_runtrig = mean(speed_runtrig,2);
    steSpd_runtrig = std(speed_runtrig,0,2)/sqrt(size(speed_runtrig,2)); % ste is reflecting trial by trial basis
else
    aveSpd_runtrig = speed_runtrig;
    steSpd_runtrig = 'no trail by trial variation, only one trial';
end
% bin speed into every 100ms and then average across trials and calculate ste
bin = length(speed)/3;
mean_speed_100ms = zeros(1,length(speed));

for x = 1:bin
    y = 3*x -2;
    mean_speed_100ms(y:y+2) = mean(speed(y:y+2));
    if abs(mean_speed_100ms(y))<2 %on average moves less than 2 units per 0.1s, consider as stationary
     mean_speed_100ms(y:y+2) = 0;
    end
end
speed_runtrig_100ms = mean_speed_100ms(frms_runTrig_mat);
avespd_runtrig_100ms = mean(speed_runtrig_100ms,2);
stespd_runtrig_100ms = std(speed_runtrig_100ms,0,2)/sqrt(size(speed_runtrig_100ms,2));

% save([image_analysis_dest sessions '_dfOvF.mat'],'ave_dfOvF_runtrig_cells',...
%     'ave_dfOvF_runtrig','ste_dfOvF_runTrig','-append');
save([behav_dest days '_behavAnalysis.mat'],'aveSpd_runtrig','steSpd_runtrig',...
    'speed_runtrig_100ms','avespd_runtrig_100ms','stespd_runtrig_100ms','-append');

%plot
% x = ((1:size(frms_runTrig_mat,1))/30)-1;
% dfOvF_Onset_fig = figure; subplot(2,1,1);
% shadedErrorBar(x,ave_dfOvF_runtrig,ste_dfOvF_runTrig,{'color',[0.8431 0.0980 0.1098]},{'Linewidth',1});hold on;
% xlim([-1.1 1.6]);ylim([0 0.8]);
% vline(0,'k'); ylabel('df/f');
% a = get(gca,'XTickLabel');
% set(gca,'XTickLabel',a,'FontSize',7);
% %title([sessions ' running triggered average']);
% subplot(2,1,2); %plot(speed_runtrig);
% %convert speed to cm/s: each unit = 2pai*r/128 (there're 128 pulses in total for each circle)
% if size(frms_runTrig_mat,2) > 1
%     shadedErrorBar(x,avespd_runtrig_100ms*2*3.1415926*7.5/128, stespd_runtrig_100ms*2*3.1415926*7.5/128,{'color',[0 0 0]},{'Linewidth',1});hold on;
% else
%     plot(x,aveSpd_runtrig*2*3.1415926*7.5/128,'k'); hold on;
% end
% xlim([-1.1 1.6]); ylim([0 15]);
% vline(0,'k');ylabel('speed(cm/s)');xlabel('time from running onset(s)');
% a = get(gca,'XTickLabel');
% set(gca,'XTickLabel',a,'FontSize',7);
% %text(0.1,min(aveSpd_runtrig*2*3.1415926*7.5/128)+2,['n = ' num2str(size(frms_runTrig_mat,2))]);
% dfOvF_Onset_fig.Units = 'centimeters';
% dfOvF_Onset_fig.Position = [1 1 7 6];
% fig_name = ['eg_session_runonset_',sessions];
% path = 'Z:\Analysis\figures\figure2_2Prun\';
% orient(dfOvF_Onset_fig,'landscape');
% print(dfOvF_Onset_fig,[path,fig_name],'-r600','-depsc');
% %savefig([image_analysis_dest sessions '_runTrigAve_dfOvF_wspeed100ms']);


% aligned with running offset ------------------------------------------------------------------------
dfOvF_runoff_mat = zeros(size(frms_runoff_mat,1),size(frms_runoff_mat,2),size(dfOvF_btm,2));
for i = 1: size(frms_runoff_mat,2)                                    % for each running trig window
    dfOvF_runoff_mat(:,i,:) = dfOvF_btm(frms_runoff_mat(:,i),:);         % df/f for all cells during that window, frame*cells 
end
ave_dfOvF_runoff_cells = squeeze(mean(dfOvF_runoff_mat,2)); % if there's no trials fullfill the criteria, then ave_dfOvF will be [] because dfOvF_runoff_mat is [] and mean [] is NaN.
ave_dfOvF_runoff = squeeze(mean(ave_dfOvF_runoff_cells,2));
ste_dfOvF_runoff = std(ave_dfOvF_runoff_cells,0,2)/sqrt(size(ave_dfOvF_runoff_cells,2)); % a bunch of zeros
% speed
speed_runoff = speed(frms_runoff_mat);
if size(frms_runoff_mat,2)>1
    aveSpd_runoff = mean(speed_runoff,2);
    steSpd_runoff = std(speed_runoff,0,2)/sqrt(size(speed_runoff,2));
else
    aveSpd_runoff = speed_runoff;
    steSpd_runoff = 'no trail by trial variation, only one trial';
end

% use bined every 100ms
speed_runoff_100ms = mean_speed_100ms(frms_runoff_mat);
avespd_runoff_100ms = mean(speed_runoff_100ms,2);
stespd_runoff_100ms = std(speed_runoff_100ms,0,2)/sqrt(size(speed_runoff_100ms,2));


% save([image_analysis_dest sessions '_dfOvF.mat'],'ave_dfOvF_runoff_cells',...
%     'ave_dfOvF_runoff','ste_dfOvF_runoff','-append');
save([behav_dest days '_behavAnalysis.mat'],'aveSpd_runoff','steSpd_runoff',...
    'frms_runTrig_mat','frms_runoff_mat','speed_runoff_100ms','avespd_runoff_100ms','stespd_runoff_100ms','-append');
%plot

% x = ((1:size(frms_runoff_mat,1))/30)-1;
% dfOvF_runOff_fig = figure; subplot(2,1,1);
% shadedErrorBar(x,ave_dfOvF_runoff,ste_dfOvF_runoff,{'color',[0.8431 0.0980 0.1098]},{'Linewidth',1});hold on;
% xlim([-1.1 1.6]);ylim([0 0.8]);
% vline(0,'k'); ylabel('df/f');
% a = get(gca,'XTickLabel');
% set(gca,'XTickLabel',a,'FontSize',7);
% %title([sessions ' running triggered average']);
% subplot(2,1,2); %plot(speed_runtrig);
% %convert speed to cm/s: each unit = 2pai*r/128 (there're 128 pulses in total for each circle)
% if size(frms_runoff_mat,2) > 1
%     shadedErrorBar(x,avespd_runoff_100ms*2*3.1415926*7.5/128, stespd_runoff_100ms*2*3.1415926*7.5/128,{'color',[0 0 0]},{'Linewidth',1});hold on;
% else
%     plot(x,aveSpd_runoff*2*3.1415926*7.5/128,'k'); hold on;
% end
% xlim([-1.1 1.6]); ylim([0 15]);
% vline(0,'k');ylabel('speed(cm/s)');xlabel('time from running offset(s)');
% a = get(gca,'XTickLabel');
% set(gca,'XTickLabel',a,'FontSize',7);
% %text(0.1,min(aveSpd_runtrig*2*3.1415926*7.5/128)+2,['n = ' num2str(size(frms_runTrig_mat,2))]);
% dfOvF_runOff_fig.Units = 'centimeters';
% dfOvF_runOff_fig.Position = [1 1 7 6];
% fig_name = ['eg_session_runoffset_',sessions];
% path = 'Z:\Analysis\figures\figure2_2Prun\';
% orient(dfOvF_runOff_fig,'landscape');
% print(dfOvF_runOff_fig,[path,fig_name],'-r600','-depsc');

%text(1.8,min(aveSpd_runtrig*2*3.1415926*7.5/128)+2,['n = ' num2str(size(frms_runoff_mat,2))]);
%savefig([image_analysis_dest sessions '_runoffset_dfOvF_wspeed']);


%% SECTION II run triggered ave and run offset ave --- df/f --- 
% REGARDLESS OF whether the animal is stationary all the time before/after running
% this part hasn't been rerunned with the new df/f, and variables are not
% saved!
%{
behav_output = load([behav_dest days '_behavAnalysis.mat']);
frames_run_cell = behav_output.frames_run_cell;
speed = double(behav_output.speed);

period = 9; % # of frames before and after running
befoRunStay = 30; %1s, # of frames that speed =0 before running
totalT = 60; %2s # of frmaes before running onset + # of frames following running onset/ # of frames before running offset + # of frames after running offset
aftRunOff = 30; %1s,# of frames that speed = 0 after running
[~,~,frms_runTrig_NC_mat,frms_runoff_NC_mat,~]= findFrames_run_no_criteria_2P(speed,...
    frames_run_cell,period,befoRunStay,totalT,aftRunOff);
% aligned with running onset ----------------------------------------------------------------------------
dfOvF_runtrig_NC_mat = zeros(size(frms_runTrig_NC_mat,1),size(frms_runTrig_NC_mat,2),size(dfOvF_btm,2));
for i = 1: size(frms_runTrig_NC_mat,2)                                    % for each running trig window
    dfOvF_runtrig_NC_mat(:,i,:) = dfOvF_btm(frms_runTrig_NC_mat(:,i),:);         % df/f for all cells during that window, frame*cells 
end
ave_dfOvF_runtrig_NC_trials = squeeze(mean(dfOvF_runtrig_NC_mat,2));
ave_dfOvF_runtrig_NC = squeeze(mean(ave_dfOvF_runtrig_NC_trials,2));
ste_dfOvF_runTrig_NC = std(ave_dfOvF_runtrig_NC_trials,0,2)/sqrt(size(ave_dfOvF_runtrig_NC_trials,2)); % a bunch of zeros
% speed
speed_runtrig_NC = speed(frms_runTrig_NC_mat);
aveSpd_runtrig_NC = mean(speed_runtrig_NC,2);
steSpd_runtrig_NC = std(speed_runtrig_NC,0,2)/sqrt(size(speed_runtrig_NC,2));
%plot
x = (1:size(frms_runTrig_NC_mat,1));
figure; subplot(2,1,1);
shadedErrorBar(x,ave_dfOvF_runtrig_NC,ste_dfOvF_runTrig_NC);hold on;
title([sessions ' running triggered average INCLUSIVE']);
vline(31,'r'); ylabel('df/f');
subplot(2,1,2); %plot(speed_runtrig);
shadedErrorBar(x,aveSpd_runtrig_NC*2*3.1415926*7.5/128, steSpd_runtrig_NC*2*3.1415926*7.5/128);hold on;
vline(31,'r');ylabel('speed(cm/s)');
text(0.1,min(aveSpd_runtrig_NC*2*3.1415926*7.5/128),['n = ' num2str(size(frms_runTrig_NC_mat,2))]);
savefig([image_analysis_dest sessions '_runTrigAve_dfOvF_wspeed_INCLUSIVE']);

% aligned with running offset ---------------------------------------------------------------------------------------
dfOvF_runoff_NC_mat = zeros(size(frms_runoff_NC_mat,1),size(frms_runoff_NC_mat,2),size(dfOvF_btm,2));
for i = 1: size(frms_runoff_NC_mat,2)                                    % for each running trig window
    dfOvF_runoff_NC_mat(:,i,:) = dfOvF_btm(frms_runoff_NC_mat(:,i),:);         % df/f for all cells during that window, frame*cells 
end
ave_dfOvF_runoff_NC_trials = squeeze(mean(dfOvF_runoff_NC_mat,2));
ave_dfOvF_runoff_NC = squeeze(mean(ave_dfOvF_runoff_NC_trials,2));
ste_dfOvF_runoff_NC = std(ave_dfOvF_runoff_NC_trials,0,2)/sqrt(size(ave_dfOvF_runoff_NC_trials,2)); % a bunch of zeros
% speed
speed_runoff_NC = speed(frms_runoff_NC_mat);
aveSpd_runoff_NC = mean(speed_runoff_NC,2);
steSpd_runoff_NC = std(speed_runoff_NC,0,2)/sqrt(size(speed_runoff_NC,2));
%plot
x = (1:size(frms_runoff_NC_mat,1));
figure; subplot(2,1,1); 
shadedErrorBar(x,ave_dfOvF_runoff_NC,ste_dfOvF_runoff_NC);hold on;
vline(31,'r'); ylabel('df/f'); xlabel('frame');
title([sessions ' aligned to running offset - average INCLUSIVE']);
subplot(2,1,2); %plot(speed_runoff);
shadedErrorBar(x,aveSpd_runoff_NC*2*3.1415926*7.5/128, steSpd_runoff_NC*2*3.1415926*7.5/128);hold on;
vline(31,'r');ylabel('speed(cm/s)'); xlabel('frame');
text(1.8,min(aveSpd_runoff_NC*2*3.1415926*7.5/128),['n = ' num2str(size(frms_runoff_NC_mat,2))]);
savefig([image_analysis_dest sessions '_runoffset_dfOvF_wspeed_INCLUSIVE']);
%}

%% df/f stationary vs. running -------probably NOT USEFUL
%{ 
behav_output = load([behav_dest days '_behavAnalysis.mat']);
frm_run_cell = behav_output.frames_run_cell;
frm_run = cell2mat(frm_run_cell);
frm_stay_cell = behav_output.frames_stay_cell;
frm_stay = cell2mat(frm_stay_cell);

avedfOvFrun = zeros(1,size(TCave,2));
for i = 1: size(TCave,2) %for each cell
avedfOvFrun(i) = mean(dfOvF(frm_run,i));
end 

avedfOvFstay = zeros(1,size(TCave,2));
for i = 1: size(TCave,2) %for each cell
avedfOvFstay(i) = mean(dfOvF(frm_stay,i));
end 

grand_dfOvFrun = mean (avedfOvFrun);
grand_dfOvFstay = mean (avedfOvFstay);

figure;
scatter(avedfOvFstay,avedfOvFrun,'filled','r'); hold on;
scatter(grand_dfOvFstay, grand_dfOvFrun,'filled','k');hold on;
xlabel('stationary'); ylabel('running');
%xlim([-0.1,0.5]); ylim([min(avedfOvFstay),max(avedfOvFstay)]);
refline(1,0);
title('mean df/f for each cell');
savefig([image_analysis_dest sessions '_meandfOvFstayVSrun']);
%}






