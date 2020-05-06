%% assign document paths and experimental sessions
clear;
sessions = {'191114_img1040','191115_img1039','191115_img1041','191115_img1042','200316_img1064_airpuff_2'};
image_analysis_base  = 'Z:\Analysis\Airpuff_analysis\imaging_analysis\';%stores the data on crash in the movingDots analysis folder
% behavior analysis results
days = {'1040-191114_1','1039-191115_1','1041-191115_1','1042-191115_1','1064-200316_2'};
color_code = {'b','r','k','c'};

%% df/f
ave_dfOvF_runairpuff_allsession = [];% frames*total number of cells in all sessions
ave_dfOvF_stayairpuff_allsession = [];
aveSpd_run_all_session = [];% frames*total number of sessions that have runtrig
aveSpd_stay_all_session = [];
avespk_runairpuff_allsession = [];
avespk_stayairpuff_allsession = [];
nPCs_total = 0;
% ntrials_run_total = 0;
% ntrials_stay_total = 0;

% 'spd_airstimAve_allrun','spd_airstimAve_stay','spk_airpuff_allrun_cells','spk_airpuff_stay_cells'
% 'dfOvFbtm_airpuff_allrun_cells','dfOvFbtm_airpuff_stay_cells'
for ii = 1:length(sessions) 
    %load things
    image_analysis_dest = [image_analysis_base, sessions{ii}, '\']; 
    dfOvF_strct = load([image_analysis_dest sessions{ii} '_dfOvF.mat']);
    dfOvFbtm_airpuff_allrun_cells = dfOvF_strct.dfOvFbtm_airpuff_allrun_cells;%frame*cells
    dfOvFbtm_airpuff_stay_cells = dfOvF_strct.dfOvFbtm_airpuff_stay_cells;
    behav_dest = ['Z:\Analysis\Airpuff_analysis\behavioral_analysis\' days{ii}];
    behav_struct = load([behav_dest '\' days{ii} '_behavAnalysis.mat']);
    spd_airstimAve_allrun = behav_struct.spd_airstimAve_allrun;
    spd_airstimAve_stay = behav_struct.spd_airstimAve_stay;
    threshold = -4;
    image_analysis_dest = [image_analysis_base, sessions{ii}, '\'];
    spk_deconv_output = load([image_analysis_dest sessions{ii},'_spk_deconvolve_threshold' num2str(threshold) '.mat']);
    spk_airpuff_allrun_cells = spk_deconv_output.spk_airpuff_allrun_cells;
    spk_airpuff_stay_cells = spk_deconv_output.spk_airpuff_stay_cells;
    
    nPCs_total = nPCs_total + size(spk_airpuff_stay_cells,2);
    %ntrials_run_total = ntrials_run_total + size(,2);
    %ntrials_stay_total = ntrials_stay_total + size(,2);
    
    ave_dfOvF_runairpuff_allsession = cat(2,ave_dfOvF_runairpuff_allsession,dfOvFbtm_airpuff_allrun_cells);%when cat matrix, the matrix with NaN will be cat. And you have matrix with NaN because for some sessions there's no running trials fullfill the runtrig/runoff criteria
    ave_dfOvF_stayairpuff_allsession = cat(2,ave_dfOvF_stayairpuff_allsession,dfOvFbtm_airpuff_stay_cells);
    aveSpd_run_all_session = cat(1,aveSpd_run_all_session,spd_airstimAve_allrun);%when cat vectors, the empty vectors won't be cat
    aveSpd_stay_all_session = cat(1,aveSpd_stay_all_session,spd_airstimAve_stay);
    avespk_runairpuff_allsession = cat(2,avespk_runairpuff_allsession,spk_airpuff_allrun_cells);
    avespk_stayairpuff_allsession = cat(2,avespk_stayairpuff_allsession,spk_airpuff_stay_cells);
end

ave_dfOvF_runairpuff_acrossessions = mean(ave_dfOvF_runairpuff_allsession,2);
ste_dfOvF_runairpuff_acrossessions = std(ave_dfOvF_runairpuff_allsession,0,2)/sqrt(size(ave_dfOvF_runairpuff_allsession,2));

ave_dfOvF_stayairpuff_acrossessions = mean(ave_dfOvF_stayairpuff_allsession,2);
ste_dfOvF_stayairpuff_acrossessions = std(ave_dfOvF_stayairpuff_allsession,0,2)/sqrt(size(ave_dfOvF_stayairpuff_allsession,2));

aveFR_runairpuff_acrosssessions = mean(avespk_runairpuff_allsession,2)*30;
steFR_runairpuff_acrosssessions = std(avespk_runairpuff_allsession,0,2)/sqrt(size(avespk_runairpuff_allsession,2))*30;

aveFR_stayairpuff_acrosssessions = mean(avespk_stayairpuff_allsession,2)*30;
steFR_stayairpuff_acrosssessions = std(avespk_stayairpuff_allsession,0,2)/sqrt(size(avespk_stayairpuff_allsession,2))*30;

aveSpd_run_acrossessions = mean(aveSpd_run_all_session,1);
steSpd_run_acrossessions = std(aveSpd_run_all_session,0,1)/sqrt(size(aveSpd_run_all_session,1));
aveSpd_stay_acrossessions = mean(aveSpd_stay_all_session,1);
steSpd_stay_acrossessions = std(aveSpd_stay_all_session,0,1)/sqrt(size(aveSpd_stay_all_session,1));

bin = length(aveSpd_run_acrossessions)/3;
mean_run_every100ms_across = zeros(1,length(aveSpd_run_acrossessions));
ste_run_every100ms_across = zeros(1,length(aveSpd_run_acrossessions));

for x1 = 1:bin
    y = 3*x1 -2;
    mean_run_every100ms_across(y:y+2) = mean(aveSpd_run_acrossessions(y:y+2));
    ste_run_every100ms_across(y:y+2) = mean(steSpd_run_acrossessions(y:y+2));
    if abs(mean_run_every100ms_across(y))<2 %on average moves less than 2 units per 0.1s, consider as stationary
     mean_run_every100ms_across(y:y+2) = 0;
     ste_run_every100ms_across(y:y+2) = 0;
    end
end


bin = length(aveSpd_stay_acrossessions)/3;
mean_stay_every100ms_across = zeros(1,length(aveSpd_stay_acrossessions));
ste_stay_every100ms_across = zeros(1,length(aveSpd_stay_acrossessions));
a = 1;
for x1 = 1:bin
    y = 3*x1 -2;
    mean_stay_every100ms_across(y:y+2) = mean(aveSpd_stay_acrossessions(y:y+2));
    ste_stay_every100ms_across(y:y+2) = mean(steSpd_stay_acrossessions(y:y+2));
    if abs(mean_stay_every100ms_across(y))<2 %on average moves less than 2 units per 0.1s, consider as stationary
     mean_stay_every100ms_across(y:y+2) = 0;
     ste_stay_every100ms_across(y:y+2) = 0;
    end
end

% save variable: df/f and speed. speed can be used when plotting spikes.

%%
%plot
fig = figure;

x = (1: size(ave_dfOvF_runairpuff_acrossessions,1))/30;
subplot(3,2,1);
shadedErrorBar(x,ave_dfOvF_runairpuff_acrossessions,ste_dfOvF_runairpuff_acrossessions,{'color',[0.1373 0.5451 0.2706]});
title('airpuff stim running'); ylabel('df/f');
xlim([0,1]);ylim([0,0.6]);
vline(15/30,'r');
subplot(3,2,3);
shadedErrorBar(x,aveFR_runairpuff_acrosssessions,steFR_runairpuff_acrosssessions,{'color',[0.1373 0.5451 0.2706]});
ylabel('firing rate');
xlim([0,1]); ylim([0 5]);
vline(15/30,'r');
subplot(3,2,5);
shadedErrorBar(x,mean_run_every100ms_across*2*3.1415926*7.5/128,ste_run_every100ms_across*2*3.1415926*7.5/128,{'color',[0.1373 0.5451 0.2706]});
xlabel('time(s)'); ylabel('speed(cm/s)'); 
xlim([0,1]);ylim([-0.5 30]);
vline(15/30,'r');

subplot(3,2,2);
shadedErrorBar(x,ave_dfOvF_stayairpuff_acrossessions,ste_dfOvF_stayairpuff_acrossessions,{'color',[0.1373 0.5451 0.2706]});
title('airpuff stim stay'); %ylabel('df/f');
xlim([0,1]);ylim([0,0.6]);
vline(15/30,'r');
subplot(3,2,4);
shadedErrorBar(x,aveFR_stayairpuff_acrosssessions,steFR_stayairpuff_acrosssessions,{'color',[0.1373 0.5451 0.2706]});
%ylabel('firing rate');
xlim([0,1]); ylim([0 5]);
vline(15/30,'r');
subplot(3,2,6);
shadedErrorBar(x,mean_stay_every100ms_across*2*3.1415926*7.5/128,mean_stay_every100ms_across*2*3.1415926*7.5/128,{'color',[0.1373 0.5451 0.2706]});
xlabel('time(s)'); %ylabel('speed(cm/s)');
xlim([0,1]);ylim([-0.5 30]); vline(15/30,'r');

title('airpuff response across sessions');
savefig(['Z:\Analysis\Airpuff_analysis\across_sessions\' 'across_sessions_puffResponse.fig']);


