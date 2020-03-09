% Continuing after "airpuff_dfOvFnFRnSpd_alingonset", 
% this script is to plot the PC activity and speed before and after airpuff onset during running,
% with speed catogerized into two categories: slow and fast. 

%% Section I: set paths and create analysis folders for each session
%define the directory and files
clear;
%NEED TO UPDATE THIS SO IT ACCESSES SPREADSHEET INSTEAD OF JUST WRITING IN THE NAMES
sessions = {'190912_img1034'};
days = {'1034-190912_1'};
%sessionID = {'1023-190923_1','','','','','','','',''};
image_analysis_base  = 'Z:\Analysis\Airpuff_analysis\imaging_analysis\';%stores the data on crash in the movingDots analysis folder

color_code = {'c','r','y','g'};

%% Section II: airpuff triggered response - response during slow running and during fast running
for ii = 1: length(sessions)
    behav_dest = ['Z:\Analysis\Airpuff_analysis\behavioral_analysis\' days{ii}];
    behav_struct = load([behav_dest '\' days{ii} '_behavAnalysis.mat']);
    %behav_struct = load([behav_dest '\' days{ii} '_first18000frames_behavAnalysis.mat']);
    
    speed = behav_struct.speed; speed = double(speed);
    spd_airstim_run = behav_struct.spd_airstim_run;
    frms_airstim_run = behav_struct.frms_airstim_run;    
    
    % categorize speed into two : slow and fast
    % slow: average speed of whole trial <= median
    % fast: average speed of whole trial > median
    avespd_trials = mean(spd_airstim_run,2);
    medianspd = median(avespd_trials);
    trials = 1: length(avespd_trials);
    
    fast_trials = trials(avespd_trials > medianspd);
    slow_trials = trials(avespd_trials <= medianspd);
    
    frms_airstim_run_fast = frms_airstim_run(fast_trials,:);
    frms_airstim_run_slow = frms_airstim_run(slow_trials,:);
    
    spd_airstim_run_fast = spd_airstim_run(fast_trials,:);
    spd_airstim_run_slow = spd_airstim_run(slow_trials,:);
    
    spd_airstimAve_run_fast = mean(spd_airstim_run_fast);
    spd_airstimSte_run_fast = std(spd_airstim_run_fast)/sqrt(size(spd_airstim_run_fast,1));
    
    spd_airstimAve_run_slow = mean(spd_airstim_run_slow);
    spd_airstimSte_run_slow = std(spd_airstim_run_slow)/sqrt(size(spd_airstim_run_slow,1));
    save([behav_dest '\' days{ii} '_behavAnalysis.mat' ],...
        'spd_airstimAve_run_fast','spd_airstimSte_run_fast','spd_airstimAve_run_slow',...
        'spd_airstimSte_run_slow','-append');
    
    % spike 
    threshold = -4;
    image_analysis_dest = [image_analysis_base, sessions{ii}, '\'];
    spk_deconv_output = load([image_analysis_dest sessions{ii},'_spk_deconvolve_threshold' num2str(threshold) '.mat']);
    spk_airpuff_run_mat = spk_deconv_output.spk_airpuff_run_mat; %trial*frame*cell
    
    spk_airpuff_fastrun_mat = spk_airpuff_run_mat(fast_trials,:,:);
    spk_airpuff_slowrun_mat = spk_airpuff_run_mat(slow_trials,:,:);
    %average across trials
    spk_airpuff_fastrun_cells = squeeze(mean(spk_airpuff_fastrun_mat));
    spk_airpuff_slowrun_cells = squeeze(mean(spk_airpuff_slowrun_mat));

    %average across cells
    FR_airpuff_fastrun_ave = mean(spk_airpuff_fastrun_cells,2)*30;
    FR_airpuff_fastrun_ste = std(spk_airpuff_fastrun_cells,0,2)*30/sqrt(size(spk_airpuff_fastrun_cells,2));
    
    FR_airpuff_slowrun_ave = mean(spk_airpuff_slowrun_cells,2)*30;
    FR_airpuff_slowrun_ste = std(spk_airpuff_slowrun_cells,0,2)*30/sqrt(size(spk_airpuff_slowrun_cells,2));
    
    %dfOvF ---- using bottom 10%
    dfOvF_strct = load([image_analysis_dest sessions{ii} '_dfOvF.mat']);
    dfOvFbtm_airpuff_run_mat = dfOvF_strct.dfOvFbtm_airpuff_run_mat;%trial*frame*cell
    
    dfOvFbtm_airpuff_fastrun_mat = dfOvFbtm_airpuff_run_mat(fast_trials,:,:);    
    dfOvFbtm_airpuff_slowrun_mat = dfOvFbtm_airpuff_run_mat(slow_trials,:,:);   
    %average across trials
    dfOvFbtm_airpuff_fastrun_cells = squeeze(mean(dfOvFbtm_airpuff_fastrun_mat));
    dfOvFbtm_airpuff_slowrun_cells = squeeze(mean(dfOvFbtm_airpuff_slowrun_mat));

    %average across cells
    dfOvFbtm_airpuff_fastrun_ave = mean(dfOvFbtm_airpuff_fastrun_cells,2);
    dfOvFbtm_airpuff_fastrun_ste = std(dfOvFbtm_airpuff_fastrun_cells,0,2)/sqrt(size(dfOvFbtm_airpuff_fastrun_cells,2));
    
    dfOvFbtm_airpuff_slowrun_ave = mean(dfOvFbtm_airpuff_slowrun_cells,2);
    dfOvFbtm_airpuff_slowrun_ste = std(dfOvFbtm_airpuff_slowrun_cells,0,2)/sqrt(size(dfOvFbtm_airpuff_slowrun_cells,2));
    
    figure;
    x = (1:length(spd_airstimAve_run_fast))/30;
    subplot(3,2,1);
    shadedErrorBar(x,FR_airpuff_slowrun_ave,FR_airpuff_slowrun_ste);
    title([sessions{ii}, '  slow']); ylabel('Firing rate'); 
    xlim([0,1]); ylim([0 4]);vline(0.5,'r');
    subplot(3,2,3);
    shadedErrorBar(x,dfOvFbtm_airpuff_slowrun_ave,dfOvFbtm_airpuff_slowrun_ste);
    ylabel('df/F');
    xlim([0,1]); ylim([0 0.4]);vline(0.5,'r');
    subplot(3,2,5);
    shadedErrorBar(x,spd_airstimAve_run_slow*2*3.1415926*7.5/128,spd_airstimSte_run_slow*2*3.1415926*7.5/128);
    ylabel('cm/s'); 
    ylim([0 25]);vline(0.5,'r');
    xlabel('time(s)');
    text(0.1,3,['n = ' num2str(size(spd_airstim_run_slow,1))]);
    
    subplot(3,2,2);
    shadedErrorBar(x,FR_airpuff_fastrun_ave,FR_airpuff_fastrun_ste);
    title('fast'); ylabel('Firing rate'); 
    xlim([0,1]);ylim([0 4]);vline(0.5,'r');
    subplot(3,2,4);
    shadedErrorBar(x,dfOvFbtm_airpuff_fastrun_ave,dfOvFbtm_airpuff_fastrun_ste);
    ylabel('df/F');
    xlim([0,1]); ylim([0 0.4]);vline(0.5,'r');
    subplot(3,2,6);
    shadedErrorBar(x,spd_airstimAve_run_fast*2*3.1415926*7.5/128,spd_airstimSte_run_fast*2*3.1415926*7.5/128);
    ylabel('cm/s'); 
    ylim([0 25]);vline(0.5,'r');
    xlabel('time(s)');
    text(0.1,3,['n = ' num2str(size(spd_airstim_run_fast,1))]);
    %save figure

    
end