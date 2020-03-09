%This script is to plot the PC activity and speed before and after reverse moving dot onset. 
%In order to be included in frms_revstim_stay and frms_revstim_running,
%the bahevioral state can't change before and after the stimulus. 
%(i.e, if mice was running before reverse moving dot, it has to be running as well after)
% this script will also create the variables that are useful:
% frms_revdot_run/frms_revdot_stay: trails*frames
% spk_revdot_stay/run_mat: frames*trials*cells. 
% need a 3D matrix because need to keep the cell by cell variation. the average is the same when do average across trials first then across cells, or average cells first then across trials. However, the std is different. 
% FR_revdotAve_run: average firing rate across cells
% FR_revdot_AveSte_run: cell by cell variation.

%% Section I: set paths and create analysis folders for each session
%define the directory and files
clear;
%NEED TO UPDATE THIS SO IT ACCESSES SPREADSHEET INSTEAD OF JUST WRITING IN THE NAMES
sessions = {'191220_img1042'};
days = {'1042-191220_1'};
%sessionID = {'1023-190923_1','','','','','','','',''};
image_analysis_base  = 'Z:\Analysis\2P_MovingDots_Analysis\imaging_analysis\';%stores the data on crash in the movingDots analysis folder
color_code = {'c','r','y','g'};


%%  sectionII: reverse moving dot triggered response-spikes
for ii = 1: length(sessions)
    %% generate frame matrix based on the criteria, calculate average speed across trials
    behav_dest = ['Z:\Analysis\2P_MovingDots_Analysis\behavioral_analysis\' days{ii}];
    behav_struct = load([behav_dest '\' days{ii} '_behavAnalysis.mat']);
    %behav_struct = load([behav_dest '\' days{ii} '_first18000frames_behavAnalysis.mat']);
    speed = behav_struct.speed; speed = double(speed);
    frames_run_cell = behav_struct.frames_run_cell;
    frames_stay_cell = behav_struct.frames_stay_cell;
    reverse = behav_struct.cReverse_vec;
    
    period = 15; %#of frames before and after reverse moving dot
    
    % generate frame matrices:how does neural activity change with reverse moving dot when the behvaior is the same?
    % frms_revstim_stay: 15 frames (0.5s)stationary before and after reverse moving dot onset(the mice has to be stationary all the time!)
    % frms_revstim_run: !!1!! frames (0.43s)running before and after reverse moving dot onset
    
    frms_revstim_stay = [];
    frms_revstim_run = [];
    for i = 1:length(reverse) 
        if reverse(i) - period < 0 || reverse(i)+period-1 > length(speed)
            continue
        elseif sum(speed((reverse(i)-period):(reverse(i)-1))==0) >= period && sum(speed(reverse(i):(reverse(i)+period-1))==0) >= period
            frms_revstim_stay = cat(1,frms_revstim_stay, reverse(i)- period : reverse(i)+period);% trail*frames
        elseif sum(speed((reverse(i)-6):(reverse(i)-1))~= 0) >= 3 && sum(speed(reverse(i):(reverse(i)+8))~= 0) >= 4
            frms_revstim_run = cat(1,frms_revstim_run, reverse(i)- period : reverse(i)+period);% trial * frames
        end
    end
  
 % generate matrix for speeds   
    spd_revstim_run = speed(frms_revstim_run);
    spd_revstim_stay = speed(frms_revstim_stay);
 % averages and stes  
    spd_revstimAve_run = mean(spd_revstim_run);
    spd_revstimAveSte_run = std(spd_revstim_run)/sqrt(size(spd_revstim_run,1));
    
    spd_revstimAve_stay = mean(spd_revstim_stay);
    spd_revstimAveSte_stay = std(spd_revstim_stay)/sqrt(size(spd_revstim_stay,1));
    
    save([behav_dest '\' days{ii} '_behavAnalysis.mat' ],...
        'frms_revstim_run','frms_revstim_stay','spd_revstim_run',...
        'spd_revstim_stay','spd_revstimAve_run','spd_revstimAve_stay',...
        'spd_revstimAveSte_run','spd_revstimAveSte_stay','-append');
    
    %% spikes -- firing rate
    threshold = -4;
    image_analysis_dest = [image_analysis_base, sessions{ii}, '\'];
    spk_deconv_output = load([image_analysis_dest sessions{ii},'_spk_deconvolve_threshold' num2str(threshold) '.mat']);
    spk_logic = spk_deconv_output.spk_logic_cl;
    %spk_logic_ave = mean(spk_logic,2);
    %FR_ave = spk_logic_ave * 30; %firing rate = spike probablity*imaging rate
    
    % generate matrix for spikes
    spk_revdot_run_mat = zeros(size(frms_revstim_run,1),size(frms_revstim_run,2),size(spk_logic,2)); %trial*frame*cell
    for i = 1: size(spk_logic,2)                                                    %for each cell
        for j = 1:size(frms_revstim_run,2)                                    %for each running trial
            spk_revdot_run_mat (:,j,i) = spk_logic(frms_revstim_run(:,j),i);%spk_revdot_run_mat: trial*frame*cell
        end
    end
    % average across trials:
    spk_revdot_run_cells = squeeze(mean(spk_revdot_run_mat)); %frame*cells
    %average across cells:
    FR_revdotAve_run = mean(spk_revdot_run_cells,2)*30; %firing rate = spike probablity*imaging rate
    FR_revdotAveSte_run = std(spk_revdot_run_cells,0,2)*30/sqrt(size(spk_revdot_run_cells,2)); % cell by cell variation
    
    % generate matrix for spikes
    spk_revdot_stay_mat = zeros(size(frms_revstim_stay,1),size(frms_revstim_stay,2),size(spk_logic,2));
    for i = 1: size(spk_logic,2)                                                    %for each cell
        for j = 1:size(frms_revstim_stay,2)                                    %for each running trial
            spk_revdot_stay_mat (:,j,i) = spk_logic(frms_revstim_stay(:,j),i);%spk_revdot_run_mat: trial*frame*cell
        end
    end
    % average across trials:
    spk_revdot_stay_cells = squeeze(mean(spk_revdot_stay_mat));
    %average across cells:
    FR_revdotAve_stay = mean(spk_revdot_stay_cells,2)*30; %firing rate = spike probablity*imaging rate
    FR_revdotAveSte_stay = std(spk_revdot_stay_cells,0,2)*30/sqrt(size(spk_revdot_stay_cells,2)); % cell by cell variation
    
    save([image_analysis_dest '\' sessions{ii}, '_spk_deconvolve_threshold' num2str(threshold) '.mat' ],...
        'spk_revdot_run_mat','spk_revdot_stay_mat','FR_revdotAve_run',...
        'FR_revdotAveSte_run','FR_revdotAve_stay','FR_revdotAveSte_stay',...
        'spk_revdot_run_cells','spk_revdot_stay_cells','-append');
 %{   
% plot only spike rate and speed
    % figure
    fig = figure;
    % if there's only one trial for running or stationary, plot that single
    % trial
    x = (1: length(spd_airstimAve_stay))/30;
    subplot(2,2,1);
    shadedErrorBar(x,FR_airstimAve_run,FR_airstimAveSte_run);
    title('airpuff stim running');  ylabel('firing rate');vline(period/30,'r');
    xlim([0,1]);ylim([0,4]);
    subplot(2,2,3);
    shadedErrorBar(x,spd_airstimAve_run*2*3.1415926*7.5/128,spd_airstimAveSte_run*2*3.1415926*7.5/128);
    xlabel('time(s)'); ylabel('speed(cm/s)'); vline(period/30,'r');
    xlim([0 1]);%ylim([0 20]);
    text(0.1,min(spd_airstimAve_run*2*3.1415926*7.5/128),['n = ' num2str(size(frms_airstim_run,1))]);
    %xlim([0,31]);
    
    subplot(2,2,2);
    shadedErrorBar(x,FR_airstimAve_stay,FR_airstimAveSte_stay);
    title('airpuff stim stay'); ylabel('firing rate');vline(period/30,'r');
    xlim([0,1]);ylim([0,4]);
    subplot(2,2,4);
    shadedErrorBar(x,spd_airstimAve_stay*2*3.1415926*7.5/128,spd_airstimAveSte_stay*2*3.1415926*7.5/128);
    xlabel('time(s)'); ylabel('speed(cm/s)'); vline(period/30,'r');
    xlim([0 1]);%ylim([-1 1])
    text(0.1,0.05,['n = ' num2str(size(frms_airstim_stay,1))]);
    %xlim([0,31]);
    
    supertitle(sessions{ii});
    
    saveas(fig, [image_analysis_dest '\' sessions{ii} '_airpuffTrigAve_FRnSpd']);    
  %}  
end



%%  sectionIII: reverse moving dot triggered response-df/F
%frm_revdot_run is from the second section of this script, so have to run
%the section section first.
for ii = 1: length(sessions)
   % load data
    image_analysis_dest = [image_analysis_base, sessions{ii}, '\']; 
    dfOvF_strct = load([image_analysis_dest sessions{ii} '_dfOvF.mat']);
    dfOvF_btm = dfOvF_strct.dfOvF_btm;
    dfOvF_stay = dfOvF_strct.dfOvF_stay;
    behav_dest = ['Z:\Analysis\2P_MovingDots_Analysis\behavioral_analysis\' days{ii}];
    behav_struct = load([behav_dest '\' days{ii} '_behavAnalysis.mat']);
    speed = behav_struct.speed; speed = double(speed);
    frms_revstim_stay = behav_struct.frms_revstim_stay;
    frms_revstim_run = behav_struct.frms_revstim_run;
    spd_revstimAve_run = behav_struct.spd_revstimAve_run;
    spd_revstimAveSte_run = behav_struct.spd_revstimAveSte_run;
    spd_revstimAve_stay = behav_struct.spd_revstimAve_stay;
    spd_revstimAveSte_stay = behav_struct.spd_revstimAveSte_stay;
 
 % use bottom 10% of fluorescence:
 % ------------------------------------------------------------------------------------------------------------------------------------------------------
    % generate matrix for dfOvF_btm during running
    dfOvFbtm_revdot_run_mat = zeros(size(frms_revstim_run,1),size(frms_revstim_run,2),size(dfOvF_btm,2)); % trial*frame*cell
    for i = 1: size(dfOvF_btm,2)                                                    %for each cell
        for j = 1:size(frms_revstim_run,1)                                    %for each running trial
            dfOvFbtm_revdot_run_mat(j,:,i) = dfOvF_btm(frms_revstim_run(j,:),i);%spk_revdot_run_mat: trial*frame*cell
        end
    end
    % average across trials:
    dfOvFbtm_revdot_run_cells = squeeze(mean(dfOvFbtm_revdot_run_mat)); %frame*cells
    %average across cells:
    dfOvFbtm_revdotAve_run = mean(dfOvFbtm_revdot_run_cells,2); %firing rate = spike probablity*imaging rate
    dfOvFbtm_revdotAveSte_run = std(dfOvFbtm_revdot_run_cells,0,2)/sqrt(size(dfOvFbtm_revdot_run_cells,2)); % cell by cell variation
    
     % generate matrix for dfOvF_btm during stationary
    dfOvFbtm_revdot_stay_mat = zeros(size(frms_revstim_stay,1),size(frms_revstim_stay,2),size(dfOvF_btm,2));
    for i = 1: size(dfOvF_btm,2)                                                    %for each cell
        for j = 1:size(frms_revstim_stay,2)                                    %for each running trial
            dfOvFbtm_revdot_stay_mat(:,j,i) = dfOvF_btm(frms_revstim_stay(:,j),i);%spk_revdot_run_mat: trial*frame*cell
        end
    end
    % average across trials:
    dfOvFbtm_revdot_stay_cells = squeeze(mean(dfOvFbtm_revdot_stay_mat)); %frame*cells
    %average across cells:
    dfOvFbtm_revdotAve_stay = mean(dfOvFbtm_revdot_stay_cells,2); 
    dfOvFbtm_revdotAveSte_stay = std(dfOvFbtm_revdot_stay_cells,0,2)/sqrt(size(dfOvFbtm_revdot_stay_cells,2)); % cell by cell variation
    
% ------------------------------------------------------------------------------------------------------------------------------------------------------  
  
% use fluorescence during stationary without airpuff
    % generate matrix for dfOvF_stay during running
    dfOvFstay_revdot_run_mat = zeros(size(frms_revstim_run,1),size(frms_revstim_run,2),size(dfOvF_stay,2));
    for i = 1: size(dfOvF_stay,2)                                                    %for each cell
        for j = 1:size(frms_revstim_run,2)                                           %for each running trial
            dfOvFstay_revdot_run_mat(:,j,i) = dfOvF_stay(frms_revstim_run(:,j),i);  %spk_airpuff_run_mat: trial*frame*cell
        end
    end
    % average across trials:
    dfOvFstay_revdot_run_cells = squeeze(mean(dfOvFstay_revdot_run_mat)); %frame*cells
    %average across cells:
    dfOvFstay_revdotAve_run = mean(dfOvFstay_revdot_run_cells,2)*30; %firing rate = spike probablity*imaging rate
    dfOvFstay_revdotAveSte_run = std(dfOvFstay_revdot_run_cells,0,2)*30/sqrt(size(dfOvFstay_revdot_run_cells,2)); % cell by cell variation
    
     % generate matrix for dfOvF_stay during stationary
    dfOvFstay_revdot_stay_mat = zeros(size(frms_revstim_stay,1),size(frms_revstim_stay,2),size(dfOvF_stay,2));
    for i = 1: size(dfOvF_stay,2)                                                    %for each cell
        for j = 1:size(frms_revstim_stay,2)                                    %for each running trial
            dfOvFstay_revdot_stay_mat(:,j,i) = dfOvF_stay(frms_revstim_stay(:,j),i);%spk_airpuff_run_mat: trial*frame*cell
        end
    end
    % average across trials:
    dfOvFstay_revdot_stay_cells = squeeze(mean(dfOvFstay_revdot_stay_mat)); %frame*cells
    %average across cells:
    dfOvFstay_revdotAve_stay = mean(dfOvFstay_revdot_stay_cells,2)*30; %firing rate = spike probablity*imaging rate
    dfOvFstay_revdotAveSte_stay = std(dfOvFstay_revdot_stay_cells,0,2)*30/sqrt(size(dfOvFstay_revdot_stay_cells,2)); % cell by cell variation
    
    save([image_analysis_dest '\' sessions{ii} '_dfOvF.mat' ],...
        'dfOvFstay_revdot_run_mat','dfOvFstay_revdot_stay_mat','dfOvFbtm_revdot_run_mat','dfOvFbtm_revdot_stay_mat',...
        'dfOvFbtm_revdot_run_cells','dfOvFbtm_revdotAve_run','dfOvFbtm_revdotAveSte_run',...
        'dfOvFbtm_revdot_stay_cells','dfOvFbtm_revdotAve_stay','dfOvFbtm_revdotAveSte_stay',...
        'dfOvFstay_revdot_run_cells','dfOvFstay_revdotAve_run','dfOvFstay_revdotAveSte_run',...
        'dfOvFstay_revdot_stay_cells','dfOvFstay_revdotAve_stay','dfOvFstay_revdotAveSte_stay','-append');
    % ------------------------------------------------------------------------------------------------------------------------------------------------------
   %{
    %plot
    period = 15; %#of frames before and after airpuff delivery
    fig = figure;
    x = (1:length(spd_airstimAve_stay))/30;
    subplot(3,2,1);
    shadedErrorBar(x,dfOvFbtm_airstimAve_run,dfOvFbtm_airstimAveSte_run);
    title('airpuff, running, df/F bottom 10%'); ylabel('df/f'); vline(period/30,'r');
    xlim([0,1]);
    subplot(3,2,3);
    shadedErrorBar(x,dfOvFstay_airstimAve_run,dfOvFstay_airstimAveSte_run);
    title('df/f F during stationary'); ylabel('df/f'); vline(period/30,'r');
    xlim([0,1]);
    subplot(3,2,5);
    shadedErrorBar(x,spd_airstimAve_run*2*3.1415926*7.5/128,spd_airstimAveSte_run*2*3.1415926*7.5/128);
    xlabel('time(s)'); ylabel('speed(cm/s)'); vline(period/30,'r');
    xlim([0 1]);
    text(0.1,min(spd_airstimAve_run*2*3.1415926*7.5/128),['n = ' num2str(size(frms_airstim_run,1))]);   
    
    subplot(3,2,2);
    shadedErrorBar(x,dfOvFbtm_airstimAve_stay,dfOvFbtm_airstimAveSte_stay);
    title('airpuff, stationary, df/F bottom 10%'); ylabel('df/f'); vline(period/30,'r');
    xlim([0,1]);
    subplot(3,2,4);
    shadedErrorBar(x,dfOvFstay_airstimAve_stay,dfOvFstay_airstimAveSte_stay);
    title('df/f F during stationary');ylabel('df/f'); vline(period/30,'r');
    xlim([0,1]);
    subplot(3,2,6);
    shadedErrorBar(x,spd_airstimAve_stay*2*3.1415926*7.5/128,spd_airstimAveSte_stay*2*3.1415926*7.5/128);
    xlabel('time(s)'); ylabel('speed(cm/s)'); vline(period/30,'r');
    xlim([0 1]);
    text(0.1,min(spd_airstimAve_stay*2*3.1415926*7.5/128),['n = ' num2str(size(frms_airstim_stay,1))]);
    
    supertitle(sessions{ii});
    saveas(fig, [image_analysis_dest '\' sessions{ii} '_airpuffTrigAve_dfOvFnSpd']);   
    %}
end

    %% plot
    % figure
    fig = figure;
    % if there's only one trial for running or stationary, plot that single
    % trial
    
    x = (1: length(spd_revstimAve_stay))/30;
    subplot(3,2,1);
    shadedErrorBar(x,dfOvFbtm_revdotAve_run,dfOvFbtm_revdotAveSte_run);
    title('reverse moving dot running'); ylabel('df/f');
    xlim([0,1]);ylim([0,0.6]);
    vline(period/30,'r');
    subplot(3,2,3);
    shadedErrorBar(x,FR_revdotAve_run,FR_revdotAveSte_run);
    ylabel('firing rate');
    xlim([0,1]); %ylim([0 4]);
    vline(period/30,'r');
    subplot(3,2,5);
    shadedErrorBar(x,spd_revstimAve_run*2*3.1415926*7.5/128,spd_revstimAveSte_run*2*3.1415926*7.5/128);
    xlabel('time(s)'); ylabel('speed(miles/hr)'); vline(period/30,'r');
    text(0.1,min(spd_revstimAve_run),['n = ' num2str(size(frms_revstim_run,1))]);
    xlim([0,1]);ylim([-0.5 30]);
    
    subplot(3,2,2);
    shadedErrorBar(x,dfOvFbtm_revdotAve_stay,dfOvFbtm_revdotAveSte_stay);
    title('reverse moving dot stim stay'); %ylabel('df/f');
    xlim([0,1]);ylim([0,0.6]);
    vline(period/30,'r');
    subplot(3,2,4);
    shadedErrorBar(x,FR_revdotAve_stay,FR_revdotAveSte_stay);
    %ylabel('firing rate');
    xlim([0,1]); %ylim([0 4]);
    vline(period/30,'r');
    subplot(3,2,6);
    shadedErrorBar(x,spd_revstimAve_stay*2*3.1415926*7.5/128,spd_revstimAveSte_stay*2*3.1415926*7.5/128);
    xlabel('time(s)'); %ylabel('speed(cm/s)'); 
    text(0.1,5,['n = ' num2str(size(frms_revstim_stay,1))]);
    xlim([0,1]);ylim([-0.5 30]); vline(period/30,'r');
    
    suptitle(sessions{ii});
    
    saveas(fig, [image_analysis_dest '\' sessions{ii} '_reverseTrigAve']);
    
