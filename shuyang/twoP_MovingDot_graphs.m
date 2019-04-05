%% 2 photon: moving dot ---- spikes
%calculate spikes(1st derivatives of fluorescence) 
%plot average spike rate for all cells in one trial, aligned with running onset and offset
%scatter plot running vs. still (each dot should be a cell)

%% assign document paths and experimental sessions
clear;
sessions = '190314_imgJ89'; 
image_analysis_base    = 'Z:\Analysis\2P_MovingDots_Analysis\imaging_analysis\'; %stores the data on crash in the movingDots analysis folder
image_analysis_dest = [image_analysis_base, sessions, '\'];

% behavior analysis results 
days = 'J89-190314_1';
behav_dest = ['Z:\Analysis\2P_MovingDots_Analysis\behavioral_analysis\' days '\'];
color_code = {'b','r','k','c'};

%% df/f
filename = dir([image_analysis_dest '*' '_TCave.mat']);
TCave = load([image_analysis_dest filename.name]);
TCave = TCave.tc_avg;
behav_output = load([behav_dest days '_behavAnalysis.mat']);
frm_stay_cell = behav_output.frames_stay_cell;
frm_stay = cell2mat(frm_stay_cell);
baseline = mean(TCave(frm_stay,:));         % baseline for each cell
dfOvF = zeros(size(TCave,1),size(TCave,2));
for i = 1:size(TCave,2)                     % for each cell
    dfOvF(:,i) = (TCave(:,i)-baseline(i))/baseline(i); %frames*cell
end


%% run triggered ave and run offset ave --- df/f
frames_run_cell = behav_output.frames_run_cell;
speed = double(behav_output.speed);
[~,~,frames_runTrigger_mat,frms_runoff_mat,~]= findFrames_runWindows_2P (speed,frames_run_cell);
% aligned with running onset ----------------------------------------------------------------------------
dfOvF_runtrig_mat = zeros(size(frames_runTrigger_mat,1),size(frames_runTrigger_mat,2),size(TCave,2));
for i = 1: size(frames_runTrigger_mat,2)                                    % for each running trig window
    dfOvF_runtrig_mat(:,i,:) = dfOvF(frames_runTrigger_mat(:,i),:);         % df/f for all cells during that window, frame*cells 
end
ave_dfOvF_runtrig_trials = squeeze(mean(dfOvF_runtrig_mat,2));
ave_dfOvF_runtrig = squeeze(mean(ave_dfOvF_runtrig_trials,2));
ste_dfOvF_runTrig = std(ave_dfOvF_runtrig,0,2)/sqrt(size(ave_dfOvF_runtrig,2)); % a bunch of zeros
% speed
speed_runtrig = speed(frames_runTrigger_mat);
aveSpd_runtrig = mean(speed_runtrig,2);
steSpd_runtrig = std(speed_runtrig,0,2)/sqrt(size(speed_runtrig,2));
%plot
x = (1:45);
figure; subplot(2,1,1);
shadedErrorBar(x,ave_dfOvF_runtrig,ste_dfOvF_runTrig);hold on;
vline(16,'r'); ylabel('df/f');
subplot(2,1,2); 
shadedErrorBar(x,aveSpd_runtrig, steSpd_runtrig);hold on;
vline(16,'r');ylabel('speed');
supertitle([sessions ' running triggered average']);
savefig([image_analysis_dest sessions '_runTrigAve_dfOvF_wspeed']);

% aligned with running offset ------------------------------------------------------------------------
dfOvF_runoff_mat = zeros(size(frms_runoff_mat,1),size(frms_runoff_mat,2),size(TCave,2));
for i = 1: size(frms_runoff_mat,2)                                    % for each running trig window
    dfOvF_runoff_mat(:,i,:) = dfOvF(frms_runoff_mat(:,i),:);         % df/f for all cells during that window, frame*cells 
end
ave_dfOvF_runoff_trials = squeeze(mean(dfOvF_runoff_mat,2));
ave_dfOvF_runoff = squeeze(mean(ave_dfOvF_runoff_trials,2));
ste_dfOvF_runoff = std(ave_dfOvF_runoff,0,2)/sqrt(size(ave_dfOvF_runoff,2)); % a bunch of zeros
% speed
speed_runoff = speed(frms_runoff_mat);
aveSpd_runoff = mean(speed_runoff,2);
steSpd_runoff = std(speed_runoff,0,2)/sqrt(size(speed_runoff,2));
%plot
x = (1:45);
figure; subplot(2,1,1);
shadedErrorBar(x,ave_dfOvF_runoff,ste_dfOvF_runoff);hold on;
vline(31,'r'); ylabel('df/f');
subplot(2,1,2); 
shadedErrorBar(x,aveSpd_runoff, steSpd_runoff);hold on;
vline(31,'r');ylabel('speed');
supertitle([sessions ' aligned to running offset -  average']);
savefig([image_analysis_dest sessions '_runTrigAve_dfOvF_wspeed']);

%% SECTION I 
% get derivatives from raw fluorescence, decide the threshold for spikes,  
% get number of spikes during stationary 
filename = dir([image_analysis_dest '*' '_TCave.mat']);
TCave = load([image_analysis_dest filename.name]);
TCave = TCave.tc_avg;
deriv = diff(TCave,1,1);
std_deriv = std(deriv);
std2 = 2*std_deriv;
std3 = 3*std_deriv;
std2_5 = 2.5*std_deriv;
behav_output = load([behav_dest days '_behavAnalysis.mat']);
frm_run_cell = behav_output.frames_run_cell;
frm_stay_cell = behav_output.frames_stay_cell;

% get the best threshold std for deciding spikes, ave FRs during stationary for each cell 
[avesWvarystd,best_thres,aveFR_allcells_stay,std_best] = best_spkthreshold_2P (deriv, frm_stay_cell,std2,std2_5,std3,std_deriv); 
aveFR_stay = mean(aveFR_allcells_stay);

%% SECTION II 
% calculate ave FRs during running for each cell and make scatter plot for
% spikes rates during running vs. still, each dot is a cell
nspk = zeros(1,size(deriv,2));                                              % num of cells
frate = zeros(size(frm_run_cell,2),size(deriv,2));                               % num of stationary windows*num of cells
for i = 1: size(frm_run_cell,2)                                                  % for each running window
    t = length(frm_run_cell{i})/30;                                              %t = nframes*1/30s, thus the unit of t is second
    der_run_temp = deriv(frm_run_cell{i},:); 
    for j = 1: size(der_run_temp,2)                                         % for each cell
        nspk(j) = sum(der_run_temp(:,j) >= std_best(j));                    % number of spikes during this window for each cell
    end
    frate(i,:) = nspk/t;                                                    % 2 dimensions: window index*cell
end
aveFR_allcells_run = mean(frate);                                           %ave firing rate for all cells across all running windows (1*cells)
aveFR_run = mean(aveFR_allcells_run);
figure;
scatter(aveFR_allcells_stay,aveFR_allcells_run,'filled','r'); hold on;
scatter(aveFR_stay, aveFR_run,'filled','k');hold on;
xlabel('stationary'); ylabel('running');
xlim([0,1.5]); ylim([0,1.5]);
refline(1,0);
title('mean firing rate per second');
savefig([image_analysis_dest sessions '_meanFRstayVSrun']);


%% SECTION III : run trigger average and run offset average--- spikes
behav_output = load([behav_dest days '_behavAnalysis.mat']);
frames_run_cell = behav_output.frames_run_cell;
speed = double(behav_output.speed);
[~,~,frames_runTrigger_mat,frms_runoff_mat,~]= findFrames_runWindows_2P (speed,frames_run_cell);
% aligned with running onset ----------------------------------------------------------------------------
spkPerFrm = zeros(size(frames_runTrigger_mat,1),size(deriv,2));             % nframes*num of cells
aveSpkTrig = zeros(size(frames_runTrigger_mat,1), size(frames_runTrigger_mat,2));
for i = 1: size(frames_runTrigger_mat,2)                                    % for each running trig window
    der_runtrig_temp = deriv(frames_runTrigger_mat(:,i),:);                 % deriv values for all cells during that window, frame*cells
    for j = 1: size(der_runtrig_temp,2)                                     % for each cell
        spkPerFrm(:,j) = der_runtrig_temp(:,j) >= std_best(j);              % binary spikes rates for each frame of each cell (0 or 1)
    end
    aveSpkTrig(:,i) = mean(spkPerFrm,2)*30;                                 % spike probablity for each frame: ave across all cells *frame rate(30)
                                                                            % frames * runtrig windows
end
aveSpkTrig_all = mean(aveSpkTrig,2);                                        % ave across all running windows
steSpkTrig_all = std(aveSpkTrig,0,2)/sqrt(size(aveSpkTrig,2));
x = (1:45);
figure; shadedErrorBar(x,aveSpkTrig_all,steSpkTrig_all);hold on;
%ylim([0,2]); hold on; 
vline(16,'r'); ylabel('firing rate');
title([sessions ' running triggered average firing rate']);
savefig([image_analysis_dest sessions '_runTrigAve_FR']);
% PLOT speed
speed_runtrig = speed(frames_runTrigger_mat);
aveSpd_runtrig = mean(speed_runtrig,2);
steSpd_runtrig = std(speed_runtrig,0,2)/sqrt(size(speed_runtrig,2));
figure; shadedErrorBar(x,aveSpd_runtrig, steSpd_runtrig);hold on;
vline(16,'r'); title([sessions ' running triggered average speed']);
savefig([image_analysis_dest sessions '_runTrigAve_speed']);

% aligned with running offset ------------------------------------------------------------------------------
spkPerFrm_off = zeros(size(frms_runoff_mat,1),size(deriv,2));             % nframes*num of cells
aveSpkoff = zeros(size(frms_runoff_mat,1), size(frms_runoff_mat,2));
for i = 1: size(frms_runoff_mat,2)                                    % for each running trig window
    der_runtrig_temp = deriv(frms_runoff_mat(:,i),:);                 % deriv values for all cells during that window, frame*cells
    for j = 1: size(der_runtrig_temp,2)                                     % for each cell
        spkPerFrm_off(:,j) = der_runtrig_temp(:,j) >= std_best(j);              % binary spikes rates for each frame of each cell (0 or 1)
    end
    aveSpkoff(:,i) = mean(spkPerFrm_off,2)*30;                                 % spike probablity for each frame: ave across all cells *frame rate(30)
                                                                            % frames * runtrig windows
end
aveSpkoff_all = mean(aveSpkoff,2);                                        % ave across all running windows
steSpkoff_all = std(aveSpkoff,0,2)/sqrt(size(aveSpkoff,2));
x = (1:45);
figure; shadedErrorBar(x,aveSpkoff_all,steSpkoff_all);hold on;
%ylim([0,2]); hold on; 
vline(31,'r'); ylabel ('firing rate');
title([sessions ' aligned to running offset']);
savefig([image_analysis_dest sessions '_runOff_FR']);
% PLOT speed
speed_runoff = speed(frms_runoff_mat);
aveSpd_runoff = mean(speed_runoff,2);
steSpd_runoff = std(speed_runoff,0,2)/sqrt(size(speed_runoff,2));
figure; shadedErrorBar(x,aveSpd_runoff, steSpd_runoff);hold on;
vline(31,'r'); title([sessions ' aligned to offset']);
savefig([image_analysis_dest sessions '_runoff_speed']);


%% plots to go back and look at the raw data

% plot first derivatives with horizontal std lines --------------------------------------------------
listfrm = [1:1000; 3001:4000; 6001:7000; 10001:11000; 15001:16000; 17001:18000; ...
  21001:22000; 24001:25000; 28001:29000];
listfrm = listfrm';
[fig_deriv] = GUI_w_hline(deriv, listfrm,std_deriv,std2,std2_5,std3);
savefig([image_analysis_dest sessions '_GUIraw_fluo_stdlines']);

% after function best_spkthreshold_2P, go back to raw fluorescence and see if the time points of spikes make sense
% plot red dots over rawTC when deriv > std_bestv for all cells and several frames---------------------------------------------------
listfrm = [1:1000; 3001:4000; 6001:7000; 10001:11000; 15001:16000; 17001:18000; ...
  21001:22000; 24001:25000; 28001:29000];
listfrm = listfrm';
[fig] = GUI(TCave, deriv, listfrm,std2_5);  !!!!!!!!!! you have to go into the function and change the super title depends on what std you give.
savefig([image_analysis_dest sessions '_GUIraw_fluo_w_spkEvent_std2_5']); !!!!!! remember to change the file name depends on the std you give

% Plot raw fluorescence for the running triggered average windows
[fig_runtrig] = GUI_w_vline(TCave, deriv, frames_runTrigger_mat,std_best); !!!!! you will get an error if you have less than 9 running windows
savefig([image_analysis_dest sessions '_GUIraw_fluo_RUNtrig']);

% -------------------------------------------NOT USEFUL--------------plot all of the first derivatives: not very useful ------------------------------------------------------------
%make a for loop to plot the shit so that the lines are not on top of each other
colord=[         0         0    1.0000
    0    0.4000         0
    1.0000         0         0
    0    0.7500    0.7500
    0.7500         0    0.7500
    0.8, 0.5, 0
    0         0    0.5
    0         0.85      0];

shift  = 0;
fig = figure;
for j = 1:size(deriv,2)
    plot(deriv(:,j)+shift,'Color',colord(mod(j-1,size(colord,1))+1,:)); hold on
    %     if ismember(i, [5,10,15,20,25,30,35,40,45,50])
    %         hline(shift);
    %     end
    shift = shift+10000;  %10000 for 2P
end
set(gca,'YTick',10000:10000:shift);
set(gca,'YTicklabel',(1:size(deriv,2)));

listfrm = [1:1000; 3001:4000; 6001:7000; 10001:11000; 15001:16000; 17001:18000; ...
  21001:22000; 24001:25000; 28001:29000];

%you're gonna get a shit ton of plots if you run this for loop, it plots first derivatives of randomly selected 10 cells in some random time periods
for i = randi(size(deriv,2),1,10)
    figure;
    for j = 1:size(listfrm,1)
        subplot(3,3,j);plot(deriv(listfrm(j,:),i));% hold on;
        hline(std_deriv(i),'r'); hold on;
        hline(std2(i),'g'); hold on;
        hline(std3(i),'k'); hold on;
        set(gca,'xticklabel',[],'yticklabel',[]);
        title(['cell' num2str(i) 'frm' num2str(listfrm(j,1)) '-' num2str(listfrm(j,1000))]);
    end
    waitforbuttonpress;
end

figure;
plot(deriv(1:1000,1)); hold on;
hline(std_deriv(1),'r'); hold on;
hline(2*std_deriv(1),'g'); hold on;
hline(3*std_deriv(1),'k'); hold on;
title(['cell1 frm 1-1000']);

