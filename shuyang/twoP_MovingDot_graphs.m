%% 2 photon: moving dot ---- spikes
%calculate spikes(1st derivatives of fluorescence) 
%plot average spike rate for all cells in one trial, aligned with running onset and offset
%scatter plot running vs. still (each dot should be a cell)

%% assign document paths and experimental sessions
clear;
sessions = '190117_img1016'; 
image_analysis_base    = 'Z:\Analysis\2P_MovingDots_Analysis\imaging_analysis\'; %stores the data on crash in the movingDots analysis folder
image_analysis_dest = [image_analysis_base, sessions, '\'];

% behavior analysis results 
days = '1016-190117_1';
behav_dest = ['Z:\Analysis\2P_MovingDots_Analysis\behavioral_analysis\' days '\'];
color_code = {'b','r','k','c'};

%% SECTION I df/f
filename = dir([image_analysis_dest 'getTC\' '*' '_TCave.mat']);
TCave = load([image_analysis_dest 'getTC\' filename.name]);
TCave = TCave.tc_avg;
behav_output = load([behav_dest days '_behavAnalysis.mat']);
speed = double(behav_output.speed);
frm_stay_cell = behav_output.frames_stay_cell;
frm_stay = cell2mat(frm_stay_cell);
baseline = mean(TCave(frm_stay,:));         % baseline for each cell
dfOvF = zeros(size(TCave,1),size(TCave,2));
for i = 1:size(TCave,2)                     % for each cell
    dfOvF(:,i) = (TCave(:,i)-baseline(i))/baseline(i); %frames*cell
end
% plot df/f over speed for the first 5000 frames (choose a period when you have both running and stationary windows)
avedfOvF = mean(dfOvF,2);
stedfOvF = std(dfOvF,0,2)/sqrt(size(dfOvF,2));
frame = (0:(length(speed)-1));
figure; n1 = 500; n2 = 9000;
[hAx,hline1,hline2] = plotyy(frame(n1:n2),speed(n1:n2),frame(n1:n2),avedfOvF(n1:n2)');
%[hAx,hline1,hline2] = plotyy(frame,speed,frame,(avedfOvF(1:length(speed))'));

%set(hline1,'color', 'b'); 
%set(hline2,'color',color_code{n});
legend('speed','df/f','Location','northeast');
ylim(hAx(1),[-30 80]);
ylim(hAx(2),[-0.4,0.4]);
title(['speed and average df/f' days]);
savefig([image_analysis_dest sessions '_avedfOvF_wspeed']);


%% SECTION II run triggered ave and run offset ave --- df/f
filename = dir([image_analysis_dest 'getTC\' '*' '_TCave.mat']);
TCave = load([image_analysis_dest 'getTC\' filename.name]);
TCave = TCave.tc_avg;
behav_output = load([behav_dest days '_behavAnalysis.mat']);
frames_run_cell = behav_output.frames_run_cell;
speed = double(behav_output.speed);
[~,~,frms_runTrig_mat,frms_runoff_mat,~]= findFrames_runWindows_2P_longer (speed,frames_run_cell);
% aligned with running onset ----------------------------------------------------------------------------
dfOvF_runtrig_mat = zeros(size(frms_runTrig_mat,1),size(frms_runTrig_mat,2),size(TCave,2));
for i = 1: size(frms_runTrig_mat,2)                                    % for each running trig window
    dfOvF_runtrig_mat(:,i,:) = dfOvF(frms_runTrig_mat(:,i),:);         % df/f for all cells during that window, frame*cells 
end
ave_dfOvF_runtrig_trials = squeeze(mean(dfOvF_runtrig_mat,2));
ave_dfOvF_runtrig = squeeze(mean(ave_dfOvF_runtrig_trials,2));
ste_dfOvF_runTrig = std(ave_dfOvF_runtrig_trials,0,2)/sqrt(size(ave_dfOvF_runtrig_trials,2)); % a bunch of zeros
% speed
speed_runtrig = speed(frms_runTrig_mat);
aveSpd_runtrig = mean(speed_runtrig,2);
steSpd_runtrig = std(speed_runtrig,0,2)/sqrt(size(speed_runtrig,2));
%plot
x = (1:size(frms_runTrig_mat,1));
figure; subplot(2,1,1);
shadedErrorBar(x,ave_dfOvF_runtrig,ste_dfOvF_runTrig);hold on;
vline(31,'r'); ylabel('df/f');
subplot(2,1,2); 
shadedErrorBar(x,aveSpd_runtrig, steSpd_runtrig);hold on;
vline(31,'r');ylabel('speed');
supertitle([sessions ' running triggered average']);
savefig([image_analysis_dest sessions '_runTrigAve_dfOvF_wspeed']);

% aligned with running offset ------------------------------------------------------------------------
dfOvF_runoff_mat = zeros(size(frms_runoff_mat,1),size(frms_runoff_mat,2),size(TCave,2));
for i = 1: size(frms_runoff_mat,2)                                    % for each running trig window
    dfOvF_runoff_mat(:,i,:) = dfOvF(frms_runoff_mat(:,i),:);         % df/f for all cells during that window, frame*cells 
end
ave_dfOvF_runoff_trials = squeeze(mean(dfOvF_runoff_mat,2));
ave_dfOvF_runoff = squeeze(mean(ave_dfOvF_runoff_trials,2));
ste_dfOvF_runoff = std(ave_dfOvF_runoff_trials,0,2)/sqrt(size(ave_dfOvF_runoff_trials,2)); % a bunch of zeros
% speed
speed_runoff = speed(frms_runoff_mat);
aveSpd_runoff = mean(speed_runoff,2);
steSpd_runoff = std(speed_runoff,0,2)/sqrt(size(speed_runoff,2));
%plot
x = (1:size(frms_runoff_mat,1));
figure; subplot(2,1,1);
shadedErrorBar(x,ave_dfOvF_runoff,ste_dfOvF_runoff);hold on;
vline(31,'r'); ylabel('df/f'); xlabel('frame');
subplot(2,1,2); 
shadedErrorBar(x,aveSpd_runoff, steSpd_runoff);hold on;
vline(31,'r');ylabel('speed'); xlabel('frame');
supertitle([sessions ' aligned to running offset -  average']);
savefig([image_analysis_dest sessions '_runoffset_dfOvF_wspeed']);


%% SECTION III  SPIKE RATES
% get derivatives from raw fluorescence
filename = dir([image_analysis_dest 'getTC\' '*' '_TCave.mat']);
TCave = load([image_analysis_dest 'getTC\' filename.name]);
TCave = TCave.tc_avg;
deriv = diff(TCave,1,1); % derivatives shoud have positive and negative values
behav_output = load([behav_dest days '_behavAnalysis.mat']);
frm_stay_cell = behav_output.frames_stay_cell;
frm_stay = cell2mat(frm_stay_cell);
%  decide the threshold for spikes,  get average number of spikes during stationary (bestFR)
[aveFRsWstds, best_thres,bestFR,std_best,spk_inx,FRstay_cells,...
    spk_bi_cellmat] = twoP_best_spkthreshold (deriv, frm_stay,TCave);
savefig([image_analysis_dest sessions '_hist_stayFR_wdiffthresh']);
save([image_analysis_dest sessions,'_spikes.mat'],'aveFRsWstds','best_thres'...
    ,'bestFR','std_best','spk_inx','spk_bi_cellmat');

%% SECTION IV FR RUNNING VS. STAYTIONARY
behav_output = load([behav_dest days '_behavAnalysis.mat']);
frm_run_cell = behav_output.frames_run_cell;
frm_run = cell2mat(frm_run_cell);
spkrun = zeros(1,size(deriv,2));
for i = 1: size(deriv,2) %for each cell
%count how many spikes are in the running periods: get the intersection
%of indices of peaks above the threshold and indices of frames for running
spkrun(i) = length(intersect(spk_inx{i},frm_run));
end 
t_run = length(frm_run)/30;
FRrun_cells = spkrun/t_run;
aveFR_run = mean(FRrun_cells);
figure;
scatter(FRstay_cells,FRrun_cells,'filled','r'); hold on;
scatter(bestFR, aveFR_run,'filled','k');hold on;
xlabel('stationary'); ylabel('running');
xlim([0.5,1.5]); ylim([0.5,1.5]);
refline(1,0);
title('mean firing rate per second');
savefig([image_analysis_dest sessions '_meanFRstayVSrun']);

%% SECTION V : run trigger average and run offset average--- spikes rates
behav_output = load([behav_dest days '_behavAnalysis.mat']);
frames_run_cell = behav_output.frames_run_cell;
speed = double(behav_output.speed);
[~,~,frms_runTrig_mat,frms_runoff_mat,~]= findFrames_runWindows_2P_longer (speed,frames_run_cell);
% aligned with running onset ----------------------------------------------------------------------------
%get a matrix of 0s and 1s, frames*windows*cells
spk_runtrig_mat = zeros(size(frms_runTrig_mat,1),size(frms_runTrig_mat,2),size(TCave,2));
for i = 1: size(TCave,2)                                                    %for each cell
    for j = 1:size(frms_runTrig_mat,2)                                    %for each running window
        spk_runtrig_mat (:,j,i) = spk_bi_cellmat(frms_runTrig_mat(:,j),i);%spk_runtrig_mat: frame*window*cell
    end
end
%calculate mean and ste
spk_runtrig_cells = squeeze(mean(spk_runtrig_mat,2))*30;                    % average across windows, *30: firing rate
spk_runtrig_grand = squeeze(mean(spk_runtrig_cells,2));                     % average across cells
ste_runtrig = std(spk_runtrig_cells,0,2)/sqrt(size(spk_runtrig_cells,2));
% speed
speed_runtrig = speed(frms_runTrig_mat);
aveSpd_runtrig = mean(speed_runtrig,2);
steSpd_runtrig = std(speed_runtrig,0,2)/sqrt(size(speed_runtrig,2));
%plot
x = (1:size(frms_runTrig_mat,1));
figure; subplot(2,1,1);
shadedErrorBar(x,spk_runtrig_grand,ste_runtrig);hold on;
ylim([0,2]); hold on; 
vline(31,'r'); ylabel('average firing rate'); xlabel('frame');
title([sessions ' aligned to running onset']);
subplot(2,1,2);shadedErrorBar(x,aveSpd_runtrig, steSpd_runtrig);hold on;
vline(31,'r'); ylabel('speed'); xlabel('frame');
%title([sessions ' running triggered average speed']);
savefig([image_analysis_dest sessions '_runTrigAve_FRwSpeed']);


% aligned with running offset ------------------------------------------------------------------------------
spk_runoff_mat = zeros(size(frms_runoff_mat,1),size(frms_runoff_mat,2),size(TCave,2));
for i = 1: size(TCave,2)                                                    %for each cell
    for j = 1:size(frms_runoff_mat,2)                                    %for each running window
        spk_runoff_mat (:,j,i) = spk_bi_cellmat(frms_runoff_mat(:,j),i);%spk_runtrig_mat: frame*window*cell
    end
end
spk_runoff_cells = squeeze(mean(spk_runoff_mat,2))*30;                    % average across windows, *30: firing rate
spk_runoff_grand = squeeze(mean(spk_runoff_cells,2));                     % average across cells
ste_runoff = std(spk_runoff_cells,0,2)/sqrt(size(spk_runoff_cells,2));

% speed
speed_runoff = speed(frms_runoff_mat);
aveSpd_runoff = mean(speed_runoff,2);
steSpd_runoff = std(speed_runoff,0,2)/sqrt(size(speed_runoff,2));

%plot
x = (1:size(frms_runoff_mat,1));
figure; subplot(2,1,1);
shadedErrorBar(x,spk_runoff_grand,ste_runoff);hold on;
ylim([0,3]); hold on; 
vline(31,'r'); ylabel('average firing rate');xlabel('frame');
title([sessions ' aligned to running offset']);

subplot(2,1,2);shadedErrorBar(x,aveSpd_runoff, steSpd_runoff);hold on;
vline(31,'r'); ylabel('speed');xlabel('frame');%title([sessions ' running triggered average speed']);

savefig([image_analysis_dest sessions '_runOffAve_FRwSpeed']);

%% plots to go back and look at the raw data

% plot first derivatives with horizontal std lines --------------------------------------------------
std_deriv = std(deriv);
std2 = 2*std_deriv;
std3 = 3*std_deriv;
std2_5 = 2.5*std_deriv;
listfrm = [1:300; 12051:12350;19001:19300;28501:28800;3301:3600; 6801:7100;  16881:17180;  ...
  21661:21960; 26001:26300 ];
listfrm = listfrm';
[fig_deriv] = GUI_w_hline(deriv, listfrm,std_deriv,std2,std2_5,std3);
savefig([image_analysis_dest sessions '_GUIderiv_stdlines']);

% after function best_spkthreshold_2P, go back to raw fluorescence(TCave) and see if the time points of spikes make sense
% plot red dots over rawTC when deriv > std_bestv for all cells and several frames---------------------------------------------------
listfrm = [1:300; 12051:12350;19001:19300;28501:28800;3301:3600; 6801:7100;  16881:17180;  ...
  21661:21960; 26001:26300 ];
listfrm = listfrm';
[fig] = GUI(TCave, deriv, listfrm,std_best);  
savefig([image_analysis_dest sessions '_GUI_TCave_w_spkEvent_stdbest']); 

% Plot raw fluorescence for the running triggered average windows -----
% this part is not working now, don't know why
[fig_runtrig] = GUI_w_vline(TCave, deriv, frms_runTrig_mat,std_best); %!!!!! you will get an error if you have less than 9 running windows
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

