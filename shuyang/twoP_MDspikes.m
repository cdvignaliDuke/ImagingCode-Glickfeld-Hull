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
frm_run = behav_output.frames_run_cell;
frm_stay = behav_output.frames_stay_cell;

% get the best threshold std for deciding spikes, ave FRs during stationary for each cell 
[avesWvarystd,best_thres,aveFR_allcells_stay,std_best] = best_spkthreshold_2P (deriv, frm_stay,std2,std2_5,std3,std_deriv); 
aveFR_stay = mean(aveFR_allcells_stay);

%% SECTION II 
% calculate ave FRs during running for each cell and make scatter plot for
% spikes rates during running vs. still, each dot is a cell
nspk = zeros(1,size(deriv,2));                                              % num of cells
frate = zeros(size(frm_run,2),size(deriv,2));                               % num of stationary windows*num of cells
for i = 1: size(frm_run,2)                                                  % for each running window
    t = length(frm_run{i})/30;                                              %t = nframes*1/30s, thus the unit of t is second
    der_run_temp = deriv(frm_run{i},:); 
    for j = 1: size(der_run_temp,2)                                         % for each cell
        nspk(j) = sum(der_run_temp(:,j) >= std_best(j));                    %number of spikes during this window for each cell
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


%% SECTION III : run trigger average and run offset average
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
ylim([0,2]); hold on; vline(16,'r');
title([sessions ' running triggered average']);
savefig([image_analysis_dest sessions '_runTrigAve']);
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
ylim([0,2]); hold on; vline(31,'r');
title([sessions ' aligned to running offset']);
savefig([image_analysis_dest sessions '_runOff']);
% PLOT speed
speed_runoff = speed(frms_runoff_mat);
aveSpd_runoff = mean(speed_runoff,2);
steSpd_runoff = std(speed_runoff,0,2)/sqrt(size(speed_runoff,2));
figure; shadedErrorBar(x,aveSpd_runoff, steSpd_runoff);hold on;
vline(31,'r'); title([sessions ' aligned to offset']);
savefig([image_analysis_dest sessions '_runoff_speed']);


%% plots when deciding the best threshold

% plot first derivatives with horizontal std lines --------------------------------------------------
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

% after function best_spkthreshold_2P, go back to raw fluorescence and see if the time points of spikes make sense
% plot red dots over rawTC when deriv > std_bestv for all cells and several frames---------------------------------------------------
listfrm = [1:1000; 3001:4000; 6001:7000; 10001:11000; 15001:16000; 17001:18000; ...
  21001:22000; 24001:25000; 28001:29000];
listfrm = listfrm';
[fig] = GUI(TCave, deriv, listfrm,std_best);
savefig([image_analysis_dest sessions '_GUIraw_fluo_w_spkEvent']);

% Plot raw fluorescence for the running triggered average windows
[fig_runtrig] = GUI_w_vline(TCave, deriv, frames_runTrigger_mat,std_best);
savefig([image_analysis_dest sessions '_GUIraw_fluo_RUNtrig']);

% plot all of the first derivatives: not very useful ------------------------------------------------------------
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
