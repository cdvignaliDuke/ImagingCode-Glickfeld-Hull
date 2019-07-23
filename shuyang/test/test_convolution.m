clear;
sessions = '190117_img1016'; 
image_analysis_base    = 'Z:\Analysis\2P_MovingDots_Analysis\imaging_analysis\'; %stores the data on crash in the movingDots analysis folder
image_analysis_dest = [image_analysis_base, sessions, '\'];

% behavior analysis results 
days = '1016-190117_1';
behav_dest = ['Z:\Analysis\2P_MovingDots_Analysis\behavioral_analysis\' days '\'];
color_code = {'b','r','k','c'};

filename = dir([image_analysis_dest 'getTC\' '*' '_TCave.mat']);
TCave = load([image_analysis_dest 'getTC\' filename.name]);
TCave = TCave.tc_avg;
deriv = diff(TCave,1,1); % derivatives shoud have positive and negative values
behav_output = load([behav_dest days '_behavAnalysis.mat']);
frm_stay_cell = behav_output.frames_stay_cell;
frm_stay = cell2mat(frm_stay_cell);
%%

% smooth the raw fluorescence: conv(vector,ones(1,n)/n, 'same')
% convolution: taking a rectangle and move it over the original signal and multiply it. 
% get rid of the spikes events within a certain period of time, there
% shouldn't be two spike events too close to each other (within however many miliseconds)

TCave2 = TCave(:,2); figure; plot(TCave2(1:200));
ConvTC2 = conv(TCave2,ones(1,5)/100, 'same'); %ones(1,a)/b, if don't want to lose too much info, a shoud be smaller than b, the smaller a is, the less info you lose
% tried a ton of numbers and doesn't seem change much, as long as all info
% looks still there.
figure; plot(ConvTC2(1:200))

%%
convTC = zeros(size(TCave,1),size(TCave,2));
for i = 1:size(TCave,2)
    convTC(:,i) = conv(TCave(:,i),ones(1,5)/100,'same');
end
derivConv = diff(convTC,1,1);
[aveFRsWstds, best_thres,bestFR,std_best,spk_inx,FRstay_cells,...
    spk_bi_cellmat] = twoP_best_spkthreshold (derivConv, frm_stay,convTC);

listfrm = [1:300; 12051:12350;19001:19300;28501:28800;3301:3600; 6801:7100;  16881:17180;  ...
  21661:21960; 26001:26300 ];
listfrm = listfrm';
[fig] = GUI(TCave,listfrm,spk_inx); 

%%
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

%%

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






