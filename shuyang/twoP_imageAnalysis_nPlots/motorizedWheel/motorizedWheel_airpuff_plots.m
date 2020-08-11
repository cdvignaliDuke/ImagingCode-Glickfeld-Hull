% average df/f during slow speed vs. high speed vs. stationary
% same thins for spikes
%dfOvF_fast: trial*frame*cell (dfOvF of all cells and frames during all fast running trials, first second of each trial is thrown away)
%dfOvF_fast_cells: dfOvF average across trials and frames, each cell has one value
%FR_fast_cells: firing rate of each cell average across trials
%runtrig ave and runoff ave for df/f and spikes


%% assign document paths and experimental sessions
clear;
sessions = '200319_img1064_airpuff_2'; 
image_analysis_base = 'Z:\Analysis\motorizedWheel_Analysis\airpuff\imaging_analysis\'; 
image_analysis_dest = [image_analysis_base, sessions, '\'];

% behavior analysis results 
days = '1064-200319_2';
behav_dest = ['Z:\Analysis\motorizedWheel_Analysis\airpuff\behavioral_analysis\' days '\'];
color_code = {'b','r','k','c'};

%% load data
behav_output = load([behav_dest days '_behavAnalysis.mat']);
airpuff_fast_mat = behav_output.airpuff_fast_mat;
airpuff_stay_mat = behav_output.airpuff_stay_mat;
airpuff_slow_mat = behav_output.airpuff_slow_mat;
runtrig_slow_mat = behav_output.runtrig_slow_mat;
runtrig_fast_mat = behav_output.runtrig_fast_mat;
runoff_slow_mat = behav_output.runoff_slow_mat;
runoff_fast_mat = behav_output.runoff_fast_mat;
spd_decrease_mat = behav_output.spd_decrease_mat;
spd_increase_mat = behav_output.spd_increase_mat;

threshold = -4;% the threshold you used in deconvolveCa function for deconvolution
spk_deconv_output = load([image_analysis_dest sessions,'_spk_deconvolve_threshold' num2str(threshold) '.mat']);
dfOvF = load([image_analysis_dest sessions '_dfOvF.mat']);
dfOvF = dfOvF.dfOvF_btm_cl;

%% df/f runtrig ave and runoffset ave and speed transit
dfOvF_runoff_fast = zeros(size(runoff_fast_mat,1),size(runoff_fast_mat,2),size(dfOvF,2)); %trial*frame*cell
dfOvF_runoff_slow = zeros(size(runoff_slow_mat,1),size(runoff_slow_mat,2),size(dfOvF,2));
dfOvF_runtrig_fast = zeros(size(runtrig_fast_mat,1),size(runtrig_fast_mat,2),size(dfOvF,2));
dfOvF_runtrig_slow = zeros(size(runtrig_slow_mat,1),size(runtrig_slow_mat,2),size(dfOvF,2));
dfOvF_spd_decrease = zeros(size(spd_decrease_mat,1),size(spd_decrease_mat,2),size(dfOvF,2));
dfOvF_spd_increase = zeros(size(spd_increase_mat,1),size(spd_increase_mat,2),size(dfOvF,2));

for c = 1:size(dfOvF,2) %for each cell
    for f = 1:size(runoff_fast_mat,1)
        dfOvF_runoff_fast(f,:,c) = dfOvF(runoff_fast_mat(f,:),c);
    end
    for g = 1:size(runoff_slow_mat,1)
        dfOvF_runoff_slow(g,:,c) = dfOvF(runoff_slow_mat(g,:),c);
    end
    for h = 1:size(runtrig_fast_mat,1)
        dfOvF_runtrig_fast(h,:,c) = dfOvF(runtrig_fast_mat(h,:),c);
    end
    for i = 1:size(runtrig_slow_mat,1)
        dfOvF_runtrig_slow(i,:,c) = dfOvF(runtrig_slow_mat(i,:),c);
    end
    for j = 1:size(spd_decrease_mat,1)
        dfOvF_spd_decrease(j,:,c) = dfOvF(spd_decrease_mat(j,:),c);
    end
    for k = 1:size(spd_increase_mat,1)
        dfOvF_spd_increase(k,:,c) = dfOvF(spd_increase_mat(k,:),c);
    end
end

%average across trials
dfOvF_runoff_fast_cells = squeeze(mean(dfOvF_runoff_fast,1)); %frame*cell
dfOvF_runoff_slow_cells = squeeze(mean(dfOvF_runoff_slow,1));
dfOvF_runtrig_fast_cells = squeeze(mean(dfOvF_runtrig_fast,1));
dfOvF_runtrig_slow_cells = squeeze(mean(dfOvF_runtrig_slow,1));
dfOvF_spd_decrease_cells = squeeze(mean(dfOvF_spd_decrease,1));
dfOvF_spd_increase_cells = squeeze(mean(dfOvF_spd_increase,1));

save([image_analysis_dest sessions '_dfOvF.mat'],'dfOvF_runoff_fast_cells',...
    'dfOvF_runoff_slow_cells','dfOvF_runtrig_fast_cells','dfOvF_runtrig_slow_cells',...
    'dfOvF_spd_decrease_cells','dfOvF_spd_increase_cells','-append');

%average across cells
dfOvF_runoff_fast_session = mean(dfOvF_runoff_fast_cells,2);
dfOvF_runoff_slow_session = mean(dfOvF_runoff_slow_cells,2);
dfOvF_runtrig_fast_session = mean(dfOvF_runtrig_fast_cells,2);
dfOvF_runtrig_slow_session = mean(dfOvF_runtrig_slow_cells,2);
dfOvF_spd_decrease_session = mean(dfOvF_spd_decrease_cells,2);
dfOvF_spd_increase_session = mean(dfOvF_spd_increase_cells,2);

ncells = size(dfOvF,2);
ste_runoff_fast_cells = std(dfOvF_runoff_fast_cells,0,2)/sqrt(ncells);
ste_runoff_slow_cells = std(dfOvF_runoff_slow_cells,0,2)/sqrt(ncells);
ste_runtrig_fast_cells = std(dfOvF_runtrig_fast_cells,0,2)/sqrt(ncells);
ste_runtrig_slow_cells = std(dfOvF_runtrig_slow_cells,0,2)/sqrt(ncells);
ste_spd_decrease_cells = std(dfOvF_spd_decrease_cells,0,2)/sqrt(ncells);
ste_spd_increase_cells = std(dfOvF_spd_increase_cells,0,2)/sqrt(ncells);

dfOvF_trig = figure;
x = (1:length(dfOvF_runoff_fast_session))/30;
subplot(2,3,1);
shadedErrorBar(x,dfOvF_runtrig_slow_session,ste_runtrig_slow_cells,{'color',[0.1373 0.5451 0.2706]});
ylabel('df/f'); ylim([0 1]); xlim([0 2.5]);
vline(1,'k');
title('stay-slow speed');
%xlim([0,1]);
subplot(2,3,4);
shadedErrorBar(x,dfOvF_runoff_slow_session,ste_runoff_slow_cells,{'color',[0.1373 0.5451 0.2706]});
ylabel('df/f'); ylim([0 1]);xlim([0 2.5]);
vline(1,'k');
title('slow speed-stay');
subplot(2,3,2);
shadedErrorBar(x,dfOvF_runtrig_fast_session,ste_runtrig_fast_cells,{'color',[0.1373 0.5451 0.2706]});
ylim([0 1]); vline(1,'k'); xlabel('time(s)');xlim([0 2.5]);
title('stay-fast speed');
subplot(2,3,5);
shadedErrorBar(x,dfOvF_runoff_fast_session,ste_runoff_fast_cells,{'color',[0.1373 0.5451 0.2706]});
ylim([0 1]); vline(1,'k'); xlabel('time(s)');xlim([0 2.5]);
title('fast speed-stay');
subplot(2,3,3);
shadedErrorBar(x,dfOvF_spd_decrease_session,ste_spd_decrease_cells,{'color',[0.1373 0.5451 0.2706]});
ylim([0 1]); xlim([0 2.5]);vline(1,'k');
title('speed decrease');
subplot(2,3,6);
shadedErrorBar(x,dfOvF_spd_increase_session,ste_spd_increase_cells,{'color',[0.1373 0.5451 0.2706]});
xlim([0 2.5]); ylim([0 1]); vline(1,'k');
title('speed increase');
supertitle(['dfOvF' sessions]);

savefig([image_analysis_dest sessions '_dfOvF_state_transit.fig']);


%% FR run trig ave and runoff ave
spk_logic = spk_deconv_output.spk_logic_cl;
%spk_logic = spk_deconv_output.spk_logic;
spk_runoff_fast = zeros(size(runoff_fast_mat,1),size(runoff_fast_mat,2),size(spk_logic,2)); %trial*frame*cell
spk_runoff_slow = zeros(size(runoff_slow_mat,1),size(runoff_slow_mat,2),size(spk_logic,2));
spk_runtrig_fast = zeros(size(runtrig_fast_mat,1),size(runtrig_fast_mat,2),size(spk_logic,2));
spk_runtrig_slow = zeros(size(runtrig_slow_mat,1),size(runtrig_slow_mat,2),size(spk_logic,2));
spk_spd_decrease = zeros(size(spd_decrease_mat,1),size(spd_decrease_mat,2),size(spk_logic,2));
spk_spd_increase = zeros(size(spd_increase_mat,1),size(spd_increase_mat,2),size(spk_logic,2));

for c = 1:size(spk_logic,2) %for each cell
    for f = 1:size(runoff_fast_mat,1)
        spk_runoff_fast(f,:,c) = spk_logic(runoff_fast_mat(f,:),c);
    end
    for g = 1:size(runoff_slow_mat,1)
        spk_runoff_slow(g,:,c) = spk_logic(runoff_slow_mat(g,:),c);
    end
    for h = 1:size(runtrig_fast_mat,1)
        spk_runtrig_fast(h,:,c) = spk_logic(runtrig_fast_mat(h,:),c);
    end
    for i = 1:size(runtrig_slow_mat,1)
        spk_runtrig_slow(i,:,c) = spk_logic(runtrig_slow_mat(i,:),c);
    end
    for j = 1:size(spd_decrease_mat,1)
        spk_spd_decrease(j,:,c) = spk_logic(spd_decrease_mat(j,:),c);
    end
    for k = 1:size(spd_increase_mat,1)
        spk_spd_increase(k,:,c) = spk_logic(spd_increase_mat(k,:),c);
    end
end

%average across trials and get average firing rate for each cell
FR_runoff_fast_cells = squeeze(mean(spk_runoff_fast,1))*30;
FR_runoff_slow_cells = squeeze(mean(spk_runoff_slow,1))*30;
FR_runtrig_fast_cells = squeeze(mean(spk_runtrig_fast,1))*30;
FR_runtrig_slow_cells = squeeze(mean(spk_runtrig_slow,1))*30;
FR_spd_decrease_cells = squeeze(mean(spk_spd_decrease,1))*30;
FR_spd_increase_cells = squeeze(mean(spk_spd_increase,1))*30;

save([image_analysis_dest sessions,'_spk_deconvolve_threshold' num2str(threshold) '.mat'],...
    'FR_runoff_fast_cells','FR_runoff_slow_cells','FR_runtrig_fast_cells',...
    'FR_runtrig_slow_cells','FR_spd_decrease_cells','FR_spd_increase_cells','-append');

%average across session
FR_runoff_fast_session = mean(FR_runoff_fast_cells,2);
FR_runoff_slow_session = mean(FR_runoff_slow_cells,2);
FR_runtrig_fast_session = mean(FR_runtrig_fast_cells,2);
FR_runtrig_slow_session = mean(FR_runtrig_slow_cells,2);
FR_spd_decrease_session = mean(FR_spd_decrease_cells,2);
FR_spd_increase_session = mean(FR_spd_increase_cells,2);

ncells = size(spk_logic,2);
steFR_runoff_fast_cells = std(FR_runoff_fast_cells,0,2)/sqrt(ncells);
steFR_runoff_slow_cells = std(FR_runoff_slow_cells,0,2)/sqrt(ncells);
steFR_runtrig_fast_cells = std(FR_runtrig_fast_cells,0,2)/sqrt(ncells);
steFR_runtrig_slow_cells = std(FR_runtrig_slow_cells,0,2)/sqrt(ncells);
steFR_spd_decrease_cells = std(FR_spd_decrease_cells,0,2)/sqrt(ncells);
steFR_spd_increase_cells = std(FR_spd_increase_cells,0,2)/sqrt(ncells);

FR_trig = figure;
x = (1:length(FR_runoff_fast_session))/30;
subplot(2,3,1);
shadedErrorBar(x,FR_runtrig_slow_session,steFR_runtrig_slow_cells,{'color',[0.1373 0.5451 0.2706]});
ylabel('firing rate'); 
xlim([0 2.5]); ylim([0-0.05 2.5]);vline(1,'k');
title('stay-slow speed');
%xlim([0,1]);
subplot(2,3,4);
shadedErrorBar(x,FR_runoff_slow_session,steFR_runoff_slow_cells,{'color',[0.1373 0.5451 0.2706]});
xlim([0 2.5]); ylim([0-0.05 2.5]);ylabel('firing rate'); 
vline(1,'k');
title('slow speed-stay');
subplot(2,3,2);
shadedErrorBar(x,FR_runtrig_fast_session,steFR_runtrig_fast_cells,{'color',[0.1373 0.5451 0.2706]});
xlabel('time(s)');xlim([0 2.5]); ylim([0-0.05 2.5]); vline(1,'k');
title('stay-fast speed');
subplot(2,3,5);
shadedErrorBar(x,FR_runoff_fast_session,steFR_runoff_fast_cells,{'color',[0.1373 0.5451 0.2706]});
xlabel('time(s)');xlim([0 2.5]); ylim([0-0.05 2.5]); vline(1,'k');
title('fast speed-stay');
subplot(2,3,3);
shadedErrorBar(x,FR_spd_decrease_session,steFR_spd_decrease_cells,{'color',[0.1373 0.5451 0.2706]});
xlim([0 2.5]); ylim([0-0.05 2.5]); vline(1,'k');
title('speed decrease');
subplot(2,3,6);
shadedErrorBar(x,FR_spd_increase_session,steFR_spd_increase_cells,{'color',[0.1373 0.5451 0.2706]});
xlim([0 2.5]); ylim([0-0.05 2.5]); vline(1,'k');
title('speed increase');
supertitle(['Firing rate' sessions]);

savefig([image_analysis_dest sessions '_FR_state_transit.fig']);

%% df/f airpuff triggered average in different behavioral states
dfOvF_airpuff_fast = zeros(size(airpuff_fast_mat,1),size(airpuff_fast_mat,2),size(dfOvF,2)); %trial*frame*cell
dfOvF_airpuff_slow = zeros(size(airpuff_slow_mat,1),size(airpuff_slow_mat,2),size(dfOvF,2));
dfOvF_airpuff_stay = zeros(size(airpuff_stay_mat,1),size(airpuff_stay_mat,2),size(dfOvF,2));

for c = 1:size(dfOvF,2) %for each cell
    for f = 1:size(airpuff_fast_mat,1)
        dfOvF_airpuff_fast(f,:,c) = dfOvF(airpuff_fast_mat(f,:),c);
    end
    for g = 1:size(airpuff_slow_mat,1)
        dfOvF_airpuff_slow(g,:,c) = dfOvF(airpuff_slow_mat(g,:),c);
    end
    for h = 1:size(airpuff_stay_mat,1)
        dfOvF_airpuff_stay(h,:,c) = dfOvF(airpuff_stay_mat(h,:),c);
    end
end

%average across trials
dfOvF_airpuff_fast_cells = squeeze(mean(dfOvF_airpuff_fast,1)); %frame*cell
dfOvF_airpuff_slow_cells = squeeze(mean(dfOvF_airpuff_slow,1));
dfOvF_airpuff_stay_cells = squeeze(mean(dfOvF_airpuff_stay,1));

%combine matrix of 2 running speeds
dfOvF_airpuff_run = cat(1,dfOvF_airpuff_fast,dfOvF_airpuff_slow);

save([image_analysis_dest sessions '_dfOvF.mat'],'dfOvF_airpuff_fast_cells',...
    'dfOvF_airpuff_slow_cells','dfOvF_airpuff_stay_cells',...
    'dfOvF_airpuff_run','dfOvF_airpuff_stay','-append');

%average across cells
dfOvF_airpuff_fast_session = mean(dfOvF_airpuff_fast_cells,2);
dfOvF_airpuff_slow_session = mean(dfOvF_airpuff_slow_cells,2);
dfOvF_airpuff_stay_session = mean(dfOvF_airpuff_stay_cells,2);

ncells = size(dfOvF,2);
ste_airpuff_fast_cells = std(dfOvF_airpuff_fast_cells,0,2)/sqrt(ncells);
ste_airpuff_slow_cells = std(dfOvF_airpuff_slow_cells,0,2)/sqrt(ncells);
ste_airpuff_stay_cells = std(dfOvF_airpuff_stay_cells,0,2)/sqrt(ncells);

dfOvFairpuff_trig = figure;
x = (1:length(dfOvF_airpuff_fast_session))/30;
subplot(1,3,1);
shadedErrorBar(x,dfOvF_airpuff_fast_session,ste_airpuff_fast_cells,{'color',[0.1373 0.5451 0.2706]});
ylabel('df/f'); ylim([0 0.7]); xlim([0 2]);
vline(1,'k');
title('fast speed');
%xlim([0,1]);
subplot(1,3,2);
shadedErrorBar(x,dfOvF_airpuff_slow_session,ste_airpuff_slow_cells,{'color',[0.1373 0.5451 0.2706]});
ylim([0 0.7]);xlim([0 2]);
xlabel('time(s)');
vline(1,'k');
title('slow speed');
subplot(1,3,3);
shadedErrorBar(x,dfOvF_airpuff_stay_session,ste_airpuff_stay_cells,{'color',[0.1373 0.5451 0.2706]});
ylim([0 0.7]); vline(1,'k'); xlim([0 2]);
title('stationary');
supertitle(['tactile stimulus' sessions]);

savefig([image_analysis_dest sessions '_dfOvF_airpuff.fig']);

%% FR airpuff triggered average in different behavioral states
spk_logic = spk_deconv_output.spk_logic_cl;
FR_airpuff_fast = zeros(size(airpuff_fast_mat,1),size(airpuff_fast_mat,2),size(spk_logic,2)); %trial*frame*cell
FR_airpuff_slow = zeros(size(airpuff_slow_mat,1),size(airpuff_slow_mat,2),size(spk_logic,2));
FR_airpuff_stay = zeros(size(airpuff_stay_mat,1),size(airpuff_stay_mat,2),size(spk_logic,2));

for c = 1:size(spk_logic,2) %for each cell
    for f = 1:size(airpuff_fast_mat,1)
        FR_airpuff_fast(f,:,c) = spk_logic(airpuff_fast_mat(f,:),c);
    end
    for g = 1:size(airpuff_slow_mat,1)
        FR_airpuff_slow(g,:,c) = spk_logic(airpuff_slow_mat(g,:),c);
    end
    for h = 1:size(airpuff_stay_mat,1)
        FR_airpuff_stay(h,:,c) = spk_logic(airpuff_stay_mat(h,:),c);
    end
end

%average across trials
FR_airpuff_fast_cells = squeeze(mean(FR_airpuff_fast,1))*30; %frame*cell
FR_airpuff_slow_cells = squeeze(mean(FR_airpuff_slow,1))*30;
FR_airpuff_stay_cells = squeeze(mean(FR_airpuff_stay,1))*30;

save([image_analysis_dest sessions,'_spk_deconvolve_threshold' num2str(threshold) '.mat'],...
    'FR_airpuff_fast_cells','FR_airpuff_slow_cells','FR_airpuff_stay_cells','-append');

%average across cells
FR_airpuff_fast_session = mean(FR_airpuff_fast_cells,2);
FR_airpuff_slow_session = mean(FR_airpuff_slow_cells,2);
FR_airpuff_stay_session = mean(FR_airpuff_stay_cells,2);

ncells = size(dfOvF,2);
steFR_airpuff_fast_cells = std(FR_airpuff_fast_cells,0,2)/sqrt(ncells);
steFR_airpuff_slow_cells = std(FR_airpuff_slow_cells,0,2)/sqrt(ncells);
steFR_airpuff_stay_cells = std(FR_airpuff_stay_cells,0,2)/sqrt(ncells);

FRairpuff_trig = figure;
x = (1:length(FR_airpuff_fast_session))/30;
subplot(1,3,1);
shadedErrorBar(x,FR_airpuff_fast_session,steFR_airpuff_fast_cells,{'color',[0.1373 0.5451 0.2706]});
ylabel('Firing Rate'); ylim([0 3]); xlim([0 2]);
vline(1,'k');
title('fast speed');
%xlim([0,1]);
subplot(1,3,2);
shadedErrorBar(x,FR_airpuff_slow_session,steFR_airpuff_slow_cells,{'color',[0.1373 0.5451 0.2706]});
ylim([0 3]);xlim([0 2]);
xlabel('time(s)');
vline(1,'k');
title('slow speed');
subplot(1,3,3);
shadedErrorBar(x,FR_airpuff_stay_session,steFR_airpuff_stay_cells,{'color',[0.1373 0.5451 0.2706]});
ylim([0 3]); vline(1,'k'); xlim([0 2]);
title('stationary');
supertitle(['tactile stimulus' sessions]);

savefig([image_analysis_dest sessions '_FR_airpuff.fig']);
