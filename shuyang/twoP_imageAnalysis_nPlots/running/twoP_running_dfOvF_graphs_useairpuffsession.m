%% 2 photon: running - dfOvF
%plot average dfOvF for all cells in one trial, aligned with running onset and offset
%ave_dfOvF_runtrig_cells: frame*cells, averaged across trials
%ave_dfOvF_runtrig:     frames*1, averaged across cells
%ste_dfOvF_runtrig:     frames*1, cell by cell variation
%aveSpd_runtrig:        1*frames, averaged across trials
%steSpd_runtrig:        1*frames if more than 1 trial, trial by trial variation.
%                       a string that says no variation if only 1 trial


%% assign document paths and experimental sessions
clear;
sessions = '191115_img1042'; 
image_analysis_base  = 'Z:\Analysis\Airpuff_analysis\imaging_analysis\'; %stores the data on crash in the movingDots analysis folder
image_analysis_dest = [image_analysis_base, sessions, '\'];

% behavior analysis results
days = '1042-191115_1';
behav_dest = ['Z:\Analysis\Airpuff_analysis\behavioral_analysis\' days '\'];
color_code = {'b','r','k','c'};

%% SECTION I df/f and speed
% load df/f bottom (F using bottom 10 percent of fluorescence values)
dfOvF_output = load([image_analysis_dest sessions '_dfOvF.mat']);
avedfOvF_btm = dfOvF_output.dfOvF_btm_cl;
behav_output = load([behav_dest days '_behavAnalysis.mat']);
speed = behav_output.speed;
% plot df/f over speed for the first 5000 frames (choose a period when you have both running and stationary windows)
frame = (0:(length(speed)-1));
% convert frames to seconds
x = frame/30;

figure; n1 = 1; n2 = 26436;
%convert speed to cm/s: each unit = 2pai*r/128 (there're 128 pulses in total for each circle)
speedcm = speed*2*3.1415926*7.5/128;
[hAx,hline1,hline2] = plotyy(x(n1:n2),speedcm(n1:n2),x(n1:n2),avedfOvF_btm(n1:n2)');
%[hAx,hline1,hline2] = plotyy(frame,speed,frame,(avedfOvF(1:length(speed))'));
%set(hline1,'color', 'b'); 
%set(hline2,'color',color_code{n});
legend('speed','df/f','Location','northwest');
ylim(hAx(1),[-15 50]);
ylim(hAx(2),[-0.5,3]);
xlabel('time (s)');
title(['speed and average df/f' days]);
savefig([image_analysis_dest sessions '_avedfOvF_wspeed']);


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
totalT = 75; %2s # total T = of frames before running onset + # of frames following running onset/ # of frames before running offset + # of frames after running offset
aftRunStay = 42; %1s,# of frames that speed = 0 after running
aftRun = 45;%# of frames to include after running ends
[~,~,frms_runTrig_mat,frms_runoff_mat,~]= findFrames_runWindows_2P(speed,...
     frames_run_cell,period,befoRunStay,totalT,aftRunStay,befoRun,aftRun);
%airpuffon = behav_output.airpuffon1;
% [~,~,frms_runTrig_mat,frms_runoff_mat,~]= findFrames_runWindows_2P_useairpuffsessions(speed,...
%     frames_run_cell,period,befoRunStay,totalT,aftRunStay,befoRun,aftRun,airpuffon);
% tried to sort out the running trials that have airpuff, and then have no running trials left at all. tried this for all 4 airpuff sessions, none left
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
%calculate speed every 100ms, but the vector length needs to be the same as
%df/f
bin = length(aveSpd_runtrig)/3;
mean_runtrig_every100ms = zeros(1,length(aveSpd_runtrig));
ste_runtrig_every100ms = zeros(1,length(aveSpd_runtrig));

for x = 1:bin
    y = 3*x -2;
    mean_runtrig_every100ms(y:y+2) = mean(aveSpd_runtrig(y:y+2));
    ste_runtrig_every100ms(y:y+2) = mean(steSpd_runtrig(y:y+2));
    if abs(mean_runtrig_every100ms(y))<2 %on average moves less than 2 units per 0.1s, consider as stationary
     mean_runtrig_every100ms(y:y+2) = 0;
     ste_runtrig_every100ms(y:y+2) = 0;
    end
end

save([image_analysis_dest sessions '_dfOvF.mat'],'ave_dfOvF_runtrig_cells',...
    'ave_dfOvF_runtrig','ste_dfOvF_runTrig','-append');
save([behav_dest days '_behavAnalysis.mat'],'aveSpd_runtrig','steSpd_runtrig',...
    'mean_runtrig_every100ms','ste_runtrig_every100ms','-append');

%plot
x = (1:size(frms_runTrig_mat,1))/30;
if size(frms_runTrig_mat,2) > 0
    figure; subplot(2,1,1);
    shadedErrorBar(x,ave_dfOvF_runtrig,ste_dfOvF_runTrig,{'color',[0.1922 0.6392 0.3294]});hold on;
    vline(befoRun/30,'k'); ylabel('df/f');
    %ylim([0 0.05]);
    title([sessions ' running triggered average']);
    subplot(2,1,2); %plot(speed_runtrig);
    %convert speed to cm/s: each unit = 2pai*r/128 (there're 128 pulses in total for each circle)
    if size(frms_runTrig_mat,2) > 1
        shadedErrorBar(x,mean_runtrig_every100ms*2*3.1415926*7.5/128, ste_runtrig_every100ms*2*3.1415926*7.5/128,{'color',[0.0196 0.4392 0.6902]});hold on;
    else
        plot(x,aveSpd_runtrig*2*3.1415926*7.5/128,'k'); hold on;
    end
    vline(befoRun/30,'k');ylabel('speed(cm/s)');xlabel('time(s)');
    text(0.1,min(aveSpd_runtrig*2*3.1415926*7.5/128)+2,['n = ' num2str(size(frms_runTrig_mat,2))]);
    %ylim([-5 20]);
    savefig([image_analysis_dest sessions '_runTrigAve_dfOvF_wspeed100ms']);
end

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

bin = length(aveSpd_runoff)/3;
mean_runoff_every100ms = zeros(1,length(aveSpd_runoff));
ste_runoff_every100ms = zeros(1,length(aveSpd_runoff));
a = 1;
for x = 1:bin
    y = 3*x -2;
    mean_runoff_every100ms(y:y+2) = mean(aveSpd_runoff(y:y+2));
    ste_runoff_every100ms(y:y+2) = mean(steSpd_runoff(y:y+2));
    if abs(mean_runoff_every100ms(y))<2 %on average moves less than 2 units per 0.1s, consider as stationary
     mean_runoff_every100ms(y:y+2) = 0;
     ste_runoff_every100ms(y:y+2) = 0;
    end
end

save([image_analysis_dest sessions '_dfOvF.mat'],'ave_dfOvF_runoff_cells',...
    'ave_dfOvF_runoff','ste_dfOvF_runoff','-append');
save([behav_dest days '_behavAnalysis.mat'],'aveSpd_runoff','steSpd_runoff',...
    'frms_runTrig_mat','frms_runoff_mat','mean_runoff_every100ms','ste_runoff_every100ms','-append');
%plot
x = (1:size(frms_runoff_mat,1))/30;
if size(frms_runoff_mat,2)>0
    figure; subplot(2,1,1);
    shadedErrorBar(x,ave_dfOvF_runoff,ste_dfOvF_runoff,{'color',[0.1922 0.6392 0.3294]});hold on;
    vline((totalT-aftRun)/30,'r'); ylabel('df/f'); %xlabel('frame');
    %ylim([0 0.05]);
    title([sessions ' aligned to running offset -  average']);
    subplot(2,1,2); %plot(speed_runoff);
    if size(frms_runoff_mat,2)>1
        shadedErrorBar(x,mean_runoff_every100ms*2*3.1415926*7.5/128, ste_runoff_every100ms*2*3.1415926*7.5/128,{'color',[0.0196 0.4392 0.6902]});hold on;
    else
        plot(x,aveSpd_runoff*2*3.1415926*7.5/128,'k');
    end
    vline((totalT-aftRun)/30,'r');ylabel('speed(cm/s)'); xlabel('time(s)');
    %ylim([-5 20]);
    text(1.8,min(aveSpd_runtrig*2*3.1415926*7.5/128)+2,['n = ' num2str(size(frms_runoff_mat,2))]);
    savefig([image_analysis_dest sessions '_runoffset_dfOvF_wspeed100ms']);
end

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






