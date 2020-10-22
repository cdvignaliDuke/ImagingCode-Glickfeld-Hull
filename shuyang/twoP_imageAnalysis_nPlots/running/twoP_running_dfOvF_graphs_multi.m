%% 2 photon: running - dfOvF
% have to plot these after deconvolution. Because deconvolution tells you
% which cells are bad, and you don't look at the fluorescence of those cells (in section one, it's TCave_cl instead of TCave).
%plot average dfOvF for all cells in one trial, aligned with running onset and offset

%% assign document paths and experimental sessions
clear;
sessions = {'190903_img1023'}; 
days = {'1023-190903_1'};
image_analysis_base    = 'Z:\Analysis\2P_MovingDots_Analysis\imaging_analysis\'; %stores the data on crash in the movingDots analysis folder
color_code = {'b','r','k','c'};

%% SECTION I df/f and speed
for i = 1:size(sessions,2)
    image_analysis_dest = [image_analysis_base, sessions{i}, '\'];
    behav_dest = ['Z:\Analysis\2P_MovingDots_Analysis\behavioral_analysis\' days{i} '\'];
   %save([image_analysis_dest sessions{i} '_deconvolution_thresh', num2str(threshold), '_TCave_cl.mat'], 'TCave_cl','threshold','badPCs');

    threshold = -4;
    filename = dir([image_analysis_dest '*' num2str(threshold) '_TCave_cl.mat']);
    TCave_cl = load([image_analysis_dest filename.name]);
    TCave_cl = TCave_cl.TCave_cl;
    behav_output = load([behav_dest days{i} '_behavAnalysis.mat']);
    speed = double(behav_output.speed);
    frm_stay_cell = behav_output.frames_stay_cell;
    frm_stay = cell2mat(frm_stay_cell);
    baseline = mean(TCave_cl(frm_stay,:));         % baseline for each cell
    dfOvF = zeros(size(TCave_cl,1),size(TCave_cl,2));
    for n = 1:size(TCave_cl,2)                     % for each cell
        dfOvF(:,n) = (TCave_cl(:,n)-baseline(n))/baseline(n); %frames*cell
    end
    % plot df/f over speed for the first 5000 frames (choose a period when you have both running and stationary windows)
    avedfOvF = mean(dfOvF,2);
    stedfOvF = std(dfOvF,0,2)/sqrt(size(dfOvF,2));
    frame = (0:(length(speed)-1));
    figure; n1 = 26000; n2 = 29000;
    [hAx,hline1,hline2] = plotyy(frame(n1:n2),speed(n1:n2),frame(n1:n2),avedfOvF(n1:n2)');
    %[hAx,hline1,hline2] = plotyy(frame,speed,frame,(avedfOvF(1:length(speed))'));
    
    %set(hline1,'color', 'b');
    %set(hline2,'color',color_code{n});
    legend('speed','df/f','Location','northeast');
    ylim(hAx(1),[-30 100]);
    ylim(hAx(2),[-0.3,0.6]);
    title(['speed and average df/f ' sessions{i}]);
    savefig([image_analysis_dest sessions{i} '_avedfOvF_wspeed']);
    save([image_analysis_dest sessions{i} '_dfOvF.mat'],'dfOvF','avedfOvF','stedfOvF');
    
end


%% SECTION II run triggered ave and run offset ave --- df/f
for i = 1: size(sessions,2)
    image_analysis_dest = [image_analysis_base, sessions{i}, '\'];
    behav_dest = ['Z:\Analysis\2P_MovingDots_Analysis\behavioral_analysis\' days{i} '\'];
    dfOvF_output = load([image_analysis_dest sessions{i} '_dfOvF.mat']);
    dfOvF = dfOvF_output.dfOvF;
    behav_output = load([behav_dest days{i} '_behavAnalysis.mat']);
    frames_run_cell = behav_output.frames_run_cell;
    speed = double(behav_output.speed);
    
    period = 9; % # of frames right before and right after running -- not useful for this section but the find frames function needs this
    befoRun = 30; %1s before running onset
    befoRunStay = 25; %1s, # of frames that speed has to be 0 before running
    totalT = 60; %2s # of frmaes before running onset + # of frames following running onset/ # of frames before running offset + # of frames after running offset
    aftRunOff = 30; %1s,# of frames that speed = 0 after running
    [~,~,frms_runTrig_mat,frms_runoff_mat,~]= findFrames_runWindows_2P(speed,...
        frames_run_cell,period,befoRunStay,totalT,aftRunOff);
    % aligned with running onset ----------------------------------------------------------------------------
    dfOvF_runtrig_mat = zeros(size(frms_runTrig_mat,1),size(frms_runTrig_mat,2),size(dfOvF,2));
    for w = 1: size(frms_runTrig_mat,2)                                    % for each running trig window
        dfOvF_runtrig_mat(:,w,:) = dfOvF(frms_runTrig_mat(:,w),:);         % df/f for all cells during that window, frame*cells
    end
    ave_dfOvF_runtrig_trials = squeeze(mean(dfOvF_runtrig_mat,2));
    ave_dfOvF_runtrig = squeeze(mean(ave_dfOvF_runtrig_trials,2));
    ste_dfOvF_runTrig = std(ave_dfOvF_runtrig_trials,0,2)/sqrt(size(ave_dfOvF_runtrig_trials,2)); % a bunch of zeros
    % speed
    speed_runtrig = speed(frms_runTrig_mat);
    aveSpd_runtrig = mean(speed_runtrig,2);
    steSpd_runtrig = std(speed_runtrig,0,2)/sqrt(size(speed_runtrig,2));
    %plot
    if size(frms_runTrig_mat,2)<=1 % if there's only one trial that fullfills the criteria, it can't plot the speed because aveSpd_runtrig/off is just a single value. so for now I made it just not do the plot if there's only one trial
        disp('!!There is only one running triggered window in fullfills the criteria, not plotting!!');
        disp(sessions{i});
    else
        x = (1:size(frms_runTrig_mat,1));
        figure; subplot(2,1,1);
        shadedErrorBar(x,ave_dfOvF_runtrig,ste_dfOvF_runTrig);hold on;
        vline(31,'r'); ylabel('df/f'); ylim([-0.1 0.5]);
        title([sessions{i} ' running triggered average']);
        subplot(2,1,2); %plot(speed_runtrig);
        shadedErrorBar(x,aveSpd_runtrig, steSpd_runtrig);hold on;
        vline(31,'r');ylabel('speed');
        savefig([image_analysis_dest sessions{i} '_runTrigAve_dfOvF_wspeed']);
    end
    % save the variables you plot
    
    % aligned with running offset ------------------------------------------------------------------------
    dfOvF_runoff_mat = zeros(size(frms_runoff_mat,1),size(frms_runoff_mat,2),size(dfOvF,2));
    for w = 1: size(frms_runoff_mat,2)                                    % for each running trig window
        dfOvF_runoff_mat(:,w,:) = dfOvF(frms_runoff_mat(:,w),:);         % df/f for all cells during that window, frame*cells
    end
    ave_dfOvF_runoff_trials = squeeze(mean(dfOvF_runoff_mat,2));
    ave_dfOvF_runoff = squeeze(mean(ave_dfOvF_runoff_trials,2));
    ste_dfOvF_runoff = std(ave_dfOvF_runoff_trials,0,2)/sqrt(size(ave_dfOvF_runoff_trials,2)); % a bunch of zeros
    % speed
    speed_runoff = speed(frms_runoff_mat);
    aveSpd_runoff = mean(speed_runoff,2);
    steSpd_runoff = std(speed_runoff,0,2)/sqrt(size(speed_runoff,2));
    %plot
    if size(frms_runoff_mat,2)<=1
        disp('!!There is only one running off window in fullfills the criteria, not plotting!!');
        disp(sessions{i});
    else
        x = (1:size(frms_runoff_mat,1));
        figure; subplot(2,1,1);
        shadedErrorBar(x,ave_dfOvF_runoff,ste_dfOvF_runoff);hold on;
        vline(31,'r'); ylabel('df/f'); xlabel('frame');
        title([sessions{i} ' aligned to running offset -  average']);
        subplot(2,1,2); %plot(speed_runoff);
        shadedErrorBar(x,aveSpd_runoff, steSpd_runoff);hold on;
        vline(31,'r');ylabel('speed'); xlabel('frame');
        savefig([image_analysis_dest sessions{i} '_runoffset_dfOvF_wspeed']);
    end
    
    %pause;
end


%% SECTION II run triggered ave and run offset ave --- df/f --- 
% REGARDLESS OF whether the animal is stationary all the time before/after running
for i = 1:size(sessions,2)
    image_analysis_dest = [image_analysis_base, sessions{i}, '\'];
    behav_dest = ['Z:\Analysis\2P_MovingDots_Analysis\behavioral_analysis\' days{i} '\'];
    dfOvF_output = load([image_analysis_dest sessions{i} '_dfOvF.mat']);
    dfOvF = dfOvF_output.dfOvF;
    behav_output = load([behav_dest days{i} '_behavAnalysis.mat']);
    frames_run_cell = behav_output.frames_run_cell;
    speed = double(behav_output.speed);
    
    period = 9; % # of frames before and after running
    befoRunStay = 30; %1s, # of frames that speed =0 before running
    totalT = 60; %2s # of frmaes before running onset + # of frames following running onset/ # of frames before running offset + # of frames after running offset
    aftRunOff = 30; %1s,# of frames that speed = 0 after running
    [~,~,frms_runTrig_NC_mat,frms_runoff_NC_mat,~]= findFrames_run_no_criteria_2P(speed,...
        frames_run_cell,period,befoRunStay,totalT,aftRunOff);
    % aligned with running onset ----------------------------------------------------------------------------
    dfOvF_runtrig_NC_mat = zeros(size(frms_runTrig_NC_mat,1),size(frms_runTrig_NC_mat,2),size(dfOvF,2));
    for w = 1: size(frms_runTrig_NC_mat,2)                                    % for each running trig window
        dfOvF_runtrig_NC_mat(:,w,:) = dfOvF(frms_runTrig_NC_mat(:,w),:);         % df/f for all cells during that window, frame*cells
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
    title([sessions{i} ' running triggered average INCLUSIVE']);
    vline(31,'r'); ylabel('df/f');
    subplot(2,1,2); %plot(speed_runtrig);
    shadedErrorBar(x,aveSpd_runtrig_NC, steSpd_runtrig_NC);hold on;
    vline(31,'r');ylabel('speed');
    savefig([image_analysis_dest sessions{i} '_runTrigAve_dfOvF_wspeed_INCLUSIVE']);
    
    % aligned with running offset ---------------------------------------------------------------------------------------
    dfOvF_runoff_NC_mat = zeros(size(frms_runoff_NC_mat,1),size(frms_runoff_NC_mat,2),size(dfOvF,2));
    for w = 1: size(frms_runoff_NC_mat,2)                                    % for each running trig window
        dfOvF_runoff_NC_mat(:,w,:) = dfOvF(frms_runoff_NC_mat(:,w),:);         % df/f for all cells during that window, frame*cells
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
    title([sessions{i} ' aligned to running offset -  average INCLUSIVE']);
    subplot(2,1,2); %plot(speed_runoff);
    shadedErrorBar(x,aveSpd_runoff_NC, steSpd_runoff_NC);hold on;
    vline(31,'r');ylabel('speed'); xlabel('frame');
    savefig([image_analysis_dest sessions{i} '_runoffset_dfOvF_wspeed_INCLUSIVE']);
   
    %pause;
end


%% df/f stationary vs. running -------probably NOT USEFUL
% % the df/f during stationary will always be very close to 0 (10 to the -16), and the df/f
% % during running will be a number at a level of 10 to the -2, so it will look wierd anyway when you plot it.
% for i = 1: size(sessions,2)
%     image_analysis_dest = [image_analysis_base, sessions{i}, '\'];
%     behav_dest = ['Z:\Analysis\2P_MovingDots_Analysis\behavioral_analysis\' days{i} '\'];
%     behav_output = load([behav_dest days{i} '_behavAnalysis.mat']);
%     frm_run_cell = behav_output.frames_run_cell;
%     frm_run = cell2mat(frm_run_cell);
%     frm_stay_cell = behav_output.frames_stay_cell;
%     frm_stay = cell2mat(frm_stay_cell);
%     dfOvF_output = load([image_analysis_dest sessions{i} '_dfOvF.mat']);
%     dfOvF = dfOvF_output.dfOvF;
%     avedfOvFrun = mean(dfOvF(frm_run,:));
%     avedfOvFstay = mean(dfOvF(frm_stay,:));
%     grand_dfOvFrun = mean (avedfOvFrun);
%     grand_dfOvFstay = mean (avedfOvFstay);
%     figure;
%     scatter(avedfOvFstay,avedfOvFrun,'filled','r'); hold on;
%     scatter(grand_dfOvFstay, grand_dfOvFrun,'filled','k');hold on;
%     xlabel('stationary'); ylabel('running');
%     xlim([-0.1,0.5]); ylim([min(avedfOvFstay),max(avedfOvFstay)]);
%     refline(1,0);
%     title('mean df/f for each cell');
%     savefig([image_analysis_dest sessions{i} '_meandfOvFstayVSrun']);
% end

