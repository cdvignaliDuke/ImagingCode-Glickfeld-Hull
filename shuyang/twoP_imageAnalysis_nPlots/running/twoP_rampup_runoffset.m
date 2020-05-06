%% SECTION - assign pathnames and datasets to be analyzed/written. 
clear;
%sessions = {'190429_img1021','190430_img1023','190507_img1024','190603_img1025'};
%days = {'1021-190429_1','1023-190430_1','1024-190507_1','1025-190603_1'};

sessions = {'190429_img1021','190430_img1023','190507_img1024','190603_img1025'};
days = {'1021-190429_1','1023-190430_1','1024-190507_1','1025-190603_1'};
%there might be more than 1 sessions on a single subject on the same day
image_analysis_base    = 'Z:\Analysis\2P_MovingDots_Analysis\imaging_analysis\'; %stores the data on crash in the movingDots analysis folder

% behavior analysis results 
color_code = {'r','g','m','y','b'};

%% 1. when does the first peak happen?
for i = 1:length(sessions)
    image_analysis_dest = [image_analysis_base, sessions{i}, '\'];
    dfOvF_output = load([image_analysis_dest sessions{i} '_dfOvF.mat']);
    avedfOvF = dfOvF_output.avedfOvF_btm_cl;
    behav_dest = ['Z:\Analysis\2P_MovingDots_Analysis\behavioral_analysis\' days{i} '\'];
    behav_output = load([behav_dest days{i} '_behavAnalysis.mat']);
    speed = behav_output.speed;
    frames_runoff_mat = behav_output.frms_runoff_mat;%frames*trials
    % generate a matrix for loooking at peak times after running offset
    offtime = 75; %2.5 seconds after running offset
    offset_mat = zeros(offtime,size(frames_runoff_mat,2)); %frames*trials
    for t = 1: size(frames_runoff_mat,2)
        offset_mat(:,t) = frames_runoff_mat(31,t):frames_runoff_mat(31,t)+offtime-1; % in frames_runoff_mat, the 31th frame in each trial is the first frame after running offset
    end
    dfOvF_offsetmat = avedfOvF(offset_mat);%frames*trials
    
    peak_loc = zeros(1,size(dfOvF_offsetmat,2));
    for t = 1:size(dfOvF_offsetmat,2) %for each trial
        [pks,locs] = findpeaks(dfOvF_offsetmat(:,t));
        %findpeaks(dfOvF_offsetmat(:,t)); %automatically plot the peaks
        %find all peaks --> if the first peak happens to be the biggest
        %one, then the first peak is the real peak.
        % if not, compare the peaks before the largest peak with the
        % largest one, if those peaks is not much smaller than the largest
        % one, then the previous peak is the first real peak.
        biggest = max(pks);
        inx = 1:length(locs);
        biggest_inx = inx(pks==biggest);
        if biggest_inx==1
            peak_loc(t) = locs(1);
        else
            for p = 1:biggest_inx
                pk_diff = biggest - pks(p);
                if pk_diff<=0.15
                    peak_loc(t) = locs(p);
                    break; %once fullfills the requirement, jump out of the for loop
                end
            end
        end
        
    end
    
    %make a binary mask, only peak frame = 1, others = 0
    
    fig_dest = [image_analysis_dest '\runoffPeak\'];
    if ~exist(fig_dest)
        mkdir(fig_dest);
    end
    
    runoff_bi = zeros(size(dfOvF_offsetmat));%frames*trials
    for t = 1:size(dfOvF_offsetmat,2)
        runoff_bi(peak_loc(t),t)=1;
    end
    x = (0.1:0.1:2.5);
    y = [0 size(dfOvF_offsetmat,2)];
    figure; imagesc(x,y,runoff_bi');
    ylabel('trial #');
    xlabel('time(s)');
    title([sessions{i}]);
    savefig([fig_dest sessions{i} '1peakTime']);
    
    save([image_analysis_dest sessions{i} '_dfOvF.mat'] ,'dfOvF_offsetmat','peak_loc','runoff_bi','-append');
end


%% 2. how many peaks does each trial have in total 
for i = 1:length(sessions)
    image_analysis_dest = [image_analysis_base, sessions{i}, '\'];
    dfOvF_output = load([image_analysis_dest sessions{i} '_dfOvF.mat']);
    dfOvF_offsetmat = dfOvF_output.dfOvF_offsetmat;
    big_pks = zeros(size(dfOvF_offsetmat)); %frames*trials
    sumpks = zeros(1,size(dfOvF_offsetmat,2));
    for t = 1:size(dfOvF_offsetmat,2) %for each trial
        [pks,locs] = findpeaks(dfOvF_offsetmat(:,t));
        ave_pks = mean(pks);
        for p = 1:length(pks)
            if pks(p)> ave_pks
                big_pks(locs(p),t) = 1;
            end
        end
        sumpks(t) = sum(big_pks(:,t)==1);
    end
    save([image_analysis_dest sessions{i} '_dfOvF.mat'] ,'big_pks','sumpks','-append');
end

%% is trial length and average speed correlated with when the first peak happens/how many peaks are there?
for i = 1:length(sessions)
    image_analysis_dest = [image_analysis_base, sessions{i}, '\'];
    dfOvF_output = load([image_analysis_dest sessions{i} '_dfOvF.mat']);
    peak_loc = dfOvF_output.peak_loc;
    sumpks = dfOvF_output.sumpks;
    behav_dest = ['Z:\Analysis\2P_MovingDots_Analysis\behavioral_analysis\' days{i} '\'];
    behav_struct = load([behav_dest '\' days{i} '_behavAnalysis.mat']);
    speed = behav_struct.speed;
    frames_runoff_include = behav_struct.frames_runoff_include;
    trial_lens = cellfun(@length,frames_runoff_include);
    speed_runoff_include = cell(size(frames_runoff_include));
    for t = 1:size(frames_runoff_include,2)
        speed_runoff_include{t} = speed(frames_runoff_include{t});
    end
    avespd_runofftrials = cellfun(@mean,speed_runoff_include);
    
    fig_dest = [image_analysis_dest '\runoffPeak\'];
    %sumpks and ave speed/trial lengths
    
    figure; scatter(avespd_runofftrials*2*3.1415926*7.5/128,sumpks,'filled','MarkerFacecolor',[0.9373 0.396 0.2824]);
    ylim([0 max(sumpks)]); xlim([0 max(avespd_runofftrials*2*3.1415926*7.5/128)]);
    ylabel('total number of df/f peaks');
    xlabel('average speed (cm/s)');
    title([sessions{i}]);
    savefig([fig_dest sessions{i} 'sum_pks_Vs_speed']);
    
    figure; scatter(trial_lens./30,sumpks,'filled','MarkerFacecolor',[0.9373 0.396 0.2824]);
    ylim([0 max(sumpks)]);xlim([0 max(trial_lens./30)]);
    ylabel('total number of df/f peaks');
    xlabel('running trial duration (s)');
    title([sessions{i}]);
    savefig([fig_dest sessions{i} 'sum_pks_Vs_duration']);
    
    figure; scatter(avespd_runofftrials*2*3.1415926*7.5/128,peak_loc,'filled','MarkerFacecolor',[0.9373 0.396 0.2824]);
    ylim([0 max(peak_loc)]); xlim([0 max(avespd_runofftrials*2*3.1415926*7.5/128)]);
    ylabel('time of first peak after running offset');
    xlabel('average speed (cm/s)');
    title([sessions{i}]);
    savefig([fig_dest sessions{i} '_1pk_time_Vs_speed']);
    
    figure; scatter(trial_lens./30,peak_loc,'filled','MarkerFacecolor',[0.9373 0.396 0.2824]);
    ylim([0 max(peak_loc)]);xlim([0 max(trial_lens./30)]);
    ylabel('time of first peak after running offset');
    xlabel('running trial duration (s)');
    title([sessions{i}]);
    savefig([fig_dest sessions{i} '_1pk_time_Vs_duration']);
    
    save([image_analysis_dest sessions{i} '_dfOvF.mat'],'avespd_runofftrials','trial_lens','-append');
    
end


