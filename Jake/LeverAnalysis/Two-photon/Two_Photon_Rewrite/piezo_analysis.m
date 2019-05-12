%Script for loading and analyzing piezo data
%to be run after WF_CRP_licking_only_overview

%this current one should merely retreive the full piezo trace and align it
%to the imaging frames
clear 
bx_source      = ['Y:\home\jake\Data\WidefieldImaging\GCaMP\behavior\'];
bx_outputs_dir = ['Y:\home\jake\Analysis\Cue_reward_pairing_analysis\BxAndAnalysisOutputs\BxOutputs\'];
CRP_fig_dir_base    = ['Y:\home\jake\Analysis\Cue_reward_pairing_analysis\CRPFigureFolder\'];
img_dir = ['Y:\home\jake\Data\2P_imaging\'];
piezo_dir = ['Y:\home\jake\Analysis\Cue_reward_pairing_analysis\2P\piezo analysis\'];
time_before_ms = 2000; %defines the window around the cue presentation which will be taken for plotting
time_after_ms = 3000;
days = {'190213_img092', '190214_img092', '190215_img092', '190216_img092', '190217_img092', '190218_img092', '190219_img092', '190220_img092', '190221_img092', '190222_img092', '190223_img092', '190224_img092', '190225_img092', '190226_img092'};
days = {'190215_img093', '190216_img093', '190217_img093', '190218_img093', '190219_img093', '190220_img093', '190222_img093', '190225_img093', '190226_img093', '190227_img093', '190228_img093', '190301_img093', '190302_img093', '190304_img093'};

gap_inx = gap_day_finder(days);

for ii = 1:length(days)
    %check for piezo data
    data_dir = [img_dir, days{ii}, '\', days{ii}(end-5:end), '\'];
    piezo_file = dir(fullfile(data_dir,['*' '.ephys']));
    if isempty(piezo_file) | ~exist([data_dir, piezo_file.name])
        continue
    end
    
    %check for existing piezo variables and skip if they already exist
    if ~exist([piezo_dir, days{ii}(end-5:end)])
        mkdir([piezo_dir, days{ii}(end-5:end)]);
    end
    
    %check for existing piezo analysis variables. If they exist - skip to plotting
    if ~exist([piezo_dir, days{ii}(end-5:end), '\', days{ii}, '.mat'])
        
        %load bx data
        days(ii)
        bx_out_dir  = [bx_outputs_dir days{ii} '_bx_outputs'];
        b_data = get_bx_data(bx_source, days{ii});  %find the correct behavior file and loads it.
        
        % --- Obtain frame info from behavior file
        [counterValues, counterTimesMs, shift_info] = counter_fixer_2P_CRP(cell2mat(cellfun(@int64,b_data.counterValues,'UniformOutput',0)), cell2mat(cellfun(@int64,b_data.counterTimesUs,'UniformOutput',0))/1000, b_data.saveTime(1:6), b_data.subjectNum);
        counterInd = counter_calculator(counterValues, counterTimesMs);
        
        %use time of first frame to align licking times to start of imaging
        lickTimes=[];
        if isfield(b_data, 'lickometerTimesUs');  %img24 and 25 use datasets which have no licking fields
            for kk = 1:length(b_data.lickometerTimesUs);
                lickTimes = [lickTimes cell2mat(b_data.lickometerTimesUs(kk))/1000];
            end
            lickTimes = double(lickTimes);%-bx_start_MWorks_time;
            
            %index licking by frame number
            licksByFrame = [];
            for i = 1:length(counterTimesMs)-1;
                licksThisFrame = sum(lickTimes>counterTimesMs(i) & lickTimes<=counterTimesMs(i+1));
                licksByFrame = [licksByFrame licksThisFrame];
            end
            lick_frames = find(licksByFrame);
        end
        
        %Collects various events during session
        react_time = double(cell2mat(b_data.reactTimesMs));
        %collect event frame nums
        IFI = mode(diff(counterTimesMs));
        cue_frame = double(cell2mat(b_data.cTargetOn));
        rew_delay_frames = round( (react_time + double(b_data.RewardDelayDurationMs)) /IFI);
        rew_frame = cue_frame + rew_delay_frames;
        
        %identify reward omission trials and unexpected reward trials
        if b_data.rewardOmissionPercent == 0 %created this if statement because in 170417_img90 there were empty cells in b_data.rewardOmissionPercent which caused the script to fail
            reward_omit_inx = [];
        elseif sum(cellfun(@isempty, b_data.tRewardOmissionTrial)) > 0
            b_data.tRewardOmissionTrial{find(cellfun(@isempty, b_data.tRewardOmissionTrial))} = int64(0);
            reward_omit_inx = find(cell2mat(b_data.tRewardOmissionTrial(1:end-1))); %exclude last trial in case it is incomplete
        else
            reward_omit_inx = find(cell2mat(b_data.tRewardOmissionTrial(1:end-1))); %exclude last trial in case it is incomplete
        end
        unexp_rew_inx = find(cell2mat(b_data.tDoNoStimulusChange(1:end-1)));
        
        %load piezo data
        piezo_data_temp = fopen([data_dir, piezo_file.name]);
        piezo_data = fread(piezo_data_temp, 'single');
        piezo_frame_ind = piezo_data([1:2:end]);
        piezo_data = piezo_data([2:2:end]);
        
        %find piezo indeces for cue onset, reward, and licking events
        cue_inds = find(ismember(piezo_frame_ind, cue_frame));
        cue_inds = [cue_inds([find(diff(cue_inds)>10)])', cue_inds(end)]-1;
        rew_inds = find(ismember(piezo_frame_ind, rew_frame));
        rew_inds = [rew_inds([find(diff(rew_inds)>10)])', rew_inds(end)]-1;
        lick_inds = find(ismember(piezo_frame_ind, lick_frames));
        lick_inds = [lick_inds([find(diff(lick_inds)>10)])', lick_inds(end)]-1;
        
        %convert time before/after to 30Hz but then take 3x +1 that number from piezo data
        piezo_post = round(time_after_ms/b_data.frameRateHz)*3+1;
        piezo_pre = round(time_before_ms/b_data.frameRateHz)*3+1;
        move_by_trial = zeros(length(cue_inds)-1, piezo_post +piezo_pre +1); %dim1=trial# dim2=ms
        %plot in increments of 11.1111ms
        x_axis_piezo = [1:size(move_by_trial,2)]*11.1111;  %this code and some subqent lines assume 90Hz piezo setting
        x_axis_piezo = x_axis_piezo- x_axis_piezo([piezo_pre+1]);
        
        for trial_num = 2:length(cue_inds)-1
            move_by_trial(trial_num,:) = piezo_data([(cue_inds(trial_num)-piezo_pre):(cue_inds(trial_num)+piezo_post)]);
        end
        
        %identify trial types
        move_by_trial_NR = move_by_trial;
        move_by_trial_NR(sort([reward_omit_inx, unexp_rew_inx]) , :) = []; %remove any unexpected rewards or reward omission trials from the rewarded trials condition.
        move_by_trial_OR = move_by_trial(reward_omit_inx, :);
        move_by_trial_UR = move_by_trial(unexp_rew_inx, :);
        
        %save variables
        save([piezo_dir, days{ii}(end-5:end), '\', days{ii}], 'move_by_trial_NR', 'move_by_trial_OR', 'move_by_trial_UR', 'x_axis_piezo', 'lick_inds', 'rew_inds', 'cue_inds', 'piezo_data', 'piezo_frame_ind', ...
            'reward_omit_inx', 'unexp_rew_inx');
    end
    
    %% analyze movement traces
    if ~exist([piezo_dir, days{ii}(end-5:end), '\', days{ii}, '_mov_thresh.mat'])
        
        %calculate a floating mean and std dev
        window_size = 9000;
        load([piezo_dir, days{ii}(end-5:end), '\', days{ii}], 'piezo_data', 'rew_inds', 'cue_inds');
        [mov_means, mov_stds] = piezo_mov_thresh(piezo_data, window_size);
        
        %find all movement onset indeces  and  
        [all_mov_onsets, mov_inds, bout_buffer] = piezo_mov_onset(piezo_data, mov_means, mov_stds);
        
        %find movement onset latency relative to cue/reward
        trial_onset_lats = piezo_onset_lat_in_trial(mov_inds, reward_omit_inx, unexp_rew_inx, rew_inds, cue_inds);
        
        %save variables
        
    end
    
    %% plot mmovement traces
    if exist([piezo_dir, days{ii}(end-5:end), '\', days{ii}, '_trial_avg.fig']) & exist([piezo_dir, days{ii}(end-5:end), '\', days{ii}, '_full_session.fig'])
        continue
    end
    
    figure;
    subplot(1,2,1);
    plot(x_axis_piezo, mean(abs(move_by_trial_NR)));
    xlabel('time (ms) from cue onset');
    ylabel('movement reading (mean)');
    title(['nomral rewarded trials: n=', num2str(size(move_by_trial_NR,1))]);
    ylim([0 0.35]); xlim([-2000 3000]);
    vline(600, 'g'); vline(0, 'k');
    
    subplot(1,2,2);
    if isempty(move_by_trial_UR)
        plot(x_axis_piezo, mean(abs(move_by_trial_OR)));
        vline(600, 'r'); vline(0, 'k');
        title(['omitted reward trials: n=', num2str(size(move_by_trial_OR,1))]);
    else
        plot(x_axis_piezo, mean(abs(move_by_trial_UR)));
        vline(600, 'g'); vline(0, '--k');
        title(['unexpected reward trials: n=', num2str(size(move_by_trial_UR,1))]);
    end
    xlabel('time (ms) from cue onset');
    ylabel('movement reading (mean)');
    ylim([0 0.35]); xlim([-2000 3000]);
    suptitle([days{ii}(1:6), ' ', days{ii}(8:end)]);
    savefig([piezo_dir, days{ii}(end-5:end), '\', days{ii}, '_trial_avg']);
    
    %plot entire movement trace with vlines for cue/reward
    figure; plot(piezo_data); vline(cue_inds, 'k'); vline(rew_inds, 'g'); hold on;
    plot(lick_inds, ones(size(lick_inds))*0.075, 'LineStyle', 'none', 'Marker', 'o', 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r');
    title([days{ii}(1:6), ' ', days{ii}(8:end)]);
    savefig([piezo_dir, days{ii}(end-5:end), '\', days{ii}, '_full_session']);
    
end


