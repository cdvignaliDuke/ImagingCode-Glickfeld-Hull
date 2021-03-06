%% Modified version of the anaylysis overview for analyzing the licking behavior of cue-reward trials. 
%SECTION ONE - assign pathnames and datasets to be analyzed/written. 
%clear;

bx_source = 'H:\home\mike\Data\Behavior\';
bx_outputs_dir = 'H:\home\mike\Analysis\BX Outputs\';
CRP_fig_dir_base = 'H:\home\mike\Analysis\CRP Figures\Sessions & Summaries\';


old_cd = cd; %save old cd so I can restore it later
time_before_ms = 2000; %defines the window around the cue presentation which will be taken for plotting
time_after_ms = 3000;

%days = {'191023_mm1505','191024_mm1505','191025_mm1505'}; 
days = {'191023_mm1506','191024_mm1506','191025_mm1506'};

%check and make sure the figure destinations exist
session_fig_dir = [CRP_fig_dir_base, days{1}(end-3:end), '_sessions\'];
sum_fig_dir = [CRP_fig_dir_base, days{1}(end-3:end), '_summary'];
if exist([CRP_fig_dir_base, days{1}(end-3:end), '_sessions'], 'file') == 0
    mkdir([CRP_fig_dir_base, days{1}(end-3:end), '_sessions']);
end
if exist([CRP_fig_dir_base, days{1}(end-3:end), '_summary'], 'file') == 0
    mkdir([CRP_fig_dir_base, days{1}(end-3:end), '_summary']);
end

gap_inx = gap_day_finder(days);

%% SECTION TWO - 
%figure; hold on;
avg_licks_post_cue = [];
avg_licks_pre_cue = [];
avg_licks_post_cue_sem = [];
avg_licks_pre_cue_sem = [];
RT_across_days = [];
RT_across_days_sem = [];
std_of_RT_across_days = [];
TFT_rates = [];
miss_rates = [];
RT_across_sessions = [];
RT_across_sessions_delay = [];
RT_across_sessions_1000ms_delay = [];
days_divider_inx = [];
days_divider_inx_delay = [];
days_divider_inx_1000ms_delay = [];
non_consecutive_inx = [];
non_consecutive_inx_delay = [];
non_consecutive_inx_1000ms_delay = [];
pre_cue_lick_window_avg = [];
pre_cue_lick_rate_sem = [];
iti_lick_window_avg = [];
iti_lick_rate_sem = [];
for ii = 1:length(days)
    days(ii)
    bx_out_dir  = [bx_outputs_dir days{ii} '_bx_outputs'];
    b_data = get_bx_data(bx_source, days{ii});  %find the correct behavior file and loads it.
    
    trial_start = round(double(cell2mat(b_data.tThisTrialStartTimeMs)));  %finds each trial's start time in free floating MWorks time
    num_trials  = length(b_data.tThisTrialStartTimeMs); %gets vectors for # of trials. Includes unimaged trials.
    
    % --- stores the time of the beginning of the first trial in MWorks time.
    bx_start_MWorks_time  = trial_start(1);
    
    %use time of first frame to align licking times to start of imaging
    lickTimes=[];
    if isfield(b_data, 'lickometerTimesUs')  %img24 and 25 use datasets which have no licking fields
        for kk = 1:length(b_data.lickometerTimesUs)
            lickTimes = [lickTimes cell2mat(b_data.lickometerTimesUs(kk))/1000];
        end
        lickTimes = double(lickTimes)-bx_start_MWorks_time;
    end
    
    %Collects various events during session
    hold_start = double(cell2mat(b_data.holdStartsMs)) - bx_start_MWorks_time;
    hold_time  = double(cell2mat(b_data.holdTimesMs));   %duration of the "lever hold" on that trial. meaningless here except to calculate cue onset
    react_time = double(cell2mat(b_data.reactTimesMs));
    req_hold   = double(cell2mat(b_data.tTotalReqHoldTimeMs));
    rnd_hold   = double(cell2mat(b_data.tRandReqHoldTimeMs));
    tot_req_hold = req_hold + rnd_hold;
    release_time = hold_start + hold_time;
    cue_presentation = release_time-react_time;   %=====================THIS MEANS CUM HISTs ARE ALIGNED TO CUE ONSET
    
    if exist(bx_out_dir)
        save(bx_out_dir, 'lickTimes', '-append');
    else
        save(bx_out_dir, 'lickTimes');
    end
    
    %identify reward omission trials and unexpected reward trials
    if b_data.rewardOmissionPercent == 0 %created this if statement because in 170417_img90 there were empty cells in b_data.rewardOmissionPercent which caused the script to fail
        reward_omit_inx = [];
    elseif sum(cell2mat(b_data.tRewardOmissionTrial)) == 1
         reward_omit_inx = [];
    elseif sum(cellfun(@isempty, b_data.tRewardOmissionTrial)) > 0
        b_data.tRewardOmissionTrial{find(cellfun(@isempty, b_data.tRewardOmissionTrial))} = int64(0);
        reward_omit_inx = find(cell2mat(b_data.tRewardOmissionTrial(1:end-1))); %exclude last trial in case it is incomplete
    else
        reward_omit_inx = find(cell2mat(b_data.tRewardOmissionTrial(1:end-1))); %exclude last trial in case it is incomplete
    end
    unexp_rew_inx = find(cell2mat(b_data.tDoNoStimulusChange(1:end-1)));
    
    %isolate the time of cue onset and divide licking into trials as such
    licks_by_trial = zeros(length(cue_presentation)-1,(time_before_ms+time_after_ms+1)); %dim1=trial# dim2=ms
    for kk = 1:length(cue_presentation)-1 %look at all trials except the last one.
        %find all the licking events Xms before and Yms after cue presentation
        licks_this_window = lickTimes(find(lickTimes>cue_presentation(kk)-time_before_ms & lickTimes<cue_presentation(kk)+time_after_ms));
        alignment_this_trial = cue_presentation(kk)-(time_before_ms+1); %subtract off this time so that the lick times are converted to numbers which will correspond to there index in licks_by_trial
        licks_this_window = licks_this_window - alignment_this_trial;
        licks_by_trial(kk, licks_this_window) = 1;
    end
    licks_by_trial_rewarded = licks_by_trial;
    licks_by_trial_rewarded(sort([reward_omit_inx, unexp_rew_inx]) , :) = []; %remove any unexpected rewards or reward omission trials from the rewarded trials condition. 
    licks_by_trial_omit = licks_by_trial(reward_omit_inx, :);
    licks_by_trial_unexp = licks_by_trial(unexp_rew_inx, :);
    
    %REWARD generate cumulative histograms of licking for rewarded trials
    avg_licks_per_ms_rew = cumsum(mean(licks_by_trial_rewarded));
    max_licks_for_norm = max(avg_licks_per_ms_rew);
    all_trials_lick_hist = cumsum(licks_by_trial_rewarded,2);
    all_trials_lick_hist = all_trials_lick_hist; %/max_licks_for_norm;
    
    %REWARD plot licking histograms for rewarded trials
    x_axis_range = [-1*time_before_ms:time_after_ms];
    if exist([session_fig_dir, days{ii}, '_rew_cum_hist.fig'], 'file')==0
        figure;
        plot(x_axis_range, avg_licks_per_ms_rew, 'r', 'LineWidth', 3);  %plot average cum lick hist for t by t analysis
        ylabel('cumulative # of licks');
        xlabel('time (ms) relative to release cue onset');
        title(['Rewarded Trials: cumulative hist of licking ', days{ii}(end-4:end), ' ' days{ii}(1:6) ' n=' num2str(size(licks_by_trial_rewarded,1))]);
        hold on;
        for kk = 1:6:size(all_trials_lick_hist,1)
            plot(x_axis_range, all_trials_lick_hist(kk,:), 'Color', [0,0,0]+(1-(kk/size(all_trials_lick_hist,1))));
        end
        vline((b_data.RewardDelayDurationMs + round(mean(react_time))), 'b');
        savefig([session_fig_dir, days{ii}, '_rew_cum_hist']);
    end
    
    if b_data.rewardOmissionPercent > 0
        %OMISSION generate cumulative histograms of licking for omission trials
        avg_licks_per_ms = cumsum(mean(licks_by_trial_omit));
        max_licks_for_norm = max(avg_licks_per_ms);
        all_trials_lick_hist = cumsum(licks_by_trial_omit,2); %/max_licks_for_norm;
        
        %OMISSION plot licking histograms for omission trials
        if exist([session_fig_dir, days{ii}, '_om_cum_hist.fig'], 'file')==0
            figure;
            plot(x_axis_range, avg_licks_per_ms, 'r', 'LineWidth', 3);  %plot average cum lick hist for t by t analysis
            ylabel('cumulative # of licks');
            xlabel('time (ms) relative to release cue onset');
            title(['Reward Omission Trials: cumulative hist of licking ', days{ii}(end-4:end), ' ', days{ii}(1:6) ' n=' num2str(size(licks_by_trial_omit,1))]);
            hold on;
            for kk = 1:size(all_trials_lick_hist,1)
                plot(x_axis_range, all_trials_lick_hist(kk,:), 'Color', [0,0,0]+(1-(kk/size(all_trials_lick_hist,1))));
            end
            savefig([session_fig_dir, days{ii}, '_om_cum_hist']);
        end
    end
    
    if b_data.rewardUnexpectPercent > 0
        %UNEXPECTED reward cumulative histograms
        avg_licks_per_ms_unexp = cumsum(mean(licks_by_trial_unexp));
        max_licks_for_norm = max(avg_licks_per_ms_unexp);
        all_trials_lick_hist = cumsum(licks_by_trial_unexp,2); %/max_licks_for_norm;
        
        %UNEXPECTED plot licking histograms for omission trials
        if exist([session_fig_dir, days{ii}, '_unexp_cum_hist.fig'], 'file')==0
            figure;
            plot(x_axis_range, avg_licks_per_ms_unexp, 'r', 'LineWidth', 3);  %plot average cum lick hist for t by t analysis
            ylabel('cumulative # of licks');
            xlabel(['reward delivered at ', num2str(b_data.RewardDelayDurationMs), 'ms']);
            title(['Unexpected Reward Trials: cumulative hist of licking ', days{ii}(end-4:end), ' ', days{ii}(1:6) ' n=' num2str(size(licks_by_trial_unexp,1))]);
            hold on;
            for kk = 1:size(all_trials_lick_hist,1)
                plot(x_axis_range, all_trials_lick_hist(kk,:), 'Color', [0,0,0]+(1-(kk/size(all_trials_lick_hist,1))));
            end
            savefig([session_fig_dir, days{ii}, '_unexp_cum_hist']);
        end
    end
    
    if exist([session_fig_dir, days{ii}, '_rew_lick_raster.fig'], 'file')==0
        %REWARDED trials lick raster plot
        figure;
        plotSpikeRaster(logical(licks_by_trial), 'PlotType', 'vertline');
        vline(time_before_ms+1, 'k');
        vline(time_before_ms+1+(b_data.RewardDelayDurationMs + round(mean(react_time))), 'b');
        ylabel('trial # (descending)');
        xlabel('time (ms) black=cue blue=reward');
        title(['Rewarded Trials: lick time raster ', days{ii}(end-4:end), ' ' days{ii}(1:6) ' n=' num2str(size(licks_by_trial_rewarded,1))]);
        savefig([session_fig_dir, days{ii}, '_rew_lick_raster']);
    end
    
    %define the windows for summing licks before and after cue
    lick_window = [time_before_ms+201:time_before_ms+700];
    pre_cue_window = [time_before_ms-700:time_before_ms-201];
    
    %determine # of licks in 500ms window before/following cue presentation for reward omission trials
    if b_data.rewardOmissionPercent >1
        avg_licks_om_trials = mean(licks_by_trial_omit); %Averages the number of licks per ms across all reward omission trials. Size is a vector with length time_before_cue+time_after_cue.
        num_om_trials = size(licks_by_trial_omit,1);
        licks_pre_window_by_trial = sum(licks_by_trial_omit(:,pre_cue_window) ,2); %sums the total # of licks within the window for each trial. Size is a vector with length = # reward omission trials.
        licks_lick_window_by_trial = sum(licks_by_trial_omit(:,lick_window) ,2);
        avg_licks_post_cue_this_day = sum(avg_licks_om_trials(lick_window));  %determine the avg number (across trials) of licks which occur throughout the analysis window for this session
        avg_licks_pre_cue_this_day = sum(avg_licks_om_trials(pre_cue_window));
        avg_licks_post_cue = [avg_licks_post_cue, avg_licks_post_cue_this_day]; %concatenate number of licks in analysis window for each session
        avg_licks_pre_cue = [avg_licks_pre_cue, avg_licks_pre_cue_this_day];
        avg_licks_post_cue_sem = [avg_licks_post_cue_sem, std(licks_lick_window_by_trial)/sqrt(num_om_trials)]; %calculate and store standard error of the mean
        avg_licks_pre_cue_sem  = [avg_licks_pre_cue_sem, std(licks_pre_window_by_trial)/sqrt(num_om_trials)];
    else 
        avg_licks_post_cue = [avg_licks_post_cue, NaN]; %this allows the training day # to be aligned to the pre/post cue lick rate for that day across animals
        avg_licks_pre_cue = [avg_licks_pre_cue, NaN];
        avg_licks_post_cue_sem = [avg_licks_post_cue_sem, NaN]; 
        avg_licks_pre_cue_sem  = [avg_licks_pre_cue_sem, NaN];
    end
    
    %determine lick rate in the middle of the iti, just before the cue, and after reward delivery. 
    avg_licks_by_trial_rewarded = mean(licks_by_trial_rewarded);
    
    %Determin the RTs for each trial this session
    RT_this_session=[];
    for kk = 1:num_trials-1 %look at each trial
         this_trial_RT = find(licks_by_trial(kk,[time_before_ms+1:end]), 1, 'first'); %look at all time points starting from cue onset for each trial, find the first lick
          RT_this_session = [RT_this_session, this_trial_RT]; %store the RTs
%         this_trial_RT = find(licks_by_trial(kk,[time_before_ms+1:end]), 3, 'first'); %look at all time points starting from cue onset for each trial, find the first lick
%         if this_trial_RT(1) > 175
%             RT_this_session = [RT_this_session, this_trial_RT(1)]; %store the RTs
%         elseif  this_trial_RT(2) > 175
%             RT_this_session = [RT_this_session, this_trial_RT(2)]; %store the RTs
%         else
%             RT_this_session = [RT_this_session, this_trial_RT(3)]; %store the RTs
%         end
    end
    RT_this_session_raw = RT_this_session;
    save(['H:\home\mike\Analysis\CRP Figures\RT Data\', days{ii}], 'RT_this_session_raw', 'time_before_ms');
%     

    %trim RT data to exclude misses and FAs  ==================================================================
    TFT_rate_this_session = length(find(RT_this_session<200))/(num_trials-1);
    TFT_rates = [TFT_rates, TFT_rate_this_session];
    if b_data.RewardDelayDurationMs == 500 && b_data.rewardDelayPercent == 100
        miss_rate_this_session = length(find(RT_this_session > 1000))/(num_trials-1);
        RT_this_session = RT_this_session((RT_this_session > 200 & RT_this_session < 1000));
    elseif b_data.RewardDelayDurationMs == 1000 && b_data.rewardDelayPercent == 100
        miss_rate_this_session = length(find(RT_this_session > 1500))/(num_trials-1);
        RT_this_session = RT_this_session((RT_this_session > 200 & RT_this_session < 1500));
    elseif b_data.RewardDelayDurationMs == 450 && b_data.rewardDelayPercent == 100
        miss_rate_this_session = length(find(RT_this_session > 1000))/(num_trials-1);
        RT_this_session = RT_this_session((RT_this_session > 200 & RT_this_session < 1000));
    else
        miss_rate_this_session = length(find(RT_this_session>500))/(num_trials-1);
        RT_this_session = RT_this_session((RT_this_session > 200 & RT_this_session < 500));
    end
    
    %Determine the FA and miss rates for this session. Store stats for summary statistics
    miss_rates = [miss_rates, miss_rate_this_session];
    RT_across_days = [RT_across_days, mean(RT_this_session)];
    RT_across_days_sem = [RT_across_days_sem, std(RT_this_session)/sqrt(size(licks_by_trial,1))];
    std_of_RT_across_days = [std_of_RT_across_days, std(RT_this_session)];
    
    %plot RT within days relative to cue onset 
    if exist([session_fig_dir, days{ii}, '_RT_plot.fig'], 'file') == 0
        figure; plot(RT_this_session);
        title(['RT for individual trials within for day ' days{ii}(1:6), days{ii}(end-4:end)]);
        xlabel('trial #');
        ylabel('RT (ms)');
        savefig([session_fig_dir, days{ii}, '_RT_plot']);
    end
    
    %store trial-by-trial RTs for across days plot 
    if b_data.rewardDelayPercent == 0
        RT_across_sessions = [RT_across_sessions, RT_this_session];
    elseif b_data.rewardDelayPercent == 100 && b_data.RewardDelayDurationMs == 500
        RT_across_sessions_delay = [RT_across_sessions_delay, RT_this_session];
    elseif b_data.rewardDelayPercent == 100 && b_data.RewardDelayDurationMs == 1000
        RT_across_sessions_1000ms_delay = [RT_across_sessions_1000ms_delay, RT_this_session];
    end
    
    %determine trial index for trials different days and non-consecutive days
    if b_data.rewardDelayPercent == 0
        days_divider_inx = [days_divider_inx, length(RT_across_sessions)+0.5];
        if ismember([ii+0.5], gap_inx)
            non_consecutive_inx = [non_consecutive_inx, length(RT_across_sessions)];
        end
    elseif b_data.rewardDelayPercent == 100 && b_data.RewardDelayDurationMs == 500
        days_divider_inx_delay = [days_divider_inx_delay, length(RT_across_sessions_delay)+0.5];
        if ismember([ii+0.5], gap_inx)
            non_consecutive_inx_delay = [non_consecutive_inx_delay, length(RT_across_sessions_delay)];
        end
        elseif b_data.rewardDelayPercent == 100 && b_data.RewardDelayDurationMs == 1000
          days_divider_inx_1000ms_delay = [days_divider_inx_1000ms_delay, length(RT_across_sessions_1000ms_delay)+0.5];
        if ismember([ii+0.5], gap_inx)
            non_consecutive_inx_1000ms_delay = [non_consecutive_inx_1000ms_delay, length(RT_across_sessions_1000ms_delay)];
        end  
    end
    
    %----bin licking by 100ms windows relative to cue delivery----
    %determine  bin size and window size
    bin_size = 100; %number of ms to bin licking. 
    trial_start = trial_start - bx_start_MWorks_time;
    min_start_to_cue = min([cue_presentation-trial_start]);
    min_cue_to_end = min([trial_start(2:end)-cue_presentation(1:end-1)]);
    if min_start_to_cue > 20000
        pre_cue_window_lick = 20000;
    else 
        pre_cue_window_lick = floor(min_start_to_cue/bin_size)*bin_size;
    end
    if min_cue_to_end > 6000
        post_cue_window_lick = 5999;
    else 
        post_cue_window_lick = floor(min_cue_to_end/bin_size)*bin_size-1;
    end
    pre_cue_window_lick = 18000;
    post_cue_window_lick = 5999;
    
    %get lick traces (1ms resolution) for rewarded and omission trials
    full_trial_licks = zeros(length(cue_presentation)-1,(pre_cue_window_lick+post_cue_window_lick+1)); %dim1=trial# dim2=ms
    for kk = 1:length(cue_presentation)-1 %look at all trials except the last one.
        %find all the licking events Xms before and Yms after cue presentation
        licks_this_window = lickTimes((lickTimes>cue_presentation(kk)-pre_cue_window_lick & lickTimes<cue_presentation(kk)+post_cue_window_lick));
        alignment_this_trial = cue_presentation(kk)-(pre_cue_window_lick+1); %subtract off this time so that the lick times are converted to numbers which will correspond to there index in licks_by_trial
        licks_this_window = licks_this_window - alignment_this_trial;
        full_trial_licks(kk, licks_this_window) = 1;
    end
    full_trial_licks_rewarded = full_trial_licks;
    full_trial_licks_rewarded(sort([reward_omit_inx, unexp_rew_inx]) , :) = []; %exclude omission and unexpected trials from rewarded trials condition
    full_trial_licks_omit = full_trial_licks(reward_omit_inx, :);
    full_trial_licks_unexp = full_trial_licks(unexp_rew_inx, :);
    
    %bin licking by 50ms bins and convert to licks/sec
    full_trial_licks_rewarded_sum = sum(full_trial_licks_rewarded,1);
    full_trial_licks_omit_sum = sum(full_trial_licks_omit,1);
    full_trial_licks_unexp_sum = sum(full_trial_licks_unexp,1);
    full_trial_licks_rewarded_bin = zeros(1,(length(full_trial_licks_rewarded_sum)/bin_size));
    full_trial_licks_omit_bin = zeros(1,(length(full_trial_licks_omit_sum)/bin_size));
    full_trial_licks_unexp_bin = zeros(1,(length(full_trial_licks_unexp_sum)/bin_size));
    cue_presentation_binned = (pre_cue_window_lick/bin_size)+1;
    iii=1;
    for kk = 1:bin_size:length(full_trial_licks_rewarded_sum)
        iii=iii+1;
        full_trial_licks_rewarded_bin(iii) = sum(full_trial_licks_rewarded_sum(kk:(kk+bin_size-1)));
        full_trial_licks_omit_bin(iii) = sum(full_trial_licks_omit_sum(kk:(kk+bin_size-1)));
        full_trial_licks_unexp_bin(iii) = sum(full_trial_licks_unexp_sum(kk:(kk+bin_size-1)));
    end
    %plot rewarded trials
    full_trial_licks_rewarded_bin = (full_trial_licks_rewarded_bin/size(full_trial_licks_rewarded,1))*(1000/bin_size); % convert to lick rate in Hz
    x_axis_bin = ([1:length(full_trial_licks_rewarded_bin)]-cue_presentation_binned)*(bin_size/1000);
    save(['H:\home\mike\Analysis\CRP Figures\Hist Data Across Animals\', days{ii}, 'rew_hist'], 'full_trial_licks_rewarded_bin');
    if exist([session_fig_dir, days{ii}, '_rew_lick_hist.fig'], 'file') == 0
        figure; bar(x_axis_bin, full_trial_licks_rewarded_bin);
        title(['Rewarded trials: lick rate per ', num2str(bin_size), 'ms', days{ii}(1:6), days{ii}(8:end), 'n=', num2str(size(full_trial_licks_rewarded,1))]);
        xlabel('time (s) relative to cue onset');
        ylabel('lick rate (Hz)');
        vline(0,'k');
        savefig([session_fig_dir, days{ii}, '_rew_lick_hist']);
    end
    %plot reward omission trials
    full_trial_licks_omit_bin = (full_trial_licks_omit_bin/size(full_trial_licks_omit,1))*(1000/bin_size); % convert to lick rate in Hz
    if b_data.rewardOmissionPercent >0
        if exist([session_fig_dir, days{ii}, '_om_lick_hist.fig'], 'file') == 0
            figure; bar(x_axis_bin, full_trial_licks_omit_bin);
            title(['Reward omission: lick rate per ', num2str(bin_size), 'ms', ' ', days{ii}(1:6), ' ', days{ii}(8:end), 'n=', num2str(size(full_trial_licks_omit,1))]);
            xlabel('time (s) relative to cue onset');
            ylabel('lick rate (Hz)');
            vline(0,'k');
            savefig([session_fig_dir, days{ii}, '_om_lick_hist']);
        end
        save(['H:\home\mike\Analysis\CRP Figures\Hist Data Across Animals\', days{ii}, 'rew_om_hist'], 'full_trial_licks_omit_bin');
    end
    %plot unexpected reward trials
    full_trial_licks_unexp_bin = (full_trial_licks_unexp_bin/size(full_trial_licks_unexp,1))*(1000/bin_size); % convert to lick rate in Hz
    if b_data.rewardUnexpectPercent >0
        if exist([session_fig_dir, days{ii}, '_unexp_lick_hist.fig'], 'file') == 0
            figure; bar(x_axis_bin, full_trial_licks_unexp_bin);
            title(['Unexpected reward: lick rate per ', num2str(bin_size), 'ms', days{ii}(1:6),' ', days{ii}(8:end), 'n=', num2str(size(full_trial_licks_unexp,1))]);
            xlabel(['Reward delivered at ' num2str(b_data.RewardDelayDurationMs), 'ms']);
            ylabel('lick rate (Hz)');
            vline(b_data.RewardDelayDurationMs,'b');
            savefig([session_fig_dir, days{ii}, '_unexp_lick_hist']);
        end
    end
%     
    %store avg lick values in 500ms window just before cue on rewarded trials and during the ITI. 
    pre_cue_lick_window = full_trial_licks_rewarded(:,[pre_cue_window_lick-500:pre_cue_window_lick]);
    pre_cue_lick_window_sum = sum(pre_cue_lick_window,2);
    pre_cue_lick_window_avg_this_session = mean(pre_cue_lick_window_sum);
    pre_cue_lick_rate_sem_this_session = std(pre_cue_lick_window_sum)/sqrt(size(pre_cue_lick_window,1));
    pre_cue_lick_window_avg = [pre_cue_lick_window_avg ,pre_cue_lick_window_avg_this_session];
    pre_cue_lick_rate_sem = [pre_cue_lick_rate_sem ,pre_cue_lick_rate_sem_this_session];
    %now do it for mid iti lick window
    mid_iti = round(pre_cue_window_lick/2);
    iti_lick_window = full_trial_licks_rewarded(:,[mid_iti-500:mid_iti]);
    iti_lick_window_sum = sum(iti_lick_window,2);
    iti_lick_window_avg_this_session = mean(iti_lick_window_sum);
    iti_lick_rate_sem_this_session = std(iti_lick_window_sum)/sqrt(size(iti_lick_window,1));
    iti_lick_window_avg = [iti_lick_window_avg , iti_lick_window_avg_this_session];
    iti_lick_rate_sem = [iti_lick_rate_sem , iti_lick_rate_sem_this_session];

   
end 

%save the pre/post cue data to be use in across animals summary
save(['H:\home\mike\Analysis\CRP Figures\Post_Cue_Lick Across Animals\', days{ii}(end-4:end), 'rew_om_post_cue_lick'], 'avg_licks_post_cue', 'avg_licks_pre_cue', 'avg_licks_post_cue_sem', 'avg_licks_pre_cue_sem'); 
avg_licks_post_cue = avg_licks_post_cue(~isnan(avg_licks_post_cue)); %this allows the training day # to be aligned to the pre/post cue lick rate for that day across animals
avg_licks_pre_cue = avg_licks_pre_cue(~isnan(avg_licks_pre_cue));
avg_licks_post_cue_sem = avg_licks_post_cue_sem(~isnan(avg_licks_post_cue_sem));
avg_licks_pre_cue_sem  = avg_licks_pre_cue_sem(~isnan(avg_licks_pre_cue_sem));

%save RT std, FA rate and miss rate
save(['H:\home\mike\Analysis\CRP Figures\RT Stats Across Animals\', days{ii}(end-4:end), 'RTstd_FA_misses'], 'std_of_RT_across_days', 'TFT_rates', 'miss_rates');

% Plot summary statistics 
if b_data.rewardOmissionPercent > 1
    figure; subplot(1,2,2); 
    bar(avg_licks_post_cue);
    hold on; errorbar(avg_licks_post_cue, avg_licks_post_cue_sem);
    title('Omission trials: avg # of licks in 500ms window post cue across days');
    xlabel('day');
    ylabel('# of licks');
    
    subplot(1,2,1); bar(avg_licks_pre_cue);
    hold on; errorbar(avg_licks_pre_cue, avg_licks_pre_cue_sem);
    title('Omission trials: avg # of licks in 500ms window before cue across days');
    xlabel('day');
    ylabel('# of licks');
    max_y_val = max([avg_licks_post_cue, avg_licks_pre_cue]);
    if max_y_val > 0
        subplot(1,2,1); ylim([0 ceil(max_y_val)]);
        subplot(1,2,2); ylim([0 ceil(max_y_val)]);
    end
    savefig([sum_fig_dir, '\', days{ii}(end-4:end), '_licks_pre_vs_post_cue']);
end

figure; subplot(2,2,1);
sgtitle(days{ii}(end-4:end));
errorbar(RT_across_days, RT_across_days_sem);
title(['RT across days']);
xlabel('day');
ylabel('RT (ms)');

subplot(2,2,2); plot(std_of_RT_across_days);
title('std of RT across days');
xlabel('day');
ylabel('standard deviation');

subplot(2,2,3); plot(TFT_rates);
title('too fast lick rates as a % of all trials by day: FA = RT<200ms');
xlabel('day');
ylabel('% false alarms');
ylim([0 1]);

subplot(2,2,4); plot(miss_rates);
title('Miss rates as a % of all trials by day: misses = RT>1000ms or 1500ms');
xlabel('day');
ylabel('% misses');
ylim([0 1]);
savefig([sum_fig_dir, '\', days{ii}(end-4:end), '_RT_RTstd_misses_TFT']);

if length(days_divider_inx) > 0
    figure; plot(RT_across_sessions);
    title(['RTs across days. No delay. ' days{ii}(end-4:end)]);
    xlabel('trial num');
    ylabel('RT(ms)');
    vline(days_divider_inx, 'k');
    if ~isempty(non_consecutive_inx)
        vline(non_consecutive_inx, 'r');
    end
    savefig([sum_fig_dir, '\', days{ii}(end-4:end), '_RT_0ms_summary']);
end

if length(RT_across_sessions_delay) > 0
    figure; plot(RT_across_sessions_delay);
    title(['RTs across days. 500ms delay. ' days{ii}(end-4:end)]);
    xlabel('trial num');
    ylabel('RT(ms)');
    vline(days_divider_inx_delay, 'k');
    if ~isempty(non_consecutive_inx_delay)
        vline(non_consecutive_inx_delay, 'r');
    end
    savefig([sum_fig_dir, '\', days{ii}(end-4:end), '_RT_500ms_summary']);
end

if length(RT_across_sessions_1000ms_delay) > 0
    figure; plot(RT_across_sessions_1000ms_delay);
    title(['RTs across days. 1000ms delay. ' days{ii}(end-4:end)]);
    xlabel('trial num');
    ylabel('RT(ms)');
    vline(days_divider_inx_1000ms_delay, 'k');
    if ~isempty(non_consecutive_inx_1000ms_delay)
        vline(non_consecutive_inx_1000ms_delay, 'r');
    end
    savefig([sum_fig_dir, '\', days{ii}(end-4:end), '_RT_1000ms_summary']);
end

%plot # of licks in iti window vs pre-cue window
figure; errorbar(iti_lick_window_avg, iti_lick_rate_sem, 'k');
hold on; 
errorbar(pre_cue_lick_window_avg, pre_cue_lick_rate_sem, 'b');
title('licking in iti (black) vs pre-cue (blue)');
ylabel('avg # of licks in 500ms window');
xlabel('day #')
savefig([sum_fig_dir, '\', days{ii}(end-4:end), '_iti_vs_pre-cue_licking']);

%save variables for across animals RT plot

out_dir = 'H:\home\mike\Analysis\CRP Figures\RT Across Animals\';
save([out_dir, days{ii}(end-4:end), 'RT_sem_across_days'], 'RT_across_days', 'RT_across_days_sem');

