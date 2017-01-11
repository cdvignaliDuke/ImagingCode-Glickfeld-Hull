%% Modified version of the anaylysis overview for analyzing the licking behavior of cue-reward trials. 
%SECTION ONE - assign pathnames and datasets to be analyzed/written. 
clear;
%days = {'161212_img73', '161213_img73', '161214_img73', '161215_img73'};
days = {'161212_img74', '161214_img74', '161215_img74', '161216_img74'};
%days = {'161212_img73', '161212_img74'};
bx_source     = ['Z:\Data\WidefieldImaging\GCaMP\behavior\'];
bx_outputs_dir = ['Z:\Analysis\Cue_reward_pairing_analysis\BxAndAnalysisOutputs\BxOutputs\'];
old_cd = cd; %save old cd so I can restore it later
time_before_ms = 5000; %defines the window around the cue presentation which will be taken for plotting
time_after_ms = 10000;
%% SECTION TWO - collect the data 
%figure; hold on;
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
    if isfield(b_data, 'lickometerTimesUs');  %img24 and 25 use datasets which have no licking fields
        for kk = 1:length(b_data.lickometerTimesUs);
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
    cue_presentation = release_time-react_time;
 
    %save(bx_out_dir, 'lickTimes', '-append');
    
    %====================
    %isolate the time of cue onset and divide licking into trials as such
    licks_by_trial = zeros(length(cue_presentation)-1,(time_before_ms+time_after_ms+1)); %dim1=trial# dim2=ms 
    for kk = 1:length(cue_presentation)-1 %look at all trials except the last one.
        %find all the licking events Xms before and Yms after cue presentation
        licks_this_window = lickTimes(find(lickTimes>cue_presentation(kk)-time_before_ms & lickTimes<cue_presentation(kk)+time_after_ms));
        alignment_this_trial = cue_presentation(kk)-(time_before_ms+1); %subtract off this time so that the lick times are converted to numbers which will correspond to there index in licks_by_trial
        licks_this_window = licks_this_window - alignment_this_trial;
        licks_by_trial(kk, licks_this_window) = 1; 
    end
   cum_hist_licks = cumsum(sum(licks_by_trial));
   cum_hist_licks = cum_hist_licks/max(cum_hist_licks); %normalize to 1
   licks_by_trial_cum = cumsum(licks_by_trial,2);  % take a cumulative sum of the licks on individual trials and normalize to 1
   licks_by_trial_cum_norm = licks_by_trial_cum/max(max(licks_by_trial_cum));
   x_axis_range = [-1*time_before_ms:time_after_ms];
   figure; 
   plot(x_axis_range, cum_hist_licks, 'r', 'LineWidth', 3);  %plot average cum lick hist for t by t analysis 
   %plot(x_axis_range, cum_hist_licks, 'Color', [0,0,0]+(1-(ii/(length(days)+3))), 'LineWidth', 3); %plotting the average lick cum hist across days 
   %plot(x_axis_range, mean(licks_by_trial_cum_norm([1:10],:)), 'k'); hold on;
   %plot(x_axis_range, mean(licks_by_trial_cum_norm([end-30:end-2],:)), 'b');
   ylabel('cumulative fraction of licks');
   xlabel('time (ms) relative to release cue onset');
   title(['cumulative histogram of lick times ', days(ii)]);
 hold on;
   for kk = 1:7:size(licks_by_trial_cum,1)
       plot(x_axis_range, licks_by_trial_cum_norm(kk,:), 'Color', [0,0,0]+(1-(kk/size(licks_by_trial_cum,1))));
   end
   
end




