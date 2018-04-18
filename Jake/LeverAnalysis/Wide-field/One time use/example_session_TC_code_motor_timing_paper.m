%script used to generate an example timecourse from 151211_img32 for the
%motor timing paper. Meant to be run after running WF_lever_plotting_TCs
%for just that one session.
pre_frames = 40;
ts = [-pre_frames:post_frames].*100;
use_ev_success = trial_outcome.success_time;
[success_roi, use_times_succ, lick_trace_succ, lick_trace_succ_10ms] = trigger_movie_by_event_licks(tc_dfoverf, frame_info, licking_data, use_ev_success, pre_frames, post_frames);

use_ev_fail = trial_outcome.early_time;
[fail_roi, use_times_fail, lick_trace_fail, lick_trace_fail_10ms] = trigger_movie_by_event_licks(tc_dfoverf, frame_info, licking_data, use_ev_fail, pre_frames, post_frames);

%get hold duration stats
succ_hold_dur = trial_outcome.succ_hold_dur;
succ_hold_std = std(succ_hold_dur);
succ_hold_sem = succ_hold_std/length(succ_hold_dur);
succ_hold_mean = mean(succ_hold_dur);
fail_hold_dur = trial_outcome.fail_hold_dur;
fail_hold_std = std(fail_hold_dur);
fail_hold_sem = fail_hold_std/length(fail_hold_dur);
fail_hold_mean = mean(fail_hold_dur);

%calculate press times
press_times_succ = trial_outcome.success_time-succ_hold_dur;
press_times_fail = trial_outcome.early_time-fail_hold_dur;
[press_roi_succ, use_times_press_succ, lick_trace_press_succ, lick_trace_press_10ms_succ] = trigger_movie_by_event_licks(tc_dfoverf, frame_info, licking_data, press_times_succ, 20, 20);
[press_roi_fail, use_times_press_fail, lick_trace_press_fail, lick_trace_press_10ms_fail] = trigger_movie_by_event_licks(tc_dfoverf, frame_info, licking_data, press_times_fail, 20, 20);

%average across ROIs
avg_success_roi = squeeze(mean(success_roi,2));
%get SEM and average across trials 
std_success = std(avg_success_roi,1);
sm_success = std_success./sqrt(size(avg_success_roi,1));
avg_success_roi = mean(avg_success_roi,1);

%baseline each curve so it passes through zero
shift = (-1)*avg_success_roi(1);
avg_success_roi = avg_success_roi+shift;

%plotting
figure; subplot(1,2,2);
errorbar(ts, avg_success_roi, sm_success, 'k'); hold on;

%average across ROIs
avg_fail_roi = squeeze(mean(fail_roi,2));
%get SEM and average across trials 
std_fail = std(avg_fail_roi,1);
sm_fail = std_fail./sqrt(size(avg_fail_roi,1));
avg_fail_roi = mean(avg_fail_roi,1);

%baseline each curve so it passes through zero
shift = (-1)*avg_fail_roi(1);
avg_fail_roi = avg_fail_roi+shift;

%plotting
errorbar(ts, avg_fail_roi, sm_fail, 'r'); hold on;
ylim([-0.02 0.1]);
xlim([-4100 2025]);
xlabel('Time from release (ms)');
ylabel('dF/F');
title('Correct: black.   Early: red.');

%plot hold duration variables
% vline(-succ_hold_mean, 'k');
% vline(-fail_hold_mean, 'r');
% plot([-succ_hold_mean-succ_hold_std, -succ_hold_mean+succ_hold_std], [0.02, 0.02], 'k');
% plot([-succ_hold_mean-succ_hold_sem, -succ_hold_mean+succ_hold_sem], [0.022, 0.022], 'k');
% plot([-fail_hold_mean-fail_hold_std, -fail_hold_mean+fail_hold_std], [0.024, 0.024], 'r');
% plot([-fail_hold_mean-fail_hold_sem, -fail_hold_mean+fail_hold_sem], [0.026, 0.026], 'r');

%average across ROIs
press_roi_avg_succ = squeeze(mean(press_roi_succ,2));
press_roi_avg_fail = squeeze(mean(press_roi_fail,2));

%get SEM and average across trials 
std_press_succ = std(press_roi_avg_succ,1);
sm_press_succ = std_press_succ./sqrt(size(press_roi_avg_succ,1));
press_roi_avg_succ = mean(press_roi_avg_succ,1);
std_press_fail = std(press_roi_avg_fail,1);
sm_press_fail = std_press_fail./sqrt(size(press_roi_avg_fail,1));
press_roi_avg_fail = mean(press_roi_avg_fail,1);
ts_press = [-20:20]*100;

%baseline each curve so it passes through zero
shift = (-1)*press_roi_avg_succ(1);
press_roi_avg_succ = press_roi_avg_succ+shift;
shift = (-1)*press_roi_avg_fail(1);
press_roi_avg_fail = press_roi_avg_fail+shift;

%plot the lever press
subplot(1,2,1);
errorbar(ts_press, press_roi_avg_succ, sm_press_succ, 'k'); hold on;
errorbar(ts_press, press_roi_avg_fail, sm_press_fail, 'r'); hold on;
xlabel('Time from press (ms)');
ylabel('dF/F');
ylim([-0.02 0.1]);
xlim([-4100 2025]);
title('Lever press');

% %plot the subtraction
% subplot(1,2,2);
% avg_sub = avg_success_roi- avg_fail_roi;
% sub_sm = sqrt(sm_fail.^2+sm_success.^2);
% errorbar(ts, avg_sub, sub_sm, 'k');
% ylim([-0.02 0.1]);
% xlim([-4100 2025]);
% xlabel('Time from release (ms)');
% ylabel('dF/F');
% title('correct-early');






