% triggers time-courses off of event times
% 1. finds frame and lever times based on events and outcomes
% 2. obtain df/f timecourse
% 3. create event triggered movies



%load frame and lever info
frame_info_dest = [dest '_frame_times.mat'];
load(frame_info_dest);
ftimes.frame_times = frame_times; clear frame_times;
b_data.input = input; clear input;

%% 1. find frame and lever times
ifi = (ftimes.frame_times(end)-ftimes.frame_times(1))/length(ftimes.frame_times);
Sampeling_rate = 1000/ifi;

if(~exist('first_frame', 'var'))
    f_frame =1;
else
    f_frame = first_frame;
end

if(~exist('last_frame', 'var'))
    l_frame =length(ftimes.frame_times);
else
    l_frame = last_frame;
end
% ---- parse behavior
holdT_min  = 500000;
[lever, frame_info, trial_outcome] = parse_behavior_for_HAD(b_data.input, ...
    f_frame, l_frame, ftimes.frame_times, holdT_min);

data_dest = [dest '_parse_behavior.mat'];
save(data_dest, 'lever', 'frame_info', 'trial_outcome', 'Sampeling_rate', 'holdT_min')

%% 2. Obtain a df/f TC from baseline times
data_tc = data_tc';
startT = round(b_data.input.counterTimesUs{1}(1)./1000);
tc_dfoverf = zeros(size(data_tc));    %this could be problematic due to the frame skipping issue
first_baseline = find(~isnan(lever.baseline_timesMs(1,:)),1, 'first');    %find the first trial / baseline_timesMs window that is not NaN
F_range = [];
for iT=2:length(lever.baseline_timesMs)-1;    %this could be problematic due to unremoved NaNs
    if ~isnan(lever.baseline_timesMs(1,iT));
        F_range = frame_info.counter(lever.baseline_timesMs(1,iT)):frame_info.counter(lever.baseline_timesMs(2,iT));
    elseif isempty(F_range)
        F_range = frame_info.counter(lever.baseline_timesMs(1,first_baseline)):frame_info.counter(lever.baseline_timesMs(2,first_baseline));
    end
    F_avg= mean(data_tc(:,F_range),2);
    t_range = frame_info.counter(cell2mat(b_data.input.tThisTrialStartTimeMs(iT))-startT):frame_info.counter(cell2mat(b_data.input.tThisTrialStartTimeMs(iT+1))-startT);
    t_df = bsxfun(@minus, double(data_tc(:,t_range)), F_avg);
    t_dfoverf = bsxfun(@rdivide, t_df, F_avg);
    tc_dfoverf(:,t_range) = t_dfoverf;
end 

%% create event triggered movies
func = @mean;
pre_release_frames = 5;
post_release_frames = 10;

ts = (-pre_release_frames:post_release_frames)*1000/round(double(Sampeling_rate));
tot_frame = pre_release_frames + post_release_frames+1;

%successes
use_ev_success = trial_outcome.success_time;
if strcmp(b_data.input.trialOutcomeCell{1}, 'success')
    use_ev_success(1) = [];
elseif strcmp(b_data.input.trialOutcomeCell{end}, 'success')
    use_ev_success(end) = [];
end

success_movie = trigger_movie_by_event(tc_dfoverf, frame_info, ...
    use_ev_success, pre_release_frames, post_release_frames);
avg_success = squeeze(func(success_movie,1));
sem_success = squeeze(std(success_movie,1)./sqrt(size(success_movie,1)));

avg_success_all = mean(avg_success,1);
sem_success_all = std(avg_success,1)./sqrt(size(avg_success,1));

%failures
use_ev_fail = trial_outcome.early_time;
if strcmp(b_data.input.trialOutcomeCell{1}, 'failure')
    use_ev_fail(1) = [];
elseif strcmp(b_data.input.trialOutcomeCell{end}, 'failure')
    use_ev_fail(end) = [];
end

% -----trigger movie by early release
fail_movie = trigger_movie_by_event(tc_dfoverf, frame_info, ...
    use_ev_fail, pre_release_frames, post_release_frames);
avg_fail = squeeze(func(fail_movie,1));
sem_fail = squeeze(std(fail_movie,1)./sqrt(size(fail_movie,1)));

%average of all ROIs
avg_fail_all = mean(avg_fail,1);
sem_fail_all = std(avg_fail,1)./sqrt(size(avg_fail,1));
tt =((-pre_release_frames:post_release_frames).*double(ifi))./1000;
figure; errorbar(tt,avg_success_all, sem_success_all,'k')
hold on;
errorbar(tt,avg_fail_all, sem_fail_all,'r')
title(['Average release: ' tc_type ' Success- black; Failure- red; n = ' num2str(size(avg_success,1)) ' cells'])
print([dest '_release_avgTCs_' tc_type '.eps'], '-depsc');
print([dest '_release_avgTCs_' tc_type '.pdf'], '-dpdf');

%average by ROI
nCells = size(data_tc,1);
z = ceil(sqrt(nCells));
figure;
avg_all = [avg_success avg_fail];
ymax = max(max(avg_all,[],2),[],1);
ymin = min(min(avg_all,[],2),[],1);
for ic = 1:nCells
    subplot(z,z,ic)
    errorbar(tt,avg_success(ic,:), sem_success(ic,:),'k')
    hold on;
    errorbar(tt,avg_fail(ic,:), sem_fail(ic,:),'r')
    ylim([ymin*1.1 ymax*1.1])
    xlim([tt(1) tt(end)])
end
suptitle(['Average release: ' tc_type ' Success- black (n = ' num2str(size(success_movie,1)) ' trials); Failure- red (n = ' num2str(size(fail_movie,1)) ' trials)'])
orient landscape
print([dest '_release_avg_allTCs_' tc_type '.eps'], '-depsc');
print([dest '_release_avg_allTCs_' tc_type '.pdf'], '-dpdf');
save([dest '_release_resp_by_outcome_' tc_type '.mat'],'fail_movie','success_movie','pre_release_frames','post_release_frames','ifi');

% ---- Trigger movie off all lever presses at trial start
pre_press_frames = 10;
post_press_frames = 10;

pressTime = NaN(1,length(lever.baseline_timesMs));
releaseTime = NaN(1,length(lever.baseline_timesMs));
for iT = 2:length(lever.baseline_timesMs)-1
    leverTimes = round((cell2mat(b_data.input.leverTimesUs(iT))-b_data.input.counterTimesUs{1}(1))./1000);
    if ~isnan(trial_outcome.ind_press_prerelease(iT))
        pressTime(1,iT) = leverTimes(trial_outcome.ind_press_prerelease(iT));
        releaseTime(1,iT) = leverTimes(trial_outcome.ind_press_prerelease(iT)+1);
    else
        pressTime(1,iT) = NaN;
        releaseTime(1,iT) = NaN;
    end
end

use_ev_press = pressTime;
use_ev_press(isnan(use_ev_press)) = [];

press_movie = trigger_movie_by_event(tc_dfoverf, frame_info, ...
    use_ev_press, pre_press_frames, post_press_frames);
avg_press = squeeze(func(press_movie,1));
sem_press = squeeze(std(press_movie,1)./sqrt(size(press_movie,1)));
avg_press_all = squeeze(func(avg_press,1));
sem_press_all = squeeze(std(avg_press,1)./sqrt(size(avg_press,1)));
figure;
tt =((-pre_press_frames:post_press_frames).*double(ifi))./1000;
errorbar(tt, avg_press_all, sem_press_all, '-k')
title(['Initiating press- ' tc_type ' all trials: n = ' num2str(size(press_movie,1))])
print([dest '_press_avgTCs_alltrials_' tc_type '.eps'], '-depsc');
print([dest '_press_avgTCs_alltrials_' tc_type '.pdf'], '-dpdf');

figure;
avg_all = [avg_press];
ymax = max(max(avg_all,[],2),[],1);
ymin = min(min(avg_all,[],2),[],1);
for ic = 1:nCells
    subplot(z,z,ic)
    errorbar(tt,avg_press(ic,:), sem_press(ic,:),'k')
    ylim([ymin*1.1 ymax*1.1])
    xlim([tt(1) tt(end)])
end
suptitle(['Initiating press- ' tc_type ' all trials: n = ' num2str(size(press_movie,1))])
orient landscape
print([dest '_press_avg_allTCs_alltrials_' tc_type '.eps'], '-depsc');
print([dest '_press_avg_allTCs_alltrials_' tc_type '.pdf'], '-dpdf');

%break up presses by hold time
holdTime = releaseTime-pressTime;
ind_200 = find(holdTime<200);
ind_500 = find(holdTime<500);
ind_both = find(ismember(ind_500,ind_200));
ind_500(ind_both) = [];
ind_long = find(holdTime>=500);
ind_both = find(ismember(ind_long,ind_500));
ind_long(ind_both) = [];

use_200_press = pressTime(ind_200);
use_500_press = pressTime(ind_500);
use_long_press = pressTime(ind_long);

press_200_movie = trigger_movie_by_event(tc_dfoverf, frame_info, ...
    use_200_press, pre_press_frames, post_press_frames);
press_500_movie = trigger_movie_by_event(tc_dfoverf, frame_info, ...
    use_500_press, pre_press_frames, post_press_frames);
press_long_movie = trigger_movie_by_event(tc_dfoverf, frame_info, ...
    use_long_press, pre_press_frames, post_press_frames);
avg_200_press = squeeze(func(press_200_movie,1));
sem_200_press = squeeze(std(press_200_movie,1)./sqrt(size(press_200_movie,1)));
avg_200_press_all = squeeze(func(avg_200_press,1));
sem_200_press_all = squeeze(std(avg_200_press,1)./sqrt(size(avg_200_press,1)));
avg_500_press = squeeze(func(press_500_movie,1));
sem_500_press = squeeze(std(press_500_movie,1)./sqrt(size(press_500_movie,1)));
avg_500_press_all = squeeze(func(avg_500_press,1));
sem_500_press_all = squeeze(std(avg_500_press,1)./sqrt(size(avg_500_press,1)));
avg_long_press = squeeze(func(press_long_movie,1));
sem_long_press = squeeze(std(press_long_movie,1)./sqrt(size(press_long_movie,1)));
avg_long_press_all = squeeze(func(avg_long_press,1));
sem_long_press_all = squeeze(std(avg_long_press,1)./sqrt(size(avg_long_press,1)));

figure;
tt =((-pre_press_frames:post_press_frames).*double(ifi))./1000;
errorbar(tt, avg_long_press_all, sem_long_press_all, '-k')
hold on;
errorbar(tt, avg_200_press_all, sem_200_press_all, '-b')
hold on;
errorbar(tt, avg_500_press_all, sem_500_press_all, '-g')
legend('hold>500ms', '200ms<hold<500ms', 'hold<200ms')
title(['Initiating press-by hold length: ' tc_type ' Short n = ' num2str(size(press_200_movie,1)) '; Mid n = '  num2str(size(press_500_movie,1)) '; Long n = '  num2str(size(press_long_movie,1))])
print([dest '_press_avgTCs_bylength_' tc_type '.eps'], '-depsc');
print([dest '_press_avgTCs_bylength_' tc_type '.pdf'], '-dpdf');

avg_all = [avg_200_press avg_500_press avg_long_press];
ymax = max(max(avg_all,[],2),[],1);
ymin = min(min(avg_all,[],2),[],1);
figure;
for ic = 1:nCells
    subplot(z,z,ic)
    errorbar(tt,avg_long_press(ic,:), sem_long_press(ic,:),'k')
    hold on
    errorbar(tt,avg_200_press(ic,:), sem_200_press(ic,:),'b')
    hold on
    errorbar(tt,avg_500_press(ic,:), sem_500_press(ic,:),'g')
    ylim([ymin*1.1 ymax*1.1])
    xlim([tt(1) tt(end)])
end
suptitle(['Initiating press-by hold length: ' tc_type ' Blue (Short) n = ' num2str(size(press_200_movie,1)) '; Green (Mid) n = '  num2str(size(press_500_movie,1)) '; Black (Long) n = '  num2str(size(press_long_movie,1))])
orient landscape
print([dest '_press_allTCs_bylength_' tc_type '.eps'], '-depsc');
print([dest '_press_allTCs_bylength_' tc_type '.pdf'], '-dpdf');

save([dest '_press_resp_by_hold_' tc_type '.mat'],'press_200_movie','press_500_movie','press_long_movie','press_movie','pre_press_frames', 'post_press_frames');

%break up presses longer than 500 ms by outcome
longHoldIx = zeros(1,length(b_data.input.trialOutcomeCell));
longHoldIx(ind_long) = 1;
successIx = strcmp(b_data.input.trialOutcomeCell,'success');
failureIx = strcmp(b_data.input.trialOutcomeCell,'failure');

longHold_success = find(and(longHoldIx,successIx));
longHold_failure = find(and(longHoldIx,failureIx));

use_press_success = pressTime(longHold_success);
use_press_failure = pressTime(longHold_failure);

press_success_movie = trigger_movie_by_event(tc_dfoverf, frame_info, ...
    use_press_success, pre_press_frames, post_press_frames);
press_failure_movie = trigger_movie_by_event(tc_dfoverf, frame_info, ...
    use_press_failure, pre_press_frames, post_press_frames);

avg_success_press = squeeze(func(press_success_movie,1));
sem_success_press = squeeze(std(press_success_movie,1)./sqrt(size(press_success_movie,1)));
avg_success_press_all = squeeze(func(avg_success_press,1));
sem_success_press_all = squeeze(std(avg_success_press,1)./sqrt(size(avg_success_press,1)));
avg_failure_press = squeeze(func(press_failure_movie,1));
sem_failure_press = squeeze(std(press_failure_movie,1)./sqrt(size(press_failure_movie,1)));
avg_failure_press_all = squeeze(func(avg_failure_press,1));
sem_failure_press_all = squeeze(std(avg_failure_press,1)./sqrt(size(avg_failure_press,1)));

figure;
tt =((-pre_press_frames:post_press_frames).*double(ifi))./1000;
errorbar(tt, avg_success_press_all, sem_success_press_all, '-k')
hold on;
errorbar(tt, avg_failure_press_all, sem_failure_press_all, '-r')
hold on;

title(['Initiating press: ' tc_type ' Success- black n = ' num2str(size(press_success_movie,1)) ' trials; Failure- red n = ' num2str(size(press_failure_movie,1)) ' trials'])
print([dest '_press_avgTCs_byoutcome_' tc_type '.eps'], '-depsc');
print([dest '_press_avgTCs_byoutcome_' tc_type '.pdf'], '-dpdf');

figure;
avg_all = [avg_success_press avg_failure_press];
ymax = max(max(avg_all,[],2),[],1);
ymin = min(min(avg_all,[],2),[],1);
for ic = 1:nCells
    subplot(z,z,ic)
    errorbar(tt,avg_success_press(ic,:), sem_success_press(ic,:),'k')
    hold on
    errorbar(tt,avg_failure_press(ic,:), sem_failure_press(ic,:),'r')
    ylim([ymin*1.1 ymax*1.1])
    xlim([tt(1) tt(end)])
end
suptitle(['Initiating press: ' tc_type ' Success- black n = ' num2str(size(press_success_movie,1)) ' trials; Failure- red n = ' num2str(size(press_failure_movie,1)) ' trials'])
orient landscape
print([dest '_press_allTCs_byoutcome_' tc_type '.eps'], '-depsc');
print([dest '_press_allTCs_byoutcome_' tc_type '.pdf'], '-dpdf');

save([dest '_press_resp_by_outcome.mat'],'press_success_movie','press_failure_movie','pre_press_frames', 'post_press_frames');

%compare single cell press and release responses
avg_all = [avg_long_press avg_success];
ymax = max(max(avg_all,[],2),[],1);
ymin = min(min(avg_all,[],2),[],1);
figure;
for ic = 1:nCells
    subplot(z,z,ic)
    errorbar(tt,avg_long_press(ic,:), sem_long_press(ic,:),'c')
    hold on
    errorbar(tt(6:end),avg_success(ic,:), sem_success(ic,:),'k')
    ylim([ymin*1.1 ymax*1.1])
    xlim([tt(1) tt(end)])
end
suptitle(['Press ' tc_type ' (cyan: n = ' num2str(size(press_long_movie,1)) ' trials); Release (black n = '  num2str(size(success_movie,1)) ' trials)'])
orient landscape
print([dest '_press_release_allTCs_' tc_type '.eps'], '-depsc');
print([dest '_press_release_allTCs_' tc_type '.pdf'], '-dpdf');
