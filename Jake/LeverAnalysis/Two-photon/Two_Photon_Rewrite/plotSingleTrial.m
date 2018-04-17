%plot single trial for omission
fig = figure; ic = 6; shift = 0;
ttp = ((-pre_release_frames:pre_release_frames).*double(ifi));
subplot(1,2,1); axis square; xlabel('time (ms)'); ylabel('df/f');
title('omission');
ni = size(omitReward_movie,1);
subplot(1,2,1); axis square
hold on
plot(ttp, squeeze(omitReward_movie(3, 17,1:pre_release_frames*2+1))+shift);
shift = shift + 0.2;
plot(ttp, squeeze(omitReward_movie(35, 18,1:pre_release_frames*2+1))+shift);
shift = shift + 0.2;
for ievent = [12,16:18,24,28,33,38]%1:size(omitReward_movie,1)
%     if mod(ievent,5) == 0
%         figure;
%         shift = 0;
%     end
        
    subplot(1,2,1); axis square
    hold on
    plot(ttp, squeeze(omitReward_movie(ievent, ic,1:pre_release_frames*2+1))+shift);
    set(gca,'XTick',[-500 0 200 400])
%     text(0, shift, num2str(ievent));
    shift = shift+0.2;
end
ylim([-1 round(shift)+1])
subplot(1,2,2); axis square; xlabel('time (ms)'); ylabel('df/f'); shift=0;
title('reward');ylim([-1 5])

for ievent = [87,61,97,103,45,109,105,130,45,52]%1:size(success_movie,1)
%     if mod(ievent,5) == 0
%         figure;
%         shift = 0;
%     end
        
    subplot(1,2,2); axis square
    hold on
    plot(ttp, squeeze(success_movie(ievent, ic,1:pre_release_frames*2+1))+shift);
    set(gca,'XTick',[-500 0 200 400])
%     text(0, shift, num2str(ievent));
    shift = shift+0.2;
end
ylim([-1 round(shift)+1])

supertitle('Single Trial for all ROIS for Img 063');
saveas(fig,['Z:\public\Motor Timing Paper\Ziye\2P_PCA_ICA\30_sessions\single_session_example_omission\', '180109_img063_singleTrial_omission_success.fig']);
print(['Z:\public\Motor Timing Paper\Ziye\2P_PCA_ICA\30_sessions\single_session_example_omission\', '180109_img063_singleTrial_omission_success.pdf'], '-dpdf');
print(['Z:\public\Motor Timing Paper\Ziye\2P_PCA_ICA\30_sessions\single_session_example_omission\', '180109_img063_singleTrial_omission_success.eps'], '-deps');

%plot single trial for omission
fig = figure; ic = 19; shift = 0;
ttp = ((-pre_release_frames:pre_release_frames).*double(ifi));
subplot(1,2,1); axis square; xlabel('time (ms)'); ylabel('df/f');
title('omission');
ni = size(omitReward_movie,1);
subplot(1,2,1); axis square
hold on
plot(ttp, squeeze(omitReward_movie(31, 1,1:pre_release_frames*2+1))+shift);
for ievent = 1:size(omitReward_movie,1)
    if mod(ievent,5) == 0
        figure;
        shift = 0;
    end
        
    subplot(1,2,1); axis square
    hold on
    plot(ttp, squeeze(omitReward_movie(ievent, ic,1:pre_release_frames*2+1))+shift);
    set(gca,'XTick',[-500 0 200 400])
    text(0, shift, num2str(ievent));
    shift = shift+0.2;
end
ylim([-1 round(shift)+1])


tt =((-pre_release_frames:post_release_frames).*double(ifi));
icell = 35;
for iT = 46%1:size(omitReward_movie,1)
    figure
    plot(tt,squeeze(omitReward_movie(iT,icell,:)))
end
xlabel('Time (ms)'); ylabel('dF/F');
xlim([-500 500]); title(['single trial time course for late early Cell', num2str(icell)]);

% new graph
fig = figure; shift  = 0;
ttp = ((-pre_release_frames:pre_release_frames).*double(ifi));
subplot(1,2,1); axis square; xlabel('time (ms)'); ylabel('df/f');
title('omission');
ni = size(omitReward_movie,1);


subplot(1,2,2); axis square; xlabel('time (ms)'); ylabel('df/f'); shift=0;
title('reward');ylim([-1 5])

% plot single trial for late early and early early for 180102_img 060
fig = figure; ievent = 1; shift = 0;
ttp = ((-pre_release_frames:pre_release_frames).*double(ifi));
subplot(1,2,1); axis square; xlabel('time (ms)'); ylabel('df/f');
title('early < 1s')
for ic = 1:size(early_fail_movie,2)
    hold on
    plot(ttp, squeeze(early_fail_movie(ievent, ic,1:pre_release_frames*2+1))*0.5+shift);
    shift = shift+0.15;
end
ylim([-0.1 round(shift)+1])
subplot(1,2,2); axis square; xlabel('time (ms)'); ylabel('df/f'); shift=0;
title('early > 3.5s')
for ic = 1:size(late_fail_movie,2)
    hold on
    plot(ttp, squeeze(late_fail_movie(ievent, ic,1:pre_release_frames*2+1))+shift);
    shift = shift+0.15;
end
ylim([-0.1 round(shift)+1])
supertitle('Single Trial for all ROIS for Img 060');
saveas(fig,['Z:\public\Motor Timing Paper\Ziye\2P_PCA_ICA\30_sessions\single_session_example_img060\', '180102_img060_singleTrial_earlyLate_early2.fig']);
print(['Z:\public\Motor Timing Paper\Ziye\2P_PCA_ICA\30_sessions\single_session_example_img060\', '180102_img060_singleTrial_earlyLate_early2.pdf'], '-dpdf');
print(['Z:\public\Motor Timing Paper\Ziye\2P_PCA_ICA\30_sessions\single_session_example_img060\', '180102_img060_singleTrial_earlyLate_early2.eps'], '-deps');

tt =((-pre_release_frames:post_release_frames).*double(ifi));
icell = 27;
for iT = 1:size(late_fail_movie,1)
    figure
    plot(tt,squeeze(late_fail_movie(iT,icell,:)))
end
xlabel('Time (ms)'); ylabel('dF/F');
xlim([-500 500]); title(['single trial time course for late early Cell', num2str(icell)]);

% plot single trial for late early and early early for 160202_img 32
early_fail_movie_temp = early_fail_movie;
early_fail_movie_temp(:,:,1:end-2) = early_fail_movie(:,:,3:end);
early_fail_movie_temp(:,:,end-2:end) = early_fail_movie(:,:,1:3);

fig = figure; ic = 13; shift = 0;
ttp = ((-pre_release_frames:pre_release_frames).*double(ifi));
subplot(1,2,1); axis square; xlabel('time (ms)'); ylabel('df/f');
title('early < 1s')
for ievent = [11,9,8,7]
    hold on
    plot(ttp, squeeze(early_fail_movie_temp(ievent, ic,1:pre_release_frames*2+1))+shift);
    shift = shift+0.5;
end
ylim([-1 round(shift)+1])

subplot(1,2,2); axis square; xlabel('time (ms)'); ylabel('df/f'); shift=0;
title('early > 3.5s')
ic = 13;
hold on
plot(ttp, squeeze(late_fail_movie(3, ic,1:pre_release_frames*2+1))+shift);
shift = shift+0.5;

ic = 10;
plot(ttp, squeeze(late_fail_movie(1, ic,1:pre_release_frames*2+1))+shift);
shift = shift+0.5;

ic = 15;
plot(ttp, squeeze(late_fail_movie(4, ic,1:pre_release_frames*2+1))+shift);
shift = shift+0.5;

ic = 7;
for ievent = [2]
    hold on
    plot(ttp, squeeze(late_fail_movie(ievent, ic,1:pre_release_frames*2+1))+shift);
    shift = shift+0.5;
end
ylim([-1 round(shift)+1])

supertitle('Single Trial for all ROIS for Img 32');
saveas(fig,['Z:\public\Motor Timing Paper\Ziye\2P_PCA_ICA\30_sessions\single_session_example_img32\', '160202_img32_singleTrial_earlyLate_early_cell.fig']);
print(['Z:\public\Motor Timing Paper\Ziye\2P_PCA_ICA\30_sessions\single_session_example_img32\', '160202_img32_singleTrial_earlyLate_early_cell.pdf'], '-dpdf');
print(['Z:\public\Motor Timing Paper\Ziye\2P_PCA_ICA\30_sessions\single_session_example_img32\', '160202_img32_singleTrial_earlyLate_early_cell.eps'], '-deps');

tt =((-pre_release_frames:post_release_frames).*double(ifi));
icell = 1;
for icell = 1:size(late_fail_movie,2)
for iT = 1:size(late_fail_movie,1)
    figure
    plot(tt,squeeze(late_fail_movie(iT,icell,:)))
    title(['single trial time course for late early Cell', num2str(icell)]);
end
end
xlabel('Time (ms)'); ylabel('dF/F');
xlim([-500 500]); title(['single trial time course for late early Cell', num2str(icell)]);


% plot single trial for late early and early early for img 59
fig = figure; ievent = 4; shift = 0;
ttp = ((-pre_release_frames:pre_release_frames).*double(ifi));
subplot(1,2,1); axis square; xlabel('time (ms)'); ylabel('df/f');
title('early < 1s')
for ic = 1:size(early_fail_movie,2)
    hold on
    plot(ttp, squeeze(early_fail_movie(ievent, ic,1:pre_release_frames*2+1))+shift);
    shift = shift+0.1;
end

subplot(1,2,2); axis square; xlabel('time (ms)'); ylabel('df/f'); shift=0;
title('early > 3.5s')
for ic = 1:size(late_fail_movie,2)
    hold on
    plot(ttp, squeeze(late_fail_movie(ievent, ic,1:pre_release_frames*2+1))+shift);
    shift = shift+0.1;
end
supertitle('Single Trial for all ROIS for Img 59');
saveas(fig,['Z:\public\Motor Timing Paper\Ziye\2P_PCA_ICA\', '161103_img59_singleTrial_earlyLate_early.fig']);
print(['Z:\public\Motor Timing Paper\Ziye\2P_PCA_ICA\', '161103_img59_singleTrial_earlyLate_early.pdf'], '-dpdf');
print(['Z:\public\Motor Timing Paper\Ziye\2P_PCA_ICA\', '161103_img59_singleTrial_earlyLate_early.eps'], '-deps');

load(['Z:\home\ziye\2P_Analysis\2P_Analysis\180109_000_img063\', '_release_movies.mat'])
%180104_000_img060 180109_000_img063
iC = 10;
tt =((-pre_release_frames:post_release_frames).*double(ifi));
for iT = 1:size(omitReward_movie,1)
    figure
    plot(tt, squeeze(omitReward_movie(iT, iC, :)))
end
xlabel('Time (ms)'); ylabel('dF/F');
title(['single trial time course for omission reward Cell', num2str(6), ' event', num2str(33)]);

load(['Z:\home\ziye\2P_Analysis\2P_Analysis\161103_000_img59\', '_release_movies.mat'])

tt =((-pre_release_frames:post_release_frames).*double(ifi));
for iT = 1:size(late_fail_movie,1)
    figure
    plot(tt,squeeze(late_fail_movie(iT,16,:)))
end
xlabel('Time (ms)'); ylabel('dF/F');
xlim([-500 500]); title(['single trial time course for late early Cell', num2str(16)]);

figure;

b_data = input;
load(['Z:\home\ziye\2P_Analysis\2P_Analysis\161103_000_img59\' '_dFOverF_TC.mat'])
tN = 11; % 21 36 39 64for fail and 11 for success
tC = 3;
tc_trial = tc_dfoverf(:,b_data.cLeverDown{tN} - 30 :b_data.cTrialEnd{tN} + 20); % +40 for plot entire early

fig = figure;subplot(1,2,2);plot((0:b_data.cTrialEnd{tN} - b_data.cLeverDown{tN} +50 )*33, tc_trial(tC, :)); % +40 for plot entire early

% subplot(1,2,2);plot((0:b_data.cTrialEnd{tN} - b_data.cLeverDown{tN})*33, tc_trial(tC, :));
if strcmp(b_data.trialOutcomeCell{tN}, 'failure')
    ylim([-0.25 0.6]); xlim([0 6000]);
    hold on; vline([31, double(b_data.cLeverUp{tN} - b_data.cLeverDown{tN}) + 30]*33, {'--k', '--k', '--k'},...
    {'press', 'release'});
    title('Early trial timecourse')
    xlabel('time(ms)'); ylabel('df/f');
    
else
    ylim([-0.25 0.6]);
    hold on; vline([31, double(b_data.cTargetOn{tN} - b_data.cLeverDown{tN})+30, double(b_data.cLeverUp{tN} - b_data.cLeverDown{tN})+30]*33, {'--k', '--k', '--k'},...
    {'press', 'cue', 'release'});
    title('Correct trial timecourse');
    xlabel('time(ms)'); ylabel('df/f');
end

saveas(fig,['Z:\public\Motor Timing Paper\Ziye\2P_PCA_ICA\', '161103_img59_singleTrial_tc.fig']);
print(['Z:\public\Motor Timing Paper\Ziye\2P_PCA_ICA\', '161103_img59_singleTrial_tc.pdf'], '-dpdf');

% plot single trial for correct and early for img 59
fig = figure; ievent = 5; shift = 0;
ttp = ((-pre_release_frames:pre_release_frames).*double(ifi));
subplot(1,2,1); axis square; xlabel('time (ms)'); ylabel('df/f');
title('Correct')
for ic = 1:size(success_movie,2)
    hold on
    plot(ttp, squeeze(success_movie(ievent, ic,1:pre_release_frames*2+1))+shift);
    shift = shift+0;
end

subplot(1,2,2); axis square; xlabel('time (ms)'); ylabel('df/f');
title('Early')
for ic = 1:size(fail_movie,2)
    hold on
    plot(ttp, squeeze(fail_movie(ievent, ic,1:pre_release_frames*2+1))+shift);
    shift = shift+0;
end
supertitle('Single Trial for all ROIS for Img 59');
saveas(fig,['Z:\public\Motor Timing Paper\Ziye\2P_PCA_ICA\', '161103_img59_singleTrial_AllTC.fig']);
print(['Z:\public\Motor Timing Paper\Ziye\2P_PCA_ICA\', '161103_img59_singleTrial_AllTC.pdf'], '-dpdf');

% plot press time

fig = figure; ievent = 3; shift = 0;
ttp = ((-pre_press_frames:post_press_frames).*double(ifi));
subplot(1,2,1); axis square; xlabel('time (ms)'); ylabel('df/f');
title('Correct')
for ic = 1:size(press_success_movie,2)
    hold on
    plot(ttp, squeeze(press_success_movie(ievent, ic,:))+shift);
    shift = shift+0.1;
end
vline(0,'--k'); vline((trial_outcome.success_ptimeEndF(ievent))*double(ifi),'--k');
xlim([-500 2500])

set(gca, 'XTick', [0, (trial_outcome.success_ptimeEndF(ievent))*double(ifi)], 'XTickLabel', {'press', 'release'});
subplot(1,2,2); axis square; xlabel('time (ms)'); ylabel('df/f');
title('Early')
ievent = 5
for ic = 1:size(press_fail_movie,2)
    hold on
    plot(ttp, squeeze(press_fail_movie(ievent, ic,:))+shift);
    shift = shift+0.1;
end
vline(0,'--k'); vline((trial_outcome.early_ptimeEndF(ievent))*double(ifi),'--k');
set(gca, 'XTick', [0, (trial_outcome.early_ptimeEndF(ievent))*double(ifi)], 'XTickLabel', {'press', 'release'});
xlim([-500 1000])
supertitle('Single Trial for all ROIS for Img 59');
saveas(fig,['Z:\public\Motor Timing Paper\Ziye\2P_PCA_ICA\', '161103_img59_singleTrial_AllTC.fig']);
print(['Z:\public\Motor Timing Paper\Ziye\2P_PCA_ICA\', '161103_img59_singleTrial_AllTC.pdf'], '-dpdf');



fig = figure;
c1 = 83; 
% for c1 = 1:100
%     figure
errorbar(tt, licks_TC_notlick(c1,:)+0.04, sem_licks_TC_notlick(c1,:), 'r');
% end
hold on; 
c2 = 117;
% for c2 = 100:120
% figure
errorbar(tt, licks_TC_lick(c2,:), sem_licks_TC_lick(c2,:),'b');
% end
xlabel('time(ms)'); ylabel('df/f');
title('single cell lick triggered avg: no response- red, response- blue');
saveas(fig, [out_base 'single_cell_LickNoLick.fig']);
print([out_base 'single_cell_LickNoLick.eps'], '-depsc');
print([out_base 'single_cell_LickNoLick.pdf'], '-dpdf');

fig = figure; ievent = 20; shift = 0;
ttp = ((-pre_release_frames:pre_release_frames).*double(ifi));
subplot(1,2,1); axis square; xlabel('time (ms)'); ylabel('df/f');
title('lick resp')
single_lick_movie = licks_movie{7}; %23
for ic = 1:size(single_lick_movie,2)
    hold on
    plot(ttp, squeeze(single_lick_movie(ievent, ic,1:pre_release_frames*2+1))+shift);
    shift = shift+0;
end
subplot(1,2,2); axis square; xlabel('time (ms)'); ylabel('df/f');
title('not lick resp')
single_lick_movie = licks_movie_notlick{11}; %3
for ic = 1:size(single_lick_movie,2)
    hold on
    plot(tt, squeeze(single_lick_movie(ievent, ic,1:pre_release_frames*2+1)));
%     shift = shift+0;
end
saveas(fig,['Z:\public\Motor Timing Paper\Ziye\2P_PCA_ICA\', 'LickNotLick_singleTrial_AllTC.fig']);
print(['Z:\public\Motor Timing Paper\Ziye\2P_PCA_ICA\', 'LickNotLick_singleTrial_AllTC.eps'], '-depsc');
print(['Z:\public\Motor Timing Paper\Ziye\2P_PCA_ICA\', 'LickNotLick_singleTrial_AllTC.pdf'], '-dpdf');

% run getTC_events
%id = 27 180110_img064 notlick resp
% load([tc_dir 'parse_behavior.mat']);
% load([tc_dir, '_dFOverF_TC.mat']);
nf = size(single_lick_movie,3); ievent = 15;
for ievent = [5 20 24] %1:size(single_lick_movie,1)
fig = figure;plot((0:nf-1)*33, squeeze(single_lick_movie(ievent,23,:))); 
hold on;b = bar((0:nf-1)*33, licks_lick(ievent, :)/30+0.15, 0.5, 'BaseValue', 0.15);
bl = b.BaseLine;
bl.Color = 'white';
xlabel('time (ms)'); ylabel('dff');
title('not lick resp cell');
saveas(fig,['Z:\public\Motor Timing Paper\Ziye\2P_PCA_ICA\', 'NotLickCell_raw', num2str(ievent),'.fig']);
print(['Z:\public\Motor Timing Paper\Ziye\2P_PCA_ICA\', 'NotLickCell_raw', num2str(ievent),'.eps'], '-depsc');
print(['Z:\public\Motor Timing Paper\Ziye\2P_PCA_ICA\', 'NotLickCell_raw', num2str(ievent),'.pdf'], '-dpdf');
end


%id = 32 180104_img060
nf = size(single_lick_movie,3); ievent = 5;
for ievent = [18 22]%1:size(single_lick_movie,1)
fig = figure;plot((0:nf-1)*33, squeeze(single_lick_movie(ievent,7,:))); 
hold on;b = bar((0:nf-1)*33, licks_lick(ievent, :)/30+0.3, 0.5, 'BaseValue', 0.3);
bl = b.BaseLine;
bl.Color = 'white';
xlabel('time (ms)'); ylabel('dff');
title('lick resp cell');
saveas(fig,['Z:\public\Motor Timing Paper\Ziye\2P_PCA_ICA\', 'LickCell_raw', num2str(ievent),'.fig']);
print(['Z:\public\Motor Timing Paper\Ziye\2P_PCA_ICA\', 'LickCell_raw', num2str(ievent),'.eps'], '-depsc');
print(['Z:\public\Motor Timing Paper\Ziye\2P_PCA_ICA\', 'LickCell_raw', num2str(ievent),'.pdf'], '-dpdf');
end

