%compare the magnitude of lick triggered responses to the magnitude of the
%difference between corrects and earlies for that trial

clear;
WF_plotting_lists_of_days;
lick_trig_dir = ['Z:\Analysis\WF Lever Analysis\licking_investigation\lick_triggered_averages\lick_trig_mags'];
TC_dir = ['Z:\Analysis\WF Lever Analysis\licking_investigation\corr_early_diff_by_ROI\'];
colors = {'k','k', 'b','b', 'g','g', 'c','c', 'm','m', 'r','r', 'y','y', 'k', 'b', 'r'}; 
peak_window = [41:80];
peak_diff = {};
peak_diff{1} = [];

%load lick triggered magnitudes
load([lick_trig_dir]);

%load tc data and determine differences
for session_num = 5:length(days)
    load([TC_dir, days{session_num}]);
    
    %find the peak values and get the difference
    peak_corr = max(success_roi_interp_avg2(peak_window, :));
    peak_early = max(fail_roi_interp_avg2(peak_window, :));
    peak_diff{session_num} = peak_corr-peak_early;
end

grand_mean_diff = [];
grand_mean_trig = [];
figure;
subplot(1,2,1);
for session_num= 5:length(days)
    for ROI_num = 1:length(peak_diff{session_num})
        plot([1,2],  [peak_diff{session_num}(ROI_num), lick_trig_mags{session_num}(ROI_num)],   strcat('-', colors{session_num},'o')); hold on;
        grand_mean_diff = [grand_mean_diff, peak_diff{session_num}(ROI_num)];
        grand_mean_trig = [grand_mean_trig, lick_trig_mags{session_num}(ROI_num)];
    end
end
plot([1,2], [mean(grand_mean_diff), mean(grand_mean_trig)],'-o', 'MarkerFaceColor', 'k');
xlim([0 3]); ylim([-0.05 0.2]);
xlabel('diff between correct and early peaks      lick triggered peak');
ylabel('df/f');
title('all ROIs');

grand_mean_diff = [];
grand_mean_trig = [];
subplot(1,2,2);
for session_num= 5:length(days)
        plot([1,2],  [mean(peak_diff{session_num}), mean(lick_trig_mags{session_num})],   strcat('-', colors{session_num},'o')); hold on;
        grand_mean_diff = [grand_mean_diff, mean(peak_diff{session_num})];
        grand_mean_trig = [grand_mean_trig, mean(lick_trig_mags{session_num})];
end
plot([1,2], [mean(grand_mean_diff), mean(grand_mean_trig)],'-o', 'MarkerFaceColor', 'k');
xlim([0 3]); ylim([-0.05 0.2]);
xlabel('diff between correct and early peaks      lick triggered peak');
ylabel('df/f');
title('ROIs averaged: one point per session');

suptitle('single licks vs difference between peak magnitude for correct and early trials');
savefig(['Z:\Analysis\WF Lever Analysis\licking_investigation\lick_trig_vs_outcome_diff']);

