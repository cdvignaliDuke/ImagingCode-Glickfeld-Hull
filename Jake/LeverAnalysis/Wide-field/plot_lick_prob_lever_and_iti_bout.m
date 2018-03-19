%plot the lick traces for corr, early, and ITI licks for each animal as a
%probability of a lick occuring on a given frame. Also plot the mean across
%animals on the same plot. 
clear; 
lever_release = 6; 
pre_frames=5;
post_frames=10;
WF_plotting_lists_of_days;
lick_TC_dir = 'Z:\Analysis\WF Lever Analysis\BxAndAnalysisOutputs\BxOutputs\'; 
ITI_lick_trace_source = ['Z:\Analysis\WF Lever Analysis\licking_investigation\lick_triggered_averages\individual lick bout TCs\'];
days = {'151021_img29', '151022_img29', [], [], '151211_img32', '151212_img32', '160129_img35', '160131_img35', '160129_img36', '160131_img36', '160314_img38', '160315_img38', '160319_img41', '160320_img41', '160606_img46', '160722_img53', '160904_img55'}; %'150718_img27', '150719_img27', '150716_img28', '150717_img28', 
colors_roi = {'b', 'r', 'g', 'k', 'c', 'm'};
colors = {'k','k', 'b','b', 'g','g', 'c','c', 'm','m', 'r','r', 'y','y', 'k', 'b', 'r'}; 
plot_symbols = {'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 's', 's', 's'};
ITI_lick_trace_all = [];
lick_trace_corr_all = [];
lick_trace_fail_all = [];

for session_num = 1:length(days)
    if isempty(days{session_num})
        continue
    end
    load([ITI_lick_trace_source, days{session_num}, '_lick_trace']);
    load([lick_TC_dir, days{session_num}, '_bx_outputs'], 'licking_data');
    lick_trace_corr = boolean(licking_data.lick_trace_succ);
    lick_trace_fail = boolean(licking_data.lick_trace_fail);
    ITI_lick_trace = boolean(ITI_lick_trace);
    corr_trial_num = size(lick_trace_corr,1);
    fail_trial_num = size(lick_trace_fail,1);
    ITI_bout_num = size(ITI_lick_trace,1);
    if corr_trial_num > 9
        if fail_trial_num > 9
            if ITI_bout_num > 9
                ITI_lick_trace_all = [ITI_lick_trace_all; mean(ITI_lick_trace)];
                lick_trace_corr_all = [lick_trace_corr_all; mean(lick_trace_corr)];
                lick_trace_fail_all = [lick_trace_fail_all; mean(lick_trace_fail)];
            end
        end
    end

end

figure; 
subplot(1,3,1);
for session_num = 1:size(ITI_lick_trace_all,1)
    plot([-pre_frames:post_frames], lick_trace_corr_all(session_num,:), 'k'); hold on;
end
errorbar([-pre_frames:post_frames],mean(lick_trace_corr_all), std(lick_trace_corr_all)/sqrt(size(lick_trace_corr_all,1)), 'm');
title('correct trials lick probability per frame');
xlabel('frames relative to lever release');
ylabel('probability of a lick on a given frame');
xlim([-5 10]); ylim([0 1.1]);

subplot(1,3,2);
for session_num = 1:size(ITI_lick_trace_all,1)
    plot([-pre_frames:post_frames], lick_trace_fail_all(session_num,:), 'k'); hold on;
end
errorbar([-pre_frames:post_frames],mean(lick_trace_fail_all), std(lick_trace_fail_all)/sqrt(size(lick_trace_fail_all,1)), 'm');
title('early trials lick probability per frame');
xlabel('frames relative to lever release');
ylabel('probability of a lick on a given frame');
xlim([-5 10]); ylim([0 1.1]);

subplot(1,3,3);
for session_num = 1:size(ITI_lick_trace_all,1)
    plot([-pre_frames:post_frames], ITI_lick_trace_all(session_num,:), 'k'); hold on;
end
errorbar([-pre_frames:post_frames],mean(ITI_lick_trace_all), std(ITI_lick_trace_all)/sqrt(size(ITI_lick_trace_all,1)), 'm');
title('ITI licks lick probability per frame');
xlabel('frames relative to first lick');
ylabel('probability of a lick on a given frame');
xlim([-5 10]); ylim([0 1.05]);





