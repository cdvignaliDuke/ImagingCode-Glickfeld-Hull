% SCRIPT RETIRED ON 10/22/18. Replaced by resp_cell_criteria



% finds and maps responsive cells for each condition
% 1. calculate average timecourses for NR Ca events
% 2. calculate variability by trial over base and resp windows, 
%    ttest for significant responses
% 3. calculate 10% rise time of F onset
% 4. define cells by response to event
% 5. plot everything

% load([dest 'ROI_TCs.mat']);
load([dest 'parse_behavior.mat']);
load([dest '_cue_movies.mat']);
load([dest '_cue_movies_lick.mat'])
nCells = size(NR_movie,2);
n = ceil(sqrt(nCells));
if nCells <((n.^2)-n) %function
    n2= n-1;
else
    n2 = n;
end

%% 1. calculate average timecourses for NR/release events
trial_cond = {'all'};
%find avg and sem across trials then across trials and cells
[avg_NR, sem_NR, avg_NR_all, sem_NR_all] = get_movie_mean_sem(NR_movie, trial_cond, 'all');
[avg_OR, sem_OR, avg_OR_all, sem_OR_all] = get_movie_mean_sem(OR_movie, trial_cond, 'all');
[avg_UR, sem_UR, avg_UR_all, sem_UR_all] = get_movie_mean_sem(UR_movie, trial_cond, 'all');

%average and sem across Trials for no lick trials
[avg_NR_nolick, sem_NR_nolick, avg_NR_all_nolick, sem_NR_all_nolick] = get_movie_mean_sem(NR_movie_nolick, trial_cond, 'nolick');
[avg_OR_nolick, sem_OR_nolick, avg_OR_all_nolick, sem_OR_all_nolick] = get_movie_mean_sem(OR_movie_nolick, trial_cond, 'nolick');
[avg_UR_nolick, sem_UR_nolick, avg_UR_all_nolick, sem_UR_all_nolick] = get_movie_mean_sem(UR_movie_nolick, trial_cond', 'nolick');

if strcmp(trial_cond, 'half')
    save([dest '_cell_TCs_half.mat'], 'avg_NR', 'sem_NR', 'avg_OR', 'sem_OR', 'avg_UR', 'sem_UR', ...
        'avg_NR_nolick', 'sem_NR_nolick', 'avg_OR_nolick', 'sem_OR_nolick', 'avg_UR_nolick', 'sem_UR_nolick');
elseif strcmp(trial_cond, 'all')
    save([dest '_cell_TCs.mat'], 'avg_NR', 'sem_NR', 'avg_OR', 'sem_OR', 'avg_UR', 'sem_UR', ...
        'avg_NR_nolick', 'sem_NR_nolick', 'avg_OR_nolick', 'sem_OR_nolick', 'avg_UR_nolick', 'sem_UR_nolick');
end

% Plot timecourse for each cell 
tt =((-pre_cue_frames:post_cue_frames).*double(ifi))./1000;
%overlay average success/failure for each ROI
fig=figure;
avg_all = [avg_NR avg_OR];
ymax = max(max(avg_all,[],2),[],1);
ymin = min(min(avg_all,[],2),[],1);
for ic = 1:nCells
    subplot_tight(n,n2,ic);
    errorbar(tt,avg_NR(ic,:), sem_NR(ic,:),'k');
    hold on;
    errorbar(tt,avg_OR(ic,:), sem_OR(ic,:),'r');
    ylim([ymin*1.1 ymax*1.1]);
    xlim([tt(1) tt(end)]);
end
supertitle([session_date ' ' mouse_ID ' timecourse of normal reward- black (n = ' num2str(size(NR_movie,1)) ' trials); Omission- red (n = ' num2str(size(OR_movie,1)) ' trials). All cells (n=', num2str(nCells), ')']);
orient landscape
saveas(fig, [dest '_NR_OR_allTCs.fig']);
print([dest '_NR_OR_allTCs.eps'], '-depsc');
print([dest 'NR_OR_allTCs.pdf'], '-dpdf');

%% 2. calculate response amplitude and variablity by trial over base and resp windows

%POSITIVE CUE responsive neurons 
ops.base_cue_buffer = 500;
ops.resp_cue_buffer = 400;
ops.effect_sign = 'pos';
ops.event_type = 'cue';
ops.ifi = ifi;
[NR_Cue_h_pos, NR_Cue_p_pos, NR_Cue_resp_cells_pos, NR_Cue_resp_avg_pos, NR_Cue_resp_sem_pos, NR_Cue_base_pos, NR_Cue_resp_pos] = findRespCell(NR_movie_nolick, pre_cue_frames, ifi, base_cue_buffer, resp_cue_buffer, effect_sign);
[OR_Cue_h_pos, OR_Cue_p_pos, OR_Cue_resp_cells_pos, OR_Cue_resp_avg_pos, OR_Cue_resp_sem_pos, OR_Cue_base_pos, OR_Cue_resp_pos] = findRespCell(OR_movie_nolick, pre_cue_frames, ifi, base_cue_buffer, resp_cue_buffer, effect_sign);
[UR_Cue_h_pos, UR_Cue_p_pos, UR_Cue_resp_cells_pos, UR_Cue_resp_avg_pos, UR_Cue_resp_sem_pos, UR_Cue_base_pos, UR_Cue_resp_pos] = findRespCell(UR_movie_nolick, pre_cue_frames, ifi, base_cue_buffer, resp_cue_buffer, effect_sign);
%omission responsive window from cue to reward 0-600 
[OR_Cue_h2_pos, OR_Cue_p2_pos, OR_Cue_resp_cells2_pos, OR_Cue_resp_avg2_pos, OR_Cue_resp_sem2_pos, OR_Cue_base2_pos, OR_Cue_resp2_pos] = findRespCell(OR_movie_nolick, pre_cue_frames, ifi, base_cue_buffer, resp_cue_buffer+200, effect_sign);

%NEGATIVE CUE responsive neurons
effect_sign = 'neg';
[NR_Cue_h_neg, NR_Cue_p_neg, NR_Cue_resp_cells_neg, NR_Cue_resp_avg_neg, NR_Cue_resp_sem_neg, NR_Cue_base_neg, NR_Cue_resp_neg] = findRespCell(NR_movie_nolick, pre_cue_frames, ifi, base_cue_buffer, resp_cue_buffer, effect_sign);
[OR_Cue_h_neg, OR_Cue_p_neg, OR_Cue_resp_cells_neg, OR_Cue_resp_avg_neg, OR_Cue_resp_sem_neg, OR_Cue_base_neg, OR_Cue_resp_neg] = findRespCell(OR_movie_nolick, pre_cue_frames, ifi, base_cue_buffer, resp_cue_buffer, effect_sign);
[UR_Cue_h_neg, UR_Cue_p_neg, UR_Cue_resp_cells_neg, UR_Cue_resp_avg_neg, UR_Cue_resp_sem_neg, UR_Cue_base_neg, UR_Cue_resp_neg] = findRespCell(UR_movie_nolick, pre_cue_frames, ifi, base_cue_buffer, resp_cue_buffer, effect_sign);
%omission responsive window from cue to reward 0-600 
[OR_Cue_h2_neg, OR_Cue_p2_neg, OR_Cue_resp_cells2_neg, OR_Cue_resp_avg2_neg, OR_Cue_resp_sem2_neg, OR_Cue_base2_neg, OR_Cue_resp2_neg] = findRespCell(OR_movie_nolick, pre_cue_frames, ifi, base_cue_buffer, resp_cue_buffer+200, effect_sign);

%POSITIVE REWARD responsive neurons
base_reward_buffer = 500;
resp_reward_buffer = 1000;
pre_rew_frames = pre_cue_frames + round((trial_outcome.normalReward(1) - trial_outcome.normalRewardCue(1))/ifi);
effect_sign = 'pos';
[NR_Rew_h_pos, NR_Rew_p_pos, NR_Rew_resp_cells_pos, NR_Rew_resp_avg_pos, NR_Rew_resp_sem_pos, NR_Rew_base_pos, NR_Rew_resp_pos] = findRespCell(NR_movie_nolick, pre_rew_frames, ifi, base_reward_buffer, resp_reward_buffer, effect_sign);
[OR_Rew_h_pos, OR_Rew_p_pos, OR_Rew_resp_cells_pos, OR_Rew_resp_avg_pos, OR_Rew_resp_sem_pos, OR_Rew_base_pos, OR_Rew_resp_pos] = findRespCell(OR_movie_nolick, pre_rew_frames, ifi, base_reward_buffer, resp_reward_buffer, effect_sign);
[UR_Rew_h_pos, UR_Rew_p_pos, UR_Rew_resp_cells_pos, UR_Rew_resp_avg_pos, UR_Rew_resp_sem_pos, UR_Rew_base_pos, UR_Rew_resp_pos] = findRespCell(UR_movie_nolick, pre_rew_frames, ifi, base_reward_buffer, resp_reward_buffer, effect_sign);

%NEGATIVE REWARD responsive neurons
effect_sign = 'neg';
[NR_Rew_h_neg, NR_Rew_p_neg, NR_Rew_resp_cells_neg, NR_Rew_resp_avg_neg, NR_Rew_resp_sem_neg, NR_Rew_base_neg, NR_Rew_resp_neg] = findRespCell(NR_movie_nolick, pre_rew_frames, ifi, base_reward_buffer, resp_reward_buffer, effect_sign);
[OR_Rew_h_neg, OR_Rew_p_neg, OR_Rew_resp_cells_neg, OR_Rew_resp_avg_neg, OR_Rew_resp_sem_neg, OR_Rew_base_neg, OR_Rew_resp_neg] = findRespCell(OR_movie_nolick, pre_rew_frames, ifi, base_reward_buffer, resp_reward_buffer, effect_sign);
[UR_Rew_h_neg, UR_Rew_p_neg, UR_Rew_resp_cells_neg, UR_Rew_resp_avg_neg, UR_Rew_resp_sem_neg, UR_Rew_base_neg, UR_Rew_resp_neg] = findRespCell(UR_movie_nolick, pre_rew_frames, ifi, base_reward_buffer, resp_reward_buffer, effect_sign);

save([dest '_cell_resp.mat'], 'NR_Cue_base_pos', 'NR_Cue_resp_pos', 'OR_Cue_base_pos', 'OR_Cue_base2_pos', 'OR_Cue_resp_pos', 'OR_Cue_resp2_pos', 'UR_Cue_base_pos', 'UR_Cue_resp_pos',...
    'NR_Rew_base_pos', 'NR_Rew_resp_pos', 'OR_Rew_base_pos', 'OR_Rew_resp_pos', 'UR_Rew_base_pos', 'UR_Rew_resp_pos', ...
	'NR_Cue_base_neg', 'NR_Cue_resp_neg', 'OR_Cue_base_neg', 'OR_Cue_base2_neg', 'OR_Cue_resp_neg', 'OR_Cue_resp2_neg', 'UR_Cue_base_neg', 'UR_Cue_resp_neg', ...
    'NR_Rew_base_neg', 'NR_Rew_resp_neg', 'OR_Rew_base_neg', 'OR_Rew_resp_neg', 'UR_Rew_base_neg', 'UR_Rew_resp_neg');

%Make combinations of responsive subtypes for plotting later
if ~isempty(OR_movie)
    allresp_cells = NR_Cue_h_pos | OR_Cue_h_pos | NR_Rew_h_pos | OR_Rew_h_pos | NR_Cue_h_neg | OR_Cue_h_neg | NR_Rew_h_neg | OR_Rew_h_neg;
    allresp_cells_pos = NR_Cue_h_pos | OR_Cue_h_pos | NR_Rew_h_pos | OR_Rew_h_pos;
    allresp_cells_neg = NR_Cue_h_neg | OR_Cue_h_neg | NR_Rew_h_neg | OR_Rew_h_neg;
    cue_cells = NR_Cue_h_pos | OR_Cue_h_pos | NR_Cue_h_neg | OR_Cue_h_neg;
    cue_cells_pos = NR_Cue_h_pos | OR_Cue_h_pos;
    cue_cells_neg = NR_Cue_h_neg | OR_Cue_h_neg;
    rew_cells = NR_Rew_h_pos | OR_Rew_h_pos | NR_Rew_h_neg | OR_Rew_h_neg;
    rew_cells_pos = NR_Rew_h_pos | OR_Rew_h_pos;
    rew_cells_neg = NR_Rew_h_neg | OR_Rew_h_neg;
else
    allresp_cells = NR_Cue_h_pos | NR_Rew_h_pos | NR_Cue_h_neg | NR_Rew_h_neg;
    allresp_cells_pos = NR_Cue_h_pos | NR_Rew_h_pos;
    allresp_cells_neg = NR_Cue_h_neg | NR_Rew_h_neg;
    cue_cells = NR_Cue_h_pos | NR_Cue_h_neg;
    cue_cells_pos = NR_Cue_h_pos;
    cue_cells_neg = NR_Cue_h_neg;
    rew_cells = NR_Rew_h_pos | NR_Rew_h_neg;
    rew_cells_pos = NR_Rew_h_pos;
    rew_cells_neg = NR_Rew_h_neg;
end

if ~isempty(UR_movie)
    allresp_cells = allresp_cells | UR_Cue_h_pos | UR_Rew_h_pos | UR_Cue_h_neg | UR_Rew_h_neg;
    allresp_cells_pos = allresp_cells_pos | UR_Cue_h_pos | UR_Rew_h_pos;
    allresp_cells_neg = allresp_cells_neg | UR_Cue_h_neg | UR_Rew_h_neg;
    cue_cells = cue_cells | UR_Cue_h_pos | UR_Cue_h_neg;
    cue_cells_pos = cue_cells_pos | UR_Cue_h_pos;
    cue_cells_neg = cue_cells_neg | UR_Cue_h_neg;
    rew_cells = rew_cells | UR_Rew_h_pos | UR_Rew_h_neg;
    rew_cells_pos = rew_cells_pos | UR_Rew_h_pos;
    rew_cells_neg = rew_cells_neg | UR_Rew_h_neg;
end

cue_cells_perc = sum(cue_cells) / nCells;
cue_cells_perc_pos = sum(cue_cells_pos) / nCells;
cue_cells_perc_neg = sum(cue_cells_neg) / nCells;
rew_cells_perc = sum(rew_cells) / nCells;
rew_cells_perc_pos = sum(rew_cells_pos) / nCells;
rew_cells_perc_neg = sum(rew_cells_neg) / nCells;

save([dest '_cell_categories.mat'], 'NR_Cue_resp_cells_pos', 'NR_Cue_resp_cells_neg', 'OR_Cue_resp_cells_pos', 'OR_Cue_resp_cells_neg', 'OR_Cue_resp_cells2_pos', 'OR_Cue_resp_cells2_neg', ...
    'UR_Cue_resp_cells_pos', 'UR_Cue_resp_cells_neg', 'NR_Rew_resp_cells_pos', 'NR_Rew_resp_cells_neg', 'OR_Rew_resp_cells_pos', 'OR_Rew_resp_cells_neg', 'NR_Rew_resp_cells_pos', 'NR_Rew_resp_cells_neg',...
    'allresp_cells', 'allresp_cells_pos', 'allresp_cells_neg', 'nCells', 'cue_cells', 'cue_cells_pos', 'cue_cells_neg',  'rew_cells', 'rew_cells_pos', 'rew_cells_neg', ...
    'cue_cells_perc', 'cue_cells_perc_pos', 'cue_cells_perc_neg', 'rew_cells_perc', 'rew_cells_perc_pos', 'rew_cells_perc_neg');
save([dest '_pvals.mat'], 'NR_Cue_p_pos', 'NR_Cue_p_neg', 'NR_Cue_h_pos', 'NR_Cue_h_neg', 'NR_Rew_p_pos', 'NR_Rew_p_neg', 'NR_Rew_h_pos', 'NR_Rew_h_neg', ...
    'OR_Cue_p_pos', 'OR_Cue_p_neg', 'OR_Cue_p2_pos', 'OR_Cue_p2_neg', 'OR_Cue_h_pos', 'OR_Cue_h_neg', 'OR_Cue_h2_pos', 'OR_Cue_h2_neg', 'OR_Rew_p_pos', 'OR_Rew_p_pos', 'OR_Rew_h_pos', 'OR_Rew_h_neg', ...
    'UR_Cue_p_pos', 'UR_Cue_p_neg', 'UR_Cue_h_pos', 'UR_Cue_h_neg', 'UR_Rew_p_pos', 'UR_Rew_p_neg', 'UR_Rew_h_pos', 'UR_Rew_h_neg');

%% 3. Find the 10% rise times of the mean TCs for each trial type.

% [NR_riseIdx] = findRisetime(avg_NR_nolick, pre_cue_frames); %need to modify so it only includes positively modulated neurons
% 
% if ~isempty(OR_movie)
%     [OR_riseIdx] = findRisetime(avg_OR_nolick, pre_cue_frames);
% else
%     OR_riseIdx = [];
% end
% 
% if ~isempty(UR_movie)
%     [UR_riseIdx] = findRisetime(avg_UR_nolick, pre_cue_frames);
% else
%     UR_riseIdx = [];
% end
% 
% save([dest '_F_risetime.mat'], 'NR_riseIdx', 'OR_riseIdx', 'UR_riseIdx');

%% another way to find peak response
% % find resp magnitude to base
% base_window = pre_cue_frames + round(350./double(ifi)) : pre_cue_frames + round(500./double(ifi)) -1; %selects a baseline window. Should rewrite this to use baseline times
% % resp_window = pre_cue_frames + round(500./double(ifi)): pre_release_frames + round(700./double(ifi));  %selecting a response window in which to analyze the TC. Selects three consecutive frames where the lever event occurs during the first fraem
% 
% poss_NR_wins = [];  %using a sliding window to determine peak response
% poss_NR_win_vals = []; 
% for i = pre_cue_frames + round(500./double(ifi)):pre_cue_frames + round(2000./double(ifi))
%     poss_NR_win_vals = [poss_NR_win_vals; mean(mean(mean(NR_movie(:,:,i-1:i+1)),3))];
%     poss_NR_wins = [poss_NR_wins, i];
% end
% resp_window = poss_NR_wins(find(poss_NR_win_vals == max(poss_NR_win_vals))); %selects window with the peak response 
% resp_window = [resp_window-1, resp_window, resp_window+1]; 
% 
% NR_base = squeeze(mean(NR_movie(:,:,base_window),3));
% NR_resp = squeeze(mean(NR_movie(:,:,resp_window),3));
% 
% poss_OR_wins = [];  %using a sliding window to determine peak response
% poss_OR_win_vals = []; 
% for i = pre_cue_frames + round(500./double(ifi)):pre_cue_frames + round(2000./double(ifi))
%     poss_OR_win_vals = [poss_OR_win_vals; mean(mean(mean(OR_movie(:,:,i-1:i+1)),3))];
%     poss_OR_wins = [poss_OR_wins, i];
% end
% resp_window = poss_OR_wins(find(poss_OR_win_vals == max(poss_OR_win_vals))); %selects window with the peak response 
% resp_window = [resp_window-1, resp_window, resp_window+1]; 
% 
% OR_base = squeeze(mean(OR_movie(:,:,base_window),3));
% OR_resp = squeeze(mean(OR_movie(:,:,resp_window),3));
% 
% poss_UR_wins = [];  %using a sliding window to determine peak response
% poss_UR_win_vals = []; 
% for i = pre_cue_frames + round(500./double(ifi)):pre_cue_frames + round(2000./double(ifi))
%     poss_UR_win_vals = [poss_UR_win_vals; mean(mean(mean(OR_movie(:,:,i-1:i+1)),3))];
%     poss_UR_wins = [poss_UR_wins, i];
% end
% resp_window = poss_UR_wins(find(poss_UR_win_vals == max(poss_UR_win_vals))); %selects window with the peak response 
% resp_window = [resp_window-1, resp_window, resp_window+1]; 
% 
% UR_base = squeeze(mean(UR_movie(:,:,base_window),3));
% UR_resp = squeeze(mean(UR_movie(:,:,resp_window),3));
% 
% save([dest '_cell_resp2.mat'], 'NR_base', 'NR_resp', 'UR_base', 'UR_resp', 'OR_base', 'OR_resp');

%% 4. plotting
%plot all no-lick trials
tt =((-pre_cue_frames:post_cue_frames).*double(ifi))./1000;
%overlay average success/failure for each ROI
figure;
avg_all_nolick = [avg_NR_nolick avg_OR_nolick avg_UR_nolick];
ymax = max(max(avg_all_nolick,[],2),[],1);
ymin = -0.1;
% ymin = min(min(avg_all_nolick,[],2),[],1);
for ic = 1:nCells
    subplot_tight(n,n2,ic);
%     errorbar(tt,avg_NR(ic,:), sem_NR(ic,:),'k');hold on;
%     errorbar(tt,avg_OR(ic,:), sem_OR(ic,:),'r');
%     errorbar(tt,avg_UR(ic,:), sem_UR(ic,:),'b');
    plot(tt,avg_NR_nolick(ic,:),'k');hold on;
    if ~isempty(OR_movie_nolick);
        plot(tt,avg_OR_nolick(ic,:),'r');
    end
    if ~isempty(UR_movie_nolick);
        plot(tt,avg_UR_nolick(ic,:),'g');
    end
    ylim([ymin*1.1 ymax*1.1]);
    xlim([-1 2]);
%     set(gca,'XTick',[)
end
suptitle([session_date ' ' mouse_ID ' Average all cells: Normal- black (n = ' num2str(size(NR_movie_nolick,1)) ' trials); Omit- red (n = ' num2str(size(OR_movie_nolick,1)) ' trials); Unexpect- green (n = ' num2str(size(UR_movie_nolick,1)) ' trials)']);
orient landscape
print([dest '_allTCs_nolick.eps'], '-depsc');
print([dest '_allTCs_nolick.pdf'], '-dpdf');

%plot all responsive no-lick trials
figure;
for ic = find(allresp_cells)
    subplot_tight(n,n2,ic);
    plot(tt,avg_NR_nolick(ic,:),'k');hold on;
    if ~isempty(OR_movie_nolick);
        plot(tt,avg_OR_nolick(ic,:),'r');
    end
    if ~isempty(UR_movie_nolick);
        plot(tt,avg_UR_nolick(ic,:),'g');
    end
    ylim([ymin*1.1 ymax*1.1]);
    xlim([-1 2]);
%     set(gca,'XTick',[)
end
suptitle([session_date ' ' mouse_ID ' Average all resp no lick: Normal- black (n = ' num2str(size(NR_movie_nolick,1)) ' trials); Omit- red (n = ' num2str(size(OR_movie_nolick,1)) ' trials); Unexpect- green (n = ' num2str(size(UR_movie_nolick,1)) ' trials). ', 'resp cells (n=', num2str(length(find(allresp_cells))),')']);
orient landscape
print([dest '_allTCs_nolick_allresp.eps'], '-depsc');
print([dest '_allTCs_nolick_allresp.pdf'], '-dpdf');

%overlay average success/failure across ROIs
h=figure; 
subplot(3,3,1)
errorbar(tt,avg_NR_all, sem_NR_all,'k')
hold on;
if ~isempty(OR_movie)
    errorbar(tt,avg_OR_all, sem_OR_all,'r');
end
if ~isempty(UR_movie)
    errorbar(tt,avg_UR_all, sem_UR_all,'g')
end
% bar(tt,mean(lick_trace_NR)/5, 'k');
% bar(tt,mean(lick_trace_OR)/5, 'r');
% if ~isempty(UR_movie)
%     bar(tt,mean(lick_trace_UR), 'g-');
% end
title(['Avg all cells (inlcudes lick trials); n = ' num2str(size(avg_NR,1))])
ylabel('dF/F')
xlabel('Time (s)')
xlim([-floor(pre_cue_frames/ifi) floor(post_cue_frames/ifi)]);

% plot NR resp cells
NR_Cue_h = NR_Cue_h_pos | NR_Cue_h_neg | NR_Rew_h_pos | NR_Rew_h_neg;
subplot(3,3,2)
errorbar(tt,mean(avg_NR(NR_Cue_h,:),1),std(avg_NR(NR_Cue_h,:),[],1)./sqrt(sum(NR_Cue_h,2)),'k')
hold on;

bar(tt,mean(lick_trace_NR)/8, 'k', 'ShowBaseLine', 'off');

title(['Avg NR trial resp cells p<0.05; (inlcudes lick trials) n = ' num2str(sum(NR_Cue_h,2))])
ylabel('dF/F')
xlabel('Time (s)')
xlim([-floor(pre_cue_frames/ifi) floor(post_cue_frames/ifi)]);

% plot OR resp cells
if ~isempty(OR_movie)
    OR_Cue_h = OR_Cue_h_pos | OR_Cue_h_neg | OR_Rew_h_pos | OR_Rew_h_neg;
subplot(3,3,3)
errorbar(tt,mean(avg_OR(OR_Cue_h,:),1),std(avg_OR(OR_Cue_h,:),[],1)./sqrt(sum(OR_Cue_h,2)),'r')
hold on;

bar(tt,mean(lick_trace_OR)/8, 'FaceColor', [1,0,0], 'EdgeColor', [1,0,0], 'ShowBaseLine', 'off');

title(['Avg OR trial all resp cells p<0.05; (inlcudes lick trials) n = ' num2str(sum(OR_Cue_h,2))])
ylabel('dF/F')
xlabel('Time (s)')
xlim([-floor(pre_cue_frames/ifi) floor(post_cue_frames/ifi)]);
end

% plot UR resp cells
if ~isempty(UR_movie)
    UR_Cue_h = UR_Cue_h_pos | UR_Cue_h_neg | UR_Rew_h_pos | UR_Rew_h_neg;
subplot(3,3,4)
errorbar(tt,mean(avg_UR(UR_Cue_h,:),1),std(avg_UR(UR_Cue_h,:),[],1)./sqrt(sum(UR_Cue_h,2)),'g')
hold on;

bar(tt,mean(lick_trace_UR)/8, 'FaceColor', [0,1,0], 'EdgeColor', [0,1,0], 'ShowBaseLine', 'off');

title(['Avg UR trial all resp cells p<0.05; n = ' num2str(sum(UR_Cue_h,2))])
ylabel('dF/F')
xlabel('Time (s)')
xlim([-floor(pre_cue_frames/ifi) floor(post_cue_frames/ifi)]);
end

% plot NR and OR resp cells
if ~isempty(OR_movie)
    NR_OR_h = OR_Cue_h | NR_Cue_h;
subplot(3,3,5)
errorbar(tt,mean(avg_NR(NR_OR_h,:),1),std(avg_NR(NR_OR_h,:),[],1)./sqrt(sum(NR_OR_h,2)),'k')
hold on;
errorbar(tt,mean(avg_OR(NR_OR_h,:),1),std(avg_OR(NR_OR_h,:),[],1)./sqrt(sum(NR_OR_h,2)),'r')
% bar(tt,mean(lick_trace_NR)/8, 'k');

title(['Avg NR and OR trial all resp cells p<0.05; (inlcudes lick trials) n = ' num2str(sum(NR_OR_h,2))])
ylabel('dF/F')
xlabel('Time (s)')
xlim([-floor(pre_cue_frames/ifi) floor(post_cue_frames/ifi)]);
end

% % plot NR and UR resp cells
% if ~isempty(UR_movie)
% subplot(3,3,6)
% errorbar(tt,mean(avg_NR(NUR_h,:),1),std(avg_NR(NUR_h,:),[],1)./sqrt(sum(NUR_h,2)),'k')
% hold on;
% errorbar(tt,mean(avg_UR(NUR_h,:),1),std(avg_UR(NUR_h,:),[],1)./sqrt(sum(NUR_h,2)),'g')
% % bar(tt,mean(lick_trace_NR)/8, 'k');
% 
% title(['Avg Normal and Unexpected Reward resp cells p<0.05; n = ' num2str(sum(NUR_h,2))])
% ylabel('dF/F')
% xlabel('Time (s)')
% xlim([-floor(pre_cue_frames/ifi) floor(post_cue_frames/ifi)]);
% end
% 
% % plot OR and UR resp cells
% if ~isempty(UR_movie) && ~isempty(OR_movie)
% subplot(3,3,7)
% errorbar(tt,mean(avg_OR(OUR_h,:),1),std(avg_OR(OUR_h,:),[],1)./sqrt(sum(OUR_h,2)),'r')
% hold on;
% errorbar(tt,mean(avg_UR(OUR_h,:),1),std(avg_UR(OUR_h,:),[],1)./sqrt(sum(OUR_h,2)),'g')
% % bar(tt,mean(lick_trace_NR)/8, 'k');
% 
% title(['Avg Omitted and Unexpected Reward resp cells p<0.05; n = ' num2str(sum(OUR_h,2))])
% ylabel('dF/F')
% xlabel('Time (s)')
% xlim([-floor(pre_cue_frames/ifi) floor(post_cue_frames/ifi)]);
% end

allresp_cells_temp = allresp_cells;
subplot(3,3,8)
errorbar(tt,mean(avg_NR(find(allresp_cells_temp),:),1),std(avg_NR(find(allresp_cells_temp),:),[],1)./sqrt(sum(allresp_cells_temp,2)),'k')
hold on;
errorbar(tt,mean(avg_OR(find(allresp_cells_temp),:),1),std(avg_OR(find(allresp_cells_temp),:),[],1)./sqrt(sum(allresp_cells_temp,2)),'r')
errorbar(tt,mean(avg_UR(find(allresp_cells_temp),:),1),std(avg_UR(find(allresp_cells_temp),:),[],1)./sqrt(sum(allresp_cells_temp,2)),'g')

bar(tt,mean(lick_trace_NR)/8, 'k', 'ShowBaseLine', 'off');
if ~isempty(OR_movie)
    bar(tt,mean(lick_trace_OR)/8, 'FaceColor', [1,0,0], 'EdgeColor', [1,0,0], 'ShowBaseLine', 'off');
end

if ~isempty(UR_movie)
    bar(tt,mean(lick_trace_UR/8), 'FaceColor', [0,1,0], 'EdgeColor', [0,1,0], 'ShowBaseLine', 'off');
end

title(['Avg all resp cells p<0.05; (inlcudes lick trials)n = ' num2str(sum(allresp_cells_temp,2))])
ylabel('dF/F')
xlabel('Time (s)')
xlim([-floor(pre_cue_frames/ifi) floor(post_cue_frames/ifi)]);

suptitle([session_date ' ' mouse_ID ' Average: Normal- black; Omit- red; Unexpect- green bars- licking data'])
print([dest '_avgTCs_pos.eps'], '-depsc');
print([dest '_avgTCs_pos.pdf'], '-dpdf');
savefig(h, [dest '_avgTCs_pos']); 
% 
% figure;
% if ~isempty(UR_movie) && ~isempty(OR_movie)
%     all_resp = [NR_resp_avg OR_resp_avg UR_resp_avg];
% elseif isempty(UR_movie) && ~isempty(OR_movie)
%     all_resp = [NR_resp_avg OR_resp_avg];
% elseif ~isempty(UR_movie) && isempty(OR_movie)
%     all_resp = [NR_resp_avg UR_resp_avg];
% elseif isempty(UR_movie) && isempty(OR_movie)
%     all_resp = [NR_resp_avg];
% end
% cmax = max(all_resp,[],2);
% cmin = min(all_resp,[],2);
% subplot(2,2,1)
% mask_final_temp = zeros(size(mask_final));
% mask_final_temp(find(mask_final>0)) = 1;
% imagesq(reshape(mask_final_temp, [sz(1) sz(2)]));
% clim([cmin cmax]);
% set(gca, 'XTickLabel', '', 'YTickLabel', '', 'Xtick', 0,'Ytick',0)
% colorbar
% 
% NR_mask = mask_final;
% for ic = 1:nCells
%     if NR_h(ic)
%         NR_mask(find(NR_mask==ic))=NR_resp_avg(1,ic);
%     else
%         NR_mask(find(NR_mask==ic))=0;
%     end
% end
% NR_mask = reshape(NR_mask,[sz(1) sz(2)]);
% subplot(2,2,2)
% imagesq(NR_mask);
% clim([cmin cmax]);
% set(gca, 'XTickLabel', '', 'YTickLabel', '', 'Xtick', 0,'Ytick',0)
% title('normal')
% colorbar
% 
% if ~isempty(OR_movie)
% OR_mask = mask_final;
% for ic = 1:nCells
%     if OR_h(ic)
%         OR_mask(find(OR_mask==ic))=OR_resp_avg(1,ic);
%     else
%         OR_mask(find(OR_mask==ic))=0;
%     end
% end
% OR_mask = reshape(OR_mask,[sz(1) sz(2)]);
% subplot(2,2,3)
% imagesq(OR_mask);
% clim([cmin cmax])
% set(gca, 'XTickLabel', '', 'YTickLabel', '', 'Xtick', 0,'Ytick',0)
% title('Omit')
% colorbar
% end
% 
% if ~isempty(UR_movie)
% UR_mask = mask_final;
% for ic = 1:nCells
%     if UR_h(ic)
%         UR_mask(find(UR_mask==ic))=UR_resp_avg(1,ic);
%     else
%         UR_mask(find(UR_mask==ic))=0;
%     end
% end
% UR_mask = reshape(UR_mask,[sz(1) sz(2)]);
% subplot(2,2,4)
% imagesq(UR_mask);
% clim([cmin cmax]);
% set(gca, 'XTickLabel', '', 'YTickLabel', '', 'Xtick', 0,'Ytick',0)
% colorbar
% title('Unexpect')
% end
% 
% suptitle([mouse_ID ' ' session_date ' cell responses']);
% print([dest '_cell_responses_FOV.eps'], '-depsc');
% print([dest '_cell_responses_FOV.pdf'], '-dpdf');

%plot amplitude response by cell
% lever = strvcat('NR', 'release','success', 'fail');
% figure;
% subplot(2,2,1)
% for ic = 1:nCells
%     plot(1:4, [NR_resp_avg(1,ic) release_resp_avg(1,ic) success_resp_avg(1,ic) fail_resp_avg(1,ic)], '-ok')
%     hold on
%     if NR_h(1,ic)
%         plot(1,NR_resp_avg(1,ic),'or')
%         hold on
%     end
%     if release_h(1,ic)
%         plot(2,release_resp_avg(1,ic),'or')
%         hold on
%     end
%     if success_h(1,ic)
%         plot(3,success_resp_avg(1,ic),'or')
%         hold on
%     end
%     if fail_h(1,ic)
%         plot(4,fail_resp_avg(1,ic),'or')
%         hold on
%     end
% end
% set(gca, 'XTick', 1:4, 'XTickLabel', lever);
% xlim([0.5 4.5])
% title('Resp amp- red: p<0.05')
% ylabel('dF/F')
% 
% subplot(2,2,2)
% errorbar(1:4, [mean(NR_resp_avg,2) mean(release_resp_avg,2) mean(success_resp_avg,2) mean(fail_resp_avg,2)], [std(NR_resp_avg,[],2) std(release_resp_avg,[],2) std(success_resp_avg,[],2) std(fail_resp_avg,[],2)]./sqrt(nCells),'ok')
% set(gca, 'XTick', 1:4, 'XTickLabel', lever);
% xlim([0.5 4.5])
% ylabel('dF/F')
% title(['Avg resp- all cells: n = ' num2str(nCells)])
% subplot(2,2,3)
% errorbar(1:4, [mean(NR_resp_avg(1,find(release_h)),2) mean(release_resp_avg(1,find(release_h)),2) mean(success_resp_avg(1,find(release_h)),2) mean(fail_resp_avg(1,find(release_h)),2)], [std(NR_resp_avg(1,find(release_h)),[],2) std(release_resp_avg(1,find(release_h)),[],2) std(success_resp_avg(1,find(release_h)),[],2) std(fail_resp_avg(1,find(release_h)),[],2)]./sqrt(sum(release_h,2)),'ok')
% set(gca, 'XTick', 1:4, 'XTickLabel', lever);
% xlim([0.5 4.5])
% ylabel('dF/F')
% title(['Avg resp- release resp cells: n = ' num2str(sum(release_h,2))])
% subplot(2,2,4)
% errorbar(1:4, [mean(NR_resp_avg(1,find(NR_h)),2) mean(release_resp_avg(1,find(NR_h)),2) mean(success_resp_avg(1,find(NR_h)),2) mean(fail_resp_avg(1,find(NR_h)),2)], [std(NR_resp_avg(1,find(NR_h)),[],2) std(release_resp_avg(1,find(NR_h)),[],2) std(success_resp_avg(1,find(NR_h)),[],2) std(fail_resp_avg(1,find(NR_h)),[],2)]./sqrt(sum(NR_h,2)),'ok')
% set(gca, 'XTick', 1:4, 'XTickLabel', lever);
% xlim([0.5 4.5])
% ylabel('dF/F')
% title(['Avg resp- NR resp cells: n = ' num2str(sum(NR_h,2))])
% suptitle([session_date ' ' mouse_ID 'avg resp amplitude'])
% print([dest '_cell_resp_amplitude.eps'], '-depsc');
% print([dest '_cell_resp_amplitude.pdf'], '-dpdf');


