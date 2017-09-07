% finds and maps responsive cells for each condition
% 1. calculate average timecourses for NR/release events
% 2. calculate 10% rise time of F onset
% 3. calculate variability by trial over base and resp windows, 
%    ttest for significant responses
% 4. define cells by response to event

% load([dest 'ROI_TCs.mat']);
load([dest 'parse_behavior.mat']);

load([dest '_cue_movies.mat'])
nCells = size(NR_movie,2);
n = ceil(sqrt(nCells));
if nCells <((n.^2)-n) %function
    n2= n-1;
else
    n2 = n;
end

%% 1. calculate average timecourses for NR/release events

%average and sem across Trials
avg_NR = squeeze(mean(NR_movie,1));
sem_NR = squeeze(std(NR_movie,1)./sqrt(size(NR_movie,1)));
%average and sem across ROIs
avg_NR_all = mean(avg_NR,1);
sem_NR_all = std(avg_NR,1)./sqrt(size(avg_NR,1));

%average and sem across Trials
avg_OR = squeeze(mean(OR_movie,1));
sem_OR = squeeze(std(OR_movie,1)./sqrt(size(OR_movie,1)));
%average and sem across ROIs
avg_OR_all = mean(avg_OR,1);
sem_OR_all = std(avg_OR,1)./sqrt(size(avg_OR,1));

%average and sem across Trials
avg_UR = squeeze(mean(UR_movie,1));
sem_UR = squeeze(std(UR_movie,1)./sqrt(size(UR_movie,1)));
%average and sem across ROIs
avg_UR_all = mean(avg_UR,1);
sem_UR_all = std(avg_UR,1)./sqrt(size(avg_UR,1));

%average and sem across Trials
avg_NR_nolick = squeeze(nanmean(NR_movie_nolick,1));
NR_movie_nolick_temp = NR_movie_nolick(~isnan(NR_movie_nolick(:,1,1)),:,:);
sem_NR_nolick = squeeze(nanstd(NR_movie_nolick,1)./sqrt(size(NR_movie_nolick_temp,1)));
%average and sem across ROIs
avg_NR_all_nolick = mean(avg_NR_nolick,1);
sem_NR_all_nolick = std(avg_NR_nolick,1)./sqrt(size(avg_NR_nolick,1));

%average and sem across Trials
avg_OR_nolick = squeeze(nanmean(OR_movie_nolick,1));
OR_movie_nolick_temp = OR_movie_nolick(~isnan(OR_movie_nolick(:,1,1)),:,:);
sem_OR_nolick = squeeze(nanstd(OR_movie_nolick,1)./sqrt(size(OR_movie_nolick_temp,1)));
%average and sem across ROIs
avg_OR_all_nolick = mean(avg_OR_nolick,1);
sem_OR_all_nolick = std(avg_OR,1)./sqrt(size(avg_OR_nolick,1));

%average and sem across Trials
avg_UR_nolick = squeeze(nanmean(UR_movie_nolick,1));
UR_movie_nolick_temp = UR_movie_nolick(~isnan(UR_movie_nolick(:,1,1)),:,:);
sem_UR_nolick = squeeze(nanstd(UR_movie_nolick,1)./sqrt(size(UR_movie_nolick_temp,1)));
%average and sem across ROIs
avg_UR_all_nolick = mean(avg_UR_nolick,1);
sem_UR_all_nolick = std(avg_UR_nolick,1)./sqrt(size(avg_UR_nolick,1));

save([dest '_cell_TCs.mat'], 'avg_NR', 'sem_NR', 'avg_OR', 'sem_OR', 'avg_UR', 'sem_UR', ...
    'avg_NR_nolick', 'sem_NR_nolick', 'avg_OR_nolick', 'sem_OR_nolick', 'avg_UR_nolick', 'sem_UR_nolick')


% % Plot timecourse for each cell 
% tt =((-pre_cue_frames:post_cue_frames).*double(ifi))./1000;
% %overlay average success/failure for each ROI
% fig=figure;
% avg_all = [avg_NR_nolick avg_OR_nolick];
% ymax = max(max(avg_all,[],2),[],1);
% ymin = min(min(avg_all,[],2),[],1);
% for ic = 1:nCells
%     subplot(n,n2,ic)
%     errorbar(tt,avg_NR_nolick(ic,:), sem_NR_nolick(ic,:),'k')
%     hold on;
%     errorbar(tt,avg_OR_nolick(ic,:), sem_OR_nolick(ic,:),'r')
%     ylim([ymin*1.1 ymax*1.1])
%     xlim([tt(1) tt(end)])
% end
% supertitle([date ' ' mouse ' timecourse of normal reward- black (n = ' num2str(size(NR_movie_nolick_temp,1)) ' trials); Omission- red (n = ' num2str(size(OR_movie_nolick_temp,1)) ' trials)'])
% orient landscape
% saveas(fig, [dest_sub '_NR_OR_allTCs.fig']);
% print([dest_sub '_NR_OR_allTCs.eps'], '-depsc');
% print([dest_sub '_NR_OR_allTCs.pdf'], '-dpdf');

%% 2. calculate 10% rise time of F onset
% [NR_riseIdx, NR_resp] = findRisetime(avg_NR, pre_cue_frames);
% 
% if ~isempty(OR_movie)
%     [OR_riseIdx, OR_resp] = findRisetime(avg_OR, pre_cue_frames);
% else
%     OR_riseIdx = [];
%     OR_resp = [];
% end
% 
% if ~isempty(UR_movie)
%     [UR_riseIdx, UR_resp] = findRisetime(avg_UR, pre_cue_frames);
% else
%     UR_riseIdx = [];
%     UR_resp = [];
% end
% 
% save([dest '_F_risetime.mat'], 'NR_riseIdx', 'OR_riseIdx', 'UR_riseIdx');

%% 3. calculate response amplitude and variablity by trial over base and resp windows

base_cue_buffer = 500;
resp_cue_buffer = 400;

[NR_Cue_h, NR_Cue_p, NR_Cue_resp_cells, NR_Cue_resp_avg, NR_Cue_resp_sem, NR_Cue_base, NR_Cue_resp, ~, ~] = findRespCell(NR_movie_nolick, pre_cue_frames, ifi, base_cue_buffer, resp_cue_buffer);

[OR_Cue_h, OR_Cue_p, OR_Cue_resp_cells, OR_Cue_resp_avg, OR_Cue_resp_sem, OR_Cue_base, OR_Cue_resp, ~, ~] = findRespCell(OR_movie_nolick, pre_cue_frames, ifi, base_cue_buffer, resp_cue_buffer);

% omission trial responsive window from cue to reward 0-600
[OR_Cue_h2, OR_Cue_p2, OR_Cue_resp_cells2, OR_Cue_resp_avg2, OR_Cue_resp_sem2, OR_Cue_base2, OR_Cue_resp2, ~, ~] = findRespCell(OR_movie_nolick, pre_cue_frames, ifi, base_cue_buffer, resp_cue_buffer+200);

[UR_Cue_h, UR_Cue_p, UR_Cue_resp_cells, UR_Cue_resp_avg, UR_Cue_resp_sem, UR_Cue_base, UR_Cue_resp, ~, ~] = findRespCell(UR_movie_nolick, pre_cue_frames, ifi, base_cue_buffer, resp_cue_buffer);

base_reward_buffer = 500;
resp_reward_buffer = 2000;
pre_rew_frames = pre_cue_frames + round((trial_outcome.normalReward(1) - trial_outcome.normalRewardCue(1))/ifi);

[NR_Rew_h, NR_Rew_p, NR_Rew_resp_cells, NR_Rew_resp_avg, NR_Rew_resp_sem, NR_Rew_base, NR_Rew_resp, ~, ~] = findRespCell(NR_movie_nolick, pre_rew_frames, ifi, base_reward_buffer, resp_reward_buffer);

[OR_Rew_h, OR_Rew_p, OR_Rew_resp_cells, OR_Rew_resp_avg, OR_Rew_resp_sem, OR_Rew_base, OR_Rew_resp, ~, ~] = findRespCell(OR_movie_nolick, pre_rew_frames, ifi, base_reward_buffer, resp_reward_buffer);

[UR_Rew_h, UR_Rew_p, UR_Rew_resp_cells, UR_Rew_resp_avg, UR_Rew_resp_sem, UR_Rew_base, UR_Rew_resp, ~, ~] = findRespCell(UR_movie_nolick, pre_rew_frames, ifi, base_reward_buffer, resp_reward_buffer);

% find cells with negative response and positive response shortly after
% reward
base_reward_buffer = 200;
resp_reward_buffer = 400;

[~, ~, ~, ~, ~, NR_Rew_pn_base, NR_Rew_pn_resp, NR_Rew_pos_h, NR_Rew_neg_h] = findRespCell(NR_movie_nolick, pre_rew_frames, ifi, base_reward_buffer, resp_reward_buffer);

[~, ~, ~, ~, ~, OR_Rew_pn_base, OR_Rew_pn_resp, OR_Rew_pos_h, OR_Rew_neg_h] = findRespCell(OR_movie_nolick, pre_rew_frames, ifi, base_reward_buffer, resp_reward_buffer);

save([dest '_cell_resp.mat'], 'NR_Cue_base', 'NR_Cue_resp', 'OR_Cue_base', 'OR_Cue_base2', 'OR_Cue_resp', 'OR_Cue_resp2', 'UR_Cue_base', 'UR_Cue_resp',...
    'NR_Rew_base', 'NR_Rew_resp', 'OR_Rew_base', 'OR_Rew_resp', 'UR_Rew_base', 'UR_Rew_resp', 'NR_Rew_pn_base', 'NR_Rew_pn_resp', 'OR_Rew_pn_base', 'OR_Rew_pn_resp');

% noresponse_cells = find((NR_h+OR_h+UR_h)==0);

if ~isempty(OR_movie)
    allresp_cells = NR_Cue_h | OR_Cue_h | NR_Rew_h | OR_Rew_h;
    cue_cells = NR_Cue_h | OR_Cue_h;
    rew_cells = NR_Rew_h | OR_Rew_h;
else
    allresp_cells = NR_Cue_h | NR_Rew_h;
    cue_cells = NR_Cue_h;
    rew_cells = NR_Rew_h;
end

if ~isempty(UR_movie)
    allresp_cells = allresp_cells | UR_Cue_h | UR_Rew_h;
    cue_cells = cue_cells | UR_Cue_h;
    rew_cells = rew_cells | UR_Rew_h;
end

cue_cells_perc = sum(cue_cells) / nCells;
rew_cells_perc = sum(rew_cells) / nCells;

NR_Rew_resp_pos_cells = find(NR_Rew_pos_h);
NR_Rew_resp_neg_cells = find(NR_Rew_neg_h);

OR_Rew_resp_pos_cells = find(OR_Rew_pos_h);
OR_Rew_resp_neg_cells = find(OR_Rew_neg_h);

save([dest '_cell_categories.mat'], 'NR_Cue_resp_cells', 'OR_Cue_resp_cells', 'OR_Cue_resp_cells2', 'UR_Cue_resp_cells', 'NR_Rew_resp_cells',...
    'OR_Rew_resp_cells', 'NR_Rew_resp_cells', 'NR_Rew_resp_pos_cells', 'NR_Rew_resp_neg_cells', 'OR_Rew_resp_pos_cells', 'OR_Rew_resp_neg_cells', ...
    'allresp_cells', 'nCells', 'cue_cells', 'rew_cells', 'cue_cells_perc', 'rew_cells_perc');

save([dest '_pvals.mat'], 'NR_Cue_p', 'NR_Cue_h', 'NR_Rew_p', 'NR_Rew_h', 'OR_Cue_p', 'OR_Cue_p2','OR_Cue_h', 'OR_Cue_h2', 'OR_Rew_p', 'OR_Rew_h', ...
    'UR_Cue_p', 'UR_Cue_h', 'UR_Rew_p', 'UR_Rew_h', 'NR_Rew_pos_h', 'NR_Rew_neg_h', 'OR_Rew_pos_h', 'OR_Rew_pos_h');

[NR_riseIdx] = findRisetime(avg_NR_nolick, pre_cue_frames);

if ~isempty(OR_movie)
    [OR_riseIdx] = findRisetime(avg_OR_nolick, pre_cue_frames);
else
    OR_riseIdx = [];
end

if ~isempty(UR_movie)
    [UR_riseIdx] = findRisetime(avg_UR_nolick, pre_cue_frames);
else
    UR_riseIdx = [];
end

save([dest '_F_risetime.mat'], 'NR_riseIdx', 'OR_riseIdx', 'UR_riseIdx');

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

%%
% %% 4. plotting
% tt =((-pre_cue_frames:post_cue_frames).*double(ifi))./1000;
% %overlay average success/failure for each ROI
% figure;
% avg_all = [avg_NR avg_OR avg_UR];
% ymax = max(max(avg_all,[],2),[],1);
% ymin = min(min(avg_all,[],2),[],1);
% for ic = 1:nCells
%     subplot(n,n2,ic)
% %     errorbar(tt,avg_NR(ic,:), sem_NR(ic,:),'k');hold on;
% %     errorbar(tt,avg_OR(ic,:), sem_OR(ic,:),'r');
% %     errorbar(tt,avg_UR(ic,:), sem_UR(ic,:),'b');
%     plot(tt,avg_NR(ic,:),'k');hold on;
%     if ~isempty(OR_movie)
%         plot(tt,avg_OR(ic,:),'r');
%     end
%     if ~isempty(UR_movie)
%         plot(tt,avg_UR(ic,:),'g');
%     end
%     ylim([ymin*1.1 ymax*1.1])
%     xlim([tt(1) tt(end)])
% end
% suptitle([date ' ' mouse ' Average cue resp: Normal- black (n = ' num2str(size(NR_movie,1)) ' trials); Omit- red (n = ' num2str(size(OR_movie,1)) ' trials); Unexpect- green (n = ' num2str(size(UR_movie,1)) ' trials)'])
% orient landscape
% print([dest '_allTCs.eps'], '-depsc');
% print([dest '_allTCs.pdf'], '-dpdf');
% 
% %overlay average success/failure across ROIs
% h=figure; 
% subplot(3,3,1)
% errorbar(tt,avg_NR_all, sem_NR_all,'k')
% hold on;
% if ~isempty(OR_movie)
%     errorbar(tt,avg_OR_all, sem_OR_all,'r');
% end
% if ~isempty(UR_movie)
%     errorbar(tt,avg_UR_all, sem_UR_all,'g')
% end
% % bar(tt,mean(lick_trace_NR)/5, 'k');
% % bar(tt,mean(lick_trace_OR)/5, 'r');
% % if ~isempty(UR_movie)
% %     bar(tt,mean(lick_trace_UR), 'g-');
% % end
% title(['Avg all cells; n = ' num2str(size(avg_NR,1))])
% ylabel('dF/F')
% xlabel('Time (s)')
% xlim([-floor(pre_cue_frames/ifi) floor(post_cue_frames/ifi)]);
% 
% % plot NR resp cells
% subplot(3,3,2)
% errorbar(tt,mean(avg_NR(NR_h,:),1),std(avg_NR(NR_h,:),[],1)./sqrt(sum(NR_h,2)),'k')
% hold on;
% 
% bar(tt,mean(lick_trace_NR)/8, 'k', 'ShowBaseLine', 'off');
% 
% title(['Avg Normal Reward resp cells p<0.05; n = ' num2str(sum(NR_h,2))])
% ylabel('dF/F')
% xlabel('Time (s)')
% xlim([-floor(pre_cue_frames/ifi) floor(post_cue_frames/ifi)]);
% 
% % plot OR resp cells
% if ~isempty(OR_movie)
% subplot(3,3,3)
% errorbar(tt,mean(avg_OR(OR_h,:),1),std(avg_OR(OR_h,:),[],1)./sqrt(sum(OR_h,2)),'r')
% hold on;
% 
% bar(tt,mean(lick_trace_OR)/8, 'FaceColor', [1,0,0], 'EdgeColor', [1,0,0], 'ShowBaseLine', 'off');
% 
% title(['Avg Omitted Reward resp cells p<0.05; n = ' num2str(sum(OR_h,2))])
% ylabel('dF/F')
% xlabel('Time (s)')
% xlim([-floor(pre_cue_frames/ifi) floor(post_cue_frames/ifi)]);
% end
% 
% % plot UR resp cells
% if ~isempty(UR_movie)
% subplot(3,3,4)
% errorbar(tt,mean(avg_UR(UR_h,:),1),std(avg_UR(UR_h,:),[],1)./sqrt(sum(UR_h,2)),'g')
% hold on;
% 
% bar(tt,mean(lick_trace_UR)/8, 'FaceColor', [0,1,0], 'EdgeColor', [0,1,0], 'ShowBaseLine', 'off');
% 
% title(['Avg Unexpected Reward resp cells p<0.05; n = ' num2str(sum(UR_h,2))])
% ylabel('dF/F')
% xlabel('Time (s)')
% xlim([-floor(pre_cue_frames/ifi) floor(post_cue_frames/ifi)]);
% end
% 
% % plot NR and OR resp cells
% if ~isempty(OR_movie)
% subplot(3,3,5)
% errorbar(tt,mean(avg_NR(NOR_h,:),1),std(avg_NR(NOR_h,:),[],1)./sqrt(sum(NOR_h,2)),'k')
% hold on;
% errorbar(tt,mean(avg_OR(NOR_h,:),1),std(avg_OR(NOR_h,:),[],1)./sqrt(sum(NOR_h,2)),'r')
% % bar(tt,mean(lick_trace_NR)/8, 'k');
% 
% title(['Avg Normal and Omitted Reward resp cells p<0.05; n = ' num2str(sum(NOR_h,2))])
% ylabel('dF/F')
% xlabel('Time (s)')
% xlim([-floor(pre_cue_frames/ifi) floor(post_cue_frames/ifi)]);
% end
% 
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
% 
% subplot(3,3,8)
% errorbar(tt,mean(avg_NR(find(allresp_cells),:),1),std(avg_NR(find(allresp_cells),:),[],1)./sqrt(sum(allresp_cells,2)),'k')
% hold on;
% errorbar(tt,mean(avg_OR(find(allresp_cells),:),1),std(avg_OR(find(allresp_cells),:),[],1)./sqrt(sum(allresp_cells,2)),'r')
% errorbar(tt,mean(avg_UR(find(allresp_cells),:),1),std(avg_UR(find(allresp_cells),:),[],1)./sqrt(sum(allresp_cells,2)),'g')
% 
% bar(tt,mean(lick_trace_NR)/8, 'k', 'ShowBaseLine', 'off');
% if ~isempty(OR_movie)
% bar(tt,mean(lick_trace_OR)/8, 'FaceColor', [1,0,0], 'EdgeColor', [1,0,0], 'ShowBaseLine', 'off');
% end
% 
% if ~isempty(UR_movie)
%     bar(tt,mean(lick_trace_UR/8), 'FaceColor', [0,1,0], 'EdgeColor', [0,1,0], 'ShowBaseLine', 'off');
% end
% 
% title(['Avg all resp cells p<0.05; n = ' num2str(sum(allresp_cells,2))])
% ylabel('dF/F')
% xlabel('Time (s)')
% xlim([-floor(pre_cue_frames/ifi) floor(post_cue_frames/ifi)]);
% 
% suptitle([date ' ' mouse ' Average cue resp: Normal- black; Omit- red; Unexpect- green bars- licking data'])
% print([dest '_avgTCs.eps'], '-depsc');
% print([dest '_avgTCs.pdf'], '-dpdf');
% savefig(h, [dest '_avgTCs']); 
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
% suptitle([mouse ' ' date ' cell responses']);
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
% suptitle([date ' ' mouse 'avg resp amplitude'])
% print([dest '_cell_resp_amplitude.eps'], '-depsc');
% print([dest '_cell_resp_amplitude.pdf'], '-dpdf');


