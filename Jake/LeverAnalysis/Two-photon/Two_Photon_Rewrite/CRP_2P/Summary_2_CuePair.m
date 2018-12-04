clear
file_info_CRP_all
out_base = fullfile('Z:\Analysis\Cue_reward_pairing_analysis\2P\PCA_ICA_recalibration\Summary_folder\TC and cell categories summary\');
data_base = 'Z:\Analysis\Cue_reward_pairing_analysis\2P\PCA_ICA_recalibration\';

outD1 = []; day = 1;
days_subset = days_1; 
for sub = 1:size(days_1,2)
    
    %collect session/mouse information
    session = days_subset{sub};
    session_date = days_subset{sub}(1:6);
    if session(end-2) == 'g'
        mouse_num = ['9', session(end-1:end)];
        mouse_ID = ['img', session(end-1:end)];
    elseif session(end-2) =='0'
        mouse_num = session(end-2:end);
        mouse_ID = ['img', session(end-2:end)];
    end
    session
    
    %set pathnames
    data_dir = fullfile('Z:\Data\2P_imaging',session, mouse_ID, '\');
    for rID = 1:2
        if exist([data_dir, mouse_ID, '_000_', runID{rID}, '.sbx'], 'file') == 2
            break
        end
    end
    if strcmp(session, '170426_img91') | strcmp(session, '170522_img89')
        rID=2;
    end
    
    %collect data and store it in a cell
    dest_sub  = fullfile(data_base, [session(1:6), '_', runID{rID}, '_', mouse_ID],'\');
    outD1 = getData2Cell2(dest_sub, outD1, sub, day);
    
    %         if ~isempty(outD1)  && length(outD1.NR_F_onset) == id
    %             outD1.NR_F_onset_RS{id} = outD1.NR_F_onset{id}(outD1.RS_cells{id});
    %             outD1.OR_F_onset_RS{id} = outD1.OR_F_onset{id}(outD1.RS_cells{id});
    %         end
end


% outPast1 = []; day = 5; doCombMask = 0;
% for id = 1:size(days_s2,2)
%     id
%     mouse = days_s2{id};
%     if doCombMask == 1
%         dest_sub  = fullfile('Z:','home','jake','Analysis','Cue_reward_pairing_analysis','2P',days_pair_folder{id},'\');
%         
%         outPast1 = getData2Cell(dest_sub, outPast1, id, day);
%         
%         outPast1.NR_F_onset_RS{id} = outPast1.NR_F_onset{id}(outD1.RS_cells{id});
%         
%         outPast1.OR_F_onset_RS{id} = outPast1.OR_F_onset{id}(outD1.RS_cells{id});
%         
%         outPast1.cue_rew_fr{id} = sum(outPast1.cue_cells{id} & outD1.rew_cells{id}) / sum(outD1.rew_cells{id});
%         
%         outPast1.rew_cue_fr{id} = sum(outPast1.rew_cells{id} & outD1.rew_cells{id}) / sum(outD1.rew_cells{id});
%         
%         % find OR resp in Day5 for omission resp cells
%         D1D5_OR_resp_cells = (outD1.OR_resp_cells{id}) & (outPast1.OR_resp_cells{id});
%         
%         dis_OR_resp_cells = (~D1D5_OR_resp_cells) & (outPast1.OR_resp_cells{id}); % OR resp cells in Day5 but not Day1, set resp to 0 for Day1
%         
%         outPast1.OR_resp_RS{id} = outPast1.OR_resp_all{id}(outPast1.OR_resp_cells{id});
%         
%         outD1.OR_resp_RS{id} = outD1.OR_resp_all{id};
%         outD1.OR_resp_RS{id}(dis_OR_resp_cells) = 0; % remove resp from non-resp omission cells 
%         outD1.OR_resp_RS{id}= outD1.OR_resp_RS{id}(outPast1.OR_resp_cells{id});
% 
%     else
%         
%         for rID  = 1:2
%             
%             dest_sub  = fullfile('Z:','home','jake','Analysis','Cue_reward_pairing_analysis','2P',[mouse(1:6), '_', runID{rID}, '_', mouse(8:end)],'\');
%             
%             
%             outPast1 = getData2Cell(dest_sub, outPast1, id, day);
%             
%             outPast1.cue_rew_fr{id} = sum(outPast1.cue_cells{id} & outD1.rew_cells{id}) / sum(outD1.rew_cells{id});
%             
%             outPast1.rew_cue_fr{id} = sum(outPast1.rew_cells{id} & outD1.rew_cells{id}) / sum(outD1.rew_cells{id});
%         end
%     end
% end

out_daysPost = [];
days_subset = days_1; 
for sub = 1:size(days_post,2)
    
    %collect session/mouse information
    session = days_subset{sub};
    session_date = days_subset{sub}(1:6);
    if session(end-2) == 'g'
        mouse_num = ['9', session(end-1:end)];
        mouse_ID = ['img', session(end-1:end)];
    elseif session(end-2) =='0'
        mouse_num = session(end-2:end);
        mouse_ID = ['img', session(end-2:end)];
    end
    session
    
    %set pathnames
    data_dir = fullfile('Z:\Data\2P_imaging',session, mouse_ID, '\');
    for rID = 1:2
        if exist([data_dir, mouse_ID, '_000_', runID{rID}, '.sbx'], 'file') == 2
            break
        end
    end
    if strcmp(session, '170426_img91') | strcmp(session, '170522_img89')
        rID=2;
    end
    
    %load and store data in cell
    dest_sub  = fullfile(data_base, [session_date, '_', runID{rID}, '_', mouse_ID],'\');
    out_daysPost = getData2Cell2(dest_sub, out_daysPost, sub, 6);
end
1
% outPN2 = []; day = 7; doCombMask = 0;
% for id = 1:size(days_1000,2)
%     mouse = days_1000{id};
%     if doCombMask == 1
%         dest_sub  = fullfile('Z:','home','jake','Analysis','Cue_reward_pairing_analysis','2P',days_pair_folder{id},'\');
%         
%         outPast1 = getData2Cell(dest_sub, outPast1, id, day);
%         
%         outPast1.NR_F_onset_RS{id} = outPast1.NR_F_onset{id}(outD1.RS_cells{id});
%         
%         outPast1.OR_F_onset_RS{id} = outPast1.OR_F_onset{id}(outD1.RS_cells{id});
%         
%         outPast1.cue_rew_fr{id} = sum(outPast1.cue_cells{id} & outD1.rew_cells{id}) / sum(outD1.rew_cells{id});
%         
%         outPast1.rew_cue_fr{id} = sum(outPast1.rew_cells{id} & outD1.rew_cells{id}) / sum(outD1.rew_cells{id});
%         
%         % find OR resp in Day5 for omission resp cells
%         D1D5_OR_resp_cells = (outD1.OR_resp_cells{id}) & (outPast1.OR_resp_cells{id});
%         
%         dis_OR_resp_cells = (~D1D5_OR_resp_cells) & (outPast1.OR_resp_cells{id}); % OR resp cells in Day5 but not Day1, set resp to 0 for Day1
%         
%         outPast1.OR_resp_RS{id} = outPast1.OR_resp_all{id}(outPast1.OR_resp_cells{id});
%         
%         outD1.OR_resp_RS{id} = outD1.OR_resp_all{id};
%         outD1.OR_resp_RS{id}(dis_OR_resp_cells) = 0; % remove resp from non-resp omission cells 
%         outD1.OR_resp_RS{id}= outD1.OR_resp_RS{id}(outPast1.OR_resp_cells{id});
% 
%     else
%         
%         for rID  = 1:2
%             
%             dest_sub  = fullfile('Z:','home','jake','Analysis','Cue_reward_pairing_analysis','2P',[mouse(1:6), '_', runID{rID}, '_', mouse(8:end)],'\');
%             
%             
%             outPN2 = getData2Cell(dest_sub, outPN2, id, day);
%             
%             
%         end
%     end
% end

day1_ncells = sum(cell2mat(outD1.ncells));
day1_respcells = sum(cell2mat(outD1.tot_resp));
day1_rewcells = sum(cell2mat(outD1.rew_cells));
day1_cuecells = sum(cell2mat(outD1.cue_cells));
day5_ncells = sum(cell2mat(outPast1.ncells));
day5_respcells = sum(cell2mat(outPast1.tot_resp));
day5_rewcells = sum(cell2mat(outPast1.rew_cells));
day5_cuecells = sum(cell2mat(outPast1.cue_cells));
day5_OR_resp_cells = sum(cell2mat(outPast1.OR_resp_cells));
% %scatter of response amplitudes
% fig=figure;
% x = [-.05:.01:.2];
% y = x;
% r = [];
% p = [];
% % col_mat = strvcat('r', 'b', 'r', 'b', 'g', 'm', 'c'); %hardcoded
col_mat = [ 0.9  0.9  0;
    1  0  1;
    0  1  1;
    0.5  0  0;
    0  1  0;
    0  0  1;
    1  0.6  1;
    0  0  0;
    1  0.8 0.4
    0  0.5 0.7
    0.5 0.4 0];
% 
col_mat_s = repmat([0.5 0.5 0.5], 11, 1);

fig = figure;
scatter_plot(days_1_mouse, outPast1.cue_rew_fr, outPast1.rew_cue_fr, col_mat);
axis square;
hold on;
x = 0:0.1:1;
plot(x,x,'-k');
xlabel('Day5 RW Cells Fraction of Day1 RW Cells');
ylabel('Day5 Cue Cells Fraction of Day1 RW Cells');
saveas(fig, [out_base 'Summary_day1-5_RWCuecell_Fraction.fig']);
print([out_base 'Summary_day1-5_RWCuecell_Fraction.eps'], '-depsc');
print([out_base 'Summary_day1-5_RWCuecell_Fraction.pdf'], '-dpdf');

fig = figure;
subplot(1,2,1)
scatter_plot(days_1_mouse, outD1.rew_cells_perc, outD1.cue_cells_perc, col_mat);
axis square;
x = 0:0.1:1;
plot(x,x,'-k');
xlabel('Reward Responsive Cells Fraction');
ylabel('Cue Responsive Cells Fraction');
title('Day1')

subplot(1,2,2)
scatter_plot(days_1_mouse, outPast1.rew_cells_perc, outPast1.cue_cells_perc, col_mat);
axis square;
x = 0:0.1:1;
plot(x,x,'-k');
xlabel('Reward Responsive Cells Fraction');
ylabel('Cue Responsive Cells Fraction');
title('Day5')

supertitle('Cell Category')
saveas(fig, [out_base 'Summary_day1-5_resp_cell.fig']);
print([out_base 'Summary_day1-5_resp_cell.eps'], '-depsc');
print([out_base 'Summary_day1-5_resp_cell.pdf'], '-dpdf');


fig = figure;
scatter_plot(days_1_mouse, outD1.NR_F_onset, outPast1.NR_F_onset, [0.5 0.5 0.5]); 
axis square;
axis([0 1500 0 1500]);
x = 0:1:1500;
hold on; plot(x,x,'k');
xlabel('day 1 latency (ms)'); ylabel('day 5 latency (ms)');
% subplot(1,2,2)
% scatter_plot(days_1_mouse, outD1.NR_F_onset_RS, outPast1.NR_F_onset_RS, col_mat_s); 
% axis square;
% axis([0 1000 0 1000]);
% x = 0:0.05:1000;
% hold on; plot(x,x,'-k');
% xlabel('day 1 latency (ms)'); ylabel('day 5 latency (ms)');
title(['Reward Responsive Cells From Day1', 'Resp Cells n = ', num2str(day1_respcells)]);

saveas(fig, [out_base 'Summary_day1-5_Fonset2Cue_resp.fig']);
print([out_base 'Summary_day1-5_Fonset2Cue_resp.eps'], '-depsc');
print([out_base 'Summary_day1-5_Fonset2Cue_resp.pdf'], '-dpdf');

fig = figure;
subplot(1,2,1)
scatter_plot(days_1_mouse, outD1.NR_F_2_bout, outPast1.NR_F_2_bout, col_mat); 
axis square;
axis([-1000 1500 -1000 1500]);
x = -1000:0.01:1500;
hold on; plot(x,x,'-k');
xlabel('day 1 latency (ms)'); ylabel('day 5 latency (ms)');
title('Reward')

subplot(1,2,2)
scatter_plot(days_1_mouse, outD1.OR_F_2_bout, outPast1.OR_F_2_bout, col_mat); 
axis square;
axis([-1000 1500 -1000 1500]);
x = -1000:0.01:1500;
hold on; plot(x,x,'-k');
xlabel('day 1 latency (ms)'); ylabel('day 5 latency (ms)');
title('Reward Omission')

supertitle(['F onset relative to lick bout onset, day1 #NR resp cell=', num2str(day1_NRcells), ' day1 #OR resp cell=', num2str(day1_ORcells), ' day5 #NR resp cell=', num2str(day5_NRcells), ' day5 #OR resp cell=', num2str(day5_ORcells)]);
saveas(fig, [out_base 'Summary_day1-5_Fonset2bout_resp.fig']);
print([out_base 'Summary_day1-5_Fonset2bout_resp.eps'], '-depsc');
print([out_base 'Summary_day1-5_Fonset2bout_resp.pdf'], '-dpdf');

% fig = figure;
% subplot(1,2,1)
% scatter_plot(days_1_mouse, outD1.NR_resp_RS, outPast1.NR_resp_RS, col_mat);
% axis square;
% axis([0 0.4 0 0.4]);
% x = 0:0.01:0.4;
% hold on; plot(x,x,'-k');
% xlabel('day 1 dF/F Peak'); ylabel('day 5 dF/F');
% title('Reward')

fig = figure;
scatter_plot(days_1_mouse, outD1.OR_resp_sc, outPast1.OR_resp_sc, [0.5 0.5 0.5]);
axis square;
axis([0 0.6 0 0.6]);
x = 0:0.1:0.6;
hold on; plot(x,x,'-k');
xlabel('day 1 dF/F Peak Amp'); ylabel('day 5 dF/F Peak Amp');
title(['Reward Omission Resp Peak Amplitude ', 'Day5 omission resp cell n = ', num2str(day1_respcells)])

saveas(fig, [out_base 'Summary_day1-5_FPeak_scatter_omission_resp.fig']);
print([out_base 'Summary_day1-FPeak_scatter_omission_resp.eps'], '-depsc');
print([out_base 'Summary_day1-FPeak_scatter_omission_resp.pdf'], '-dpdf');

fig = figure;
subplot(2,2,1)
[h_NR1, p_NR1] = scatter_plot(days_1_mouse, outD1.NR_plot_ones, outD1.NR_F_onset,  col_mat_s);
hold on;
scatter_plot(days_1_mouse, outPast1.NR_plot_ones, outPast1.NR_F_onset,  col_mat_s);
xlim([0 3]);
xt={''; 'day 1' ; 'day 5' ; ''};
set(gca,'xticklabel',xt);
ylabel('dF/F to cue(ms)');
title(['Reward p=', num2str(p_NR1)]);

subplot(2,2,2)
[h_NR2, p_NR2] = scatter_plot(days_1_mouse, outD1.NR_plot_ones, outD1.NR_F_2_bout,  col_mat_s);
hold on;
scatter_plot(days_1_mouse, outPast1.NR_plot_ones, outPast1.NR_F_2_bout,  col_mat_s);
xlim([0 3]);
set(gca,'xticklabel',xt);
ylabel('dF/F to lick bout(ms)');
title(['Reward p=', num2str(p_NR2)]);

subplot(2,2,3)
[h_OR1, p_OR1] = scatter_plot(days_1_mouse, outD1.OR_plot_ones, outD1.OR_F_onset,  col_mat_s);
hold on;
outPast1.OR_F_onset(2) = [];
outPast1.OR_plot_ones(2) = [];
scatter_plot(days_1_mouse([1,3:end]), outPast1.OR_plot_ones, outPast1.OR_F_onset,  col_mat_s);
xlim([0 3]);
set(gca,'xticklabel',xt);
ylabel('dF/F to cue(ms)');
title(['Reward Omission p=', num2str(p_OR1)]);

subplot(2,2,4)
[h_OR2, p_OR2] = scatter_plot(days_1_mouse, outD1.OR_plot_ones, outD1.OR_F_2_bout,  col_mat_s);
hold on;
outPast1.OR_F_2_bout(2) = [];
scatter_plot(days_1_mouse([1,3:end]), outPast1.OR_plot_ones, outPast1.OR_F_2_bout,  col_mat_s);
xlim([0 3]);
set(gca,'xticklabel',xt);
ylabel('dF/F to lick bout(ms)');
title(['Reward Omission p=', num2str(p_OR2)]);

saveas(fig, [out_base 'Summary_day1-5_FonsetLatency_scatter_cells_resp.fig']);
print([out_base 'Summary_day1-5_FonsetLatency_scatter_cells_resp.eps'], '-depsc');
print([out_base 'Summary_day1-5_FonsetLatency_scatter_cells_resp.pdf'], '-dpdf');

dest_sub = 'Z:\home\jake\Analysis\Cue_reward_pairing_analysis\2P\180108_000_img070\';
load([dest_sub '_cue_movies.mat']);

tth = [-pre_cue_frames:post_cue_frames].*double(min((ifi)));
OR_TC_prevR = cell2mat(out_daysPost.OR_TC_prevR');

OR_TC_prevNR = out_daysPost.OR_TC_prevNR;
OR_TC_prevNR = cell2mat(OR_TC_prevNR(~cellfun('isempty',OR_TC_prevNR))');

fig = figure;
hold on
shadedErrorBar(tth, nanmean(OR_TC_prevR,1), nanstd(OR_TC_prevR,[],1)./sqrt(size(OR_TC_prevR,1)), 'b');
shadedErrorBar(tth, nanmean(OR_TC_prevNR,1)+0.005, nanstd(OR_TC_prevNR,[],1)./sqrt(size(OR_TC_prevNR,1)), 'r');
xlabel('Time (ms)')
ylabel('dff')
% title('Omission after NoReward: red after Reward: blue');
title(['Omission after NoReward: red n = ', num2str(size(OR_TC_prevNR,1)), ' after Reward: blue n = ', num2str(size(OR_TC_prevR,1))]);
saveas(fig, [out_base 'Summary_TC_OR.fig']);
print([out_base 'Summary_TC_OR.eps'], '-depsc');
print([out_base 'Summary_TC_OR.pdf'], '-dpdf');

trialLen_prevR = cell2mat(out_daysPost.trialLen_prevR);
all_resp = cell2mat(out_daysPost.all_resp_RS')';
trialLen_prevR(isnan(all_resp)) = [];
all_resp(isnan(all_resp)) = [];
[LinearCoeff, fit] = polyfit(trialLen_prevR, all_resp, 1);
Corrfit = polyval(LinearCoeff, trialLen_prevR);

fig = figure;
scatter_plot(days_1_mouse, out_daysPost.trialLen_prevR, out_daysPost.all_resp_RS', [0.5 0.5 0.5]);
hold on
plot(trialLen_prevR, Corrfit, 'color', [0.5, 0.5,0.5], 'Linewidth', 1.5);
axis square;

xlabel('Time (s)'); ylabel('dF/F Peak Amp');
title(['Reward Omission Resp Peak Amplitude and Trial Length for trials with previous reward'])

saveas(fig, [out_base 'Summary_PostLearning_dff_trialLength.fig']);
print([out_base 'Summary_PostLearning_dff_trialLength.eps'], '-depsc');
print([out_base 'Summary_PostLearning_dff_trialLength.pdf'], '-dpdf');

% hold on
% plot(x,y,'-k')
% xlim([-.05 .3]);
% ylim([-.05 .3]);
% xlabel('Correct dF/F')
% ylabel('Early dF/F')
% hold on
% vline(0,'--k')
% hline(0,'--k')
% title(['All cells'])


% 
% subplot(3,3,1)
% [h_prall, p_prall] = scatter_plot(mouseID, release_resp_all, press_resp_all, col_mat_s);
% hold on 
% plot(x,y,'-k')
% xlim([-.05 .3]);
% ylim([-.05 .3]);
% xlabel('Release dF/F')
% ylabel('Press dF/F')
% hold on
% vline(0,'--k')
% hline(0,'--k')
% title(['All cells- p = ' num2str(p_prall)])
% 
% subplot(3,3,2)
% [h_prRS, p_prRS]=scatter_plot(mouseID, release_resp_RS, press_resp_RS, col_mat_s);
% hold on 
% plot(x,y,'-k')
% xlim([-.05 .3]);
% ylim([-.05 .3]);
% xlabel('Release dF/F')
% ylabel('Press dF/F')
% 
% hold on
% vline(0,'--k')
% hline(0,'--k')
% title(['Resp cells- p = ' num2str(p_prRS)])
% subplot(3,3,3)
% [h_prRL, p_prRL] = scatter_plot(mouseID, release_resp_RL, press_resp_RL, col_mat_s);
% hold on 
% plot(x,y,'-k')
% xlim([-.05 .3]);
% ylim([-.05 .3]);
% xlabel('Release dF/F')
% ylabel('Press dF/F')
% hold on
% vline(0,'--k')
% hline(0,'--k')
% title(['Release resp cells- p = ' num2str(p_prRL)])
% 
% subplot(3,3,4)
% [h_sfall, p_sfall] = scatter_plot(mouseID, success_resp_all, fail_resp_all, col_mat_s);
% hold on 
% xlim([-.025 .3]);
% ylim([-.005 .25]);
% plot(x,y,'-k')
% xlabel('Success dF/F')
% ylabel('Fail dF/F')
% 
% hold on
% vline(0,'--k')
% hline(0,'--k')
% title(['All cells- p = ' num2str(chop(p_sfall,2))])
% subplot(3,3,5)
% [h_sfRS, p_sfRS] = scatter_plot(mouseID, success_resp_RS, fail_resp_RS, col_mat_s);
% hold on 
% plot(x,y,'-k')
% xlim([-.05 .3]);
% ylim([-.05 .3]);
% xlabel('Success dF/F')
% ylabel('Fail dF/F')
% 
% hold on
% vline(0,'--k')
% hline(0,'--k')
% title(['Resp cells- p = ' num2str(chop(p_sfRS,2))])
% 
% subplot(3,3,6)
% [h_sfRL, p_sfRL] = scatter_plot(mouseID, success_resp_RL, fail_resp_RL, col_mat_s);
% hold on 
% plot(x,y,'-k')
% xlim([-.05 .3]);
% ylim([-.05 .3]);
% xlabel('Success dF/F')
% ylabel('Fail dF/F')
% 
% hold on
% vline(0,'--k')
% hline(0,'--k')
% title(['Release resp cells- p = ' num2str(chop(p_sfRL,2))])
% 
% subplot(3,3,7)
% [h_stall, p_sfall] = scatter_plot(mouseID, success_resp_all, tooFast_resp_all, col_mat_s);
% hold on 
% plot(x,y,'-k')
% xlim([-.05 .3]);
% ylim([-.05 .3]);
% xlabel('Success dF/F')
% ylabel('TooFast dF/F')
% 
% hold on
% vline(0,'--k')
% hline(0,'--k')
% title(['All cells- p = ' num2str(chop(p_sfall,2))])
% 
% subplot(3,3,8)
% [h_sfRS, p_sfRS] = scatter_plot(mouseID, success_resp_RS, tooFast_resp_RS, col_mat_s);
% hold on 
% plot(x,y,'-k')
% xlim([-.05 .3]);
% ylim([-.05 .3]);
% xlabel('Success dF/F')
% ylabel('TooFast dF/F')
% 
% hold on
% vline(0,'--k')
% hline(0,'--k')
% title(['Resp cells- p = ' num2str(chop(p_sfRS,2))])
% 
% subplot(3,3,9)
% [h_stRL, p_stRL] = scatter_plot(mouseID, success_resp_RL, tooFast_resp_RL, col_mat_s);
% hold on 
% plot(x,y,'-k')
% xlim([-.05 .3]);
% ylim([-.05 .3]);
% xlabel('Success dF/F')
% ylabel('TooFast dF/F')
% 
% hold on
% vline(0,'--k')
% hline(0,'--k')
% title(['Release resp cells- p = ' num2str(chop(p_stRL,2))])
% 
% supertitle(['Summary of cell response amplitudes'])
% 
% saveas(fig, [out_base 'Summary_cell_response_amp_scatter.fig']);
% print([out_base 'Summary_cell_response_amp_scatter.eps'], '-depsc');
% print([out_base 'Summary_cell_response_amp_scatter.pdf'], '-dpdf');


% %avg scatter of amplitudes by mouse
% fig=figure;
% x = [-.05:.01:.25];
% y = x;
% % col_mat = strvcat('r', 'b', 'r', 'b', 'g', 'm', 'c');   %HARDCODED
% subplot(3,3,1)
% scatter_plot(mouseID, release_resp_all, press_resp_all, col_mat);
% 
% hold on
% plot(x,y,'-k')
% xlim([-.025 .3]);
% ylim([-.025 .2]);
% xlabel('Release dF/F')
% ylabel('Press dF/F')
% hold on
% vline(0,'--k')
% hline(0,'--k')
% title(['All cells'])
% subplot(3,3,2)
% scatter_plot(mouseID, release_resp_RS, press_resp_RS, col_mat);
% hold on
% plot(x,y,'-k')
% xlim([-.05 .3]);
% ylim([-.05 .3]);
% xlabel('Release dF/F')
% ylabel('Press dF/F')
% hold on
% vline(0,'--k')
% hline(0,'--k')
% title(['Resp cells'])
% subplot(3,3,3)
% [h_prRS, p_prRS]=scatter_plot(mouseID, release_resp_RS, press_resp_RS, col_mat_s);
% hold on 
% plot(x,y,'-k')
% xlim([-.05 .3]);
% ylim([-.05 .3]);
% xlabel('Release dF/F')
% ylabel('Press dF/F')
% title(['Resp cells'])
% % subplot(3,3,3)
% % scatter_plot(mouseID, release_resp_RL, press_resp_RL, col_mat);
% % hold on
% % plot(x,y,'-k')
% % xlim([-.05 .3]);
% % ylim([-.05 .3]);
% % xlabel('Release dF/F')
% % ylabel('Press dF/F')
% % hold on
% % vline(0,'--k')
% % hline(0,'--k')
% % title(['Release resp cells'])
% subplot(3,3,4)
% scatter_plot(mouseID, success_resp_all, fail_resp_all, col_mat);
% hold on
% plot(x,y,'-k')
% xlim([-.05 .3]);
% ylim([-.05 .3]);
% xlabel('Correct dF/F')
% ylabel('Early dF/F')
% hold on
% vline(0,'--k')
% hline(0,'--k')
% title(['All cells'])
% % subplot(3,3,5)
% % scatter_plot(mouseID, success_resp_RS, fail_resp_RS, col_mat);
% % hold on
% % plot(x,y,'-k')
% % xlim([-.05 .3]);
% % ylim([-.05 .3]);
% % xlabel('Success dF/F')
% % ylabel('Fail dF/F')
% % hold on
% % vline(0,'--k')
% % hline(0,'--k')
% % title(['Resp cells'])
% subplot(3,3,5)
% scatter_plot(mouseID, success_resp_RL, fail_resp_RL, col_mat);
% hold on
% plot(x,y,'-k')
% xlim([-.05 .3]);
% ylim([-.05 .3]);
% xlabel('Correct dF/F')
% ylabel('Early dF/F')
% hold on
% vline(0,'--k')
% hline(0,'--k')
% title(['Release resp cells'])
% 
% subplot(3,3,6)
% [h_sfRL, p_sfRL] = scatter_plot(mouseID, success_resp_RL, fail_resp_RL, col_mat_s);
% hold on 
% plot(x,y,'-k')
% xlim([-.05 .3]);
% ylim([-.05 .3]);
% xlabel('Correct dF/F')
% ylabel('Early dF/F')
% hold on
% vline(0,'--k')
% hline(0,'--k')
% title(['Release resp cells'])
% % title(['Release resp cells- p = ' num2str(chop(p_sfRL,2))])
% 
% subplot(3,3,7)
% scatter_plot(mouseID, success_resp_all, tooFast_resp_all, col_mat);
% hold on
% plot(x,y,'-k')
% xlim([-.05 .3]);
% ylim([-.05 .3]);
% xlabel('Correct dF/F')
% ylabel('TooFast Correct dF/F')
% hold on
% vline(0,'--k')
% hline(0,'--k')
% title(['All cells'])
% 
% % subplot(3,3,8)
% % scatter_plot(mouseID, success_resp_RS, tooFast_resp_RS, col_mat);
% % hold on
% % plot(x,y,'-k')
% % xlim([-.05 .3]);
% % ylim([-.05 .3]);
% % xlabel('Success dF/F')
% % ylabel('Toofast Correct dF/F')
% % hold on
% % vline(0,'--k')
% % hline(0,'--k')
% % title(['Resp cells'])
% 
% subplot(3,3,8)
% scatter_plot(mouseID, success_resp_RL, tooFast_resp_RL, col_mat);
% hold on
% plot(x,y,'-k')
% xlim([-.05 .3]);
% ylim([-.05 .3]);
% xlabel('Correct dF/F')
% ylabel('Toofast Correct dF/F')
% hold on
% vline(0,'--k')
% hline(0,'--k')
% title(['Release resp cells'])
% 
% subplot(3,3,9)
% [h_stRL, p_stRL] = scatter_plot(mouseID, success_resp_RL, tooFast_resp_RL, col_mat_s);
% hold on 
% plot(x,y,'-k')
% xlim([-.05 .3]);
% ylim([-.05 .3]);
% xlabel('Correct dF/F')
% ylabel('TooFast Correct dF/F')
% hold on
% vline(0,'--k')
% hline(0,'--k')
% title(['Release resp cells'])
% % title(['Release resp cells- p = ' num2str(chop(p_stRL,2))])
% 
% supertitle(['Summary of cell response amplitudes'])
% saveas(fig, [out_base 'Summary_avg_response_amp_scatter.fig']);
% print([out_base 'Summary_avg_response_amp_scatter.eps'], '-depsc');
% print([out_base 'Summary_avg_response_amp_scatter.pdf'], '-dpdf');

% % plot avg response amp scatter for only fail/success and tooFast/success
% fig = figure;
% subplot(1,2,1)
% scatter_plot(mouseID, success_resp_RS, fail_resp_RS, col_mat);
% hold on
% plot(x,y,'-k')
% xlim([-.05 .3]);
% ylim([-.05 .3]);
% xlabel('Success dF/F')
% ylabel('Fail dF/F')
% hold on
% vline(0,'--k')
% hline(0,'--k')
% title(['Resp cells'])
% subplot(1,2,2)
% scatter_plot(mouseID, success_resp_RL, tooFast_resp_RL, col_mat);
% hold on
% plot(x,y,'-k')
% xlim([-.05 .3]);
% ylim([-.05 .3]);
% xlabel('Success dF/F')
% ylabel('Toofast Correct dF/F')
% hold on
% vline(0,'--k')
% hline(0,'--k')
% title(['Release resp cells'])
% 
% supertitle(['Summary of cell response amplitudes'])
% saveas(fig, [out_base 'Summary_avg_response_amp_scatter_2.fig']);
% print([out_base 'Summary_avg_response_amp_scatter_2.eps'], '-depsc');
% print([out_base 'Summary_avg_response_amp_scatter_2.pdf'], '-dpdf');
% %average timecourse for expts

% fig=figure;
% for id = 1:size(mouseID,2)
%     subplot(6,6,id)
%     tt =((-pre_frames{id}:post_frames{id}).*double(TC_ifi{id}))./1000;
%     shadedErrorBar(tt, success_TC_mean{id},success_TC_sem{id}, 'k');
%     hold on
%     shadedErrorBar(tt, fail_TC_mean{id},fail_TC_sem{id}, 'r');
%     hold on
%     shadedErrorBar(tt, press_TC_mean{id},press_TC_sem{id}, 'c');
%     xlim([-pre_frames{id}.*double(TC_ifi{id})./1000 post_frames{id}.*double(TC_ifi{id})./1000])
%     xlabel('Time (ms)')
%     ylabel('dF/F')
%     title([date{id} ' ' mouseID{id}])
% end
% supertitle(['Summary of all cell timecourses'])
% saveas(fig, [out_base 'Summary_allexptTCs_allcells.fig']);
% print([out_base 'Summary_allexptTCs_allcells.eps'], '-depsc');
% print([out_base 'Summary_allexptTCs_allcells.pdf'], '-dpdf');

% fig=figure;
% for id = 1:size(mouseID,2)
%     subplot(6,6,id)
%     tt =((-pre_frames{id}:post_frames{id}).*double(TC_ifi{id}))./1000;
%     shadedErrorBar(tt, success_TC_RS_mean{id},success_TC_RS_sem{id}, 'k');
%     hold on
%     shadedErrorBar(tt, fail_TC_RS_mean{id},fail_TC_RS_sem{id}, 'r');
%     hold on
%     shadedErrorBar(tt, press_TC_RS_mean{id},press_TC_RS_sem{id}, 'c');
%     xlim([-pre_frames{id}.*double(TC_ifi{id})./1000 post_frames{id}.*double(TC_ifi{id})./1000])
%     xlabel('Time (ms)')
%     ylabel('dF/F')
%     title([date{id} ' ' mouseID{id}])
% end
% supertitle(['Summary of responsive cell timecourses'])
% saveas(fig, [out_base 'Summary_allexptTCs_respcells.fig']);
% print([out_base 'Summary_allexptTCs_respcells.eps'], '-depsc');
% print([out_base 'Summary_allexptTCs_respcells.pdf'], '-dpdf');


%% commented for now until a decision is made on how to average across experiments with different acquisition rates
%averaging across all cells- specific to different acquisition rates

frame_size = cell2mat(cellfun(@size, outD1.NR_TC, 'UniformOutput', 0));
max_frame = max(frame_size(2:2:end));

NR_TC_all = interp_frame(outD1.NR_TC, max_frame);
OR_TC_all    = interp_frame(outD1.OR_TC, max_frame);
% UR_TC_all   = cell2mat(outD1.UR_TC);

NR_TC_all_RS = interp_frame(outD1.NR_TC_RS, max_frame);
OR_TC_all_RS = interp_frame(outD1.OR_TC_RS, max_frame);

% UR_TC_all_RS = interp_frame(outD1.UR_TC_RS, max_frame);

total_cellD1 = sum(cell2mat(outD1.ncells));
total_respD1 = sum(cell2mat(outD1.tot_resp));

tt =(-max(cell2mat(outD1.pre_frames)):max(cell2mat(outD1.post_frames))).*double(min(cell2mat(outD1.ifi)))./1000;

fig = figure;
subplot(2,2,1)
shadedErrorBar(tt, mean(NR_TC_all,1), std(NR_TC_all,[],1)./sqrt(size(NR_TC_all,1)), 'k');
hold on;
shadedErrorBar(tt, mean(OR_TC_all,1), std(OR_TC_all,[],1)./sqrt(size(OR_TC_all,1)), 'r');
% hold on
% shadedErrorBar(tt, mean(UR_TC_all,1), std(UR_TC_all,[],1)./sqrt(size(UR_TC_all,1)), 'g');
title(['Day1 Timecourses of average for all cells- n = ' num2str(total_cellD1)])
xlabel('Time (ms)')
ylabel('dF/F')
axis([-2 2 -0.05 0.15]);

subplot(2,2,3)
shadedErrorBar(tt, mean(NR_TC_all - OR_TC_all,1), std(NR_TC_all - OR_TC_all,[],1)./sqrt(size(NR_TC_all,1)), 'b');
title(['Day1 Subtraction of Rewarded and Omission'])
xlabel('Time (ms)')
ylabel('dF/F')
axis([-2 2 -0.05 0.15]);

% plot day after 1
NR_TC_all = interp_frame(outPast1.NR_TC, max_frame);
OR_TC_all = interp_frame(outPast1.OR_TC, max_frame);
% UR_TC_all_RS = cell2mat(outPast1.UR_TC_RS);

total_cellP1 = sum(cell2mat(outPast1.ncells));

subplot(2,2,2)
shadedErrorBar(tt, mean(NR_TC_all,1), std(NR_TC_all,[],1)./sqrt(size(NR_TC_all,1)), 'k');
hold on;
shadedErrorBar(tt, mean(OR_TC_all,1), std(OR_TC_all,[],1)./sqrt(size(OR_TC_all,1)), 'r');
% hold on
% shadedErrorBar(tt, mean(UR_TC_all,1), std(UR_TC_all,[],1)./sqrt(size(UR_TC_all,1)), 'g');
title(['After Day1 Timecourses of average for all cells- n = ' num2str(total_cellP1)])
xlabel('Time (ms)')
ylabel('dF/F')
axis([-2 2 -0.05 0.15]);

subplot(2,2,4)
shadedErrorBar(tt, mean(NR_TC_all - OR_TC_all,1), std(NR_TC_all - OR_TC_all,[],1)./sqrt(size(NR_TC_all,1)), 'b');
hold on;
title(['After Day1 Subtraction of Rewarded and Omission'])
xlabel('Time (ms)')
ylabel('dF/F')
axis([-2 2 -0.05 0.15]);

supertitle(['All cells across experiments Black- Normal; Red- Omitted; Green- Unexpected']);
saveas(fig, [out_base 'Summary_avgexpt_TCs_all.fig']);
print([out_base 'Summary_avgCellexpt_TCs_all.eps'], '-depsc');
print([out_base 'Summary_avgCellexpt_TCs_all.pdf'], '-dpdf');


fig = figure;
subplot(2,2,1)
shadedErrorBar(tt, mean(NR_TC_all_RS,1), std(NR_TC_all_RS,[],1)./sqrt(size(NR_TC_all_RS,1)), 'k');
hold on;
shadedErrorBar(tt, mean(OR_TC_all_RS,1), std(OR_TC_all_RS,[],1)./sqrt(size(OR_TC_all_RS,1)), 'r');
% hold on
% shadedErrorBar(tt, mean(UR_TC_all_RS,1), std(UR_TC_all_RS,[],1)./sqrt(size(UR_TC_all_RS,1)), 'g');
title(['Day1 Timecourses of average for responsive cells- n = ' num2str(total_respD1)])
xlabel('Time (ms)')
ylabel('dF/F')
axis([-2 2 -0.05 0.15]);

subplot(2,2,3)
shadedErrorBar(tt, mean(NR_TC_all_RS - OR_TC_all_RS,1), std(NR_TC_all_RS - OR_TC_all_RS,[],1)./sqrt(size(NR_TC_all_RS,1)), 'b');
title(['Day1 Subtraction of Rewarded and Omission'])
xlabel('Time (ms)')
ylabel('dF/F')
axis([-2 2 -0.05 0.15]);

% plot day after 1
NR_TC_all_RS = interp_frame(outPast1.NR_TC_RS, max_frame);
OR_TC_all_RS = interp_frame(outPast1.OR_TC_RS, max_frame);
% UR_TC_all_RS = cell2mat(outPast1.UR_TC_RS);

total_respP1 = sum(cell2mat(outPast1.tot_resp));

subplot(2,2,2)
shadedErrorBar(tt, mean(NR_TC_all_RS,1), std(NR_TC_all_RS,[],1)./sqrt(size(NR_TC_all_RS,1)), 'k');
hold on;
shadedErrorBar(tt, mean(OR_TC_all_RS,1), std(OR_TC_all_RS,[],1)./sqrt(size(OR_TC_all_RS,1)), 'r');
% hold on
% shadedErrorBar(tt, mean(UR_TC_all_RS,1), std(UR_TC_all_RS,[],1)./sqrt(size(UR_TC_all_RS,1)), 'g');
title(['After Day1 Timecourses of average for responsive cells- n = ' num2str(total_respP1)])
xlabel('Time (ms)')
ylabel('dF/F')
axis([-2 2 -0.05 0.15]);

subplot(2,2,4)
shadedErrorBar(tt, mean(NR_TC_all_RS - OR_TC_all_RS,1), std(NR_TC_all_RS - OR_TC_all_RS,[],1)./sqrt(size(NR_TC_all_RS,1)), 'b');
hold on;
title(['After Day1 Subtraction of Rewarded and Omission'])
xlabel('Time (ms)')
ylabel('dF/F')
axis([-2 2 -0.05 0.15]);

supertitle(['All responsive cells across experiments Black- Normal; Red- Omitted; Green- Unexpected']);
saveas(fig, [out_base 'Summary_avgexpt_TCs_resp.fig']);
print([out_base 'Summary_avgCellexpt_TCs_resp.eps'], '-depsc');
print([out_base 'Summary_avgCellexpt_TCs_resp.pdf'], '-dpdf');

NR_TC_RewPos_P1 = interp_frame(outPast1.NR_TC_RewPos, max_frame);
NR_TC_RewNeg_P1 = interp_frame(outPast1.NR_TC_RewNeg, max_frame);

NR_TC_NL_RewPos_P1 = interp_frame(outPast1.NR_TC_nolick_RewPos, max_frame);
NR_TC_NL_RewNeg_P1 = interp_frame(outPast1.NR_TC_nolick_RewNeg, max_frame);

OR_TC_RewPos_P1 = interp_frame(outPast1.OR_TC_RewPos, max_frame);
OR_TC_RewNeg_P1 = interp_frame(outPast1.OR_TC_RewNeg, max_frame);

OR_TC_NL_RewPos_P1 = interp_frame(outPast1.OR_TC_nolick_RewPos, max_frame);
OR_TC_NL_RewNeg_P1 = interp_frame(outPast1.OR_TC_nolick_RewNeg, max_frame);

NR_TC_RewPos_PN2 = interp_frame(outPN2.NR_TC_RewPos, max_frame);
NR_TC_RewNeg_PN2 = interp_frame(outPN2.NR_TC_RewNeg, max_frame);

NR_TC_NL_RewPos_PN2 = interp_frame(outPN2.NR_TC_nolick_RewPos, max_frame);
NR_TC_NL_RewNeg_PN2 = interp_frame(outPN2.NR_TC_nolick_RewNeg, max_frame);

OR_TC_RewPos_PN2 = interp_frame(outPN2.OR_TC_RewPos, max_frame);
OR_TC_RewNeg_PN2 = interp_frame(outPN2.OR_TC_RewNeg, max_frame);

OR_TC_NL_RewPos_PN2 = interp_frame(outPN2.OR_TC_nolick_RewPos, max_frame);
OR_TC_NL_RewNeg_PN2 = interp_frame(outPN2.OR_TC_nolick_RewNeg, max_frame);

NR_TC_ORcellNeg_P1 = interp_frame(outPast1.NR_TC_ORcellNeg, max_frame);
NR_TC_ORcellNeg_PN2 = interp_frame(outPN2.NR_TC_ORcellNeg, max_frame);

NR_TC_NL_ORcellNeg_P1 = interp_frame(outPast1.NR_TC_nolick_ORcellNeg, max_frame);
NR_TC_NL_ORcellNeg_PN2 = interp_frame(outPN2.NR_TC_nolick_ORcellNeg, max_frame);

fig = figure;
subplot(2,2,1)
shadedErrorBar(tt, mean(NR_TC_RewPos_P1,1), std(NR_TC_RewPos_P1,[],1)./sqrt(size(NR_TC_RewPos_P1,1)), 'g');
hold on;
shadedErrorBar(tt, mean(NR_TC_RewNeg_P1,1), std(NR_TC_RewNeg_P1,[],1)./sqrt(size(NR_TC_RewNeg_P1,1)), 'r');
shadedErrorBar(tt, mean(NR_TC_ORcellNeg_P1,1), std(NR_TC_ORcellNeg_P1,[],1)./sqrt(size(NR_TC_ORcellNeg_P1,1)), 'm');
% hold on
% shadedErrorBar(tt, mean(UR_TC_all_RS,1), std(UR_TC_all_RS,[],1)./sqrt(size(UR_TC_all_RS,1)), 'g');
title(['DayN Timecourses of average Normal Reward'])
xlabel('Time (ms)')
ylabel('dF/F')
axis([-2 2 -0.15 0.25]);

subplot(2,2,2)
shadedErrorBar(tt, mean(OR_TC_RewPos_P1,1), std(OR_TC_RewPos_P1,[],1)./sqrt(size(OR_TC_RewPos_P1,1)), 'g');
hold on;
shadedErrorBar(tt, mean(OR_TC_RewNeg_P1,1), std(OR_TC_RewNeg_P1,[],1)./sqrt(size(OR_TC_RewNeg_P1,1)), 'r');
% hold on
% shadedErrorBar(tt, mean(UR_TC_all_RS,1), std(UR_TC_all_RS,[],1)./sqrt(size(UR_TC_all_RS,1)), 'g');
title(['DayN Timecourses of average Reward Omission'])
xlabel('Time (ms)')
ylabel('dF/F')
axis([-2 2 -0.15 0.25]);

subplot(2,2,3)
shadedErrorBar(tt, mean(NR_TC_RewPos_PN2,1), std(NR_TC_RewPos_PN2,[],1)./sqrt(size(NR_TC_RewPos_PN2,1)), 'g');
hold on;
shadedErrorBar(tt, mean(NR_TC_RewNeg_PN2,1), std(NR_TC_RewNeg_PN2,[],1)./sqrt(size(NR_TC_RewNeg_PN2,1)), 'r');
shadedErrorBar(tt, mean(NR_TC_ORcellNeg_PN2,1), std(NR_TC_ORcellNeg_PN2,[],1)./sqrt(size(NR_TC_ORcellNeg_PN2,1)), 'm');
% hold on
% shadedErrorBar(tt, mean(UR_TC_all_RS,1), std(UR_TC_all_RS,[],1)./sqrt(size(UR_TC_all_RS,1)), 'g');
title(['1000ms Delay Timecourses of average Normal Reward'])
xlabel('Time (ms)')
ylabel('dF/F')
axis([-2 2 -0.1 0.15]);

subplot(2,2,4)
shadedErrorBar(tt, mean(OR_TC_RewPos_PN2,1), std(OR_TC_RewPos_PN2,[],1)./sqrt(size(OR_TC_RewPos_PN2,1)), 'g');
hold on;
shadedErrorBar(tt, mean(OR_TC_RewNeg_PN2,1), std(OR_TC_RewNeg_PN2,[],1)./sqrt(size(OR_TC_RewNeg_PN2,1)), 'r');
% hold on
% shadedErrorBar(tt, mean(UR_TC_all_RS,1), std(UR_TC_all_RS,[],1)./sqrt(size(UR_TC_all_RS,1)), 'g');
title(['1000ms Delay Timecourses of average Reward Omission'])
xlabel('Time (ms)')
ylabel('dF/F')
axis([-2 2 -0.1 0.15]);

supertitle('Time Course, Green: positive; Red: Negative; Magneta: Omission Negative');
saveas(fig, [out_base 'Summary_avgexpt_TCs_withlick_RewPN.fig']);
print([out_base 'Summary_avgCellexpt_TCs_withlick_RewPN.eps'], '-depsc');
print([out_base 'Summary_avgCellexpt_TCs_withlick_RewPN.pdf'], '-dpdf');


fig = figure;
subplot(2,2,1)
shadedErrorBar(tt, mean(NR_TC_NL_RewPos_P1,1), std(NR_TC_NL_RewPos_P1,[],1)./sqrt(size(NR_TC_NL_RewPos_P1,1)), 'g');
hold on;
shadedErrorBar(tt, mean(NR_TC_NL_RewNeg_P1,1), std(NR_TC_NL_RewNeg_P1,[],1)./sqrt(size(NR_TC_NL_RewNeg_P1,1)), 'r');
shadedErrorBar(tt, mean(NR_TC_NL_ORcellNeg_P1,1), std(NR_TC_NL_ORcellNeg_P1,[],1)./sqrt(size(NR_TC_NL_ORcellNeg_P1,1)), 'm');
% hold on
% shadedErrorBar(tt, mean(UR_TC_all_RS,1), std(UR_TC_all_RS,[],1)./sqrt(size(UR_TC_all_RS,1)), 'g');
title(['DayN Timecourses of average Normal Reward'])
xlabel('Time (ms)')
ylabel('dF/F')
axis([-2 2 -0.15 0.25]);

subplot(2,2,2)
shadedErrorBar(tt, mean(OR_TC_NL_RewPos_P1,1), std(OR_TC_NL_RewPos_P1,[],1)./sqrt(size(OR_TC_NL_RewPos_P1,1)), 'g');
hold on;
shadedErrorBar(tt, mean(OR_TC_NL_RewNeg_P1,1), std(OR_TC_NL_RewNeg_P1,[],1)./sqrt(size(OR_TC_NL_RewNeg_P1,1)), 'r');
% hold on
% shadedErrorBar(tt, mean(UR_TC_all_RS,1), std(UR_TC_all_RS,[],1)./sqrt(size(UR_TC_all_RS,1)), 'g');
title(['DayN Timecourses of average Reward Omission'])
xlabel('Time (ms)')
ylabel('dF/F')
axis([-2 2 -0.15 0.25]);

subplot(2,2,3)
shadedErrorBar(tt, mean(NR_TC_NL_RewPos_PN2,1), std(NR_TC_NL_RewPos_PN2,[],1)./sqrt(size(NR_TC_NL_RewPos_PN2,1)), 'g');
hold on;
shadedErrorBar(tt, mean(NR_TC_NL_RewNeg_PN2,1), std(NR_TC_NL_RewNeg_PN2,[],1)./sqrt(size(NR_TC_NL_RewNeg_PN2,1)), 'r');
shadedErrorBar(tt, mean(NR_TC_NL_ORcellNeg_PN2,1), std(NR_TC_NL_ORcellNeg_PN2,[],1)./sqrt(size(NR_TC_NL_ORcellNeg_PN2,1)), 'm');
% hold on
% shadedErrorBar(tt, mean(UR_TC_all_RS,1), std(UR_TC_all_RS,[],1)./sqrt(size(UR_TC_all_RS,1)), 'g');
title(['1000ms Delay Timecourses of average Normal Reward'])
xlabel('Time (ms)')
ylabel('dF/F')
axis([-2 2 -0.1 0.15]);

subplot(2,2,4)
shadedErrorBar(tt, mean(OR_TC_NL_RewPos_PN2,1), std(OR_TC_NL_RewPos_PN2,[],1)./sqrt(size(OR_TC_NL_RewPos_PN2,1)), 'g');
hold on;
shadedErrorBar(tt, mean(OR_TC_NL_RewNeg_PN2,1), std(OR_TC_NL_RewNeg_PN2,[],1)./sqrt(size(OR_TC_NL_RewNeg_PN2,1)), 'r');
% hold on
% shadedErrorBar(tt, mean(UR_TC_all_RS,1), std(UR_TC_all_RS,[],1)./sqrt(size(UR_TC_all_RS,1)), 'g');
title(['1000ms Delay Timecourses of average Reward Omission'])
xlabel('Time (ms)')
ylabel('dF/F')
axis([-2 2 -0.1 0.15]);

supertitle('Time Course, Green: positive; Red: Negative; Magneta: Omission Negative');
saveas(fig, [out_base 'Summary_avgexpt_TCs_nolick_RewPN.fig']);
print([out_base 'Summary_avgCellexpt_TCs_nolick_RewPN.eps'], '-depsc');
print([out_base 'Summary_avgCellexpt_TCs_nolick_RewPN.pdf'], '-dpdf');

%% Avg for each animal first 
frame_size = cell2mat(cellfun(@size, outPast1.NR_TC_RS_mean, 'UniformOutput', 0));
max_frame = max(frame_size(2:2:end));

NR_TC_all_mean_RS = interp_frame(outD1.NR_TC_RS_mean, max_frame);
OR_TC_all_mean_RS = interp_frame(outD1.OR_TC_RS_mean, max_frame);
% UR_TC_all_mean_RS = interp_frame(outD1.UR_TC_RS_mean, max_frame);

fig = figure;
subplot(2,2,1)
shadedErrorBar(tt, mean(NR_TC_all_mean_RS,1), std(NR_TC_all_mean_RS,[],1)./sqrt(size(NR_TC_all_mean_RS,1)), 'k');
hold on;
shadedErrorBar(tt, mean(OR_TC_all_mean_RS,1), std(OR_TC_all_mean_RS,[],1)./sqrt(size(OR_TC_all_mean_RS,1)), 'r');
% hold on
% shadedErrorBar(tt, mean(UR_TC_all_mean_RS,1), std(UR_TC_all_mean_RS,[],1)./sqrt(size(UR_TC_all_mean_RS,1)), 'g');
title(['Day1 Timecourses of average for responsive cells- n = ' num2str(total_respD1)])
xlabel('Time (ms)')
ylabel('dF/F')
axis([-2 2 -0.05 0.15]);

subplot(2,2,3)
shadedErrorBar(tt, mean(NR_TC_all_mean_RS - OR_TC_all_mean_RS,1), std(NR_TC_all_mean_RS - OR_TC_all_mean_RS,[],1)./sqrt(size(NR_TC_all_mean_RS,1)), 'b');
title(['Day1 Subtraction of Rewarded and Omission'])
xlabel('Time (ms)')
ylabel('dF/F')
axis([-2 2 -0.05 0.15]);

NR_TC_all_mean_RS = interp_frame(outPast1.NR_TC_RS_mean, max_frame);
OR_TC_all_mean_RS = interp_frame(outPast1.OR_TC_RS_mean, max_frame);
% UR_TC_all_mean_RS = cell2mat(outPast1.UR_TC_RS_mean);


subplot(2,2,2)
shadedErrorBar(tt, mean(NR_TC_all_mean_RS,1), std(NR_TC_all_mean_RS,[],1)./sqrt(size(NR_TC_all_mean_RS,1)), 'k');
hold on;
shadedErrorBar(tt, mean(OR_TC_all_mean_RS,1), std(OR_TC_all_mean_RS,[],1)./sqrt(size(OR_TC_all_mean_RS,1)), 'r');
% hold on
% shadedErrorBar(tt, mean(UR_TC_all_mean_RS,1), std(UR_TC_all_mean_RS,[],1)./sqrt(size(UR_TC_all_mean_RS,1)), 'g');
title(['After Day1 Timecourses of average for responsive cells- n = ' num2str(total_respP1)])
xlabel('Time (ms)')
ylabel('dF/F')
axis([-2 2 -0.05 0.15]);

subplot(2,2,4)
shadedErrorBar(tt, mean(NR_TC_all_mean_RS - OR_TC_all_mean_RS,1), std(NR_TC_all_mean_RS - OR_TC_all_mean_RS,[],1)./sqrt(size(NR_TC_all_mean_RS,1)), 'b');
hold on;
title(['After Day1 Subtraction of Rewarded and Omission'])
xlabel('Time (ms)')
ylabel('dF/F')
axis([-2 2 -0.05 0.15]);

supertitle(['Animal averaged across experiments Black- Normal; Red- Omitted; Green- Unexpected']);
saveas(fig, [out_base 'Summary_avgAnimalexpt_TCs.fig']);
print([out_base 'Summary_avgAnimalexpt_TCs.eps'], '-depsc');
print([out_base 'Summary_avgAnimalexpt_TCs.pdf'], '-dpdf');

save([out_base 'cell_count.mat'], 'total_cellD1', 'total_respD1', 'total_cellP1', 'total_respP1');