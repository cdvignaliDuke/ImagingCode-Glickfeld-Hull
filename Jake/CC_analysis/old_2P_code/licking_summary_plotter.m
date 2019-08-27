% licking_summary_plotter
%script for plotting licking PSTHs aligned to cue

function output_fig = licking_summary_plotter(normalR_licks_D1, omitR_licks_D1, ...
    normalR_licks_P1, omitR_licks_P1, ...
    normalR_licks_PN1, unexpR_licks_PN1, ...
    normalR_licks_PN2, omitR_licks_PN2, ...
    tth_lick, lick_cond, save_switch);

%make matrices from cell histograms
hist_mat_d1_NR = cell2mat(transpose(normalR_licks_D1)).*30;
hist_mat_d1_OR = cell2mat(transpose(omitR_licks_D1)).*30;
%dayN
hist_mat_dN_NR = cell2mat(transpose(normalR_licks_P1)).*30;
hist_mat_dN_OR = cell2mat(transpose(omitR_licks_P1)).*30;
% day UR
hist_mat_dUR_NR = cell2mat(transpose(normalR_licks_PN1)).*30;
hist_mat_dUR_UR = cell2mat(transpose(unexpR_licks_PN1)).*30;
%1000ms delay
hist_mat_d1000_NR = cell2mat(transpose(normalR_licks_PN2)).*30;
hist_mat_d1000_OR = cell2mat(transpose(omitR_licks_PN2)).*30;
    
%Determine titles and axes for each subplot
[cell_countNR,~] = cellfun(@size, normalR_licks_D1);         
[cell_countOR,~] = cellfun(@size, omitR_licks_D1); 
subplot_1_title = ['Day1: ', num2str(size(hist_mat_d1_NR,1)), ' trials. NR: ', num2str(length(cell_countNR>1)), 'sessions.  OR: ', num2str(length(cell_countOR>1)), ' sessions'];
[cell_countNR,~] = cellfun(@size, normalR_licks_P1);         
[cell_countOR,~] = cellfun(@size, omitR_licks_P1); 
subplot_2_title = ['PL: ', num2str(size(hist_mat_dN_NR,1)), ' trials. NR: ', num2str(length(cell_countNR>1)), 'sessions.  OR: ', num2str(length(cell_countOR>1)), ' sessions'];
[cell_countNR,~] = cellfun(@size, normalR_licks_PN1);         
[cell_countOR,~] = cellfun(@size, unexpR_licks_PN1); 
subplot_3_title = ['Day1: ', num2str(size(hist_mat_dUR_NR,1)), ' trial. NR: ', num2str(length(cell_countNR>1)), 'sessions.  OR: ', num2str(length(cell_countOR>1)), ' sessions'];
[cell_countNR,~] = cellfun(@size, normalR_licks_PN2);         
[cell_countOR,~] = cellfun(@size, omitR_licks_PN2); 
subplot_4_title = ['1000ms: ', num2str(size(hist_mat_d1000_NR,1)), ' trials. NR: ', num2str(length(cell_countNR>1)), 'sessions.  OR: ', num2str(length(cell_countOR>1)), ' sessions'];

all_ylim = [];
x_axis = [-1000 2500];

%% plot the mean across all trials
output_fig=figure('rend', 'painters', 'pos', [50 150 1200 550]);
%SUBPLOT 1
subplot(2,2,1)
shadedErrorBar(tth_lick, nanmean(hist_mat_d1_NR,1), nanstd(hist_mat_d1_NR,[],1)./sqrt(size(hist_mat_d1_NR,1)), 'k'); hold on;
shadedErrorBar(tth_lick, nanmean(hist_mat_d1_OR,1), nanstd(hist_mat_d1_OR,[],1)./sqrt(size(hist_mat_d1_OR,1)), 'r');
vline(0,'c');       
vline(600, 'b'); 
if lick_cond ==1
    vline(500,'k')
end
xlim(x_axis);
xlabel('Time (ms)');       
ylabel('Lick rate (Hz)');
title(subplot_1_title);
all_ylim =[all_ylim; ylim];

%SUBPLOT 2
subplot(2,2,2)
shadedErrorBar(tth_lick, nanmean(hist_mat_dN_NR,1), nanstd(hist_mat_dN_NR,[],1)./sqrt(size(hist_mat_dN_NR,1)), 'k');  hold on;
shadedErrorBar(tth_lick, nanmean(hist_mat_dN_OR,1), nanstd(hist_mat_dN_OR,[],1)./sqrt(size(hist_mat_dN_OR,1)), 'r');
vline(0,'c'); 
vline(600, 'b'); 
if lick_cond ==1
    vline(500,'k')
end
xlim(x_axis);
xlabel('Time (ms)');
ylabel('Lick rate (Hz)');
title(subplot_2_title);
all_ylim =[all_ylim; ylim];

%SUBPLOT 3
subplot(2,2,3)
shadedErrorBar(tth_lick, nanmean(hist_mat_dUR_NR,1), nanstd(hist_mat_dUR_NR,[],1)./sqrt(size(hist_mat_dUR_NR,1)), 'k');    hold on;
shadedErrorBar(tth_lick, nanmean(hist_mat_dUR_UR,1), nanstd(hist_mat_dUR_UR,[],1)./sqrt(size(hist_mat_dUR_UR,1)), 'g'); 
vline(0, '--c');      
vline(600,'b');
if lick_cond ==1
    vline(500,'k')
end
xlim(x_axis);  
xlabel('Time (ms)')
ylabel('Lick rate (Hz)')
title(subplot_3_title);
all_ylim =[all_ylim; ylim];

%SUBPLOT 4
subplot(2,2,4)
shadedErrorBar(tth_lick, nanmean(hist_mat_d1000_NR,1), nanstd(hist_mat_d1000_NR,[],1)./sqrt(size(hist_mat_d1000_NR,1)), 'k'); hold on;
shadedErrorBar(tth_lick, nanmean(hist_mat_d1000_OR,1), nanstd(hist_mat_d1000_OR,[],1)./sqrt(size(hist_mat_d1000_OR,1)), 'r');
xlim(x_axis);
vline(0, 'c'); 
vline(600,'--b'); 
vline(1100, 'b');
if lick_cond ==1
    vline(500,'k')
end
xlabel('Time (ms)')
ylabel('Lick rate (Hz)')
title(subplot_4_title); 
all_ylim =[all_ylim; ylim];

%determine x and y lims
all_ylim_min = min(min(all_ylim));
all_ylim_max = max(max(all_ylim));
subplot(2,2,1); ylim([all_ylim_min all_ylim_max]);
subplot(2,2,2); ylim([all_ylim_min all_ylim_max]);
subplot(2,2,3); ylim([all_ylim_min all_ylim_max]);
subplot(2,2,4); ylim([all_ylim_min all_ylim_max]);

%determine suptitle    lick_cond, cell_subset
if strcmp(lick_cond, 'all')
    sup_sec_2 = 'all trials';
elseif strcmp(lick_cond, 'no-lick')
    sup_sec_2 = 'no-lick condition';
end
suptitle(['lick rate ' sup_sec_2, ':']);
%save figure
if ~isempty(save_switch)
    savefig([save_switch, 'all days ', lick_cond, ' trails'])
end

%% plot the mean across trials for each session
%plot Day 1 NR and OR
color_leg = {'r', 'k', 'b', 'm', 'c', 'g', '--r', '--k', '--b', '--m', '--c', '--g', ':r', ':k', ':b', ':m', ':c', ':g'};
color_leg = {'b', 'b', 'b', 'b', 'b', 'b', 'k', 'k', 'k', 'k', 'k', 'k', 'k', 'k', 'k', 'k', 'k', 'k'};
figure('rend', 'painters', 'pos', [50 150 1200 550]);
for sess_num = 1:length(normalR_licks_D1)
    subplot(2,1,1); hold on;
    plot(tth_lick, nanmean(normalR_licks_D1{sess_num})*30, color_leg{sess_num});
    subplot(2,1,2); hold on;
    plot(tth_lick, nanmean(omitR_licks_D1{sess_num})*30, color_leg{sess_num});
end
subplot(2,1,1); title('normal reward trials.'); ylabel('lick rate (Hz)'); xlabel('time (ms)');
vline(0, 'c'); vline(600, 'b'); ylim([0 15]); xlim([-500 1500]);
subplot(2,1,2); title('omitted reward trials.'); ylabel('lick rate (Hz)'); xlabel('time (ms)');
vline(0,'c'); ylim([0 15]); xlim([-500 1500]);
suptitle(['Day 1: Session averages for each animal: ', lick_cond, ' trials']);
%save figure
if ~isempty(save_switch)
    savefig([save_switch, 'day1 ', lick_cond, ' trails ', num2str(sum(~cellfun(@isempty, normalR_licks_D1))), ' sessions'])
end

%POST LEARNING----------------------------------------------------------------
figure('rend', 'painters', 'pos', [50 150 1200 550]);
for sess_num = 1:length(normalR_licks_P1)
    subplot(2,1,1); hold on;
    plot(tth_lick, nanmean(normalR_licks_P1{sess_num})*30, color_leg{sess_num});
    subplot(2,1,2); hold on;
    plot(tth_lick, nanmean(omitR_licks_P1{sess_num})*30, color_leg{sess_num});
end
subplot(2,1,1); title('normal reward trials.'); ylabel('lick rate (Hz)'); xlabel('time (ms)');
vline(0, 'c'); vline(600, 'b'); ylim([0 15]); xlim([-500 1500]);
subplot(2,1,2); title('omitted reward trials.'); ylabel('lick rate (Hz)'); xlabel('time (ms)');
vline(0,'c'); ylim([0 15]); xlim([-500 1500]);
suptitle(['Day N+1: Session averages for each animal: ', lick_cond, ' trials']);
%save figure
if ~isempty(save_switch)
    savefig([save_switch, 'dayN+1 ', lick_cond, ' trails ', num2str(sum(~cellfun(@isempty, normalR_licks_P1))), ' sessions'])
end

% UNEXPECTED REWARD--------------------------------------------------
figure('rend', 'painters', 'pos', [50 150 1200 550]);
for sess_num = 1:length(normalR_licks_PN1)
    subplot(2,1,1); hold on;
    plot(tth_lick, nanmean(normalR_licks_PN1{sess_num})*30, color_leg{sess_num});
    subplot(2,1,2); hold on;
    plot(tth_lick, nanmean(unexpR_licks_PN1{sess_num})*30, color_leg{sess_num});
end
subplot(2,1,1); title('normal reward trials.'); ylabel('lick rate (Hz)'); xlabel('time (ms)');
vline(0, 'c'); vline(600, 'b'); ylim([0 15]); xlim([-500 1500]);
subplot(2,1,2); title('unexpected reward trials.'); ylabel('lick rate (Hz)'); xlabel('time (ms)');
vline(600,'b'); ylim([0 15]); xlim([-500 1500]);
suptitle(['Day N+2: Session averages for each animal: ', lick_cond, ' trials']);
%save figure
if ~isempty(save_switch)
    savefig([save_switch, 'dayN+2 ', lick_cond, ' trails ', num2str(sum(~cellfun(@isempty, normalR_licks_PN1))), ' sessions'])
end

% 1000ms DELAY
figure('rend', 'painters', 'pos', [50 150 1200 550]);
for sess_num = 1:length(normalR_licks_PN2)
    subplot(2,1,1); hold on;
    plot(tth_lick, nanmean(normalR_licks_PN2{sess_num})*30, color_leg{sess_num});
    subplot(2,1,2); hold on;
    plot(tth_lick, nanmean(omitR_licks_PN2{sess_num})*30, color_leg{sess_num});
end
subplot(2,1,1); title('normal reward trials.'); ylabel('lick rate (Hz)'); xlabel('time (ms)');
vline(0, 'c'); vline(600, '--b'); vline(1100, '--b'); ylim([0 15]); xlim([-500 1500]);
subplot(2,1,2); title('omitted reward trials.'); ylabel('lick rate (Hz)'); xlabel('time (ms)');
vline(0,'c'); ylim([0 15]); xlim([-500 1500]);
suptitle(['Day N+3: Session averages for each animal: ', lick_cond, ' trials']);
%save figure
if ~isempty(save_switch)
    savefig([save_switch, 'dayN+3 ', lick_cond, ' trails ', num2str(sum(~cellfun(@isempty, normalR_licks_PN2))), ' sessions'])
end

%% plot each trial for a given session 


    return
    
    