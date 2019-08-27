% PSTH plotting function inputs:
% histgrams for both conditions on all four imaging days in cell format
%    1) no-lick vs all trials
%    2) all cells, responsive cells, resp subtypes, etc
%    3) save data or not

function output = PSTH_summary_plotter(hist_cell_d1_NR, hist_cell_d1_OR, ...
    hist_cell_dN_NR, hist_cell_dN_OR, ...
    hist_cell_dUR_NR, hist_cell_dUR_UR, ...
    hist_cell_d1000_NR, hist_cell_d1000_OR, ...
    tth, lick_cond, cell_subset, save_switch);

%make matrices from cell histograms
hist_mat_d1_NR = cell2mat(hist_cell_d1_NR)';
hist_mat_d1_OR = cell2mat(hist_cell_d1_OR)';
%dayN
hist_mat_dN_NR = cell2mat(hist_cell_dN_NR)';
hist_mat_dN_OR = cell2mat(hist_cell_dN_OR)';
% day UR
hist_mat_dUR_NR = cell2mat(hist_cell_dUR_NR)';
hist_mat_dUR_UR = cell2mat(hist_cell_dUR_UR)';
%1000ms delay
hist_mat_d1000_NR = cell2mat(hist_cell_d1000_NR)';
hist_mat_d1000_OR = cell2mat(hist_cell_d1000_OR)';

%Determine titles and axes for each subplot
[~,cell_countNR] = cellfun(@size, hist_cell_d1_NR);         
[~,cell_countOR] = cellfun(@size, hist_cell_d1_OR); 
subplot_1_title = ['Day1: ', num2str(size(hist_mat_d1_NR,1)), 'neurons. NR/OR: ', num2str(sum(~cellfun(@isempty, hist_cell_d1_NR))), ' sessions'];
[~,cell_countNR] = cellfun(@size, hist_cell_dN_NR);         
[~,cell_countOR] = cellfun(@size, hist_cell_dN_OR); 
subplot_2_title = ['PL: ', num2str(size(hist_mat_dN_NR,1)), 'neurons. NR/OR: ', num2str(sum(~cellfun(@isempty, hist_cell_dN_NR))), ' sessions'];
[~,cell_countNR] = cellfun(@size, hist_cell_dUR_NR);         
[~,cell_countOR] = cellfun(@size, hist_cell_dUR_UR); 
subplot_3_title = ['Day1: ', num2str(size(hist_mat_dUR_NR,1)), 'neurons. NR/UR: ', num2str(sum(~cellfun(@isempty, hist_cell_dUR_NR))), ' sessions'];
[~,cell_countNR] = cellfun(@size, hist_cell_d1000_NR);         
[~,cell_countOR] = cellfun(@size, hist_cell_d1000_OR); 
subplot_4_title = ['1000ms: ', num2str(size(hist_mat_d1000_NR,1)), 'neurons. NR/OR: ', num2str(sum(~cellfun(@isempty, hist_cell_d1000_NR))), ' sessions'];

%determine xlim
x_axis = [-1000 1950];
all_ylim = [];

output=figure('rend', 'painters', 'pos', [50 150 1200 550]);
%SUBPLOT 1
subplot(2,2,1)
shadedErrorBar(tth, nanmean(hist_mat_d1_OR,1), nanstd(hist_mat_d1_OR,[],1)./sqrt(size(hist_mat_d1_OR,1)), 'r'); hold on;
shadedErrorBar(tth, nanmean(hist_mat_d1_NR,1), nanstd(hist_mat_d1_NR,[],1)./sqrt(size(hist_mat_d1_NR,1)), 'k'); 
vline(0,'c');       
vline(580, 'b'); 
if lick_cond ==1
    vline(500,'k')
end
xlim(x_axis);
xlabel('Time (ms)');       
ylabel('Firing rate (Hz)');
title(subplot_1_title);
all_ylim =[all_ylim; ylim];

%SUBPLOT 2
subplot(2,2,2)
shadedErrorBar(tth, nanmean(hist_mat_dN_OR,1), nanstd(hist_mat_dN_OR,[],1)./sqrt(size(hist_mat_dN_OR,1)), 'r'); hold on;
shadedErrorBar(tth, nanmean(hist_mat_dN_NR,1), nanstd(hist_mat_dN_NR,[],1)./sqrt(size(hist_mat_dN_NR,1)), 'k');  
vline(0,'c'); 
vline(580, 'b'); 
if lick_cond ==1
    vline(500,'k')
end
xlim(x_axis);
xlabel('Time (ms)');
ylabel('Firing rate (Hz)');
title(subplot_2_title);
all_ylim =[all_ylim; ylim];

%SUBPLOT 3
subplot(2,2,3)
shadedErrorBar(tth, nanmean(hist_mat_dUR_UR,1), nanstd(hist_mat_dUR_UR,[],1)./sqrt(size(hist_mat_dUR_UR,1)), 'g'); hold on;
shadedErrorBar(tth, nanmean(hist_mat_dUR_NR,1), nanstd(hist_mat_dUR_NR,[],1)./sqrt(size(hist_mat_dUR_NR,1)), 'k');   
vline(0, '--c');      
vline(580,'b');
if lick_cond ==1
    vline(500,'k')
end
xlim(x_axis);  
xlabel('Time (ms)')
ylabel('Firing rate (Hz)')
title(subplot_3_title);
all_ylim =[all_ylim; ylim];

%SUBPLOT 4
subplot(2,2,4);
shadedErrorBar(tth, nanmean(hist_mat_d1000_OR,1), nanstd(hist_mat_d1000_OR,[],1)./sqrt(size(hist_mat_d1000_OR,1)), 'r');hold on;
shadedErrorBar(tth, nanmean(hist_mat_d1000_NR,1), nanstd(hist_mat_d1000_NR,[],1)./sqrt(size(hist_mat_d1000_NR,1)), 'k'); 
xlim(x_axis);
ylim([0 1.5]); 
vline(0, 'c'); 
vline(580,'--b'); 
vline(1080, 'b');
if lick_cond ==1
    vline(500,'k')
end
xlabel('Time (ms)')
ylabel('Firing rate (Hz)')
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
suptitle([cell_subset, ' neurons: ', sup_sec_2, ':']);

%if plotting ITI lick/bout responsive neurons then do not plot UR or 1000ms days 
if strcmp(cell_subset, 'ITI_lick_resp_cells') | strcmp(cell_subset, 'ITI_bout_resp_cells'); 
    subplot(2,2,3); plot(1,1);
    subplot(2,2,4); plot(1,1);
end

%save figure
if ~isempty(save_switch)
    savefig([save_switch, cell_subset, '_cells_', sup_sec_2, '_', num2str(length(hist_cell_d1_NR)), '_sessions_ex_trials'])
end

return