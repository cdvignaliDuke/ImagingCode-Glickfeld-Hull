%plot % of failed trials with a lick vs lick triggered average magnitude 
clear;

%define some variables and pathnames
WF_plotting_lists_of_days;
early_TC_dir = 'Z:\Analysis\WF Lever Analysis\licking_investigation\early_trials_lick_v_no_lick\no_licks_-2_to_2\';
corr_TC_dir = 'Z:\Analysis\WF Lever Analysis\licking_investigation\correct_trials_lick_v_no_lick\no_licks_-2_to_2\';
struct_dir = 'Z:\Analysis\WF Lever Analysis\StructuresForPlotting\';
output_dir = 'Z:\Analysis\WF Lever Analysis\licking_investigation\percent licks on earlies scatterplots\';
LTA_dir = 'Z:\Analysis\WF Lever Analysis\licking_investigation\lick_triggered_averages\lick_trig_mags_bout';
colors = {'k','k', 'b','b', 'g','g', 'c','c', 'm','m', 'r','r', 'y','y', 'k', 'b', 'r'}; 
plot_symbols = {'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 's', 's', 's'};

%generate some empty cell arrays to be filled in
percent_early_with_lick = {};
percent_early_with_lick{1} = [];

%load LTA data
load(LTA_dir);

%load tc data and determine differences
for session_num = [[5:length(days)]]
    %load indeces and struct
    load([early_TC_dir, days{session_num}, '_trial_indeces']);
    load([struct_dir, days{session_num}]);
    
    %get percent of early/corr trials with lick
    percent_early_with_lick{session_num} = length(lick_fail_trials)/(length(no_lick_fail_trials) + length(lick_fail_trials));
    
end

figure;
for session_num = [[5:length(days)]]
        plot(mean(lick_trig_mags{session_num}), percent_early_with_lick{session_num}*100, strcat(plot_symbols{session_num}, colors{session_num})); hold on; 
end
xlabel('df/f');
ylabel('Percen of earlies with lick');
title('magnitude of lick triggered average vs the percent of early trials with a lick');
ylim([0 100]);