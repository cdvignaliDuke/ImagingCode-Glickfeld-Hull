% Script for isolating NR, OR, or both responsive dendrites and
% comparing their relative fractions within a session. 
clear 
file_info_CRP;
behav_dir = 'Z:\Data\2P_imaging\behavior\';
TC_dir = 'Z:\Analysis\Cue_reward_pairing_analysis\2P\';
comparison_days = [1,2];

tot_num_cells_day1 = [];
tot_num_cells_post = [];
rew_only_resp_cells_post = [];
rew_only_resp_cells_day1 = [];
for ii = 1:length(days_1)
    %load TC and cell category variables
    days_to_load = days_1;
    load_TC_variables;
    
    %Further segregate neuron IDs according to responsivity to cue or reward
    tot_num_cells_day1 = [tot_num_cells_day1, length(allresp_cells)];
    rew_only_resp_cells_temp = NR_resp_cells(~ismember(NR_resp_cells, OR_resp_cells));
    rew_only_resp_cells_day1 = [rew_only_resp_cells_day1, length(rew_only_resp_cells_temp)]; 
    %use the derivative of NR traces to ID cue vs cue and reward resp
    %cells. Start with cells that are resp to NR and OR cells 
    avg_NR_diff = diff(avg_NR,1,2);
    figure; 
    for iii = 6:10%size(avg_NR_diff,1);
        plot(avg_NR_diff(iii,:)); hold on;
    end
    std_traces = std(avg_NR_diff([6:10],:),[],2)
end

for ii = 1:length(days_post)
    %load TC and cell category variables
    days_to_load = days_post;
    load_TC_variables;
    
    %Further segregate neuron IDs according to responsivity to cue or reward
    tot_num_cells_post = [tot_num_cells_post, length(allresp_cells)];
    rew_only_resp_cells_temp = NR_resp_cells(~ismember(NR_resp_cells, OR_resp_cells));
    rew_only_resp_cells_post = [rew_only_resp_cells_post, length(rew_only_resp_cells_temp)];
end

%plot some scatterplots
rew_resp_fraction_day1 = rew_only_resp_cells_day1./tot_num_cells_day1;
rew_resp_fraction_post = rew_only_resp_cells_post./tot_num_cells_post;
figure; 
c = [1:6];
colors_scatter = [1,0,1; 0,1,1; 1,0,0; 0,1,0; 0,0,1; 0,0,0]; %{'r', 'b', 'k', 'm', 'c', 'g'}
scatter(rew_resp_fraction_day1, rew_resp_fraction_post, [], colors_scatter, 'filled'); hold on;
days_used = {'img90',  'img91', 'img92', 'img93', 'img89', 'img94'};
%legend(days_used);
x = 0:.1:1;
y=x; 
plot(x,y);
xlabel('day 1');
ylabel('post learning');
title('500ms fraction of dendrites responsive to reward only');



