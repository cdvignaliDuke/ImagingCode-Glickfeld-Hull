% script for finding the peak F on a trial by trial basis for two-photon
% cue-reward pairing datasets. It uses the same techniques as the scatter
% plot code for the wide-field lever analysis.

clear
file_info_CRP;
behav_dir = 'Z:\Data\2P_imaging\behavior\';
TC_dir = 'Z:\Analysis\Cue_reward_pairing_analysis\2P\';
NR_peak_times_avg_pop = [];
OR_peak_times_avg_pop = [];
NR_peak_vals_avg_pop = [];
OR_peak_vals_avg_pop = [];
cell_to_animal_ind = []; 

%modifiables
days_to_load = days_1;
session_type = '500ms day 1';

for ii = 1:length(days_to_load);  %session indicator   post learning datasets: 23,26,30,45,47,52   day1 datasets: 22,25,29,42,43,51
    %load TC and cell category variables
    load_TC_variables;
    
    %allocate memory and define variables
    cells_to_extract = NR_resp_cells; %NR_resp_cells;   %OR_resp_cells;   [1:size(allresp_cells,2)];   MODIFIABLE
    cell_type_extracted = 'NR responsive cells';                                     %MODIFIABLE
    ifi = double(ifi);
    peak_window = (pre_cue_frames+1):(pre_cue_frames+round(1500/ifi));
    NR_peak_times = nan(1, length(cells_to_extract));
    OR_peak_times = nan(1, length(cells_to_extract));
    NR_peak_vals = nan(1, length(cells_to_extract));
    OR_peak_vals = nan(1, length(cells_to_extract));
    
    %main for loop extracts peak times and latencies from each cell's avg TC
    for cell_inx = 1:length(cells_to_extract)
        this_peak_val = max(avg_NR(cells_to_extract(cell_inx), peak_window));
        NR_peak_vals(cell_inx) = this_peak_val;
        NR_peak_times(cell_inx) = find(avg_NR(cells_to_extract(cell_inx),:)==this_peak_val, 1, 'first');
        
        this_peak_val = max(avg_OR(cells_to_extract(cell_inx), peak_window));
        OR_peak_vals(cell_inx) = this_peak_val;
        OR_peak_times(cell_inx) = find(avg_OR(cells_to_extract(cell_inx),:)==this_peak_val, 1, 'first');
    end
    
    %convert to ms relative to cue
    NR_peak_times = (NR_peak_times-pre_cue_frames)*ifi;
    OR_peak_times = (OR_peak_times-pre_cue_frames)*ifi;
    
    %get mean and sem of peak times and vals across dendrites
    NR_peak_times_avg = mean(NR_peak_times);
    OR_peak_times_avg = mean(OR_peak_times);
    NR_peak_vals_avg = mean(NR_peak_vals);
    OR_peak_vals_avg = mean(OR_peak_vals);
    NR_peak_times_sem = std(NR_peak_times)/sqrt(length(NR_peak_times));
    OR_peak_times_sem = std(OR_peak_times)/sqrt(length(OR_peak_times));
    NR_peak_vals_sem = std(NR_peak_vals)/sqrt(length(NR_peak_times));
    OR_peak_vals_sem = std(OR_peak_vals)/sqrt(length(OR_peak_times));
    
    %scatter plot of peak times and mags across trials for each neuron
%     figure; 
%     subplot(1,2,1);
%     scatter(NR_peak_times, OR_peak_times); hold on;
%     x_plot = [0:1500]; plot(x_plot,x_plot);
%     plot(NR_peak_times_avg, OR_peak_times_avg,  'r', 'Marker', 'o', 'MarkerFaceColor', 'r');
%     title(['Latency to peak df/f for rewarded and reward omission trials relative to cue onset: ', session_date, ' ', session_mouseID]);
%     xlabel('Rewarded trials (ms)');
%     ylabel('Reward omission trials (ms)');
%     
%     subplot(1,2,2);
%     scatter(NR_peak_vals, OR_peak_vals); hold on;
%     x_plot = [-0.2:1.5]; plot(x_plot,x_plot);
%     plot(NR_peak_vals_avg, OR_peak_vals_avg,  'r', 'Marker', 'o', 'MarkerFaceColor', 'r');
%     title(['magnitude of peak df/f for rewarded and reward omission trials: ', session_date, ' ', session_mouseID]);
%     xlabel('Rewarded trials (df/f)');
%     ylabel('Reward omission trials (df/f)');
%     if min([NR_peak_vals, OR_peak_vals])<0;
%         min_val = min([NR_peak_vals, OR_peak_vals])*1.1;
%     else 
%         min_val = 0;
%     end
%     ylim([min_val max([NR_peak_vals, OR_peak_vals])*1.1]);
%     xlim([min_val max([NR_peak_vals, OR_peak_vals])*1.1]);
%     suptitle([session_type, ' ', cell_type_extracted]);
%     
%     savefig(['Z:\Analysis\Cue_reward_pairing_analysis\2P\CRP_scatterplots\', session_date, '_', session_mouseID, ' ', session_type, ' ', cell_type_extracted]);
    
    %store stats for summary fig
    NR_peak_times_avg_pop = [NR_peak_times_avg_pop, NR_peak_times];
    OR_peak_times_avg_pop = [OR_peak_times_avg_pop, OR_peak_times];
    NR_peak_vals_avg_pop = [NR_peak_vals_avg_pop, NR_peak_vals];
    OR_peak_vals_avg_pop = [OR_peak_vals_avg_pop, OR_peak_vals];
    
    %generate an index which can match each cell to the animal it came from
    cell_to_animal_ind = [cell_to_animal_ind, repmat(ii, 1, length(cells_to_extract))];
end

%plot summary scatter
NR_peak_times_sem_pop = std(NR_peak_times_avg_pop)/sqrt(length(NR_peak_times_avg_pop));
OR_peak_times_sem_pop = std(OR_peak_times_avg_pop)/sqrt(length(OR_peak_times_avg_pop));
NR_peak_vals_sem_pop = std(NR_peak_vals_avg_pop)/sqrt(length(NR_peak_vals_avg_pop));
OR_peak_vals_sem_pop = std(OR_peak_vals_avg_pop)/sqrt(length(OR_peak_vals_avg_pop));

figure;
subplot(1,2,1);
scatter(NR_peak_times_avg_pop, OR_peak_times_avg_pop); hold on;
plot(mean(NR_peak_times_avg_pop), mean(OR_peak_times_avg_pop),  'r', 'Marker', 'o', 'MarkerFaceColor', 'r');
x_plot = [0:1500]; plot(x_plot,x_plot, 'k');
title(['Population data for all animals: peak latency relative to cue n= ', num2str(length(NR_peak_times_avg_pop))]);
xlabel('Rewarded trials (ms)');
ylabel('Reward omission trials (ms)');
ylim([0 1500]);
xlim([0 1500]);
errorbarxy(mean(NR_peak_times_avg_pop), mean(OR_peak_times_avg_pop), NR_peak_times_sem_pop, OR_peak_times_sem_pop);

subplot(1,2,2);
scatter(NR_peak_vals_avg_pop, OR_peak_vals_avg_pop); hold on;
plot(mean(NR_peak_vals_avg_pop), mean(OR_peak_vals_avg_pop),  'r', 'Marker', 'o', 'MarkerFaceColor', 'r'); 
x_plot = [0:2]; plot(x_plot,x_plot, 'k');
title(['Population data for all animals: peak magnitude n = ',  num2str(length(OR_peak_times_avg_pop))]);
xlabel('Rewarded trials (df/f)');
ylabel('Reward omission trials (df/f)');
suptitle(['All animals ', session_type, ' ', cell_type_extracted]);
if min([NR_peak_vals_avg_pop, OR_peak_vals_avg_pop])<0;
    min_val = min([NR_peak_vals_avg_pop, OR_peak_vals_avg_pop])*1.1;
else
    min_val = 0;
end
ylim([min_val max([NR_peak_vals_avg_pop, OR_peak_vals_avg_pop])*1.1]);
xlim([min_val max([NR_peak_vals_avg_pop, OR_peak_vals_avg_pop])*1.1]);
errorbarxy(mean(NR_peak_vals_avg_pop), mean(OR_peak_vals_avg_pop), NR_peak_vals_sem_pop, OR_peak_vals_sem_pop);
savefig(['Z:\Analysis\Cue_reward_pairing_analysis\2P\CRP_scatterplots\', 'All animals ', session_type, ' ', cell_type_extracted]);

% NR_peak_times_avg_anim = [];
% OR_peak_times_avg_anim = [];
% NR_peak_vals_avg_anim = [];
% OR_peak_vals_avg_anim = [];
% for ii = 1:length(days_to_load)
%     cells_from_this_animal = find(cell_to_animal_ind==ii);
%     NR_peak_times_avg_anim = [NR_peak_times_avg_anim, mean(NR_peak_times_avg_pop(cells_from_this_animal))];
%     OR_peak_times_avg_anim = [OR_peak_times_avg_anim, mean(OR_peak_times_avg_pop(cells_from_this_animal))];
%     NR_peak_vals_avg_anim = [NR_peak_vals_avg_anim, mean(NR_peak_vals_avg_pop(cells_from_this_animal))];
%     OR_peak_vals_avg_anim = [OR_peak_vals_avg_anim, mean(OR_peak_vals_avg_pop(cells_from_this_animal))];
%     
%     NR_peak_times_sem_anim = [NR_peak_times_sem_anim, std(NR_peak_times_avg_pop(cells_from_this_animal))/sqrt(length(cells_from_this_animal))];
%     OR_peak_times_sem_anim = [OR_peak_times_sem_anim, std(OR_peak_times_avg_pop(cells_from_this_animal))/sqrt(length(cells_from_this_animal))];
%     NR_peak_vals_sem_anim = [NR_peak_vals_sem_anim, std(NR_peak_vals_avg_pop(cells_from_this_animal))/sqrt(length(cells_from_this_animal))];
%     OR_peak_vals_sem_anim = [OR_peak_vals_sem_anim, std(OR_peak_vals_avg_pop(cells_from_this_animal))/sqrt(length(cells_from_this_animal))];
% end

% x_axis = [1,2];
% for ii = 1:length(days_to_load)
%     subplot(2,2,3)
%     plot(x_axis, )
%     
%     
% end






   %GOES INSIDE FORLOOP. old code for determining peak times and vales on a
   %trial by trial basis
%     %main for loop which extracts peak times and magnitudes for each cell for each trial
%     for cell_inx = 1:length(cells_to_extract)  
%         num_NR = size(NR_movie,1);
%         num_OR = size(OR_movie,1);
%         for trial_no = 1:num_NR
%             if mean(NR_movie(trial_no, cells_to_extract(cell_inx), peak_window))== 0;
%                 continue
%             end
%             this_peak_val = max(NR_movie(trial_no, cells_to_extract(cell_inx), peak_window));
%             NR_peak_vals(trial_no, cell_inx) = this_peak_val;
%             NR_peak_times(trial_no, cell_inx) = find(NR_movie(trial_no, cells_to_extract(cell_inx),:)==this_peak_val, 1, 'first');
%         end
%         for trial_no = 1:num_OR
%             if mean(OR_movie(trial_no, cells_to_extract(cell_inx), peak_window))== 0;
%                 continue
%             end
%             this_peak_val = max(OR_movie(trial_no, cells_to_extract(cell_inx), peak_window));
%             OR_peak_vals(trial_no, cell_inx) = this_peak_val;
%             OR_peak_times(trial_no, cell_inx) = find(OR_movie(trial_no, cells_to_extract(cell_inx),:)==this_peak_val, 1, 'first');
%         end
%     end