%script for investigating licking onset relative to F onset. 

clear 
file_info_CRP; clear days89 days90 days91 days92 days93 days94;
behav_dir = 'Z:\Data\2P_imaging\behavior\';
TC_dir = 'Z:\Analysis\Cue_reward_pairing_analysis\2P\';
NR_peak_vals = [];
NR_peak_times = [];
trial_to_animal_inx = [];

for ii = 1:length(days_post)
    %load TC and cell category variables
    days_to_load = days_post;
    load_TC_variables;
    peak_window = (pre_cue_frames+1):(pre_cue_frames+round(1500/ifi));
    
    %takes the mean across cells for each neuron which was responsice on NR trials 
    NRr_avg_across_cells = squeeze(mean(NR_movie(:,NR_resp_cells,:),2));
    
    %for loop for determining peak time for each trial
    for trial_inx = 1:size(NRr_avg_across_cells,1)
        this_peak_val = max(NRr_avg_across_cells(trial_inx, peak_window));
        NR_peak_vals = [NR_peak_vals, this_peak_val];
        NR_peak_times = [NR_peak_times, find(NRr_avg_across_cells(cells_to_extract(trial_inx),:)==this_peak_val, 1, 'first')];
    end
    
    %generates index where the location matches to a value in NR_peak_vals and NR_peak_times and the value tells you which animal it belongs to
    trial_to_animal_inx = [trial_to_animal_inx, repmat(ii, [1,trial_inx])];
    
    %find the lick bout onset for each trial 
    
end



