%plot the magnitude of the df/f response vs the time of the first lick
%relative to lever release
clear;
WF_plotting_lists_of_days;
bx_source     = ['Z:\Data\WidefieldImaging\GCaMP\behavior\'];
image_source_base  = ['Z:\Data\WidefieldImaging\GCaMP\']; %location of permanently stored image files for retreiving meta data
image_dest_base    = ['Z:\Analysis\WF Lever Analysis\BxAndAnalysisOutputs\']; %stores the data on crash in the lever analysis folder
bx_outputs_dir = ['Z:\Analysis\WF Lever Analysis\BxAndAnalysisOutputs\BxOutputs\'];
TC_dir = ['Z:\Analysis\WF Lever Analysis\LeverSummaryFolder\'];
colors = [1,0,0; 0,1,0; 0,0,1; 0.5,0.5,0.5; 1,0,1; 1,1,0; 0,1,1]; %sets up the color scheme for plotting multiple ROIs with errorbar     
old_cd = cd; %save old cd so I can restore it later
analysis_window = [3:9];
x_axis = [-4:3];
release_frame = 6;
peak_window = [4:8];
ROI_colors = {'r', 'b', 'k', 'm', 'g', 'c', 'y'};

%% SECTION TWO  CORR TRIALS
for ii = 1:length(days);
    %load bx and corr trials TCs
    b_data = get_bx_data(bx_source, days{ii});  %find the correct behavior file and loads it.
    bx_out_dir  = [bx_outputs_dir days{ii} '_bx_outputs'];
    load(bx_out_dir);
    load([TC_dir, days{ii}, '_success']);
    lick_trace_corr = licking_data.lick_trace_succ; 
    
    %plot failed trials with vs without licks in window
    offset = 0;
    x_axis = [-5:10];
    
    %determine time of first lick relative to lever release
    first_lick_frame=NaN(1, size(lick_trace_corr,1));
    for trial_num = 1:size(lick_trace_corr,1)
        if isempty(find(lick_trace_corr(trial_num, analysis_window), 1, 'first'));
            first_lick_frame(trial_num) = analysis_window(end) -release_frame +2;
        else
            first_lick_frame(trial_num) = find(lick_trace_corr(trial_num, analysis_window), 1, 'first') - 3;
        end
    end
    
    figure; hold on;
    for ROI_num = 1:size(success_roi,2)
        %baseline each trial 
        success_this_ROI = squeeze(success_roi(:, ROI_num, :));
        shift = mean(success_this_ROI(:,[1:3]),2);
        shift = repmat(shift,1,16);
        success_this_ROI = success_this_ROI-shift;
        success_this_ROI_sem = squeeze(std(success_this_ROI, 1))'/sqrt(size(success_this_ROI,1));
        
        %determine the peak df/f for each trial
        peak_vals = max(success_this_ROI(:,peak_window), [],  2)';
        plot(first_lick_frame, peak_vals, 'o', 'MarkerEdgeColor', ROI_colors{ROI_num});
        
        %plot the mean values for each ROI
        all_slots = unique(first_lick_frame);
        this_ROI_mean = [];
        for time_slot = 1:length(all_slots)
            this_ind = find(first_lick_frame==all_slots(time_slot));
            this_ROI_mean = [this_ROI_mean, mean(peak_vals(this_ind))];
        end
        plot(all_slots, this_ROI_mean, ROI_colors{ROI_num});
        
    end
    axis tight;
    title([days(ii)]);
    xlabel('frame num relative to lever release');
    ylabel(['peak df/f magnitude in frames ' num2str(peak_window-release_frame)]);
    xlim([-4 6]);
    savefig(['Z:\Analysis\WF Lever Analysis\licking_investigation\df_vs_lick_scatters\', days{ii}]);
end










