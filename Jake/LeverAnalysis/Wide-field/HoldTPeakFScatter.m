%make scatter plots of the hold time duration vs the peak df/f for corr or
%early trials (trial by trial)
% also vs the peak df/f of corr-early (averaged across trials)
%adding in portions to normalize df/f to the mean of that session and then
%plot take a subset of trials from each session and plot across sessions

clear;
WF_plotting_lists_of_days;
struct_dir = 'Z:\Analysis\WF Lever Analysis\StructuresForPlotting\';
TC_dir = 'Z:\Analysis\WF Lever Analysis\LeverSummaryFolder\';
bx_data_dir = 'Z:\Data\WidefieldImaging\GCaMP\behavior\';
output_dir = 'Z:\Analysis\WF Lever Analysis\hold time vs peak f\';
colors = {'k','k', 'b','b', 'g','g', 'c','c', 'm','m', 'r','r', 'y','y', 'k', 'b', 'r'}; 
plot_symbols = {'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 's', 's', 's'};
colors_roi = {'b', 'r', 'g', 'k', 'c', 'm'};
peak_window = [4:8];
peak_diff{1} = [];
peak_corr_long_cell{1} = [];
peak_corr_short_cell{1} = [];
peak_early_long_cell{1} = [];
peak_early_short_cell{1} = [];
corr_hold_dur_cell{1} = []; 
early_hold_dur_cell{1} = [];

for session_num = 1:length(days)
    %load data
    bx_data = get_bx_data(bx_data_dir, days{session_num}); %use the bx data file to determine max Random hold duration 1500 vs 5000
    load([TC_dir, days{session_num}, '_success']);
    load([TC_dir, days{session_num}, '_fail']);
    load([struct_dir, days{session_num}]);
    corr_hold_dur = curr_struct.corr_hold_dur;
    early_hold_dur = curr_struct.early_hold_dur;
    
    %determine max rand hold
    if bx_data.randReqHoldMaxMs > 4000
        long_hold = 1;
    elseif bx_data.randReqHoldMaxMs < 2000
        long_hold = 0;
    else
        disp('error in max rand hold duration value'); 
        pause
    end
      
    %remove non-LS ROIs
    success_roi = success_roi(:,[LS_ROIs{session_num}],:);
    fail_roi = fail_roi(:,[LS_ROIs{session_num}],:);
    
    %find the peak values and get the difference
    peak_corr = max(success_roi(:,:,peak_window), [], 3);
    peak_early = max(fail_roi(:,:,peak_window), [], 3);
    peak_diff{session_num} = mean(peak_corr)-mean(peak_early);
    
    %plot each trials hold time vs peak df/f 
    figure;  
    subplot(1,2,1);
    for ROI_num = 1:size(success_roi,2)
        scatter(corr_hold_dur, peak_corr(:,ROI_num), colors_roi{ROI_num}); hold on;
        model_1 = fitlm(corr_hold_dur, peak_corr(:,ROI_num));
        plot(model_1)
    end
    ylabel('correct trials peak df/f'); xlabel('hold duration (ms)');
    
    subplot(1,2,2);
    for ROI_num = 1:size(fail_roi,2)
        scatter(early_hold_dur, peak_early(:,ROI_num), colors_roi{ROI_num}); hold on;
    end
    ylabel('early trials peak df/f'); xlabel('hold duration (ms)');
    savefig([output_dir, days{session_num}]);
    
    %collapse across ROIs
    peak_corr = mean(peak_corr,2);
    peak_early = mean(peak_early,2);
    
    %normalize df/f values to the mean
    peak_corr_mean = mean(peak_corr);
    peak_early_mean = mean(peak_early);
    peak_corr_norm = peak_corr/peak_corr_mean;
    peak_early_norm = peak_early/peak_early_mean;
    
    %store normalized df/f for corr/early in long_dur or short_dur
    if long_hold == 1
        peak_corr_long_cell{session_num} = peak_corr_norm;
        peak_early_long_cell{session_num} = peak_early_norm;
    else
        peak_corr_short_cell{session_num} = peak_corr_norm;
        peak_early_short_cell{session_num} = peak_early_norm;
    end
    
    %store hold dur corr/early no need for seperate long or short. Will
    %just use indeces to extract the right ones
    corr_hold_dur_cell{session_num} = corr_hold_dur; 
    early_hold_dur_cell{session_num} = early_hold_dur;
end

%determine the minimum number of trials for corr/early long/short
corr_long_min = cellfun(@length, peak_corr_long_cell);
corr_short_min = cellfun(@length, peak_corr_short_cell);
early_long_min = cellfun(@length, peak_early_long_cell);
early_short_min = cellfun(@length, peak_early_short_cell);
%
corr_long_min(corr_long_min==0) = [];
corr_short_min(corr_short_min==0) = [];
early_long_min(early_long_min==0) = [];
early_short_min(early_short_min==0) = [];
%
corr_long_min = min(corr_long_min);
corr_short_min = min(corr_short_min);
early_long_min = min(early_long_min);
early_short_min = min(early_short_min);

%plot 
%corr long
figure; 
subplot(1,2,1);
peak_corr_long_all = [];
hold_corr_long_all = [];
for session_num = 1:length(days)
    if ~isempty(peak_corr_long_cell{session_num})
        %take the same number of trials from each session
        rand_idx = sort(randsample(length(peak_corr_long_cell{session_num}), corr_long_min));
        scatter(corr_hold_dur_cell{session_num}([rand_idx]),   peak_corr_long_cell{session_num}([rand_idx]), strcat(plot_symbols{session_num}, colors{session_num})); hold on;
        
        %store points to fit a line
        peak_corr_long_all = [peak_corr_long_all, peak_corr_long_cell{session_num}([rand_idx])'];
        hold_corr_long_all = [hold_corr_long_all, corr_hold_dur_cell{session_num}([rand_idx])];
    end
end
corr_long_model = fitlm(hold_corr_long_all, peak_corr_long_all);
plot(corr_long_model);
title(['Correct trials R^2=', num2str(corr_long_model.Rsquared.Ordinary)]);
xlabel('hold duration (ms)');
ylabel('peak df/f');
hline(1);
ylim([-4 5]);

%early long
subplot(1,2,2);
peak_early_long_all = [];
hold_early_long_all = [];
for session_num = 1:length(days)
    if ~isempty(peak_early_long_cell{session_num})
        %take the same number of trials from each session
        rand_idx = sort(randsample(length(peak_early_long_cell{session_num}), early_long_min));
        scatter(early_hold_dur_cell{session_num}([rand_idx]),   peak_early_long_cell{session_num}([rand_idx]), strcat(plot_symbols{session_num}, colors{session_num})); hold on;
        
        %store points to fit a line
        peak_early_long_all = [peak_early_long_all, peak_early_long_cell{session_num}([rand_idx])'];
        hold_early_long_all = [hold_early_long_all, early_hold_dur_cell{session_num}([rand_idx])];
    end
end
early_long_model = fitlm(hold_early_long_all, peak_early_long_all);
plot(early_long_model);
title(['Early trials R^2=', num2str(early_long_model.Rsquared.Ordinary)]);
xlabel('hold duration (ms)');
ylabel('peak df/f');
hline(1);
ylim([-4 5]);
suptitle('hold duration vs peak df/f across animals: long duration random hold');
savefig([output_dir, 'across_animals_long_dur_hold']);


%plot short duration random hold sessions
%corr short
figure;
subplot(1,2,1);
peak_corr_short_all = [];
hold_corr_short_all = [];
for session_num = 1:length(peak_corr_short_cell)
    if ~isempty(peak_corr_short_cell{session_num})
        %take the same number of trials from each session
        rand_idx = sort(randsample(length(peak_corr_short_cell{session_num}), corr_short_min));
        scatter(corr_hold_dur_cell{session_num}([rand_idx]),   peak_corr_short_cell{session_num}([rand_idx]), strcat(plot_symbols{session_num}, colors{session_num})); hold on;
        
        %store points to fit a line
        peak_corr_short_all = [peak_corr_short_all, peak_corr_short_cell{session_num}([rand_idx])'];
        hold_corr_short_all = [hold_corr_short_all, corr_hold_dur_cell{session_num}([rand_idx])];
    end
end
corr_short_model = fitlm(hold_corr_short_all, peak_corr_short_all);
plot(corr_short_model);
title(['Correct trials R^2=', num2str(corr_short_model.Rsquared.Ordinary)]);
xlabel('hold duration (ms)');
ylabel('peak df/f');
hline(1);
ylim([-4 5]);

%early short
subplot(1,2,2);
peak_early_short_all = [];
hold_early_short_all = [];
for session_num = 1:length(peak_early_short_cell)
    if ~isempty(peak_early_short_cell{session_num})
        %take the same number of trials from each session
        rand_idx = sort(randsample(length(peak_early_short_cell{session_num}), early_short_min));
        scatter(early_hold_dur_cell{session_num}([rand_idx]),   peak_early_short_cell{session_num}([rand_idx]), strcat(plot_symbols{session_num}, colors{session_num})); hold on;
        
        %store points to fit a line
        peak_early_short_all = [peak_early_short_all, peak_early_short_cell{session_num}([rand_idx])'];
        hold_early_short_all = [hold_early_short_all, early_hold_dur_cell{session_num}([rand_idx])];
    end
end
early_short_model = fitlm(hold_early_short_all, peak_early_short_all);
plot(early_short_model);
title(['Early trials R^2=', num2str(early_short_model.Rsquared.Ordinary)]);
xlabel('hold duration (ms)');
ylabel('peak df/f');
hline(1);
ylim([-4 5]);
suptitle('hold duration vs peak df/f across animals: short duration random hold');
savefig([output_dir, 'across_animals_short_dur_hold']);








