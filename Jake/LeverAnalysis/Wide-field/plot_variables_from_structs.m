%script for loading and plotting a bunch of variables
%list of current days being analyzed
clear;
WF_plotting_lists_of_days;
struct_dir = 'Z:\Analysis\WF Lever Analysis\StructuresForPlotting\'; 
xls_dir = 'Z:\Data\WidefieldImaging\GCaMP\WF_exp_spreadsheet';
struct_files = dir(struct_dir);
meta_struct = struct('one', []);
color_palette = colormap;
color_palette = color_palette([round(linspace(1,64,length(days_chrono_order)))],:);

%days = {'160209_img36', '151222_img32', '151019_img30', '160725_img53', '160905_img55'};  %no lever controls 161109_img61 161109_img59   Do not meet criterion 160208_img35

for ii = 1:2
    if length(struct_files(ii).name) == 1 | length(struct_files(ii).name) == 2;
        struct_files(ii) = [];
    end
end
for ii = 1:2
    if length(struct_files(ii).name) == 1 | length(struct_files(ii).name) == 2;
        struct_files(ii) = [];
    end
end
for ii = 1:length(struct_files)
    if struct_files(ii).name(1:3) == 'Old';
        struct_files(ii) = [];
        break
    end
end
for ii = 1:length(struct_files)
    if struct_files(ii).name(1:3) == 'cur';
        struct_files(ii) = [];
        break
    end
end
for ii = 1:length(struct_files);
    load([struct_dir, struct_files(ii).name]);
    meta_struct.(['i' struct_files(ii).name(1:12)]) = load([struct_dir, struct_files(ii).name]);
end
meta_struct = rmfield(meta_struct, 'one');
x = [-1000:10:1000];
y = x;

%% peak magnitudes CORRECT vs EARLY
corr_magnitude_mean = [];
early_magnitude_mean = [];
corr_magnitude_sem = [];
early_magnitude_sem = [];
for ii = 1:length(struct_files)
    corr_magnitude_mean = [corr_magnitude_mean, meta_struct.(['i', struct_files(ii).name(1:12)]).curr_struct.corr_magnitude_mean];
    early_magnitude_mean = [early_magnitude_mean, meta_struct.(['i', struct_files(ii).name(1:12)]).curr_struct.early_magnitude_mean];
    corr_magnitude_sem = [corr_magnitude_sem, meta_struct.(['i', struct_files(ii).name(1:12)]).curr_struct.corr_magnitude_std / sqrt(size(meta_struct.(['i', struct_files(ii).name(1:12)]).curr_struct.corr_TCs,2))];
    early_magnitude_sem = [early_magnitude_sem, meta_struct.(['i', struct_files(ii).name(1:12)]).curr_struct.early_magnitude_std / sqrt(size(meta_struct.(['i', struct_files(ii).name(1:12)]).curr_struct.early_TCs,2))];
end

%peak magnitude df/f correct vs early
figure; hold on;
for ii = 1:length(corr_magnitude_mean)
    plot(corr_magnitude_mean(ii), early_magnitude_mean(ii),'Marker', 'o', 'Linestyle', 'none', 'Color', color_palette(ii,:), 'MarkerFaceColor', color_palette(ii,:)); hold on;
    errorbarxy(corr_magnitude_mean(ii), early_magnitude_mean(ii), corr_magnitude_sem(ii), early_magnitude_sem(ii)); hold on;
end
hold on; plot(x,y, 'Color', 'k', 'Linestyle', '--'); vline(0,'k'); hline(0,'k')
xlim([0, 0.4]);
ylim([0, 0.4]);
title('mean peak df/f: correct vs early');
xlabel('correct trials mean peak magnitude (df/f)');
ylabel('early trials mean peak magnitude (df/f)');
legend(days_chrono_order);

%% peak magnitudes CORRECT vs tooFast correct
corr_magnitude_mean = [];
tooFast_magnitude_mean = [];
corr_magnitude_sem = [];
tooFast_magnitude_sem = [];
days_used = [];
for ii = 1:length(struct_files)
    if isfield(meta_struct.(['i', struct_files(ii).name(1:12)]).curr_struct, 'tooFast_magnitude_mean') & size(meta_struct.(['i', struct_files(ii).name(1:12)]).curr_struct.tooFast_magnitude,2) > 5 
    corr_magnitude_mean = [corr_magnitude_mean, meta_struct.(['i', struct_files(ii).name(1:12)]).curr_struct.corr_magnitude_mean];
    tooFast_magnitude_mean = [tooFast_magnitude_mean, meta_struct.(['i', struct_files(ii).name(1:12)]).curr_struct.tooFast_magnitude_mean];
    corr_magnitude_sem = [corr_magnitude_sem, meta_struct.(['i', struct_files(ii).name(1:12)]).curr_struct.corr_magnitude_std / sqrt(size(meta_struct.(['i', struct_files(ii).name(1:12)]).curr_struct.corr_TCs,2))];
    tooFast_magnitude_sem = [tooFast_magnitude_sem, meta_struct.(['i', struct_files(ii).name(1:12)]).curr_struct.tooFast_magnitude_std / sqrt(size(meta_struct.(['i', struct_files(ii).name(1:12)]).curr_struct.tooFast_TCs,2))];
    days_used = [days_used, ii];
    end
end

%peak magnitude df/f correct vs tooFaast
figure; hold on;
for ii = 1: length(corr_magnitude_mean)
    plot(corr_magnitude_mean(ii), tooFast_magnitude_mean(ii),'Marker', 'o', 'Linestyle', 'none', 'Color', color_palette(days_used(ii),:), 'MarkerFaceColor', color_palette(ii,:)); hold on;
end
legend(days_chrono_order(days_used));
for ii = 1: length(corr_magnitude_mean)
    errorbarxy(corr_magnitude_mean(ii), tooFast_magnitude_mean(ii), corr_magnitude_sem(ii), tooFast_magnitude_sem(ii)); hold on;
end
hold on; plot(x,y, 'Color', 'k', 'Linestyle', '--');
xlim([0, 0.4]);
ylim([0, 0.4]);
title('mean peak df/f: correct vs tooFast corrects');
xlabel('correct trials mean peak magnitude (df/f)');
ylabel('tooFast trials mean peak magnitude (df/f)');

%% mean lick rate over 0.5s post release   VS   peak df/f
interval = 5; %interval overwhich we will look at the lick rate (in frame numbers at 10Hz)
lever_release = 6; %frame on which lever release occurs
bx_outputs_dir = 'Z:\Analysis\WF Lever Analysis\BxAndAnalysisOutputs\BxOutputs\';
corr_magnitude_mean = [];
early_magnitude_mean = [];
corr_magnitude_sem = [];
early_magnitude_sem = [];
corr_lick_rate_mean = [];
early_lick_rate_mean = [];
corr_lick_rate_sem = [];
early_lick_rate_sem = [];
for ii = 1:length(struct_files)
    corr_magnitude_mean = [corr_magnitude_mean, meta_struct.(['i', struct_files(ii).name(1:12)]).curr_struct.corr_magnitude_mean];
    early_magnitude_mean = [early_magnitude_mean, meta_struct.(['i', struct_files(ii).name(1:12)]).curr_struct.early_magnitude_mean];
    corr_magnitude_sem = [corr_magnitude_sem, meta_struct.(['i', struct_files(ii).name(1:12)]).curr_struct.corr_magnitude_std / sqrt(size(meta_struct.(['i', struct_files(ii).name(1:12)]).curr_struct.corr_TCs,2))];
    early_magnitude_sem = [early_magnitude_sem, meta_struct.(['i', struct_files(ii).name(1:12)]).curr_struct.early_magnitude_std / sqrt(size(meta_struct.(['i', struct_files(ii).name(1:12)]).curr_struct.early_TCs,2))];
    load([bx_outputs_dir, struct_files(ii).name(1:12), '_bx_outputs'], 'licking_data');
    corr_tbyt_lick_rate = mean(licking_data.lick_trace_succ(:,[lever_release+1:lever_release+interval]),2)*10;   %find the lick rate for the first 0.5s after lever release for each trial
    early_tbyt_lick_rate = mean(licking_data.lick_trace_fail(:,[lever_release+1:lever_release+interval]),2)*10;
    corr_lick_rate_mean = [corr_lick_rate_mean, mean(corr_tbyt_lick_rate)];
    early_lick_rate_mean = [early_lick_rate_mean, mean(early_tbyt_lick_rate)];
    corr_lick_rate_sem = [corr_lick_rate_sem, std(corr_tbyt_lick_rate)/sqrt(length(corr_tbyt_lick_rate))];
    early_lick_rate_sem = [early_lick_rate_sem, std(early_lick_rate_mean)/sqrt(length(early_lick_rate_mean))];
end

figure; hold on; 
for ii=3:length(corr_magnitude_mean)
    plot(corr_magnitude_mean(ii), corr_lick_rate_mean(ii),'Marker', '*', 'Linestyle', 'none', 'Color', 'k');
    plot(early_magnitude_mean(ii), early_lick_rate_mean(ii), 'Marker', '*', 'Linestyle', 'none', 'Color', 'r');
end
for ii=3:length(corr_magnitude_mean)
    errorbar(corr_magnitude_mean(ii), corr_lick_rate_mean(ii), corr_lick_rate_sem(ii), 'Color', 'k');
    errorbar(early_magnitude_mean(ii), early_lick_rate_mean(ii), early_lick_rate_sem(ii),'Color', 'r')
end
for ii=3:length(corr_magnitude_mean)
    plot([corr_magnitude_mean(ii)-corr_magnitude_sem(ii), corr_magnitude_mean(ii)+corr_magnitude_sem(ii)],[corr_lick_rate_mean(ii), corr_lick_rate_mean(ii)], 'Color', 'k')
    plot([early_magnitude_mean(ii)-early_magnitude_sem(ii), early_magnitude_mean(ii)+early_magnitude_sem(ii)],[early_lick_rate_mean(ii), early_lick_rate_mean(ii)], 'Color', 'r')

end
title(['peak df/f  vs  lick rate for the 0.5s after lever release']);
xlabel('peak df/f');
ylabel('lick rate (Hz) in 0.5ms post lever release');

% Lick rate ratio (early/corr)  vs  df/f ratio (early/corr)
lick_rate_ratio = early_lick_rate_mean./corr_lick_rate_mean;
corr_early_dfoverf_ratio = [];
for ii = 1:length(struct_files)
      corr_early_dfoverf_ratio = [corr_early_dfoverf_ratio, meta_struct.(['i', struct_files(ii).name(1:12)]).curr_struct.corr_early_dfoverf_ratio];
end

figure; hold on; 
for ii = 3:length(lick_rate_ratio)
    plot(1/corr_early_dfoverf_ratio(ii),lick_rate_ratio(ii), 'Marker', '*', 'Linestyle', 'none', 'Color', 'k')
end
title('peak df/f ratio (early/correct) vs lick rate ratio first 0.5s (early/corr)');
xlabel('peak df/f ratio (early/correct)')
ylabel('lick rate ratio (early/correct)')
xlim([-0.1 1]);
ylim([-0.1 1]);


%% get onset latency std for Corr v Early 
corr_onset_latency_std = [];
early_onset_latency_std = [];
corr_cue_aligned_onset_latency_std = [];
for ii = 1:length(struct_files)
    corr_onset_latency_std = [corr_onset_latency_std, meta_struct.(['i', struct_files(ii).name(1:12)]).curr_struct.corr_onset_latency_std];
    early_onset_latency_std = [early_onset_latency_std, meta_struct.(['i', struct_files(ii).name(1:12)]).curr_struct.early_onset_latency_std];
    corr_cue_aligned_onset_latency_std = [corr_cue_aligned_onset_latency_std, meta_struct.(['i', struct_files(ii).name(1:12)]).curr_struct.corr_cue_aligned_onset_latency_std];
end

%Plot onset latency std for Corr v early
figure; hold on;
for ii = 1: length(corr_onset_latency_std)
    plot(corr_onset_latency_std(ii), early_onset_latency_std(ii),'Marker', 'o', 'Linestyle', 'none', 'Color', color_palette(ii,:), 'MarkerFaceColor', color_palette(ii,:));
end
hold on; plot(x,y, 'Color', 'k', 'Linestyle', '--');
xlim([60, 260]);
ylim([60, 260]);
title('onset latency standard deviations: correct vs early');
xlabel('onset latency standard deviations for correct trials (ms)');
ylabel('onset latency standard deviations for early trials (ms)');
legend(days_chrono_order);

%Plot onset latency std for Corr v CueCorr
figure; 
subplot(1,2,1); hold on;
for ii = 1:length(corr_onset_latency_std)
    plot(corr_onset_latency_std(ii), corr_cue_aligned_onset_latency_std(ii),'Marker', 'o', 'Linestyle', 'none', 'Color', color_palette(ii,:), 'MarkerFaceColor', color_palette(ii,:));
end
hold on; plot(x,y, 'Color', 'k', 'Linestyle', '--');
xlim([50, 250]);
ylim([50, 250]);
hline(0); vline(0);
title('onset latency standard deviations: correct vs cue aligned');
xlabel('onset latency standard deviations for correct trials (ms)');
ylabel('onset latency standard deviations for cue aligned correct trials (ms)');
legend(days_chrono_order);
subplot(1,2,2);
bar([1,2],[mean(corr_onset_latency_std), mean(corr_cue_aligned_onset_latency_std)]);
errorbar([1,2],[mean(corr_onset_latency_std), mean(corr_cue_aligned_onset_latency_std)], [std(corr_onset_latency_std), std(corr_cue_aligned_onset_latency_std)]);
%[(std(corr_onset_latency_std)/sqrt(length(corr_onset_latency_std))), (std(corr_cue_aligned_onset_latency_std)/sqrt(length(corr_cue_aligned_onset_latency_std)))])

%%  %get mean onset latency for Corr v Early 
corr_onset_latency = [];
early_onset_latency = [];
corr_cue_aligned_onset_latency = [];
for ii = 1:length(struct_files)
    corr_onset_latency = [corr_onset_latency, meta_struct.(['i', struct_files(ii).name(1:12)]).curr_struct.corr_onset_latency_mean];
    early_onset_latency = [early_onset_latency, meta_struct.(['i', struct_files(ii).name(1:12)]).curr_struct.early_onset_latency_mean];
    corr_cue_aligned_onset_latency = [corr_cue_aligned_onset_latency, meta_struct.(['i', struct_files(ii).name(1:12)]).curr_struct.corr_cue_aligned_onset_latency_mean];
end

figure; hold on;
for ii = 1:length(corr_onset_latency)
    plot(corr_onset_latency(ii), early_onset_latency(ii),'Marker', '*', 'Linestyle', 'none', 'Color', color_palette(ii,:));
end
hold on; plot(x,y, 'Color', 'k', 'Linestyle', '--');
xlim([-250, 100]);
ylim([-250, 100]);
title('mean onset latencies: correct vs early');
xlabel('mean onset latencies for correct trials (ms)');
ylabel('mean onset latencies for early trials (ms)');
legend(days_chrono_order);

%% latency to peak for corr vs cue_corr
corr_peak_latency_std = [];
early_peak_latency_std = [];
corr_cue_aligned_peak_latency_std = [];
for ii = 1:length(struct_files)
    corr_peak_latency_std = [corr_peak_latency_std, meta_struct.(['i', struct_files(ii).name(1:12)]).curr_struct.corr_onset_latency_std];
    corr_cue_aligned_peak_latency_std = [corr_cue_aligned_peak_latency_std, meta_struct.(['i', struct_files(ii).name(1:12)]).curr_struct.corr_cue_aligned_onset_latency_std];
end

figure; subplot(1,2,1);
hold on;
for ii = 1:length(corr_peak_latency_std)
    plot(corr_peak_latency_std(ii), corr_cue_aligned_peak_latency_std(ii),'Marker', 'o', 'Linestyle', 'none', 'Color', color_palette(ii,:), 'MarkerFaceColor', color_palette(ii,:));
end
hold on; plot(x,y, 'Color', 'k', 'Linestyle', '--');
xlim([-350, 350]);
ylim([-350, 350]);
title('peak df/f latency std: correct vs cue aligned');
xlabel('Lever aligned corrects latency to peak std (ms)');
ylabel('Cue aligned corrects latency to peak std (ms)');
legend(days_chrono_order);
subplot(1,2,2);
bar([1,2],[mean(corr_peak_latency_std), mean(corr_cue_aligned_peak_latency_std)]);
errorbar([1,2],[mean(corr_peak_latency_std), mean(corr_cue_aligned_peak_latency_std)], [(std(corr_peak_latency_std)/sqrt(length(corr_peak_latency_std))), (std(corr_cue_aligned_peak_latency_std)/sqrt(length(corr_cue_aligned_peak_latency_std)))])
%-------------------------------------
%% corr/earlyF ratio versus   peak percent correct
%EXCLUDING 28 and 29 due to now %corr data
corrEarlyRatio = [];
peak_percent_corr = []; 
percent_corr = []; 
recent_corr = []; 
training_day_num = [];
for ii = 1:length(struct_files)
    corrEarlyRatio    = [corrEarlyRatio, meta_struct.(['i', struct_files(ii).name(1:12)]).curr_struct.corr_early_dfoverf_ratio];
    peak_percent_corr = [peak_percent_corr, meta_struct.(['i', struct_files(ii).name(1:12)]).curr_struct.peak_percent_corr];
    percent_corr      = [percent_corr, meta_struct.(['i', struct_files(ii).name(1:12)]).curr_struct.percent_correct];
    recent_corr       = [recent_corr, meta_struct.(['i', struct_files(ii).name(1:12)]).curr_struct.recent_percent_corr];
    training_day_num  = [training_day_num, meta_struct.(['i', struct_files(ii).name(1:12)]).curr_struct.training_day_num];
end

%plot ratio as a function of percent correct
figure; hold on;
for ii = 1:length(corrEarlyRatio)
    plot(corrEarlyRatio(ii), percent_corr(ii),'Marker', '*', 'Linestyle', 'none', 'Color', color_palette(ii,:));
end
xlim([0, 3]);
ylim([0, 100]);
title('fluorescence response ratio and % correct');
xlabel('ratio of fluorescence response correct/early');
ylabel('Percent of trials correct on the day of imaging');
legend(days_chrono_order(4:end));

%plot rstruct_diratio as a function of peak percent correct
figure; hold on;
for ii = 1:length(corrEarlyRatio)
    plot(corrEarlyRatio(ii), peak_percent_corr(ii),'Marker', '*', 'Linestyle', 'none', 'Color', color_palette(ii,:));
end
xlim([0, 3]);
ylim([0, 100]);
title('fluorescence response ratio and  peak % correct');
xlabel('ratio of fluorescence response correct/early');
ylabel('peak % correct throughout training history');
legend(days_chrono_order(4:end));

%Plot ratio vs total training days
figure; hold on;
for ii = 1:length(corrEarlyRatio)
    plot(corrEarlyRatio(ii), training_day_num(ii),'Marker', '*', 'Linestyle', 'none', 'Color', color_palette(ii,:));
end
xlim([1, 3]);
ylim([15, 55]);
title('fluorescence response ratio and training day number');
xlabel('ratio of fluorescence response correct/early');
ylabel('Number of training days before imaging');

%Plot ratio vs recent history 
figure; hold on;
for ii = 1:length(corrEarlyRatio)
    plot(corrEarlyRatio(ii), recent_corr(ii),'Marker', '*', 'Linestyle', 'none', 'Color', color_palette(ii,:));
end
xlim([1, 3]);
ylim([0, 100]);
title('fluorescence response ratio and recent percent correct');
xlabel('ratio of fluorescence response correct/early');
ylabel('percent correct for previous three days');

%% variance in peak time, onset, latency and FRatio as a function of reaction times
corr_react_times = [];
corr_onset_latency = []; 
corr_magnitude = [];
corr_peak_times = [];
corr_hold_dur = []; 
early_magnitude = [];
early_hold_dur = [];
corr_react_times    = [corr_react_times, meta_struct.(['i', struct_files(ii).name(1:12)]).curr_struct.corr_react_times];
corr_onset_latency  = [corr_onset_latency, meta_struct.(['i', struct_files(ii).name(1:12)]).curr_struct.corr_onset_latency];
corr_magnitude      = [corr_magnitude, meta_struct.(['i', struct_files(ii).name(1:12)]).curr_struct.corr_magnitude];
corr_peak_times     = [corr_peak_times, meta_struct.(['i', struct_files(ii).name(1:12)]).curr_struct.corr_peak_latency_release];
corr_hold_dur       = [corr_hold_dur, meta_struct.(['i', struct_files(ii).name(1:12)]).curr_struct.corr_hold_dur];
early_magnitude     = [early_magnitude, meta_struct.(['i', struct_files(ii).name(1:12)]).curr_struct.early_magnitude];
early_hold_dur      = [early_hold_dur, meta_struct.(['i', struct_files(ii).name(1:12)]).curr_struct.early_hold_dur];

%=========================================================================================
%correct trials: reaction times vs onset latency
%figure; hold on;
    for ii = 1:length(days_chrono_order);
        figure;
        corr_react_times    = meta_struct.(['i', struct_files(ii).name(1:12)]).curr_struct.corr_react_times;
        corr_onset_latency  = meta_struct.(['i', struct_files(ii).name(1:12)]).curr_struct.corr_onset_latency;
        %sample_trial_nums = round(linspace(1,length(corr_react_times), 25));
        %plot(corr_react_times(sample_trial_nums), corr_onset_latency(sample_trial_nums), 'Marker', 'o', 'Linestyle', 'none', 'Color', color_palette(ii,:));
        plot(corr_react_times, corr_onset_latency, 'Marker', 'o', 'Linestyle', 'none');
        xlabel('reaction times for correct trials (ms)')
        ylabel('onset latency of F response for correct trials (relative to release)');
        title([days_chrono_order(ii), 'correct trials reaction times vs onset latency']);
        xlim([100 1000]);
        ylim([-600 500])
    end
    xlabel('reaction times for correct trials (ms)')
    ylabel('onset latency of F response for correct trials (relative to release)');
    title('correct trials reaction times vs onset latency');
    
%=========================================================================================
% corrects: reaction times versus F magnitude on corrects
figure; hold on;
    for ii = 1:length(days_chrono_order);
        corr_react_times    = meta_struct.(['i', struct_files(ii).name(1:12)]).curr_struct.corr_react_times;
        corr_magnitude  = meta_struct.(['i', struct_files(ii).name(1:12)]).curr_struct.corr_magnitude;
        sample_trial_nums = round(linspace(1,length(corr_react_times), 25));
        plot(corr_react_times(sample_trial_nums), corr_magnitude(sample_trial_nums), 'Marker', 'o', 'Linestyle', 'none', 'Color', color_palette(ii,:));
    end
xlabel('reaction times for correct trials (ms)')
ylabel('magnitude of F response for correct trials (relative to release)');
title('correct trials reaction times vs magnitude of F response');

%corrects: reaction times versus peak times
figure; hold on;
    for ii = 1:length(days_chrono_order);
        corr_react_times    = meta_struct.(['i', struct_files(ii).name(1:12)]).curr_struct.corr_react_times;
        corr_peak_times  = meta_struct.(['i', struct_files(ii).name(1:12)]).curr_struct.corr_peak_latency_release;
        sample_trial_nums = round(linspace(1,length(corr_react_times), 25));
        plot(corr_react_times(sample_trial_nums), corr_peak_times(sample_trial_nums), 'Marker', 'o', 'Linestyle', 'none', 'Color', color_palette(ii,:));
    end
xlabel('reaction times for correct trials (ms)')
ylabel('latency to peak F for correct trials (relative to release)');
title('correct trials reaction times vs latency to peak');

%CUE aligned corrects: reaction times vs lat to peak
figure; hold on;
    for ii = 1:length(days_chrono_order);
        corr_react_times    = meta_struct.(['i', struct_files(ii).name(1:12)]).curr_struct.corr_react_times;
        corr_cue_aligned_peak_latency  = meta_struct.(['i', struct_files(ii).name(1:12)]).curr_struct.corr_cue_aligned_peak_latency;
        sample_trial_nums = round(linspace(1,length(corr_react_times), 25));
        plot(corr_react_times(sample_trial_nums), corr_cue_aligned_peak_latency(sample_trial_nums), 'Marker', 'o', 'Linestyle', 'none', 'Color', color_palette(ii,:));
    end
xlabel('reaction times for correct trials (ms)')
ylabel('latency to peak F for correct trials (relative to CUE)');
title('Cue aligned correct trials reaction times vs latency to peak');
ylim([50 750]);

%CUE aligned corrects: reaction times vs onset latency
figure; hold on;
    for ii = 1:length(days_chrono_order);
        corr_react_times    = meta_struct.(['i', struct_files(ii).name(1:12)]).curr_struct.corr_react_times;
        corr_cue_aligned_onset_latency  = meta_struct.(['i', struct_files(ii).name(1:12)]).curr_struct.corr_cue_aligned_onset_latency;
        sample_trial_nums = round(linspace(1,length(corr_react_times), 25));
        plot(corr_react_times(sample_trial_nums), corr_cue_aligned_onset_latency(sample_trial_nums), 'Marker', 'o', 'Linestyle', 'none', 'Color', color_palette(ii,:));
    end
xlabel('reaction times for correct trials (ms)')
ylabel('onset latency of F response for correct trials (relative to CUE)');
title('Cue aligned correct trials reaction times vs onset latency');
ylim([-600 800])

%corrects: reaction times versus rise times
figure; hold on;
    for ii = 1:length(days_chrono_order);
        corr_react_times    = meta_struct.(['i', struct_files(ii).name(1:12)]).curr_struct.corr_react_times;
        corr_rise_times  = meta_struct.(['i', struct_files(ii).name(1:12)]).curr_struct.corr_rise_times;
        sample_trial_nums = round(linspace(1,length(corr_react_times), 25));
        plot(corr_react_times(sample_trial_nums), corr_rise_times(sample_trial_nums), 'Marker', 'o', 'Linestyle', 'none', 'Color', color_palette(ii,:));
    end
xlabel('reaction times for correct trials (ms)')
ylabel('rise times correct trials');
title('correct trials reaction times vs rise times');
xlim([100 900]);
ylim([0 400]);

%correct trials: hold duration and magnitude
figure; hold on;
    for ii = 1:length(days_chrono_order);
        corr_hold_dur    = meta_struct.(['i', struct_files(ii).name(1:12)]).curr_struct.corr_hold_dur;
        corr_magnitude  = meta_struct.(['i', struct_files(ii).name(1:12)]).curr_struct.corr_magnitude;
        sample_trial_nums = round(linspace(1,length(corr_hold_dur), 25));
        plot(corr_hold_dur(sample_trial_nums), corr_magnitude(sample_trial_nums), 'Marker', 'o', 'Linestyle', 'none', 'Color', color_palette(ii,:));
    end
xlabel('hold duration for correct trials (ms)')
ylabel('magnitude of F response for correct trials (relative to release)');
title('correct trials hold duration vs magnitude of F response');

%% rise time for corr v early 
figure; hold on;
    for ii = 1:length(days_chrono_order);
        corr_rise_times    = meta_struct.(['i', struct_files(ii).name(1:12)]).curr_struct.corr_rise_times;
        early_rise_times  = meta_struct.(['i', struct_files(ii).name(1:12)]).curr_struct.early_rise_times;
        sample_trial_nums = round(linspace(1,length(corr_rise_times), 17));
        early_sample_trial_nums = round(linspace(1,length(early_rise_times), 17));
        plot(corr_rise_times(sample_trial_nums), early_rise_times(early_sample_trial_nums), 'Marker', 'o', 'Linestyle', 'none', 'Color', color_palette(ii,:));
    end
    xlim([0 400]);
    ylim([0 400]);
    xlabel('correct trials rise times');
    ylabel('early trials rise times');
    title('rise times for correct vs early trials');
    plot(x,y,'Color','k');
    
%% mean rise times corr v early
corr_rise_time_mean = [];
early_rise_time_mean = [];
for ii = 1:length(struct_files)
    corr_rise_time_mean = [corr_rise_time_mean, meta_struct.(['i', struct_files(ii).name(1:12)]).curr_struct.corr_rise_time_mean];
    early_rise_time_mean = [early_rise_time_mean, meta_struct.(['i', struct_files(ii).name(1:12)]).curr_struct.early_rise_time_mean];
end

%Plot mean rise times for Corr v early
figure; hold on;
for ii = 1: length(corr_rise_time_mean)
    plot(corr_rise_time_mean(ii), early_rise_time_mean(ii),'Marker', 'o', 'Linestyle', 'none', 'Color', color_palette(ii,:));
end
hold on; plot(x,y, 'Color', 'k', 'Linestyle', '--');
xlim([0, 400]);
ylim([0, 400]);
title('Mean rise times Correct vs Early');
xlabel('Correct mean rise times');
ylabel('Early mean rise times');
legend(days_chrono_order);

%% find and plot mean react time for each animal (and sem. Then calculate a
%grand mean and 
corr_react_time_mean = [];
corr_react_time_std = [];
for ii = 1:length(struct_files)
    corr_react_time_mean = [corr_react_time_mean, meta_struct.(['i', struct_files(ii).name(1:12)]).curr_struct.corr_react_time_mean];
    corr_react_time_std = [corr_react_time_std, meta_struct.(['i', struct_files(ii).name(1:12)]).curr_struct.corr_react_time_std];
end
%plot mean react times. 
figure; hold on; 
for ii = 1: length(corr_react_time_mean)
    bar([1:length(corr_react_time_mean)], corr_react_time_mean);
    errorbar([1:length(corr_react_time_mean)], corr_react_time_mean, corr_react_time_std, 'Linestyle', 'none');
end
mean(corr_react_time_mean);
mean(corr_react_time_std);
