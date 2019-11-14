%script to calculate spike rate and variance for each cell. 

%% PART ONE - run this segment first in order to collect all the necessary spike related variables for aav and transgenic mice

%set initialization variables
id = 2; %look at post learning condition first since it likely has most homogeneous behavior. 
aav_mouse = 0;
new_mice=0;
doCollectVars = 1;

%set directories
if doCollectVars
output_dir = 'Y:\home\jake\Analysis\Cue_reward_pairing_analysis\CC_analysis_2P\spike rate and CV comparison';
if aav_mouse 
    CRP_expt_list_all;
    exp_nums = [1:10];
    data_dir = 'Y:\home\lindsey\Analysis\2P\Jake';
else
    if new_mice
        CRP_OT_expt_list_LS_jh;
        exp_nums = [1:length(expt(id).mouse)];
        data_dir =  'Y:\home\jake\Analysis\Cue_reward_pairing_analysis\CC_analysis_2P';
    else
        CRP_expt_list_all;
        exp_nums = [11:15];
        data_dir = 'Y:\home\lindsey\Analysis\2P\Jake';
    end
end
aav_list = {};
transgen_list = {};

%extract spike info data from each dataset
for iexp =  exp_nums; 
    mouse = strtrim(expt(id).mouse(iexp,:));
    date = expt(id).date(iexp,:);
    run = expt(id).run(iexp,:);
    fprintf([date ' ' mouse '\n'])
    img_fn = [date '_' mouse];
    if aav_mouse
        aav_list{end+1} = img_fn;
    else 
        transgen_list{end+1} = img_fn;
    end
    
    %load data  _targetAlign
    load(fullfile(data_dir, img_fn, [img_fn '_targetAlign.mat']), 'targetAlign_events', 'frameRateHz', 'prewin_frames', 'postwin_frames');
    
    %calculate variables
    num_trials = size(targetAlign_events,3);
    format short;
    IFI=round((1000/frameRateHz),1); %interframe interval in ms
    times_included =  [((prewin_frames+postwin_frames+1)*IFI)-2000:((prewin_frames+postwin_frames+1)*IFI)]; [1:((prewin_frames+postwin_frames+1)*IFI)];%  %[1:1500];%
    tot_time = length(times_included)*num_trials; %total time measured for each cell
    
    %allocate memory 
    animal_ISI_mean = [];
    animal_ISI_std = [];
    animal_ISI_CV = [];
    animal_ISI_rate = [];
    
    %remove NaNs from sweeps
    targetAlign_events(isnan(targetAlign_events)) = 0;
    
    for this_cell = 1:size(targetAlign_events,2)
        
        %find the spike times for each trial for this cell
        this_cell_events = squeeze(targetAlign_events(:,this_cell,:)); %dim1=frame# dim2=trial#
        max_spikes = max(sum(this_cell_events,2)); %largets # of spikes in any sweep this cell
        this_cell_spikeT = NaN(max_spikes, num_trials);
        for this_trial = 1:num_trials
            this_cell_spikeT([1:sum(this_cell_events(:,this_trial))], this_trial) =  find(this_cell_events(:,this_trial)==1)*IFI;
        end
        
        %calculate rate and CV of spike times
        [cell_ISI_mean, cell_ISI_std, cell_ISI_CV, ~] = calc_spike_var(times_included,this_cell_spikeT);
        cell_ISI_rate = (nansum(nansum(this_cell_events)) / tot_time ) * 1000;
        
        animal_ISI_mean = [animal_ISI_mean, cell_ISI_mean];
        animal_ISI_std = [animal_ISI_std, cell_ISI_std];
        animal_ISI_CV = [animal_ISI_CV, cell_ISI_CV];
        animal_ISI_rate = [animal_ISI_rate, cell_ISI_rate];
        
        %save variables for this animal
        save(fullfile(output_dir, [img_fn '_spikeInfo.mat']), 'animal_ISI_mean', 'animal_ISI_std', 'animal_ISI_CV', 'animal_ISI_rate');
    end
end

%save list of sessions collected
if aav_mouse
    save(fullfile(output_dir, ['id', num2str(id), '_aav_list.mat']), 'aav_list');
else
    save(fullfile(output_dir, ['id', num2str(id), '_transgen_list.mat']), 'transgen_list');
end
end

%% PART TWO - run this segment after having saved all the variables from part 1 
%makes comparisons between aav and transgenic mice

%set directory and load lists
output_dir = 'Y:\home\jake\Analysis\Cue_reward_pairing_analysis\CC_analysis_2P\spike rate and CV comparison';
id = 2; 
load(fullfile(output_dir, ['id', num2str(id), '_aav_list.mat']));
load(fullfile(output_dir, ['id', num2str(id), '_transgen_list.mat']));

%allocate memory
aav_ISI_mean = [];
aav_ISI_std = [];
aav_ISI_CV = [];
aav_ISI_rate = [];
transgen_ISI_mean = [];
transgen_ISI_std = [];
transgen_ISI_CV = [];
transgen_ISI_rate = [];

%load ISI variables
for this_expt = 1:length(aav_list)
    img_fn = aav_list{this_expt};
    load(fullfile(output_dir, [img_fn '_spikeInfo.mat']));
    aav_ISI_mean = [aav_ISI_mean, animal_ISI_mean];
    aav_ISI_std = [aav_ISI_std, animal_ISI_std];
    aav_ISI_CV = [aav_ISI_CV, animal_ISI_CV];
    aav_ISI_rate = [aav_ISI_rate, animal_ISI_rate];
end
for this_expt = 1:length(transgen_list)
    img_fn = transgen_list{this_expt};
    load(fullfile(output_dir, [img_fn '_spikeInfo.mat']));
    transgen_ISI_mean = [transgen_ISI_mean, animal_ISI_mean];
    transgen_ISI_std = [transgen_ISI_std, animal_ISI_std];
    transgen_ISI_CV = [transgen_ISI_CV, animal_ISI_CV];
    transgen_ISI_rate = [transgen_ISI_rate, animal_ISI_rate];
end

%plot comparisons of ISI variables
figure;
set(gcf, 'Position', [100,100,1200,600]);

subplot(1,4,1);
bar([mean(aav_ISI_mean), mean(transgen_ISI_mean)],'b'); hold on;
errorbar([1,2],   [mean(aav_ISI_mean), mean(transgen_ISI_mean)],    [std(aav_ISI_mean)/(sqrt(length(aav_ISI_mean))), std(transgen_ISI_mean)/(sqrt(length(transgen_ISI_mean)))],   'c')
ylabel('Inter-Spike Interval');
xlabel('AAV   Ai148');
[hh, pp] = ttest2([aav_ISI_mean], [transgen_ISI_mean]);
title(['ttest: sig=', num2str(hh), '  p=', num2str(pp)]);
ylim([0 1200]);

subplot(1,4,2);
bar([mean(aav_ISI_std), mean(transgen_ISI_std)],'b'); hold on;
errorbar([1,2],   [mean(aav_ISI_std), mean(transgen_ISI_std)],    [std(aav_ISI_std)/(sqrt(length(aav_ISI_std))), std(transgen_ISI_std)/(sqrt(length(transgen_ISI_std)))],   'c')
ylabel('std of Inter-Spike Interval');
xlabel('AAV   Ai148');
[hh, pp] = ttest2([aav_ISI_std], [transgen_ISI_std]);
title(['ttest: sig=', num2str(hh), '  p=', num2str(pp)]);
ylim([0 1200]);

subplot(1,4,3);
bar([mean(aav_ISI_CV), mean(transgen_ISI_CV)],'b'); hold on;
errorbar([1,2],   [mean(aav_ISI_CV), mean(transgen_ISI_CV)],    [std(aav_ISI_CV)/(sqrt(length(aav_ISI_CV))), std(transgen_ISI_CV)/(sqrt(length(transgen_ISI_CV)))],   'c')
ylabel('CV of ISI');
xlabel('AAV   Ai148');
[hh, pp] = ttest2([aav_ISI_CV], [transgen_ISI_CV]);
title(['ttest: sig=', num2str(hh), '  p=', num2str(pp)]);

subplot(1,4,4);
bar([mean(aav_ISI_rate), mean(transgen_ISI_rate)],'b'); hold on;
errorbar([1,2],   [mean(aav_ISI_rate), mean(transgen_ISI_rate)],    [std(aav_ISI_rate)/(sqrt(length(aav_ISI_rate))), std(transgen_ISI_rate)/(sqrt(length(transgen_ISI_rate)))],   'c')
ylabel('Spike rate');
xlabel('AAV   Ai148');
[hh, pp] = ttest2([aav_ISI_rate], [transgen_ISI_rate]);
title(['ttest: sig=', num2str(hh), '  p=', num2str(pp)]);

suptitle(['aav: ', num2str(length(aav_ISI_mean)), 'neurons, ', num2str(length(aav_list)), 'animals.   Ai148: ', num2str(length(transgen_ISI_mean)), 'neurons, ', num2str(length(transgen_list)), 'animals.']);



