%% WF_CRP_licking_summaries
%across animal RT summary
clear;
WF_CRP_list_of_days;
data_dir = 'Z:\Analysis\Cue_reward_pairing_analysis\CRPFigureFolder\across_animals_RT\';
avg_RTs_all = cell(size(all_bx_animals));
avg_sem_all = cell(size(all_bx_animals));
post_cue_lick_all = cell(size(all_bx_animals));
pre_cue_lick_all = cell(size(all_bx_animals));
post_cue_sem_all = cell(size(all_bx_animals));
pre_cue_sem_all = cell(size(all_bx_animals));
miss_rates_all = cell(size(all_bx_animals));
RTstd_all = cell(size(all_bx_animals));
%rew_om_hist_all

figure;
hold on;
x_axis = 1:1:max(days_to_post_learning);
for ii = 1:length(all_bx_animals)
    load([data_dir, all_bx_animals{ii}, 'RT_sem_across_days'], 'RT_across_days_sem', 'RT_across_days');
    avg_sem_all{ii} = RT_across_days_sem; 
    avg_RTs_all{ii} =  RT_across_days;
    load(['Z:\Analysis\Cue_reward_pairing_analysis\CRPFigureFolder\across_animals_post_cue_lick\', all_bx_animals{ii}, 'rew_om_post_cue_lick']);
    post_cue_lick_all{ii} = avg_licks_post_cue;
    pre_cue_lick_all{ii} = avg_licks_pre_cue;
    post_cue_sem_all{ii} = avg_licks_post_cue_sem;
    pre_cue_sem_all{ii} = avg_licks_pre_cue_sem;
    
    load(['Z:\Analysis\Cue_reward_pairing_analysis\CRPFigureFolder\across_animals_RT_stats\', all_bx_animals{ii}, 'RTstd_FA_misses']);
    miss_rates_all{ii} = miss_rates;
    RTstd_all{ii} = std_of_RT_across_days;
    
    plot(x_axis(1:days_to_post_learning(ii)), RT_across_days(1:days_to_post_learning(ii)), 'k');
    ylabel('RT relative to cue onset');
    xlabel('training day #');
    title('mean reaction time across days for all animals');
    hline(600, 'r');
end


%% make a grand mean RT across days and across animals ending on the second day of imaging 
curr_day_mean = [];
grand_RT_mean = [];
grand_RT_sem = [];
for day_num = 1:max(days_to_post_learning);
    for ii = 1:length(all_bx_animals)
        if day_num <= days_to_post_learning(ii)
            curr_day_mean = [curr_day_mean, avg_RTs_all{ii}(day_num)];
        end
    end
    curr_day_sem = std(curr_day_mean)/sqrt(length(curr_day_mean));
    curr_day_mean = mean(curr_day_mean);
    grand_RT_mean = [grand_RT_mean , curr_day_mean];
    grand_RT_sem = [grand_RT_sem , curr_day_sem];
end

figure;
errorbar(x_axis, grand_RT_mean, grand_RT_sem);
ylabel('RT relative to cue onset. Reward at 600ms');
xlabel('training day #');
title('mean reaction time across days for all animals: includes days from 1 to post learning of 500ms delay');
hline(600);

%% check on RTs on post learning day and the day before that in order to generate criterion for bx. 

day_1_RT=[];
last_training_day_RT=[];
PL_RT = [];
for ii = 1:length(all_bx_animals)
    day_1_RT = [day_1_RT, avg_RTs_all{ii}(1)];
    last_training_day_RT = [last_training_day_RT, avg_RTs_all{ii}(days_to_post_learning(ii)-1)];
    PL_RT = [PL_RT, avg_RTs_all{ii}(days_to_post_learning(ii))];
end

day_1_RT_mean = mean(day_1_RT);
last_training_day_RT_mean = mean(last_training_day_RT);
PL_RT_mean = mean(PL_RT);

day_1_RT_sem = std(day_1_RT)/sqrt(length(day_1_RT));
last_training_day_RT_sem = std(last_training_day_RT)/sqrt(length(last_training_day_RT));
PL_RT_sem = std(PL_RT)/sqrt(length(PL_RT));

figure;
bar([1:3], [day_1_RT_mean, last_training_day_RT_mean, PL_RT_mean]); hold on;
errorbar([1:3], [day_1_RT_mean, last_training_day_RT_mean, PL_RT_mean], [day_1_RT_sem, last_training_day_RT_sem, PL_RT_sem], 'LineStyle', 'none'); 
title('mean RT across animals for day 1, last day of bx training, and second day of imaging');
ylabel('RT');
hline(600);
round(max(last_training_day_RT))

%% plot the post cue lick rate on omission trials before and after learning across animals 

day_1_pre_cue=[];
PL_pre_cue = [];
day_1_post_cue=[];
PL_post_cue = [];

day_1_pre_cue_sem=[];
PL_pre_cue_sem = [];
day_1_post_cue_sem=[];
PL_post_cue_sem = [];
for ii = 1:length(all_bx_animals)
    if ~isnan(pre_cue_lick_all{ii}(1)) && ~isnan(pre_cue_lick_all{ii}(days_to_post_learning(ii)));
        day_1_pre_cue = [day_1_pre_cue, pre_cue_lick_all{ii}(1)];
        PL_pre_cue = [PL_pre_cue, pre_cue_lick_all{ii}(days_to_post_learning(ii))];
        day_1_post_cue = [day_1_post_cue, post_cue_lick_all{ii}(1)];
        PL_post_cue = [PL_post_cue, post_cue_lick_all{ii}(days_to_post_learning(ii))];
    else
        day_1_pre_cue = [day_1_pre_cue, NaN];
        PL_pre_cue = [PL_pre_cue, NaN];
        day_1_post_cue = [day_1_post_cue, NaN];
        PL_post_cue = [PL_post_cue, NaN];
    end
end


day_1_pre_cue_mean = nanmean(day_1_pre_cue*2);
day_1_pre_cue_sem = nanstd(day_1_pre_cue*2)/sqrt(sum(~isnan(day_1_pre_cue)));

PL_pre_cue_mean = nanmean(PL_pre_cue*2);
PL_pre_cue_sem = nanstd(PL_pre_cue*2)/sqrt(sum(~isnan(PL_pre_cue)));

day_1_post_cue_mean = nanmean(day_1_post_cue*2);
day_1_post_cue_sem = nanstd(day_1_post_cue*2)/sqrt(sum(~isnan(day_1_post_cue)));

PL_post_cue_mean = nanmean(PL_post_cue*2);
PL_post_cue_sem = nanstd(PL_post_cue*2)/sqrt(sum(~isnan(PL_post_cue)));

figure;  hold on;
x_axis= [1:1:4];
bar(x_axis, [day_1_pre_cue_mean, PL_pre_cue_mean, day_1_post_cue_mean, PL_post_cue_mean]);
errorbar(x_axis, [day_1_pre_cue_mean, PL_pre_cue_mean, day_1_post_cue_mean, PL_post_cue_mean], [day_1_pre_cue_sem, PL_pre_cue_sem, day_1_post_cue_sem, PL_post_cue_sem], '.');
ylabel('lick rate Hz');
xlabel('day1 pre-cue lick rate,   day1 post-cue lick rate,   post-learning pre-cue lick rate,   post-learning post-cue lick rate');

%% plot the RTstd and miss rate across animals across days
%collect the data
day_1_MR=[];
PL_MR = [];
day_1_RTstd=[];
PL_RTstd = [];
for ii = 1:length(all_bx_animals)
        day_1_MR = [day_1_MR, miss_rates_all{ii}(1)];
        PL_MR = [PL_MR, miss_rates_all{ii}(days_to_post_learning(ii))];
        
        day_1_RTstd = [day_1_RTstd, RTstd_all{ii}(1)];
        PL_RTstd = [PL_RTstd, RTstd_all{ii}(days_to_post_learning(ii))];
end
 
%plot miss rates across days and before vs after learning
day_1_MR_mean = mean(day_1_MR);
day_1_MR_sem = std(day_1_MR)/sqrt(length(day_1_MR));
PL_MR__mean = mean(PL_MR);
PL_MR_sem = std(PL_MR)/sqrt(length(PL_MR));
figure; 
subplot(1,3,3); hold on;
x_axis= [1:1:2];
bar(x_axis, [day_1_MR_mean, PL_MR__mean]);
errorbar(x_axis, [day_1_MR_mean, PL_MR__mean], [day_1_MR_sem, PL_MR_sem], '.');
ylabel('miss rate as a percent of all trials');
xlabel('Left: day 1       Right: post-learning');
title('miss rate before and after learning');
suptitle('miss rate for 500ms delay training days (1000ms+ RT)');
%plot all miss rates
for day_num =1:max(days_to_post_learning)
    this_day_miss_rate = [];
    for ii = 1:length(miss_rates_all)
    if day_num <= days_to_post_learning(ii);
        this_day_miss_rate = [this_day_miss_rate, miss_rates_all{ii}(day_num)];
    end
    end
    [miss_rates_mean(day_num), miss_rates_sem(day_num)] = get_mean_and_sem(this_day_miss_rate);
end
subplot(1,3,[1:2]);
x_axis2 = 1:max(days_to_post_learning);
errorbar(x_axis2, miss_rates_mean, miss_rates_sem);
ylabel('miss rate'); xlabel('training day number'); title('miss rate as a percent of all trials over days');


%plot RTstd across days and before vs after learning
day_1_RTstd_mean = mean(day_1_RTstd);
day_1_RTstd_sem = std(day_1_RTstd)/sqrt(length(day_1_RTstd));
PL_RTstd_mean = mean(PL_RTstd);
PL_RTstd_sem = std(PL_RTstd)/sqrt(length(PL_RTstd));
figure;
subplot(1,3,3); hold on;
bar(x_axis, [day_1_RTstd_mean, PL_RTstd_mean]);
errorbar(x_axis, [day_1_RTstd_mean, PL_RTstd_mean], [day_1_RTstd_sem, PL_RTstd_sem], '.');
ylabel('RT std dev');
xlabel('Left: day 1       Right: post-learning');
title('RT standard deviation before and after learning');
suptitle('RT st.dev for 500ms delay training days');
%plot all RTstd
for day_num =1:max(days_to_post_learning)
    this_day_RTstd = [];
    for ii = 1:length(RTstd_all)
    if day_num <= days_to_post_learning(ii);
        this_day_RTstd = [this_day_RTstd, RTstd_all{ii}(day_num)];
    end
    end
    [RTstd_mean(day_num), RTstd_sem(day_num)] = get_mean_and_sem(this_day_RTstd);
end
subplot(1,3,[1:2]);
x_axis2 = 1:max(days_to_post_learning);
errorbar(x_axis2, RTstd_mean, RTstd_sem);
ylabel('RT standard deviation'); xlabel('training day number'); title('mean st. dev of the RT across animals all days');












