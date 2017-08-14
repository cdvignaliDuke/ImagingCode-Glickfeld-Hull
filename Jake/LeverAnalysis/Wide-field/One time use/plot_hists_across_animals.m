
%WF_CRP_list_of_days;
clear;
days_1 = {'170321_img86',      [],           '170408_img88', '170513_img89', '170416_img90',  '170417_img91', '170420_img92', '170510_img93', '170524_img94', '170605_img95', '170612_img96', '170628_img98', '170705_img99'};  %'170408_img87' first day of imaging but animal quit licking part way through
days_post = {'170325_img86', '170419_img87', '170415_img88', '170519_img89', '170426_img90',  '170425_img91', '170428_img92', '170518_img93', '170529_img94', '170614_img95', '170622_img96', '170704_img98', '170710_img99'}; %'170326_img86',  are also 10%omission post learning
data_dir = 'Z:\Analysis\Cue_reward_pairing_analysis\CRPFigureFolder\across_animals_hist_data\';

useable_days_inx = [];
lick_hist_OR_all = [];
lick_hist_OR_N_all = [];
for ii = 1:length(days_1)
    if ii == 3
        continue
    else
        if ~isempty(days_1{ii})
            if exist([data_dir, days_1{ii}, 'rew_hist.mat'], 'file')
                load([data_dir, days_1{ii}, 'rew_hist']);
                lick_hist_OR_all = [lick_hist_OR_all; full_trial_licks_rewarded_bin];
                useable_days_inx = [useable_days_inx, ii];
            end
        end
    end
end
[lick_hist_OR_all_mean, lick_hist_OR_all_sem] = get_mean_and_sem(lick_hist_OR_all);
num_animals = size(lick_hist_OR_all, 1);
x_axis = ([1:length(lick_hist_OR_all_mean)]-180)*100;
figure; subplot(2,1,1);
bar(x_axis, lick_hist_OR_all_mean); hold on;
errorbar(x_axis, lick_hist_OR_all_mean, lick_hist_OR_all_sem, '.', 'CapSize', 1);
title('Rewarded trials: Day 1');
xlabel('time in ms relative to cue onset');
ylabel('lick rate (Hz)');
ylim([0 10]); xlim([-18100 6200]);
vline(600, 'g')

for ii = 1:length(days_post(useable_days_inx))
    if exist([data_dir, days_post{ii}, 'rew_hist.mat'], 'file')
            load([data_dir, days_post{ii}, 'rew_hist']);
            lick_hist_OR_N_all = [lick_hist_OR_N_all; full_trial_licks_rewarded_bin];
    end
end
[lick_hist_OR_N_all_mean,lick_hist_OR_N_all_sem] = get_mean_and_sem(lick_hist_OR_N_all);
num_animals_N = size(lick_hist_OR_N_all, 1);
assert(size(lick_hist_OR_N_all,1) == size(lick_hist_OR_all, 1));
subplot(2,1,2); bar(x_axis, lick_hist_OR_N_all_mean); hold on;
errorbar(x_axis, lick_hist_OR_N_all_mean, lick_hist_OR_N_all_sem, '.',  'CapSize', 1);
title('Rewarded trials: day N');
xlabel('time in ms relative to cue onset');
ylabel('lick rate (Hz)');
ylim([0 10]); xlim([-18100 6200]);
vline(600, 'g');



%img86 img89 img90 img91 img92 img93 img94 img95 img96 img98 img99
%%animals used in histogram OR trials





% days_img86 = {'170321_img86', '170322_img86', '170323_img86', '170325_img86', '170326_img86', '170327_img86', '170330_img86'};
% days_img87 = {'170408_img87', '170410_img87', '170411_img87', '170412_img87', '170413_img87', '170414_img87', '170415_img87', '170417_img87', '170418_img87', '170419_img87', '170420_img87', '170422_img87', '170423_img87', '170424_img87', '170425_img87', '170426_img87'};
% days_img88 = {'170408_img88', '170410_img88', '170411_img88', '170412_img88', '170413_img88', '170414_img88', '170415_img88', '170417_img88', '170418_img88', '170419_img88', '170420_img88', '170421_img88', '170422_img88', '170423_img88', '170424_img88', '170425_img88', '170426_img88'};
% days_img89 = {'170513_img89', '170515_img89', '170516_img89', '170517_img89', '170518_img89', '170519_img89', '170522_img89', '170523_img89', '170524_img89', '170525_img89', '170526_img89', '170527_img89', '170529_img89', '170530_img89', '170531_img89'};
% days_img90 = {'170416_img90', '170417_img90', '170418_img90', '170419_img90', '170422_img90', '170423_img90', '170424_img90', '170425_img90', '170426_img90', '170427_img90', '170428_img90', '170429_img90', '170501_img90', '170502_img90', '170503_img90', '170504_img90', '170505_img90', '170506_img90', '170508_img90', '170509_img90', '170510_img90', '170511_img90', '170512_img90', '170513_img90', '170515_img90'};
% days_img91 = {'170417_img91', '170418_img91', '170419_img91', '170420_img91', '170422_img91', '170423_img91', '170424_img91', '170425_img91', '170426_img91', '170427_img91', '170428_img91', '170429_img91', '170501_img91', '170502_img91', '170503_img91', '170504_img91'};
% days_img92 = {'170420_img92', '170422_img92', '170423_img92', '170424_img92', '170425_img92', '170426_img92', '170427_img92', '170428_img92', '170429_img92', '170501_img92', '170503_img92', '170504_img92', '170505_img92', '170506_img92', '170509_img92'};
% days_img93 = {'170510_img93', '170511_img93', '170512_img93', '170513_img93', '170515_img93', '170516_img93', '170517_img93', '170518_img93', '170519_img93', '170520_img93', '170522_img93', '170523_img93', '170524_img93', '170525_img93'};
% days_img94 = {'170524_img94', '170525_img94', '170526_img94', '170527_img94', '170529_img94', '170530_img94', '170531_img94', '170601_img94', '170602_img94', '170604_img94', '170605_img94', '170606_img94'};
% days_img95 = {'170605_img95', '170606_img95', '170607_img95', '170608_img95', '170609_img95', '170610_img95', '170611_img95', '170612_img95', '170613_img95', '170614_img95'};
% days_img96 = {'170612_img96', '170613_img96', '170614_img96', '170615_img96', '170620_img96', '170621_img96', '170622_img96', '170623_img96', '170624_img96'};
% days_img98 = {'170628_img98', '170629_img98', '170701_img98', '170703_img98', '170704_img98', '170705_img98', '170706_img98'};
% days_img99 = {'170705_img99', '170706_img99', '170707_img99', '170708_img99', '170710_img99', '170711_img99', '170712_img99'};

