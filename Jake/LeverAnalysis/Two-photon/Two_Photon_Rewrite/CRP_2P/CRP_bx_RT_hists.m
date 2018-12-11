% 
% clear;
% days = {'170416_img90',  '170417_img91', '170420_img92', '170510_img93', '170513_img89', '170524_img94', '170911_img036', '170921_img044', ...
%     '171113_img050', '171227_img067', '180104_img070', '180322_img077', '180507_img081', '180425_img084', '180509_img085'}; 
% WF_CRP_licking_only_overview
% clear;
% days = {'170426_img90',  '170425_img91', '170428_img92', '170518_img93', '170519_img89', '170529_img94', '170926_img036', '170926_img044', ...
%     '171122_img050', '180104_img067', '180108_img070',  '180403_img077', '180518_img081', '180502_img084', '180517_img085'};
% WF_CRP_licking_only_overview
% clear;
% days ={'170506_img90', '170426_img91', '170429_img92', '170519_img93', '170522_img89', '170530_img94'};
% WF_CRP_licking_only_overview
% clear;
% days = {'170508_img90', '170427_img91', '170501_img92', '170520_img93', '170523_img89', '170531_img94',  ...
%      '180404_img077', '180522_img081', '180504_img084', '180519_img085'};
% WF_CRP_licking_only_overview

days = {'170426_img90',  '170425_img91', '170428_img92', '170518_img93', '170519_img89', '170529_img94', '170926_img036', '170926_img044', ...
    '171122_img050', '180104_img067', '180108_img070',  '180403_img077', '180518_img081', '180502_img084', '180517_img085'};
bx_source      = ['Z:\Data\WidefieldImaging\GCaMP\behavior\'];
bx_outputs_dir = ['Z:\Analysis\Cue_reward_pairing_analysis\BxAndAnalysisOutputs\BxOutputs\'];
CRP_fig_dir_base    = ['Z:\Analysis\Cue_reward_pairing_analysis\CRPFigureFolder\'];
old_cd = cd; %save old cd so I can restore it later
time_before_ms = 2000; %defines the window around the cue presentation which will be taken for plotting
time_after_ms = 3000;

RT_raw_all = [];
figure;
for ii = 1:length(days)
    days(ii)
    load(['Z:\Analysis\Cue_reward_pairing_analysis\CRPFigureFolder\RT_data\', days{ii}, '.mat']);
    subplot(4,4,ii);
    hist(RT_this_session_raw(RT_this_session_raw<1001), 40);
    title([days(ii)]);
    if ii > 12
        xlabel('ms relative to cue onset (25ms bins)')
    end
    ylabel('# licks per bin');
    vline(600, 'k');
    xlim([0 1000]);
    
    RT_raw_all = [RT_raw_all, RT_this_session_raw];
end
suptitle('post-learning first lick time relative to cue onset');

figure;
hist(RT_raw_all(RT_raw_all<1001), 40);
title(['all sessions n=', num2str(length(days)), ' post-learning: first lick time relative to cue onset']);
xlabel('ms relative to cue onset (25ms bins)')
ylabel('# licks per bin');
vline(600, 'k');







