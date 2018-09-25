%% assign pathnames and datasets to be analyzed/written for moving dot experiments
clear;
%NEED TO UPDATE THIS SO IT ACCESSES SPREADSHEET INSTEAD OF JUST WRITING IN THE NAMES
sessions = {'180417_img1005_1','180419_img1005_1','180419_img1007_1','180423_img1005_1',...
    '180424_img1008_1','180425_img1008_1','180428_img1008_1','180429_img1008_1','180430_img1005_1',...
    '180430_img1007_1','180430_img1008_1','180430_img1010_1','180505_img1007_1','180505_img1008_1','180505_img1010_1'};

days = {'1005-180417_1','1005-180419_1','1007-180419_1','1005-180423_1','1008-180424_1','1008-180425_1','1008-180428_1',...
    '1008-180429_1','1005-180430_1','1007-180430_1','1008-180430_1','1010-180430_1','1007-180505_1','1008-180505_1',...
    '1010-180505_1'};
image_dest_base  = ['Z:\Analysis\WF_MovingDots_Analysis\BxAndAnalysisOutputs\']; %stores the data on crash in the movingDots analysis folder
% behavior analysis results 
color_code = {'b','r','k'};

%% SECTION I - reverse trig ave across sessions
ave_dfOvF_rev_all = []; 
ACS_dest = [image_dest_base 'acrossSessions\'];
for i = 1: length(sessions)
    image_dest = [image_dest_base sessions{i} '\' sessions{i}];
    img_anal = load([image_dest '_imgAnalysis.mat' ]);
    ave_dfOvF_rev = img_anal.ave_dfOvF_rev;
    ave_dfOvF_rev_all = cat(1,ave_dfOvF_rev_all,ave_dfOvF_rev');
    
end
%acs: across sessions
ave_dfOvF_rev_ACS = mean(ave_dfOvF_rev_all);
ste_dfOvF_rev_ACS = std(ave_dfOvF_rev_all)/sqrt(length(ave_dfOvF_rev_all));

dfOvF_reverseTrigger_ACS = figure;
errorbar(ave_dfOvF_rev_ACS,ste_dfOvF_rev_ACS,'linewidth', 1.5); hold on;
xlim([0 21]);
%ylim([-0.05 0.05]);
vline(6, 'r');vline(16, 'r');
ylabel('df/f'); 

supertitle(['reverse triggered average across sessions']); 
saveas(dfOvF_reverseTrigger_ACS, [ ACS_dest 'dfOvF_reverseTrigger_ACS']);
save([ ACS_dest 'ACSanalysis.mat' ],'ave_dfOvF_rev_ACS','ste_dfOvF_rev_ACS','-append' );



%% SECTION II - reverse trig ave during run across sessions
ave_dfOvF_revRun_all = []; ave_speed_revRun_all = []; 
ACS_dest = [image_dest_base 'acrossSessions\'];
for i = 1: length(sessions)
    image_dest = [image_dest_base sessions{i} '\' sessions{i}];
    img_anal = load([image_dest '_imgAnalysis.mat' ]);
    ave_dfOvF_revRun = img_anal.ave_dfOvF_revRun;
    ave_dfOvF_revRun_all = cat(1,ave_dfOvF_revRun_all,ave_dfOvF_revRun');
    ave_speed_revRun = img_anal.ave_speed_revRun;
    ave_speed_revRun_all = cat(1,ave_speed_revRun_all,ave_speed_revRun );
end
%acs: across sessions
ave_dfOvF_revRun_ACS = mean(ave_dfOvF_revRun_all);
ste_dfOvF_revRun_ACS = std(ave_dfOvF_revRun_all)/sqrt(length(ave_dfOvF_revRun_all));
ave_speed_revRun_ACS = mean(ave_speed_revRun_all);
ste_speed_revRun_ACS = std(ave_speed_revRun_all)/sqrt(length(ave_speed_revRun_all));

dfOvF_revRun_ACS = figure;
subplot(2,1,1);hold on;
errorbar(ave_dfOvF_revRun_ACS,ste_dfOvF_revRun_ACS,'linewidth', 1.5); hold on;
xlim([0 21]);
%ylim([-0.05 0.05]);
vline(6, 'r');vline(16, 'r');
ylabel('df/f'); 

subplot(2,1,2);hold on;
errorbar(ave_speed_revRun_ACS,ste_speed_revRun_ACS,'linewidth', 1.5); hold on;
xlabel('frames');
ylabel('speed');
vline(6, 'r');vline(16, 'r');
xlim([0 21]);

supertitle(['reverse triggered average across sessions during running']); 
saveas(dfOvF_revRun_ACS, [ ACS_dest 'dfOvF_revRun_ACS']);
save([ ACS_dest 'ACSanalysis.mat' ],'ave_dfOvF_revRun_ACS','ste_dfOvF_revRun_ACS',...
   'ave_speed_revRun_ACS','ste_speed_revRun_ACS', '-append' );


%% SECTION III reverse trig ave during stay across sessions
ave_dfOvF_revStay_all = []; ave_speed_revStay_all = []; 
ACS_dest = [image_dest_base 'acrossSessions\'];
for i = 1: length(sessions)
    image_dest = [image_dest_base sessions{i} '\' sessions{i}];
    img_anal = load([image_dest '_imgAnalysis.mat' ]);
    ave_dfOvF_revStay = img_anal.ave_dfOvF_revStay;
    ave_dfOvF_revStay_all = cat(1,ave_dfOvF_revStay_all,ave_dfOvF_revStay');
    ave_speed_revStay = img_anal.ave_speed_revStay;
    ave_speed_revStay_all = cat(1,ave_speed_revStay_all,ave_speed_revStay );
end
%acs: across sessions
ave_dfOvF_revStay_ACS = mean(ave_dfOvF_revStay_all);
ste_dfOvF_revStay_ACS = std(ave_dfOvF_revStay_all)/sqrt(length(ave_dfOvF_revStay_all));
ave_speed_revStay_ACS = mean(ave_speed_revStay_all);
ste_speed_revStay_ACS = std(ave_speed_revStay_all)/sqrt(length(ave_speed_revStay_all));

dfOvF_revStay_ACS = figure;
subplot(2,1,1);hold on;
errorbar(ave_dfOvF_revStay_ACS,ste_dfOvF_revStay_ACS,'linewidth', 1.5); hold on;
%ylim([-0.05 0.05]);
vline(6, 'r');vline(16, 'r');
xlim([0 21]);ylabel('df/f'); 

subplot(2,1,2);hold on;
errorbar(ave_speed_revStay_ACS,ste_speed_revStay_ACS,'linewidth', 1.5); hold on;
xlabel('frames');
ylabel('speed');
vline(6, 'r');vline(16, 'r');
xlim([0 21]);

supertitle(['reverse triggered average across sessions during stationary']); 
saveas(dfOvF_revStay_ACS, [ ACS_dest 'dfOvF_revStay_ACS']);
save([ ACS_dest 'ACSanalysis.mat' ],'ave_dfOvF_revStay_ACS','ste_dfOvF_revStay_ACS',...
   'ave_speed_revStay_ACS','ste_speed_revStay_ACS', '-append' );



