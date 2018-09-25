%% SECTION II - reverse trig ave during run across sessions
ave_dfOvF_revRun_all = []; ave_speed_revRun_all = []; 
ACS_dest = [image_dest_base 'acrossSessions\'];
for i = 1: length(sessions)
    image_dest = [image_dest_base sessions{i} '\' sessions{i}];
    img_anal = load([image_dest '_imgAnalysis.mat' ]);
    ave_dfOvF_revRun = img_anal.ave_dfOvF_revRun;
    ave_dfOvF_revRun_all = cat(1,ave_dfOvF_revRun_all,ave_dfOvF_revRun);
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
%xlim([-5 10]);
%ylim([-0.05 0.05]);
vline(4, 'r','running start');
ylabel('df/f'); 

subplot(2,1,2);hold on;
errorbar(ave_speed_revRun_ACS,ste_speed_revRun_ACS,'linewidth', 1.5); hold on;
xlabel('frames');
ylabel('speed');
vline(4, 'r','running start');
%xlim([-5 10]);

supertitle(['run triggered average across sessions']); 
saveas(dfOvF_revRun_ACS, [ ACS_dest 'dfOvF_revRun_ACS']);
save([ ACS_dest 'ACSanalysis.mat' ],'ave_dfOvF_revRun_ACS','ste_dfOvF_revRun_ACS',...
   'ave_speed_revRun_ACS','ste_speed_revRun_ACS', '-append' );