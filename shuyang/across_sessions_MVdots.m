%% assign pathnames and datasets to be analyzed/written for moving dot experiments
clear;
%NEED TO UPDATE THIS SO IT ACCESSES SPREADSHEET INSTEAD OF JUST WRITING IN THE NAMES
sessions = {'180414_img1005_1','180414_img1007_1','180414_img1008_1','180417_img1008_1',...
    '180419_img1008_1','180423_img1010_1'}; 
days = {'1005-180414_1','1007-180414_1','1008-180414_1','1008-180417_1','1008-180419_1','1010-180423_1'};
image_dest_base  = ['Z:\Analysis\WF_MovingDots_Analysis\BxAndAnalysisOutputs\']; %stores the data on crash in the movingDots analysis folder
% behavior analysis results 
color_code = {'b','r','k'};

%% SECTION I - draw df/f vs. speed across sessions 
ACS_dest = [image_dest_base 'acrossSessions\'];

dfOvF_spd_all = {}; isp_all = {};
for i = 1: length(sessions)
    image_dest = [image_dest_base sessions{i} '\' sessions{i}];
    img_anal = load([image_dest '_imgAnalysis.mat' ]);
    isp = img_anal.isp; isp = double(isp);
    dfOvF_spdmean = img_anal.dfOvF_spdmean;
    dfOvF_spd_all = cat(2,dfOvF_spd_all, dfOvF_spdmean); %problem: don't know how to put all dfOvF_spd together. this cell can't be cell2matted later.
    isp_all = cat(2,isp_all, isp);
end
%generate matrix for speeds
isp_temp = cell2mat(cellfun(@size,isp_all, 'UniformOutput',0));
isp_length = max(isp_temp);
isp_mat = nan(size(isp_all,2), isp_length);
for n = 1: size(isp_mat,1)
    temp_dfOvF = isp_all{n};
    isp_mat(n,1:size(temp_dfOvF,2)) = temp_dfOvF;
end

dfOvF_spd_all_mat = []; isp_all_mat = [];
for f = 1: size(dfOvF_spd_all,2)
    temp_dfOvF = dfOvF_spd_all{f}; 
    dfOvF_spd_all_mat = [dfOvF_spd_all_mat; temp_dfOvF(:)];
    %make the matrix of speed exactly the same as df/f
    [isp_mat1,isp_mat2] = size(temp_dfOvF);
    temp_isp = isp_mat(f,:);
    temp_isp(isnan(temp_isp)) = []; %delate NaN
    temp_isp2 = repmat(temp_isp,isp_mat1,1);
    isp_all_mat = [isp_all_mat; temp_isp2(:)];
end

spd_all_plotx = sort(unique(isp_all_mat));
dfOvF_all_ploty = zeros(length(spd_all_plotx),1);
dfOvF_errbar = zeros(length(spd_all_plotx),1);
for j = 1: length(spd_all_plotx)
    c = find(isp_all_mat == spd_all_plotx(j));
    dfOvF_all_ploty(j) = mean(dfOvF_spd_all_mat(c));
    dfOvF_errbar(j) = std(dfOvF_spd_all_mat(c))/sqrt(length(c));
end
dfOvF_vs_spd_ACS = figure; 
errorbar(spd_all_plotx,dfOvF_all_ploty,dfOvF_errbar,'.','LineStyle','-','linewidth', 1.25,'MarkerSize',20); hold on;
xlabel ('speed');
ylabel('ave df/f');
title(['df/f vs. speed across sessions']); 
saveas(dfOvF_vs_spd_ACS, [ ACS_dest 'dfOvF_vs_spd_ACS']);
save([ ACS_dest 'ACSanalysis.mat' ],'spd_all_plotx','dfOvF_all_ploty','dfOvF_errbar','-append' );



%% SECTION II - df/f for stay&run across sessions 
ave_dfOvF_behav_all = [];
ACS_dest = [image_dest_base 'acrossSessions\'];
for i = 1: length(sessions)
    image_dest = [image_dest_base sessions{i} '\' sessions{i}];
    img_anal = load([image_dest '_imgAnalysis.mat' ]);
    ave_dfOvF_behav = img_anal.ave_dfOvF_behav;
    ave_dfOvF_behav_all = cat(1,ave_dfOvF_behav_all,ave_dfOvF_behav ); %problem: don't know how to put all dfOvF_spd together. this cell can't be cell2matted later.
end
%acs: across sessions
ave_dfOvF_behav_ACS = mean(ave_dfOvF_behav_all);
ste_dfOvF_behav_ACS = std(ave_dfOvF_behav_all)/sqrt(length(ave_dfOvF_behav_all));

x = [1,2];
dfOvF_behavStates_ACS = figure;
errorbar(x,ave_dfOvF_behav_ACS,ste_dfOvF_behav_ACS,'.','LineStyle','-','linewidth', 1.25,'MarkerSize',20); hold on;
xlabel ('behavioral state');xlim([0.5 2.5]);
set(gca,'XTick',x,'XTicklabel',{'stationary','run'});
ylabel('df/f');
ylim([-0.12 0.02]);
title(['df/f for each beahvioral state across sessions']); 
saveas(dfOvF_behavStates_ACS, [ ACS_dest 'dfOvF_behavStates_ACS']);
save([ ACS_dest 'ACSanalysis.mat' ],'ave_dfOvF_behav_ACS','ste_dfOvF_behav_ACS' );


%% SECTION III - df/f for beforeRunAfter across sessions 
ave_dfOvF_befoRunAft_all = [];
ACS_dest = [image_dest_base 'acrossSessions\'];
for i = 1: length(sessions)
    image_dest = [image_dest_base sessions{i} '\' sessions{i}];
    img_anal = load([image_dest '_imgAnalysis.mat' ]);
    ave_befoRunaft = img_anal.ave_befoRunaft;
    ave_dfOvF_befoRunAft_all = cat(1,ave_dfOvF_befoRunAft_all,ave_befoRunaft ); %problem: don't know how to put all dfOvF_spd together. this cell can't be cell2matted later.
end
%acs: across sessions
ave_dfOvF_befoRunAft_ACS = mean(ave_dfOvF_befoRunAft_all);
ste_dfOvF_befoRunAft_ACS = std(ave_dfOvF_befoRunAft_all)/sqrt(length(ave_dfOvF_befoRunAft_all));

x = [1,2,3];
dfOvF_befoRunAft_ACS = figure;
errorbar(x,ave_dfOvF_befoRunAft_ACS,ste_dfOvF_befoRunAft_ACS,'.','LineStyle','-','linewidth', 1.25,'MarkerSize',20); hold on;
plot(x,ave_dfOvF_befoRunAft_ACS,'.','Color','b','MarkerSize',20,'LineStyle','none');

xlabel ('behavioral state');xlim([0.5 3.5]);
set(gca,'XTick',x,'XTicklabel',{'before','run','after'});
ylabel('df/f');
title(['df/f around running across sessions']); 
saveas(dfOvF_befoRunAft_ACS, [ ACS_dest 'dfOvF_befoRunAft_ACS']);
save([ ACS_dest 'ACSanalysis.mat' ],'ave_dfOvF_befoRunAft_ACS','ste_dfOvF_befoRunAft_ACS','-append' );


%% SECTION IV - df/f and speed 300ms across sessions
dfOvF_bl300ms_all = []; speed_bl300ms_all = [];
ACS_dest = [image_dest_base 'acrossSessions\'];
period = 9;
for i = 1: length(sessions)
    image_dest = [image_dest_base sessions{i} '\' sessions{i}];
    img_anal = load([image_dest '_imgAnalysis.mat' ]);
    speed_plot300ms = img_anal.speed_plot300ms;
    dfOvF_plot300ms = img_anal.dfOvF_plot300ms;
    dfOvF_bgPart300ms = dfOvF_plot300ms(:,(1: 1+period)); %first point is before running
    dfOvF_lstPart300ms = dfOvF_plot300ms(:,(end-period:end));%last point is after running
    dfOvF_blPart300ms = [dfOvF_bgPart300ms,dfOvF_lstPart300ms];%bl: begin and last
    speed_bgPart300ms = speed_plot300ms(:,(1: 1+period));
    speed_lstPart300ms = speed_plot300ms(:,(end-period:end));
    speed_blPart300ms = [speed_bgPart300ms,speed_lstPart300ms];%bl: begin and last
    
    dfOvF_bl300ms_all = cat(1,dfOvF_bl300ms_all,dfOvF_blPart300ms );
    speed_bl300ms_all = cat(1,speed_bl300ms_all,speed_blPart300ms );
    
end

ave_dfOvF_bl300ms = mean(dfOvF_bl300ms_all);
ave_speed_bl300ms = mean(speed_bl300ms_all);
ste_dfOvF_bl300ms = std(dfOvF_bl300ms_all)/sqrt(length(dfOvF_bl300ms_all));
ste_speed_bl300ms = std(speed_bl300ms_all)/sqrt(length(speed_bl300ms_all));

x = (1: 1: 20);
dfOvF_run300ms_ACS = figure;
subplot(2,1,1);hold on;
errorbar(x,ave_dfOvF_bl300ms,ste_dfOvF_bl300ms,'.','LineStyle','-','linewidth', 1.25,'MarkerSize',20); hold on;
%xlim([-5 10]);
ylim([-0.1 0]);
ylabel('df/f'); 

subplot(2,1,2);hold on;
errorbar(x,ave_speed_bl300ms,ste_speed_bl300ms,'.','LineStyle','-','linewidth', 1.25,'MarkerSize',20); hold on;
xlabel('frames');
ylabel('speed');
%xlim([-5 10]);

supertitle(['run 300ms average across sessions']); 
saveas(dfOvF_run300ms_ACS, [ ACS_dest 'dfOvF_run300ms_ACS']);
save([ ACS_dest 'ACSanalysis.mat' ],'ave_dfOvF_bl300ms','ave_speed_bl300ms',...
   'ste_dfOvF_bl300ms','ste_speed_bl300ms', '-append' );


%% SECTION V - run trig ave across sessions
ave_dfOvF_runTrigger_all = []; ave_speed_runTrigger_all = []; 
ACS_dest = [image_dest_base 'acrossSessions\'];
for i = 1: length(sessions)
    image_dest = [image_dest_base sessions{i} '\' sessions{i}];
    img_anal = load([image_dest '_imgAnalysis.mat' ]);
    ave_dfOvF_runTrigger = img_anal.ave_dfOvF_runTrigger;
    if size(ave_dfOvF_runTrigger,1) > 1
        ave_dfOvF_runTrigger = ave_dfOvF_runTrigger';
    end
    ave_dfOvF_runTrigger_all = cat(1,ave_dfOvF_runTrigger_all,ave_dfOvF_runTrigger);
    ave_speed_runTrigger = img_anal.ave_speed_runTrigger;
    ave_speed_runTrigger_all = cat(1,ave_speed_runTrigger_all,ave_speed_runTrigger );
end
%acs: across sessions
ave_dfOvF_runTrigger_ACS = mean(ave_dfOvF_runTrigger_all);
ste_dfOvF_runTrigger_ACS = std(ave_dfOvF_runTrigger_all)/sqrt(length(ave_dfOvF_runTrigger_all));
ave_speed_runTrigger_ACS = mean(ave_speed_runTrigger_all);
ste_speed_runTrigger_ACS = std(ave_speed_runTrigger_all)/sqrt(length(ave_speed_runTrigger_all));

x = (1:15);
dfOvF_runTrigger_ACS = figure;
subplot(2,1,1);hold on;
errorbar(x,ave_dfOvF_runTrigger_ACS,ste_dfOvF_runTrigger_ACS,'.','LineStyle','-','linewidth', 1.25,'MarkerSize',20); hold on;
%xlim([-5 10]);
%ylim([-0.05 0.05]);
vline(6, 'r','running start');
ylabel('df/f'); 

subplot(2,1,2);hold on;
errorbar(x,ave_speed_runTrigger_ACS,ste_speed_runTrigger_ACS,'.','LineStyle','-','linewidth', 1.25,'MarkerSize',20); hold on;
xlabel('frames');
ylabel('speed');
vline(6, 'r','running start');
%xlim([-5 10]);

supertitle(['run triggered average across sessions']); 
saveas(dfOvF_runTrigger_ACS, [ ACS_dest 'dfOvF_runTrigger_ACS']);
save([ ACS_dest 'ACSanalysis.mat' ],'ave_dfOvF_runTrigger_ACS','ste_dfOvF_runTrigger_ACS',...
   'ave_speed_runTrigger_ACS','ste_speed_runTrigger_ACS', '-append' );


