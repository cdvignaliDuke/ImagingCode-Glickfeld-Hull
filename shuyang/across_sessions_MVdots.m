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
dfOvF_spd_all = {}; isp_all = {};
for i = 1: length(sessions)
    image_dest = [image_dest_base sessions{i} '\' sessions{i}];
    img_anal = load([image_dest '_imgAnalysis.mat' ]);
    isp = img_anal.isp; isp = double(isp);
    dfOvF_spdmean = img_anal.dfOvF_spdmean;
    dfOvF_spd_all = cat(2,dfOvF_spd_all, dfOvF_spdmean); %problem: don't know how to put all dfOvF_spd together. this cell can't be cell2matted later.
    isp_all = cat(2,isp_all, isp);
end

for f = 1: size(dfOvF_spd_all,2)
    %dfOvF_spd_all(f,:) = 

end
% ----------------------------------------same thing if can access dfOvF_spd_all
isp = unique(speed);
    dfOvF_spd = []; dfOvF_spdmean = []; dfOvF_spdste = [];
    for k = 1:length(isp)
        dfOvF_spd = dfOvF(:,speed==isp(k));
        dfOvF_spdmean(:,k)  = mean(dfOvF_spd,2);
        dfOvF_spdste(:,k) = std(dfOvF_spd,[],2)/sqrt(length(dfOvF_spd));
    end
    isp_plot =  repmat(isp,size(dfOvF_spdmean,1),1);
    dfOvF_ave_vs_spd = figure;
    errorbar(isp_plot',dfOvF_spdmean',dfOvF_spdste','linewidth', 1.5);
    %if do errorbar(y,ste y), it can do multiple lines at once. But if do
    %errorbar(x,y,ste y), size of x and y must match(if your y contains n lines, even though the x for all lines are the same, x must have n lines too. 
    xlabel ('speed');
    ylabel('ave df/f');
    title(['df/f vs. speed',sessions{i}]); legend;
    saveas(dfOvF_ave_vs_spd, [image_dest, '_dfOvf_vs_speed']);
    save([image_dest '_imgAnalysis.mat' ],'isp', 'dfOvF_spdmean');
    
    
    
%---------------------------------------------prob. don't need this
%generate matrix for speeds
isp_temp = cell2mat(cellfun(@size,isp_all, 'UniformOutput',0));
isp_length = max(isp_temp);
isp_mat = nan(size(isp_all,2), isp_length);
for n = 1: size(isp_mat,1)
    temp = isp_all{n};
    isp_mat(n,1:size(temp,2)) = temp;
end
ave_isp_all = mean(isp_mat);

%generate matrix for df/f_speedmean
dfOvF_spd_temp = cell2mat(cellfun(@size,dfOvF_spd_all, 'UniformOutput',0));
dfOvF_spd_length = max(dfOvF_spd_temp);
dfOvF_spd_mat = nan(size(dfOvF_spd_all,2), dfOvF_spd_length);
for n = 1: size(dfOvF_spd_mat,1)
    temp = dfOvF_spd_all{n};
    dfOvF_spd_mat(n,1:size(temp,2)) = temp;
end


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
errorbar(ave_dfOvF_behav_ACS,ste_dfOvF_behav_ACS,'linewidth', 2);
xlabel ('behavioral state');xlim([0.5 2.5]);
set(gca,'XTick',x,'XTicklabel',{'stationary','run'});
ylabel('df/f');
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
errorbar(ave_dfOvF_befoRunAft_ACS,ste_dfOvF_befoRunAft_ACS,'linewidth', 2);
xlabel ('behavioral state');xlim([0.5 3.5]);
set(gca,'XTick',x,'XTicklabel',{'before','run','after'});
ylabel('df/f');
title(['df/f around running across sessions']); 
saveas(dfOvF_befoRunAft_ACS, [ ACS_dest 'dfOvF_befoRunAft_ACS']);
save([ ACS_dest 'ACSanalysis.mat' ],'ave_dfOvF_befoRunAft_ACS','ste_dfOvF_befoRunAft_ACS','-append' );

%% SECTION IV - df/f and speed 300ms across sessions


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

dfOvF_runTrigger_ACS = figure;
subplot(2,1,1);hold on;
errorbar(ave_dfOvF_runTrigger_ACS,ste_dfOvF_runTrigger_ACS,'linewidth', 1.5); hold on;
%xlim([-5 10]);
%ylim([-0.05 0.05]);
vline(4, 'r','running start');
ylabel('df/f'); 

subplot(2,1,2);hold on;
errorbar(ave_speed_runTrigger_ACS,ste_speed_runTrigger_ACS,'linewidth', 1.5); hold on;
xlabel('frames');
ylabel('speed');
vline(4, 'r','running start');
%xlim([-5 10]);

supertitle(['run triggered average across sessions']); 
saveas(dfOvF_runTrigger_ACS, [ ACS_dest 'dfOvF_runTrigger_ACS']);
save([ ACS_dest 'ACSanalysis.mat' ],'ave_dfOvF_runTrigger_ACS','ste_dfOvF_runTrigger_ACS',...
   'ave_speed_runTrigger_ACS','ste_speed_runTrigger_ACS', '-append' );


