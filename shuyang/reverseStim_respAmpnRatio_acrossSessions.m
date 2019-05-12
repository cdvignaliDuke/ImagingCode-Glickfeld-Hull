% response amplitude and ratio for each session, distribution of a scores
% scatter plot of mean resp amplitude of all session during stationary and running
% scatter for resp ratio of all session during stationary and running
%%
clear;
sessions = {'180417_img1005_1','180419_img1005_1','180419_img1007_1','180423_img1005_1','180424_img1008_1',...
    '180425_img1008_1','180428_img1008_1','180429_img1008_1','180430_img1005_1','180430_img1007_1',...
    '180430_img1008_1','180430_img1010_1','180505_img1007_1','180505_img1008_1','180505_img1010_1'}; 
image_dest_base    = ['Z:\Analysis\WF_MovingDots_Analysis\BxAndAnalysisOutputs\']; %stores the data on crash in the movingDots analysis folder
overall_dest = ['Z:\Analysis\WF_MovingDots_Analysis\BxAndAnalysisOutputs\acrossSessions\'];

color_code = {'b','r','k'};

%% calculate response amplitude and response ratio for each session
for i = 1:length(sessions)
   %% find amplitudes after reverse and calculate response of each reverse.

    % load things
    image_dest = [image_dest_base sessions{i} '\' sessions{i}];
    dfOvF_rev_struct = load([image_dest, '_dfOvF_rev.mat']);
    dfOvF_rev_stay = dfOvF_rev_struct.dfOvF_rev_stay;
    dfOvF_rev_run = dfOvF_rev_struct.dfOvF_rev_run;
   
   %find df/f 500ms before reverse stim, and 100ms to 700ms after stim
   %dfOvF_after_rev : 3-7s after reverse (9~13)
   %dfOvF_befor_rev : -5~-1s(1~5)
   dfOvF_stayAftRev = dfOvF_rev_stay(:,(9:13),:);
   dfOvF_stayBefoRev = dfOvF_rev_stay(:,(1:5),:);
   dfOvF_runAftRev = dfOvF_rev_run(:,(9:13),:);
   dfOvF_runBefoRev = dfOvF_rev_run(:,(1:5),:);
   
   % calculate average and std of df/f before reverse
   dfOvF_stayBefoRev_ave = squeeze(mean(dfOvF_stayBefoRev,2));
   dfOvF_runBefoRev_ave  = squeeze(mean(dfOvF_runBefoRev,2));
   dfOvF_stayBefoRev_std = squeeze(std(dfOvF_stayBefoRev,0,2));
   dfOvF_runBefoRev_std  = squeeze(std(dfOvF_runBefoRev,0,2));
   
   % find the amplitudes in each window during reverse, calculate diff between amp and ave before reverse,and ste of the diff of this session
   dfOvF_stayAmp = [];
   response_stay = [];
   for n = 1:size(dfOvF_stayAftRev,1)
       dfOvF_stayAmp(n,:) = max(dfOvF_stayAftRev(n,:,:));
       response_stay (n,:) = dfOvF_stayAmp(n,:) - dfOvF_stayBefoRev_ave(n,:);
   end
   stayRespAve = mean(response_stay);
   stayRespSte = std(response_stay)/sqrt(size(response_stay,1));
   %stayAmp Ave and ste for all session
   stayRespAve_all{i} = stayRespAve;
   stayRespSte_all{i} = stayRespSte;
   
   dfOvF_runAmp = [];
   response_run = [];
   for n = 1:size(dfOvF_runAftRev,1)
       dfOvF_runAmp(n,:) = max(dfOvF_runAftRev(n,:,:));
       response_run (n,:) = dfOvF_runAmp(n,:) - dfOvF_runBefoRev_ave(n,:);       
   end
   runRespAve = mean(response_run);
   runRespSte = std(response_run)/sqrt(size(response_run,1));
  %runAmp Ave and ste for all session
   runRespAve_all{i} = runRespAve;
   runRespSte_all{i} = runRespSte;
   
   
   %% find amplitude with a z socre > 2, draw distribution of z score, 
   % find the amps that are bigger than befoRevAve + 2 std and color code when draw df/f for individual reverses plot
   z = 2;
   
   %STAY%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   inx_stayamp_z = zeros(size(dfOvF_stayAmp,1),size(dfOvF_stayAmp,2));
   Ampzscore_stay = zeros(size(dfOvF_stayAmp,1),size(dfOvF_stayAmp,2));
   for roi = 1:size(dfOvF_stayAmp,2)
       dfOvF_stayAmp_roi = dfOvF_stayAmp(:,roi);
       Ampzscore_stay(:,roi) = (dfOvF_stayAmp_roi - dfOvF_stayBefoRev_ave(:,roi))./dfOvF_stayBefoRev_std(:,roi);
       inx_stayamp_z(:,roi) = Ampzscore_stay(:,roi) > z;
       set(0,'DefaultFigureVisible','off');
       zscore_distfig(roi) = figure;
       scatter((1:size(Ampzscore_stay,1)),Ampzscore_stay(:,roi),'filled'); hold on;
       refline(0,2);
       xlabel('reverse trials'); ylabel('z score');
       title(['STAY reverse stimuli z score ROI', num2str(roi) ' ' sessions{i}]);
       saveas(zscore_distfig(roi),[image_dest '_reverseTrigAve_stayZscore' num2str(roi)]);
       %ROIvals{roi}   = dfOvF_stayAmp_roi(Ampzscore > z);
   end
   stay_respRatio = sum(inx_stayamp_z == 1)./size(inx_stayamp_z,1);
   stay_respRatio_all{i} = stay_respRatio;
   
   t = 1:21;
   for n = 1: size(dfOvF_rev_stay,3)
       %make every line start at the same height
       startMean_stay = mean(dfOvF_rev_stay(:,1,n));
       diff = dfOvF_rev_stay(:,1,n)-startMean_stay;
       diff_array = repmat(diff,1,size(dfOvF_rev_stay,2));
       dfOvF_rev_stay_plot = dfOvF_rev_stay(:,:,n) - diff_array;
       dfOvF_rev_stay_plot = dfOvF_rev_stay_plot';
       indi_revStay_fig(n) = figure;
       for p = 1: size(dfOvF_rev_stay_plot,2)
           if inx_stayamp_z(p,n) == 1
               plot(t,dfOvF_rev_stay_plot(:,p),'color','r');hold on;
           else
               plot(t,dfOvF_rev_stay_plot(:,p),'color','b'); hold on;
           end
       end
       xlabel('frames'); ylabel('df/f');
       xlim([1 21]); %ylim([-0.2 0.15]);
       vline(6, 'k');vline(16, 'k');
       title(['STAY reverse stimuli ROI', num2str(n) ' ' sessions{i}]);
       saveas(indi_revStay_fig(n),[image_dest '_reverseTrigAve_stayROI' num2str(n)]);
   end
   
   %RUN%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   inx_runamp_z = zeros(size(dfOvF_runAmp,1),size(dfOvF_runAmp,2));
   Ampzscore_run = zeros(size(dfOvF_runAmp,1),size(dfOvF_runAmp,2));
   for roi = 1:size(dfOvF_runAmp,2)
       dfOvF_runAmp_roi = dfOvF_runAmp(:,roi);
       Ampzscore_run(:,roi) = (dfOvF_runAmp_roi - dfOvF_runBefoRev_ave(:,roi))./dfOvF_runBefoRev_std(:,roi);
       inx_runamp_z(:,roi) = Ampzscore_run(:,roi) > z;
       zscore_distfig(roi) = figure;
       scatter((1:size(Ampzscore_run,1)),Ampzscore_run(:,roi),'filled'); hold on;
       refline(0,2);
       xlabel('reverse trials'); ylabel('z score');
       title(['RUN reverse stimuli z score ROI', num2str(roi) ' ' sessions{i}]);      
       saveas(zscore_distfig(roi),[image_dest '_reverseTrigAve_runZscore' num2str(roi)]);
   end
   run_respRatio = sum(inx_runamp_z == 1)./size(inx_runamp_z,1);
   run_respRatio_all{i} = run_respRatio;
   
   for n = 1: size(dfOvF_rev_run,3)
       %make every line start at the same height
       startMean_run = mean(dfOvF_rev_run(:,1,n));
       diff = dfOvF_rev_run(:,1,n)-startMean_run;
       diff_array = repmat(diff,1,size(dfOvF_rev_run,2));
       dfOvF_rev_run_plot = dfOvF_rev_run(:,:,n) - diff_array;
       dfOvF_rev_run_plot = dfOvF_rev_run_plot';
       indi_revRun_fig(n) = figure;
       for p = 1: size(dfOvF_rev_run_plot,2)
           if inx_runamp_z(p,n) == 1
               plot(t,dfOvF_rev_run_plot(:,p),'color','r');hold on;
           else
               plot(t,dfOvF_rev_run_plot(:,p),'color','b'); hold on;
           end
       end
       xlabel('frames'); ylabel('df/f');
       xlim([1 21]);%ylim([-0.2 0.15]); 
       vline(6, 'k');vline(16, 'k');
       title(['Run reverse stimuli ROI', num2str(n) ' ' sessions{i}]);
       saveas(indi_revRun_fig(n),[image_dest '_reverseTrigAve_runROI' num2str(n)]);
   end
   [stay_H,stay_P] = ttest(dfOvF_stayAmp,dfOvF_stayBefoRev_ave);
   %for n = 1:size(dfOvF_stayBefoRev,1)
       %baseMax(n,:) = max(dfOvF_stayBefoRev(n,:,:));
   %end 

 
end
save([overall_dest 'reverse.mat' ],'stayRespAve_all', 'stayRespSte_all','runRespAve_all','runRespSte_all','stay_respRatio_all','run_respRatio_all');


%% plot mean amp of run vs. stay
Response_struct = load([overall_dest, 'reverse.mat']);
stayRespAve_all = Response_struct.stayRespAve_all;
stayRespSte_all = Response_struct.stayRespSte_all;
runRespAve_all = Response_struct.runRespAve_all;
runRespSte_all = Response_struct.runRespSte_all;

stayRespAves = cell2mat(stayRespAve_all);
stayRespStes = cell2mat(stayRespSte_all);
runRespAves = cell2mat(runRespAve_all);
runRespStes = cell2mat(runRespSte_all);

stay_overall_ave = mean(stayRespAves);
stay_overall_ste = std(stayRespAves)/sqrt(length(stayRespAves));
run_overall_ave = mean(runRespAves);
run_overall_ste = std(runRespAves)/sqrt(length(runRespAves));

set(0,'DefaultFigureVisible','on');
respAllFig = figure; clf;
errorbar(stayRespAves,runRespAves,stayRespStes,'o','horizontal','Color','m','Marker','.','MarkerSize',20); hold on; 
errorbar(stayRespAves,runRespAves,runRespStes,'o','Color','m');hold on;
errorbar(stay_overall_ave,run_overall_ave,run_overall_ste,'o','Color','k','Marker','.','MarkerSize',20);hold on;
errorbar(stay_overall_ave,run_overall_ave,stay_overall_ste,'o','horizontal','Color','k');
xlabel('stay'); ylabel('run');
title('mean response amplitude')
refline(1,0);
saveas(respAllFig,[overall_dest 'stay_run_ave_resp']);

%% plot percentage of response of run vs. stay
Response_struct = load([overall_dest, 'reverse.mat']);
stay_respRatio_all = Response_struct.stay_respRatio_all;
run_respRatio_all = Response_struct.run_respRatio_all;

stay_respRatios = cell2mat(stay_respRatio_all);
run_respRatios = cell2mat(run_respRatio_all);

stay_aveRespRatio = mean(stay_respRatios);
stay_steRespRatio = std(stay_respRatios)/sqrt(length(stay_respRatios));
run_aveRespRatio = mean(run_respRatios);
run_steRespRatio = std(run_respRatios)/sqrt(length(run_respRatios));

ratioAllFig = figure; clf;
scatter(stay_respRatios,run_respRatios,'filled','r');
hold on;
errorbar(stay_aveRespRatio,run_aveRespRatio,run_steRespRatio,'o','Color','k','Marker','.','MarkerSize',20);hold on;
errorbar(stay_aveRespRatio,run_aveRespRatio,stay_steRespRatio,'o','horizontal','Color','k');
%xlim([0 0.95]); %ylim([0 0.95]);
xlabel('stay'); ylabel('run');
refline(1,0);
title('mean response ratio')
saveas(ratioAllFig,[overall_dest 'stay_run_respRatios']);


