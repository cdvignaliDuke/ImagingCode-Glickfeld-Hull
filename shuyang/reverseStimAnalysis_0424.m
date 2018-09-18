%% SECTION ONE - assign pathnames and datasets to be analyzed/written. 
clear;
sessions = {'180505_img1010_1'};
days = {'1010-180505_1'};
%sessions = {'180413_img1002_1','180417_img1001_1','180419_img1002_1','180423_img1001_1',...
%    '180425_img1001_1'};
%days = {'1002-180413_1','1001-180417_1','1002-180419_1','1001-180423_1','1001-180425_1'};
image_dest_base    = ['Z:\Analysis\WF_MovingDots_Analysis\BxAndAnalysisOutputs\']; %stores the data on crash in the movingDots analysis folder
% behavior analysis results 
color_code = {'b','r','k','m'};

%% load
for i = 1:length(sessions)
    image_dest = [image_dest_base sessions{i} '\' sessions{i}];
    behav_dest = ['Z:\Analysis\WF_MovingDots_Analysis\behavioral_analysis\' days{i}];
    behav_output = load([behav_dest '\' days{i} '_behavAnalysis.mat']);
    speed = behav_output.speed;
    cReverse = behav_output.cReverse_vec;
    dfOvF_struct = load([image_dest, '_dfOvF_staybase.mat']);
    dfOvF = dfOvF_struct.dfOvF_staybase;
end

%% determine if reverse stimuli does anything to df/f. ave df/f before, during, and after reverse
for i = 1:length(sessions)
    befo = 5;
    aft = 15;
    dfOvF_rev = [];
    plot_x = -befo : aft;
    t = -5:1:15;    
    dfOvF = dfOvF';
    %If there's no reverse stimuli during the experiment, it will just not run through this for loop b/c length cReverse is empty.
    for f = 1:length(cReverse)
        dfOvF_rev(end+1,:,:) = dfOvF(cReverse(f)-befo:cReverse(f)+aft,:);
    end
   
    % save dfOvF_revfor later analysis
    save([image_dest '_dfOvF_rev.mat' ],'dfOvF_rev');

    %plot average
    ave_dfOvF_rev = mean(dfOvF_rev,1);
    ste_dfOvF_rev = std(dfOvF_rev,0,1)/sqrt(size(dfOvF_rev,1));
    if size(ave_dfOvF_rev,3) > 1
        ave_dfOvF_rev = squeeze(ave_dfOvF_rev);
        ste_dfOvF_rev = squeeze(ste_dfOvF_rev);
    else
        ave_dfOvF_rev = ave_dfOvF_rev';
        ste_dfOvF_rev = ste_dfOvF_rev';
    end
    
    ave_rev_fig = figure;
    errorbar(ave_dfOvF_rev,ste_dfOvF_rev); hold on;
    xlim([1 21]);
    ylim([-0.3 0.1]);
    xlabel('frames');
    vline(6, 'k');
    vline(16, 'k');
    ylabel('df/f');
    title(['df/f before and after reverse stimuli', days{i}]);
    saveas(ave_rev_fig, [image_dest '_reverseTrigAve']);
    
    % Plot individual reverse windows for each ROI  
   % for n = 1: size(dfOvF_rev,3)
    %    indi_rev_fig = figure;
     %   plot(t,dfOvF_rev(:,:,n)');
      %  xlabel('frames'); ylabel('df/f');
       % title(['df/f before and after reverse stimuli ROI', num2str(n)]);
       % saveas(indi_rev_fig,[image_dest '_reverseTrigROI' num2str(n)]);
   % end

end


%% determine if reverse stim does anything differently during running vs. stay
for i = 1:length(sessions)
    speed_rev_stay = [];
    speed_rev_run = [];
    ealier = 5;
    later = 15;
    plot_x = -ealier : later;
    
    % find speed before and after reverse, calculate average speed across reverse windows.
    for v = 1: length(cReverse)
        if sum(speed(cReverse(v)-ealier:cReverse(v)) == 0) == 1+ealier && cReverse(v)+later <= length(speed)
            speed_rev_stay = cat(1, speed_rev_stay,speed(cReverse(v)-ealier:cReverse(v)+later));
        elseif cReverse(v)+later <= length(speed)
            speed_rev_run = cat(1,speed_rev_run, speed(cReverse(v)-ealier:cReverse(v)+later));
        else
            continue
        end
    end
    ave_speed_revStay = mean(speed_rev_stay);
    ste_speed_revStay = std(double(speed_rev_stay),0,1)/sqrt(length(speed_rev_stay));
    ave_speed_revRun = mean(speed_rev_run);
    ste_speed_revRun = std(double(speed_rev_run),0,1)/sqrt(length(speed_rev_run));

    % find df/f before and after every reverse window
    dfOvF = dfOvF';
    dfOvF_rev_stay = [];
    dfOvF_rev_run  = [];
    for f = 1:length(cReverse)
        if sum(speed(cReverse(f)-ealier:cReverse(f)) == 0) == 1+ealier && cReverse(f)+later <= length(speed)
            dfOvF_rev_stay(end+1,:,:) = dfOvF(cReverse(f)-ealier:cReverse(f)+later,:);
        elseif cReverse(f)+later <= length(speed)
            dfOvF_rev_run(end+1,:,:) = dfOvF(cReverse(f)-ealier:cReverse(f)+later,:);
        else
            continue
        end
    end
    
     % save dfOvF_rev_stay and run for later analysis
    save([image_dest '_dfOvF_rev.mat' ],'dfOvF_rev_stay', 'dfOvF_rev_run', '-append');
    
    ave_dfOvF_revStay = mean(dfOvF_rev_stay,1);
    ste_dfOvF_revStay = std(dfOvF_rev_stay,0,1)/sqrt(size(dfOvF_rev_stay,2)); 
    ave_dfOvF_revRun = mean(dfOvF_rev_run,1);
    ste_dfOvF_revRun = std(dfOvF_rev_run,0,1)/sqrt(size(dfOvF_rev_run,2)); 
    
    if size(ave_dfOvF_revStay,3) > 1
        ave_dfOvF_revStay = squeeze(ave_dfOvF_revStay);
        ste_dfOvF_revStay = squeeze(ste_dfOvF_revStay);
        ave_dfOvF_revRun = squeeze(ave_dfOvF_revRun);
        ste_dfOvF_revRun = squeeze(ste_dfOvF_revRun);
    else
        ave_dfOvF_revStay = ave_dfOvF_revStay';
        ste_dfOvF_revStay = ste_dfOvF_revStay';
        ave_dfOvF_revRun = ave_dfOvF_revRun';
        ste_dfOvF_revRun = ste_dfOvF_revRun';
    end
    
    % plot STAY / RUN Averages
    behav_ave_fig = figure;clf
    %t = -5:1:15;
    subplot(2,2,1);
    errorbar(ave_dfOvF_revStay,ste_dfOvF_revStay); hold on;
    xlim([1 21]); ylim([-0.3 0.1]);
    vline(6, 'k');vline(16, 'k');
    ylabel('df/f'); title('stay');
    
    subplot(2,2,2);
    errorbar(ave_dfOvF_revRun,ste_dfOvF_revRun); hold on;
    xlim([1 21]); ylim([-0.3 0.1]);
    vline(6, 'k');vline(16, 'k');
    title('run'); 

    subplot(2,2,3);
    errorbar(ave_speed_revStay,ste_speed_revStay); hold on;
    xlim([1 21]); ylim([0 60]);
    vline(6, 'k');vline(16, 'k');
    ylabel('speed'); xlabel('frames');
    
    subplot(2,2,4);
    errorbar(ave_speed_revRun,ste_speed_revRun); hold on;
    xlim([1 21]); ylim([0 60]);
    vline(6, 'k');vline(16, 'k');
    xlabel('frames');
    
    supertitle(['df/f before and after reverse stimuli', days{i}]);

    saveas(behav_ave_fig,[image_dest '_reverseTrigAve_behav']);
        
    % STAY / RUN Individual reverses, make all of the lines start at the same point
    %for n = 1: size(dfOvF_rev_stay,3)
     %   startMean_stay = mean(dfOvF_rev_stay(:,1,n));
     %   diff = dfOvF_rev_stay(:,1,n)-startMean_stay;
     %   diff_array = repmat(diff,1,size(dfOvF_rev_stay,2));
     %   dfOvF_rev_stay_plot = dfOvF_rev_stay(:,:,n) - diff_array;
     %   indi_revStay_fig = figure;
     %   plot(t,dfOvF_rev_stay_plot(:,:,n)');
     %   xlabel('frames'); ylabel('df/f');
     %   ylim([-0.2 0.6])
     %   title(['STAY reverse stimuli ROI', num2str(n)]);
     %   saveas(indi_revStay_fig,[image_dest '_reverseTrigAve_stayROI' num2str(n)]);

   %end
    
   % for n = 1: size(dfOvF_rev_stay,3)
   %     indi_revStay_fig = figure;
   %     plot(t,dfOvF_rev_stay(:,:,n)');
   %     xlabel('frames'); ylabel('df/f');
   %     title(['STAY reverse stimuli ROI', num2str(n)]);
   %     saveas(indi_revStay_fig,[image_dest '_reverseTrigAve_stayROI' num2str(n)]);
   % end
    
   % for n = 1: size(dfOvF_rev_run,3)
   %     indi_revRun_fig = figure;
   %     plot(t,dfOvF_rev_run(:,:,n)');
   %     xlabel('frames'); ylabel('df/f');
   %     title(['RUN reverse stimuli ROI', num2str(n)]);
   %     saveas(indi_revRun_fig,[image_dest '_reverseTrigAve_runROI' num2str(n)]);
   % end
   
end