runTriggerDura = 15;
for ii = 1: length(sessions)
    behav_dest = ['Z:\Analysis\WF_MovingDots_Analysis\behavioral_analysis\' days{ii}];
    behav_struct = load([behav_dest '\' days{ii} '_behavAnalysis.mat']);
    speed = behav_struct.speed; speed = double(speed);
    image_dest = [image_dest_base sessions{ii} '\' sessions{ii}];
    dfOvF_strct = load([image_dest, '_dfOvF_staybase.mat']);
    dfOvF = dfOvF_strct.dfOvF_staybase;
    frames_runTrigger_mat = behav_struct.frames_runTrigger_mat;

    runTrigSpeed = (speed(frames_runTrigger_mat))';
    ave_speed_runTrigger = mean(runTrigSpeed);
    ste_speed_runTrigger = std(runTrigSpeed,0,1)/sqrt(size(runTrigSpeed,1));
    
    % generate matrixes for df/f and plot----------------------------------
    dfOvF_runTrigger = zeros(size(frames_runTrigger_mat,2),size(frames_runTrigger_mat,1),size(dfOvF,1));
    x = -5:runTriggerDura-6;
    for n = 1:size(dfOvF,1)
        temp = dfOvF(n,:);
        dfOvF_runTrigger(:,:,n) = (temp(frames_runTrigger_mat))';
        set(0,'DefaultFigureVisible','off'); %doesn't show figures for the single sessions on the screen
        windows_fig(n) = figure; 
        %make every line start at the same height
        startMean = mean(dfOvF_runTrigger(:,1,n));
        diff = dfOvF_runTrigger(:,1,n)-startMean;
        diff_array = repmat(diff,1, size(dfOvF_runTrigger,2));
        dfOvF_runTrigger_plot = dfOvF_runTrigger(:,:,n) - diff_array;
        dfOvF_runTrigger_plot = dfOvF_runTrigger_plot';
        plot(x,dfOvF_runTrigger_plot);
        xlabel('frames');
        ylabel('df/f');
        set(gca,'XTick',x);
        vline(-1, 'k','running start');
        title(['df/f for each running window ROI' num2str(n)]); 
        saveas(windows_fig(n), [image_dest '_runTrigSessions_ROI' num2str(n)]);
    end
    
    ave_dfOvF_runTrigger = squeeze(mean(dfOvF_runTrigger,1));
    ste_dfOvF_runTrigger = squeeze(std(dfOvF_runTrigger,0,1)/sqrt(size(dfOvF_runTrigger,1)));
    
    set(0,'DefaultFigureVisible','on'); 
    mean_fig = figure; 
    subplot(2,1,1);hold on;
    errorbar(ave_dfOvF_runTrigger,ste_dfOvF_runTrigger,'linewidth', 1.5); hold on;
    %xlim([-5 10]);
    %ylim([-0.05 0.05]);
    vline(4, 'r','running start');
    ylabel('df/f'); legend;
    
    subplot(2,1,2);hold on;
    errorbar(ave_speed_runTrigger,ste_speed_runTrigger,'linewidth', 1.5); hold on;
    xlabel('frames');
    ylabel('speed');
    vline(4, 'r','running start');
    %xlim([-5 10]);
    
    supertitle(['run triggered average',sessions{ii}] );
    saveas(mean_fig, [image_dest '_runTrigAve']);
    save([image_dest '_imgAnalysis.mat' ],'ave_dfOvF_runTrigger', 'ave_speed_runTrigger', '-append');
    
end

