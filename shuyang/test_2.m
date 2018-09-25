for ii = 1:length(sessions)
    image_dest = [image_dest_base sessions{ii} '\' sessions{ii}];
    behav_dest = ['Z:\Analysis\WF_MovingDots_Analysis\behavioral_analysis\' days{ii}];
    dfOvF_strct = load([image_dest, '_dfOvF_staybase.mat']);
    dfOvF = dfOvF_strct.dfOvF_staybase;
    behav_struct = load([behav_dest '\' days{ii} '_behavAnalysis.mat']);
    frames_run_cell = behav_struct.frames_run_cell; 
    frames_run = cell2mat(frames_run_cell);
    frames_befo_run = behav_struct.frames_befo_run_cell;
    frames_befo_run = cell2mat(frames_befo_run);
    frames_aft_run = behav_struct.frames_aft_run_cell;
    frames_aft_run = cell2mat(frames_aft_run);
    speed = behav_struct.speed; speed = double(speed);

    
    % calculate ave of dfOvF and put before, during run, after run together-----------------------------------------------
    dfOvF_befoRun = dfOvF(:,frames_befo_run);
    % this give you dfOvF of first, second, and third column in frames_befor_run.
    %ROI1 are now placed in a row(1-78:first column of frames_befo_run; 79-157: second column; 158-234: third column)
    dfOvF_aftRun = dfOvF(:,frames_aft_run);
    dfOvF_run = dfOvF(:,frames_run);
    
    ave_dfOvF_befoRun = mean(dfOvF_befoRun');
    %tranpose it to each column is an ROI later
    ste_dfOvF_befoRun = std(dfOvF_befoRun')/sqrt(length(dfOvF_befoRun'));
    ave_dfOvF_run = mean(dfOvF_run');
    ste_dfOvF_run = std(dfOvF_run')/sqrt(length(dfOvF_run'));
    ave_dfOvF_aftRun = mean(dfOvF_aftRun');
    ste_dfOvF_aftRun = std(dfOvF_aftRun')/sqrt(length(dfOvF_aftRun'));
    ave_befoRunaft = [ave_dfOvF_befoRun',ave_dfOvF_run',ave_dfOvF_aftRun'];
    ste_befoRunaft = [ste_dfOvF_befoRun',ste_dfOvF_run',ste_dfOvF_aftRun'];
    
    %do the same thing for speed-------------------------------------------
    speed_befoRun = speed(frames_befo_run);
    speed_aftRun = speed(frames_aft_run);
   
    ave_speed_befoRun = mean(speed_befoRun);
    %tranpose it to each column is an ROI later
    ste_speed_befoRun = std(speed_befoRun)/sqrt(length(speed_befoRun));
    ave_speed_aftRun = mean(speed_aftRun);
    ste_speed_aftRun = std(speed_aftRun)/sqrt(length(speed_aftRun));
    
    %----------------------------------------------------------------------
    dfOvF_run_buffer = figure;
    %plot(ave_befoRunaft','.-');
    %hold on;
    errorbar(ave_befoRunaft',ste_befoRunaft','linewidth', 1.5);hold on;
    x = [1,2,3];
    %xlabel ('');
    set(gca,'XTick',x,'XTicklabel',{'right before','run','right after'});
    ylabel('df/f');
    title(['df/f right before and after running',sessions{ii}]); legend;
    saveas(dfOvF_run_buffer, [image_dest '_dfOvF_runVsSurround']);
    save([behav_dest '\' sessionID{ii} '_' num2str(ii) '_behavAnalysis.mat' ],...
        'ave_dfOvF_befoRun','ave_dfOvF_aftRun','ste_dfOvF_befoRun','ste_dfOvF_aftRun',...
        'ave_speed_befoRun','ste_speed_befoRun','ave_speed_aftRun','ste_speed_aftRun','-append');
    save([image_dest '_imgAnalysis.mat' ],'ave_befoRunaft','-append');
end