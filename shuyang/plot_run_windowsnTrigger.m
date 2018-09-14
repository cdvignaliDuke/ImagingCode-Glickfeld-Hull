%% SECTION - assign pathnames and datasets to be analyzed/written. 
clear;
%NEED TO UPDATE THIS SO IT ACCESSES SPREADSHEET INSTEAD OF JUST WRITING IN THE NAMES
sessions = {'180414_img1008_1'}; 
days = '1008-180414_1';
sessionID = '1008-180414';% this  variable name is confusing, this session ID is just tha date and the subject#, 
%there might be more than 1 sessions on a single subject on the same day
%bx_source     = ['Z:\Data\Behv_MovingDots\behavior_raw'];
%image_source_base  = ['Z:\Data\WF imaging\']; %location of permanently stored image files for retreiving meta data
image_dest_base    = ['Z:\Analysis\WF_MovingDots_Analysis\BxAndAnalysisOutputs\']; %stores the data on crash in the movingDots analysis folder
% behavior analysis results 
behav_dest = ['Z:\Analysis\WF_MovingDots_Analysis\behavioral_analysis\' days];
color_code = {'c','r','y','g'};


%% SECTION I: generate matrix of frames needed for plotting and load data needed for plotting
for ii = 1: length(sessions)
    image_dest = [image_dest_base sessions{ii} '\' sessions{ii}];
    dfOvF_strct = load([image_dest, '_dfOvF_staybase.mat']);
    dfOvF = dfOvF_strct.dfOvF_staybase;
    behav_struct = load([behav_dest '\' days '_behavAnalysis.mat']);
    frames_run_cell = behav_struct.frames_run_cell;
    speed = behav_struct.speed;
    [frames_befo_run,frames_aft_run,frames_runTrigger,frames_run_mat] = findFrames_runWindows (speed,frames_run_cell);
    
     save([behav_dest '\' sessionID '_' num2str(ii) '_behavAnalysis.mat' ],...
        'frames_befo_run','frames_aft_run','frames_runTrigger','frames_run_mat', '-append');
end


%% SECTION II:draw scatter plot for df/f of right before and right after every running window
for ii = 1:length(sessions)
    % calculate ave of dfOvF and put before, during run, after run together-----------------------------------------------
    dfOvF_befoRun = dfOvF(:,frames_befo_run);
    % this give you dfOvF of first, second, and third column in frames_befor_run.
    %ROI1 are now placed in a row(1-78:first column of frames_befo_run; 79-157: second column; 158-234: third column)
    dfOvF_aftRun = dfOvF(:,frames_aft_run);
    dfOvF_run = dfOvF(:,cell2mat(frames_run_cell));
    
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
    speed_befoRun = speed(:,frames_befo_run);
    speed_aftRun = speed(:,frames_aft_run);
   
    ave_speed_befoRun = mean(speed_befoRun);
    %tranpose it to each column is an ROI later
    ste_speed_befoRun = std(speed_befoRun)/sqrt(length(speed_befoRun));
    ave_speed_aftRun = mean(speed_aftRun);
    ste_speed_aftRun = std(speed_aftRun)/sqrt(length(speed_aftRun));
    
    %----------------------------------------------------------------------
    dfOvF_run_buffer = figure;
    plot(ave_befoRunaft','.-');
    hold on;
    errorbar(ave_befoRunaft',ste_befoRunaft');hold on;
    x = [1,2,3];
    %xlabel ('');
    set(gca,'XTick',x,'XTicklabel',{'right before','run','right after'});
    ylabel('df/f');
    title(['df/f right before and after running',sessions{ii}]);
    saveas(dfOvF_run_buffer, [image_dest '_dfOvF_runVsSurround']);
    save([behav_dest '\' sessionID '_' num2str(ii) '_behavAnalysis.mat' ],...
        'ave_dfOvF_befoRun','ave_dfOvF_aftRun','ste_dfOvF_befoRun','ste_dfOvF_aftRun',...
        'ave_speed_befoRun','ste_speed_befoRun','ave_speed_aftRun','ste_speed_aftRun','-append');
end


%% SECTION III: df/f before, after and during running, every 300ms.
for ii = 1: length(sessions)
    %frames_run_buffer_mat = behav_output.frames_run_buffer_mat;
    ave_dfOvF_befoRun = behav_struct.ave_dfOvF_befoRun;
    ave_dfOvF_aftRun = behav_struct.ave_dfOvF_aftRun;
    ste_dfOvF_befoRun = behav_struct.ste_dfOvF_befoRun;
    ste_dfOvF_aftRun = behav_struct.ste_dfOvF_aftRun;
    ave_speed_befoRun = behav_struct.ave_speed_befoRun;
    ste_speed_befoRun = behav_struct.ste_speed_befoRun;
    ave_speed_aftRun = behav_struct.ave_speed_aftRun;
    ste_speed_aftRun = behav_struct.ste_speed_aftRun;
    
    % calculate average speed of every 300ms ------------------------------
    speed_run_mat = nan(size(frames_run_mat,1), size(frames_run_mat,2));
    for v = 1:size(speed_run_mat,1)
        temp = frames_run_mat(v,~isnan(frames_run_mat(v,:)));
        speed_run_mat(v,1:size(temp,2)) = speed(temp);
    end
    speed_run_mat(speed_run_mat(v)==0)=nan;
    mean_columns = mean(speed_run_mat,'omitnan');
    bin = floor(length(mean_columns)/3);
    mean_run_every300ms = zeros(1,bin);
    ste_run_every300ms = zeros(1,bin);
    for x = 1:bin
        y = 3*x -2;
        mean_run_every300ms(x) = mean(mean_columns(y:y+2));
        ste_run_every300ms(x) = std(mean_columns(y:y+2))/sqrt(3);
    end
   % if length(mean_columns) can't be divided by 3 in total, you're throwing away the remainers.
   
    % calculate average fluorescence of every 300ms------------------------
    dfOvF_staybase_mat = nan(size(frames_run_mat,1), size(frames_run_mat,2));
    dfOvF_run_cell = {};
    mean_dfOvF_300ms = []; ste_dfOvF_300ms = [];
    for q = 1:size(dfOvF,1)
        for n = 1: size(dfOvF_staybase_mat,1)
            temp = frames_run_mat(n,~isnan(frames_run_mat(n,:)));
            dfOvF_run_cell{q}(n,1:size(temp,2)) = dfOvF(q,temp);
        end
        dfOvF_run_cell{q}(dfOvF_run_cell{q}==0)=nan;
        mean_columns = mean(dfOvF_run_cell{q},'omitnan');
        %calculate average for every 300ms
        bin = floor(length(mean_columns)/3);
        mean_every300ms = zeros(1,bin);
        ste_every300ms = zeros(1,bin);
        for x = 1:bin
            y = 3*x -2;
            mean_every300ms(x) = mean(mean_columns(y:y+2));
            ste_every300ms(x) = std(mean_columns(y:y+2))/sqrt(3);
        end
       mean_dfOvF_300ms(q,:) = mean_every300ms;
       ste_dfOvF_300ms(q,:) = ste_every300ms;
    end
        
    % plot ROIs and speed--------------------------------------------------
    dfOvF_plot = [ave_dfOvF_befoRun',mean_dfOvF_300ms,ave_dfOvF_aftRun'];
    ste_dfOvF_plot = [ste_dfOvF_befoRun',ste_dfOvF_300ms,ste_dfOvF_aftRun'];
    speed_plot = [ave_speed_befoRun,mean_run_every300ms,ave_speed_aftRun];
    ste_speed_plot = [ste_speed_befoRun,ste_run_every300ms,ste_speed_aftRun];
    
    dfOvF_run_300ms = figure;
    subplot(2,1,1);
    plot(dfOvF_plot', '.-');hold on;
    errorbar(dfOvF_plot',ste_dfOvF_plot');
    ylabel ('df/f');

    subplot(2,1,2);
    plot(speed_plot, '.-'); hold on;
    errorbar(speed_plot, ste_speed_plot);
    ylabel('speed');
    %xticks([ ])
    %xticklabels({'','','','','',''});
    supertitle(['df/f before, during, and after running',sessions{ii}]);
    saveas(dfOvF_run_300ms, [image_dest, '_dfOvF_runWindows']);
     
    %set(gca,'XTick',x,'XTicklabel',{'before run','1','','','','5','','','','',...
        %'10','','','','','15','','','','','20','','after run'});
    %'','','','25','','','', '','30',...
    % '','','','','35','','','','','40','','','','','45','','','','','50','','',...
    %   '','','55','','','','','60','','','','','65','','','','','70','','','','',...
    %  '75','','','','','80','','','','','85','','','','','90','','','','','','',...
    %  '','','','','','','','','','','','','','','','','',
 
    % legend(handle,'ROI1', 'ROI2');
    
end


%%  SECTION IV: run trigger ave
runTriggerDura = 15;
for ii = 1: length(sessions)
    runTrigSpeed = (speed(frames_runTrigger))';
    ave_speed_runTrigger = mean(runTrigSpeed);
    ste_speed_runTrigger = std(runTrigSpeed,0,1)/sqrt(size(runTrigSpeed,1));
    
    % generate matrixes for df/f and plot----------------------------------
    dfOvF_runTrigger = zeros(size(frames_runTrigger,2),size(frames_runTrigger,1),size(dfOvF,1));
    x = -5:runTriggerDura-6;
    for n = 1:size(dfOvF,1)
        temp = dfOvF(n,:);
        dfOvF_runTrigger(:,:,n) = (temp(frames_runTrigger))';
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
    
    mean_fig = figure;
    subplot(2,1,1);hold on;
    errorbar(ave_dfOvF_runTrigger,ste_dfOvF_runTrigger); hold on;
    %xlim([-5 10]);
    %ylim([-0.05 0.05]);
    vline(4, 'r','running start');
    ylabel('df/f');
    
    subplot(2,1,2);hold on;
    errorbar(ave_speed_runTrigger,ste_speed_runTrigger); hold on;
    xlabel('frames');
    ylabel('speed');
    vline(4, 'r','running start');
    %xlim([-5 10]);
    
    supertitle('run triggered average' );
    saveas(mean_fig, [image_dest '_runTrigAve']);
    
end
