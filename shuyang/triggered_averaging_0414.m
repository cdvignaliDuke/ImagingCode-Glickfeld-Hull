%% SECTION ONE - assign pathnames and datasets to be analyzed/written. 
clear;
%NEED TO UPDATE THIS SO IT ACCESSES SPREADSHEET INSTEAD OF JUST WRITING IN THE NAMES
sessions = {'180423_img1010_1'}; 
days = '1010-180423_1';
bx_source     = ['Z:\Data\Behv_MovingDots\behavior_raw'];
image_source_base  = ['Z:\Data\WF imaging\']; %location of permanently stored image files for retreiving meta data
image_dest_base    = ['Z:\Analysis\WF_MovingDots_Analysis\BxAndAnalysisOutputs\']; %stores the data on crash in the movingDots analysis folder
% behavior analysis results 
behav_dest = ['Z:\Analysis\WF_MovingDots_Analysis\behavioral_analysis\' days];
color_code = {'b','r','k','c'};

%% run trigger averaging
runTriggerDura = 15;
for ii = 1: length(sessions)
    %% load things
    image_dest = [image_dest_base sessions{ii} '\' sessions{ii}];
    dfOvF_strct_staybase = load([image_dest, '_dfOvF_staybase.mat']);
    dfOvF_staybase = dfOvF_strct_staybase.dfOvF_staybase;
    behav_output = load([behav_dest '\' days '_behavAnalysis.mat']);
    speed = behav_output.speed;
    frames_runTrigger = behav_output.frames_runTrigger;
    %% find speed
    tempSpeed = speed(frames_runTrigger);
    speed_runTrigger = reshape(tempSpeed, runTriggerDura,length(frames_runTrigger)/runTriggerDura)';

    %% generate matrixes for df/f and plot
     
    ave_speed_runTrigger = mean(speed_runTrigger,1);
    speed_runTrigger = double(speed_runTrigger);
    ste_speed_runTrigger = std(speed_runTrigger,0,1)/sqrt(size(speed_runTrigger,1));
    ave_dfOvF_runTrigger = [];
    ste_dfOvF_runTrigger = [];
    x = -5:runTriggerDura-6;

    for n = 1:size(dfOvF_staybase,1)
        tempdfOvf = dfOvF_staybase(n,frames_runTrigger);
        dfOvF_runTrigger = reshape(tempdfOvf, runTriggerDura,length(frames_runTrigger)/runTriggerDura)';
        ave_dfOvF_runTrigger = [ave_dfOvF_runTrigger;mean(dfOvF_runTrigger,1)];
        ste_dfOvF_runTrigger = [ste_dfOvF_runTrigger;std(dfOvF_runTrigger,0,1)/sqrt(size(dfOvF_runTrigger,1))];
        
        windows_fig(n) = figure;
        %make every line start at the same height
       startMean = mean(dfOvF_runTrigger(:,1));
       diff = dfOvF_runTrigger(:,1)-startMean;
       diff_array = repmat(diff,1,size(dfOvF_runTrigger,2));
       dfOvF_runTrigger_plot = dfOvF_runTrigger - diff_array;
       dfOvF_runTrigger_plot = dfOvF_runTrigger_plot';
       plot(x,dfOvF_runTrigger_plot);
        xlabel('frames');
        ylabel('df/f');
        set(gca,'XTick',x);
        vline(-1, 'k','running start');
        title(['df/f for each running window ROI' num2str(n)]);
        saveas(windows_fig(n), [image_dest '_runTrigSessions_ROI' num2str(n)]);
    end
    
    mean_fig = figure;
    for n = 1: size(ave_dfOvF_runTrigger,1)
        subplot(2,1,1);hold on;
        scatter(x,ave_dfOvF_runTrigger(n,:),color_code{n}); hold on;
        errorbar(x,ave_dfOvF_runTrigger(n,:),ste_dfOvF_runTrigger(n,:),'color',color_code{n}); hold on;
        xlim([-5 10]);
        %ylim([-0.05 0.05]);
        ylabel('df/f'); 
    end
    subplot(2,1,2);hold on;
    scatter(x,ave_speed_runTrigger,'m');hold on;
    errorbar(x,ave_speed_runTrigger,ste_speed_runTrigger,'color', 'm'); hold on;
    xlabel('frames');
    ylabel('speed');
    xlim([-5 10]);
    
    supertitle('run triggered average' );
    subplot(2,1,1);
    vline(0, 'r','running start');
    subplot(2,1,2);
    vline(0, 'r','running start');
    saveas(mean_fig, [image_dest '_runTrigAve']);
    
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% % draw rawF run triggered average
    %rawF_struct = load([image_dest, '_raw_F.mat']);
    %rawF = rawF_struct.data_tc;
%     for n = 1:size(rawF,1)
%        tempF = rawF(n,frames_runTrigger);
%        rawF_runTrigger = reshape(tempF, runTriggerDura,length(frames_runTrigger)/runTriggerDura)';
%        x = -5:runTriggerDura-6;
%
%        mean_fig{n} = figure;
%        ave_rawF_runTrigger = mean(rawF_runTrigger,1);
%        ave_speed_runTrigger = mean(speed_runTrigger,1);
%        x = -5:runTriggerDura-6;
%        subplot(2,1,1);
%        plot(x,ave_rawF_runTrigger);
%        subplot(2,1,2);
%        plot(x,ave_speed_runTrigger);
%
%        windows_fig{n} = figure;
%        for f = 1: size(rawF_runTrigger,1)
%            hold on;
%            plot(x, rawF_runTrigger(f,:));
%        end
%    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
