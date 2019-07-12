%plotting fluorescence data for wide field imaging, moving dots behavior.
% df/f &speed
% df/f vs. unique speed values
% ave df/f during stationary and running
%% SECTION - assign pathnames and datasets to be analyzed/written. 
clear;
%NEED TO UPDATE THIS SO IT ACCESSES SPREADSHEET INSTEAD OF JUST WRITING IN THE NAMES
sessions = {'190409_img1018_1'};
days = {'1018-190409_1'};
sessionID = {'1018-190409'};
% this  variable name is confusing, this session ID is just tha date and the subject#, 
%there might be more than 1 sessions on a single subject on the same day
%bx_source     = ['Z:\Data\Behv_MovingDots\behavior_raw'];
%image_source_base  = ['Z:\Data\WF imaging\']; %location of permanently stored image files for retreiving meta data
image_dest_base = ['Z:\Analysis\WF_MovingDots_Analysis\BxAndAnalysisOutputs\']; %stores the data on crash in the movingDots analysis folder
% behavior analysis results 
color_code = {'r','g','m','y'};

%% load variables
%for i = 1:length(sessions)
   % image_dest = [image_dest_base sessions{i} '\' sessions{i}];
   % behav_struct = load([behav_dest '\' days '_behavAnalysis.mat']);
   % speed = behav_struct.speed;
   % cReverse_vec = behav_struct.cReverse_vec;
   % dfOvF_struct = load([image_dest, '_dfOvF_staybase.mat']);
   % dfOvF = dfOvF_struct.dfOvF_staybase; 
%end

%% SECTION I: draw df/f &. speed
for i = 1:length(sessions)
    image_dest = [image_dest_base sessions{i} '\' sessions{i}];
    behav_dest = ['Z:\Analysis\WF_MovingDots_Analysis\behavioral_analysis\' days{i}];
    behav_struct = load([behav_dest '\' days{i} '_behavAnalysis.mat']);
    speed = behav_struct.speed;
    cReverse_vec = behav_struct.cReverse_vec;
    dfOvF_struct = load([image_dest, '_dfOvF_staybase.mat']);
    dfOvF = dfOvF_struct.dfOvF_staybase;
    
    dfOvF_speed_fig = figure;
    for n = 1:size(dfOvF,1)
        dfOvF_plot = dfOvF(n,1:length(speed));
        time = (0:(length(speed)-1));
        hold on;
        %plot only part of the session so that the plot is not too crowded
        %[hAx,hline1,hline2(n)] = plotyy(time(1000:2000),speed(1000:2000),time(1000:2000),dfOvF_plot(1000:2000));
        [hAx,hline1,hline2(n)] = plotyy(time,speed,time,dfOvF_plot);
        set(hline1,'color', 'b'); 
        set(hline2(n),'color',color_code{n});
        %legend;%('speed','ROI1','ROI2','Location','northeast');
        %ylim(hAx(2),[-1.5 0.5]);
    end 
    if isempty(cReverse_vec) == 0
    vline(cReverse_vec, 'r');
    end
    xlabel ('frames');
    ylabel(hAx(1),'speed');
    ylabel(hAx(2),'df/f');
    ylim(hAx(1),[0,max(speed)]);
    title(['df/f and speed',sessions{i}])
    saveas(dfOvF_speed_fig, [image_dest, '_dfOvF_and_speed']);
    
end


%% SECTION II : draw mean df/f vs. every speed value
for i = 1: length(sessions)
    image_dest = [image_dest_base sessions{i} '\' sessions{i}];
    behav_dest = ['Z:\Analysis\WF_MovingDots_Analysis\behavioral_analysis\' days{i}];
    behav_struct = load([behav_dest '\' days{i} '_behavAnalysis.mat']);
    speed = behav_struct.speed;
    cReverse_vec = behav_struct.cReverse_vec;
    dfOvF_struct = load([image_dest, '_dfOvF_staybase.mat']);
    dfOvF = dfOvF_struct.dfOvF_staybase;
    
    isp = unique(speed);
    dfOvF_spd = []; dfOvF_spdmean = []; dfOvF_spdste = [];
    for k = 1:length(isp)
        dfOvF_spd = dfOvF(:,speed==isp(k));
        dfOvF_spdmean(:,k)  = mean(dfOvF_spd,2);
        dfOvF_spdste(:,k) = std(dfOvF_spd,[],2)/sqrt(length(dfOvF_spd));
    end
    isp_plot =  repmat(isp,size(dfOvF_spdmean,1),1);
    dfOvF_ave_vs_spd = figure;
    errorbar(isp_plot',dfOvF_spdmean',dfOvF_spdste','.','LineStyle','-','linewidth', 1.25,'MarkerSize',20);
    legend('ROI1','ROI2'); hold on;
    %if do errorbar(y,ste y), it can do multiple lines at once. But if do
    %errorbar(x,y,ste y), size of x and y must match(if your y contains n lines, even though the x for all lines are the same, x must have n lines too. 
    xlabel ('speed');
    ylabel('ave df/f');
    title(['df/f vs. speed',sessions{i}]); 
    saveas(dfOvF_ave_vs_spd, [image_dest, '_dfOvf_vs_speed']);
    save([image_dest '_imgAnalysis.mat' ],'isp', 'dfOvF_spdmean');

end

%% SECTION III: draw ave df/f for stay and run 
for ii = 1: length(sessions)
    image_dest = [image_dest_base sessions{ii} '\' sessions{ii}];
    behav_dest = ['Z:\Analysis\WF_MovingDots_Analysis\behavioral_analysis\' days{ii}];
    behav_struct = load([behav_dest '\' days{ii} '_behavAnalysis.mat']);
    speed = behav_struct.speed;
    cReverse_vec = behav_struct.cReverse_vec;
    dfOvF_struct = load([image_dest, '_dfOvF_staybase.mat']);
    dfOvF = dfOvF_struct.dfOvF_staybase;
    
    %data_tc is the fluorescence data for this session, and it will be a
    %n*num of frames double, n=number of ROIs.
    stay = behav_struct.frames_stay_cell;
   % bf  = behav_struct.frames_bf_cell;
    run = behav_struct.frames_run_cell;
    ave_dfOvF_stay = []; ste_dfOvF_stay = [];
   % ave_dfOvF_bf = []; ste_dfOvF_bf = []; 
    ave_dfOvF_run = []; ste_dfOvF_run = [];
    
    for n = 1:size(dfOvF,1)
        dfOvF_stay = dfOvF(n,cell2mat(stay));
       % dfOvF_bf = dfOvF(n,cell2mat(bf));
        dfOvF_run = dfOvF(n,cell2mat(run));
        ave_dfOvF_stay(n,:) = mean(dfOvF_stay);
        ste_dfOvF_stay(n,:) = std(dfOvF_stay)/sqrt(length(dfOvF_stay));
        %ave_dfOvF_bf(n,:) = mean(dfOvF_bf);
       % ste_dfOvF_bf(n,:) = std(dfOvF_bf)/sqrt(length(dfOvF_bf));
        ave_dfOvF_run(n,:) = mean(dfOvF_run);
        ste_dfOvF_run(n,:) = std(dfOvF_run)/sqrt(length(dfOvF_run)); 
    end
    ave_dfOvF_behav = [ave_dfOvF_stay,ave_dfOvF_run];
    ste_dfOvF_behav = [ste_dfOvF_stay,ste_dfOvF_run];
    
    x = [1,2]; x_plot = repmat(x,size(dfOvF,1),1);
    dfOvF_behavStates = figure;
    errorbar(x_plot',ave_dfOvF_behav',ste_dfOvF_behav','.','LineStyle','-','linewidth', 1.25,'MarkerSize',20);
    xlabel ('behavioral state');xlim([0.5 2.5]);
    set(gca,'XTick',x,'XTicklabel',{'stationary','run'});
    ylabel('df/f');
    title(['df/f for each beahvioral state',sessions{ii}]); legend('ROI1','ROI2');
    saveas(dfOvF_behavStates, [image_dest '_dfOvF_behavStates_2']);
    save([image_dest '_imgAnalysis.mat' ],'ave_dfOvF_behav', '-append');
end

