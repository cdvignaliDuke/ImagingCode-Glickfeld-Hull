%% SECTION - assign pathnames and datasets to be analyzed/written. 
clear;
%NEED TO UPDATE THIS SO IT ACCESSES SPREADSHEET INSTEAD OF JUST WRITING IN THE NAMES
sessions = {'180414_img1005_1','180414_img1007_1','180414_img1008_1','180417_img1008_1','180419_img1008_1',...
    '180423_img1010_1','180417_img1005_1','180419_img1005_1','180419_img1007_1','180423_img1005_1',...
    '180424_img1008_1','180425_img1008_1','180428_img1008_1','180429_img1008_1','180430_img1005_1',...
    '180430_img1007_1','180430_img1008_1','180430_img1010_1','180505_img1007_1','180505_img1008_1','180505_img1010_1'}; 

days = {'1005-180414_1','1007-180414_1','1008-180414_1','1008-180417_1','1008-180419_1','1010-180423_1',...
    '1005-180417_1','1005-180419_1','1007-180419_1','1005-180423_1','1008-180424_1','1008-180425_1','1008-180428_1',...
    '1008-180429_1','1005-180430_1','1007-180430_1','1008-180430_1','1010-180430_1','1007-180505_1','1008-180505_1',...
    '1010-180505_1'};
sessionID = {'1005-180414','1007-180414','1008-180414','1008-180417','1008-180419','1010-180423',...
    '1005-180417','1005-180419','1007-180419','1005-180423','1008-180424','1008-180425','1008-180428',...
    '1008-180429','1005-180430','1007-180430','1008-180430','1010-180430','1007-180505','1008-180505','1010-180505'};
% this  variable name is confusing, this session ID is just tha date and the subject#, 
%there might be more than 1 sessions on a single subject on the same day
%bx_source     = ['Z:\Data\Behv_MovingDots\behavior_raw'];
%image_source_base  = ['Z:\Data\WF imaging\']; %location of permanently stored image files for retreiving meta data
image_dest_base = ['Z:\Analysis\WF_MovingDots_Analysis\BxAndAnalysisOutputs\']; %stores the data on crash in the movingDots analysis folder
% behavior analysis results 
color_code = {'c','r','y','g'};

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
    behav_struct = load([behav_dest '\' days '_behavAnalysis.mat']);
    speed = behav_struct.speed;
    cReverse_vec = behav_struct.cReverse_vec;
    dfOvF_struct = load([image_dest, '_dfOvF_staybase.mat']);
    dfOvF = dfOvF_struct.dfOvF_staybase;
    
    dfOvF_speed_fig = figure;
    for n = 1:size(dfOvF,1)
        dfOvF_plot = dfOvF(n,1:length(speed));
        time = (0:(length(speed)-1));
        hold on;
        [hAx,hline1,hline2(n)] = plotyy(time,speed,time,dfOvF_plot);
        set(hline1,'color', 'b');
        set(hline2(n),'color',color_code{n});
        ylim(hAx(2),[-1.5 0.5]);
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
    errorbar(isp_plot',dfOvF_spdmean',dfOvF_spdste','linewidth', 1.5);
    %if do errorbar(y,ste y), it can do multiple lines at once. But if do
    %errorbar(x,y,ste y), size of x and y must match(if your y contains n lines, even though the x for all lines are the same, x must have n lines too. 
    xlabel ('speed');
    ylabel('ave df/f');
    title(['df/f vs. speed',sessions{i}]); legend;
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
    errorbar(x_plot',ave_dfOvF_behav',ste_dfOvF_behav','linewidth', 2);
    xlabel ('behavioral state');xlim([0.5 2.5]);
    set(gca,'XTick',x,'XTicklabel',{'stationary','run'});
    ylabel('df/f');
    title(['df/f for each beahvioral state',sessions{ii}]); legend;
    saveas(dfOvF_behavStates, [image_dest '_dfOvF_behavStates_2']);
    save([image_dest '_imgAnalysis.mat' ],'ave_dfOvF_behav', '-append');
end
