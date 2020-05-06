%plotting fluorescence data for wide field imaging, moving dots behavior/running.
% need: df/f,output from behav_analysis_movingDots_WF .
% df/f &speed
% df/f vs. unique speed values
% ave df/f during stationary and running
%% SECTION - assign pathnames and datasets to be analyzed/written. 
clear;
%NEED TO UPDATE THIS SO IT ACCESSES SPREADSHEET INSTEAD OF JUST WRITING IN THE NAMES
%sessions = {'190617_img1021_1','190617_img1023_1','190617_img1024_1',...
%    '190617_img1027_2','190618_img1025_1','200321_img1042_1','200321_img1049_1',...
%    '200321_img1063_1','200321_img1064_1'}; 
%days = {'1021-190617_1','1023-190617_1','1024-190617_1',...
%    '1027-190617_1','1025-190618_1','1042-200321_1','1049-200321_1',...
%    '1063-200321_1','1064-200321_1'};
sessions = {'190618_img1025_1'};
days = {'1025-190618_1'};
%there might be more than 1 sessions on a single subject on the same day
image_dest_base = 'Z:\Analysis\WF_MovingDots_Analysis\BxAndAnalysisOutputs\'; %stores the data on crash in the movingDots analysis folder
% behavior analysis results 
color_code = {'r','g','m','y','b'};


%% SECTION I: draw df/f &. speed
for i = 1:length(sessions)
    image_dest = [image_dest_base sessions{i} '\' sessions{i}];
    behav_dest = ['Z:\Analysis\WF_MovingDots_Analysis\behavioral_analysis\' days{i}];
    behav_struct = load([behav_dest '\' days{i} '_behavAnalysis.mat']);
    speed = behav_struct.speed;
    cReverse_vec = behav_struct.cReverse_vec;
    dfOvF_struct = load([image_dest, '_dfOvF_btmbaseline.mat']);
    dfOvF = dfOvF_struct.dfOvF_btmbase;
    
    dfOvF_speed_fig = figure;
    for n = 1:4 %size(dfOvF,1)
        dfOvF_plot = dfOvF(n,1:length(speed));
        time = (0:(length(speed)-1));
        hold on;
        %plot only part of the session so that the plot is not too crowded
        n1 = 1100;
        n2 = 1200;
        x = (n1:n2);
        %x = (1:(n2-n1+1))/10; %make x axis to second instead of frame
        [hAx,hline1,hline2(n)] = plotyy(x,speed(n1:n2)*2*3.1415926*7.5/128,x,dfOvF_plot(n1:n2));
        %[hAx,hline1,hline2(n)] = plotyy(time,speed,time,dfOvF_plot);
        set(hline1,'color', 'k');
        set(hline2(n),'color',color_code{n});
        %ylim(hAx(2),[0 0.5]);  
    end
    
    if isempty(cReverse_vec) == 0
    vline(cReverse_vec, 'r');
    end
    xlabel ('time(s)');
    ylabel(hAx(1),'speed(cm/s)');
    ylabel(hAx(2),'df/f');
    %xlim([0 max(x)]);
    %ylim(hAx(2),[0,0.5]);
    %ylim(hAx(1),[-10 25]);
    %ylim(hAx(1),[min(speed),max(speed)]);
    title(['df/f and speed',sessions{i}]);
    % legend('speed','ROI1','ROI2','ROI3','ROI4');
    saveas(dfOvF_speed_fig, [image_dest, '_dfOvF_and_speed']);
    
end


%% SECTION II : draw mean df/f vs. every speed value
for i = 1: length(sessions)
    image_dest = [image_dest_base sessions{i} '\' sessions{i}];
    behav_dest = ['Z:\Analysis\WF_MovingDots_Analysis\behavioral_analysis\' days{i}];
    behav_struct = load([behav_dest '\' days{i} '_behavAnalysis.mat']);
    speed = behav_struct.speed;
    dfOvF_struct = load([image_dest, '_dfOvF_btmbaseline.mat']);
    dfOvF = dfOvF_struct.dfOvF_btmbase;
    
    isp = unique(speed);
    dfOvF_spd = []; dfOvF_spdmean = []; dfOvF_spdste = [];
    for k = 1:length(isp)
        dfOvF_spd = dfOvF(:,speed==isp(k));
        dfOvF_spdmean(:,k)  = mean(dfOvF_spd,2); % average across frames, ROI*number of unique speeds
        dfOvF_spdste(:,k) = std(dfOvF_spd,[],2)/sqrt(length(dfOvF_spd));
    end
    isp_plot =  repmat(isp,size(dfOvF_spdmean,1),1);
    dfOvF_ave_vs_spd = figure;
    errorbar(isp_plot'*2*3.1415926*7.5/128,dfOvF_spdmean',dfOvF_spdste','.','LineStyle','-','linewidth', 1.25,'MarkerSize',20);
    legend('ROI1','ROI2','ROI3','ROI4','ROI5'); hold on;
    %if do errorbar(y,ste y), it can do multiple lines at once. But if do
    %errorbar(x,y,ste y), size of x and y must match(if your y contains n lines, even though the x for all lines are the same, x must have n lines too. 
    xlabel ('speed(cm/s)');
    ylabel('ave df/f');
    xlim([0 max(isp)*2*3.1415926*7.5/128]);
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
    dfOvF_struct = load([image_dest, '_dfOvF_btmbaseline.mat']);
    dfOvF = dfOvF_struct.dfOvF_btmbase;
    
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
    xlabel ('behavioral state');xlim([0.5 2.5]); axis square;
    set(gca,'XTick',x,'XTicklabel',{'stationary','run'});
    ylabel('df/f');
    title(['df/f for each beahvioral state',sessions{ii}]); legend('ROI1','ROI2','ROI3','ROI4','ROI5');
    saveas(dfOvF_behavStates, [image_dest '_dfOvF_behavStates']);
    save([image_dest '_imgAnalysis.mat' ],'ave_dfOvF_behav', '-append');
end

