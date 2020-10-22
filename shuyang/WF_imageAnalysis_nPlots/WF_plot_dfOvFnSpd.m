%plotting fluorescence data for wide field imaging, moving dots behavior/running.
% need: df/f,output from behav_analysis_movingDots_WF .
% df/f &speed
% df/f vs. unique speed values
% ave df/f during stationary and running
%% SECTION - assign pathnames and datasets to be analyzed/written. 
clear;
sessions = {'190617_img1021_1','190617_img1023_1','190617_img1024_1',...
   '190617_img1027_2','190618_img1025_1','200321_img1042_1','200321_img1049_1',...
   '200321_img1063_1','200321_img1064_1'}; 
days = {'1021-190617_1','1023-190617_1','1024-190617_1',...
   '1027-190617_1','1025-190618_1','1042-200321_1','1049-200321_1',...
   '1063-200321_1','1064-200321_1'};
%sessions = {'190617_img1021_1'};
%days = {'1021-190617_1'};
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
    
    n1 = 6400;
    n2 = 6800;
    dfOvF_plot = dfOvF(:,(n1:n2));
    dfOvF_plot = dfOvF_plot';
    
    dfOvF_speed_fig = figure;
    x = ((n1:n2)/10)-(n1/10);
    yyaxis right;
    p = plot(x,dfOvF_plot(:,1),x,dfOvF_plot(:,2),x,dfOvF_plot(:,3),x,dfOvF_plot(:,4));
    p(1).Color = [0.4000 0.7608 0.6471]; p(1).LineStyle  = '-'; p(1).LineWidth = 1;
    p(2).Color = [0.9882 0.5529 0.3843]; p(2).LineStyle  = '-'; p(2).LineWidth = 1;
    p(3).Color = [0.5529 0.6275 0.7961]; p(3).LineStyle  = '-'; p(3).LineWidth = 1;
    p(4).Color = [0.9059 0.5412 0.7647]; p(4).LineStyle  = '-'; p(4).LineWidth = 1;
    %,'-','color',[0.4000 0.7608 0.6471],[0.9882 0.5529 0.3843],[0.5529 0.6275 0.7961]);
    %colororder();
    legend({'ROI1','ROI2','ROI3','ROI4'},'FontSize',5);legend('boxoff'); 
    ylabel('df/F'); hold on;
    ylim([-0.2 0.7]);hold on;
    yyaxis left;
    plot(x,speed(n1:n2),'color','k','LineWidth',1);
    ylabel('speed(cm/s)','Color','k'); ylim([-25 60]); 
    set(gca,'YColor','k');
    xlabel('time(s)');
    a = get(gca,'XTickLabel');
    set(gca,'XTickLabel',a,'FontSize',8);
    dfOvF_speed_fig.Units = 'centimeters';
    dfOvF_speed_fig.Position = [1 3 8 4];
    fig_name = ['eg_session_dfOvFnSpd_',sessions{i}];
    path = 'Z:\Analysis\figures\figure1_WF\';
    orient(dfOvF_speed_fig,'landscape')
    print(dfOvF_speed_fig,[path,fig_name],'-r600','-depsc');
    %saveas(dfOvF_speed_fig, [path, fig_name]);
    
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

