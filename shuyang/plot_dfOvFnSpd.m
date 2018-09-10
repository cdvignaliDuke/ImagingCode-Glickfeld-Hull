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

%% SECTION I: draw df/f &. speed
for i = 1:length(sessions)
    image_dest = [image_dest_base sessions{i} '\' sessions{i}];
    behav_struct = load([behav_dest '\' days '_behavAnalysis.mat']);
    speed = behav_struct.speed;
    cReverse_vec = behav_struct.cReverse_vec;
    dfOvF_struct = load([image_dest, '_dfOvF_staybase.mat']);
    dfOvF = dfOvF_struct.dfOvF_staybase;
    dfOvF_speed_fig = figure;
    for n = 1:size(dfOvF,1);
        dfOvF_plot = dfOvF(n,1:length(speed));
        time = (0:(length(speed)-1));
        hold on;
        [hAx,hline1,hline2(n)] = plotyy(time,speed,time,dfOvF_plot);
        set(hline1,'color', 'b');
        set(hline2(n),'color',color_code{n});
        ylim(hAx(2),[-1.5 0.5]);
    end 
    if isempty(cReverse_vec) == 0;
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
for i = length(sessions)
    image_dest = [image_dest_base sessions{i} '\' sessions{i}];
    behav_struct = load([behav_dest '\' days '_behavAnalysis.mat']);
    speed = behav_struct.speed;
    dfOvF_struct = load([image_dest, '_dfOvF_staybase.mat']);
    dfOvF = dfOvF_struct.dfOvF_staybase;
    isp = unique(speed);
    for k = 1:length(isp)
        dfOvF_spd = dfOvF(:,speed==isp(k));
        dfOvF_spdmean(:,k)  = mean(dfOvF_spd,2);
        dfOvF_spdste(:,k) = std(dfOvF_spd,[],2)/sqrt(length(dfOvF_spd));
    end
    dfOvF_ave_vs_spd = figure;
    for n = 1:size(dfOvF,1);
        scatter(isp,dfOvF_spdmean(n,:), 'MarkerFaceColor', color_code{n});
        hold on;
        errorbar(isp,dfOvF_spdmean(n,:),dfOvF_spdste(n,:), 'Color', color_code{n});
        
    end
    xlabel ('speed');
    ylabel('ave df/f');
    title(['df/f vs. speed',sessions{i}]);
    saveas(dfOvF_ave_vs_spd, [image_dest, '_dfOvf_vs_speed']);
end

%% SECTION III: draw ave df/f for stay, bf, and run 
for ii = 1: length(sessions)
    image_dest = [image_dest_base sessions{ii} '\' sessions{ii}];
    dfOvF_strct = load([image_dest, '_dfOvF_staybase.mat']);
    dfOvF = dfOvF_strct.dfOvF_staybase;
    %data_tc is the fluorescence data for this session, and it will be a
    %n*num of frames double, n=number of ROIs.
    behav_struct = load([behav_dest '\' days '_behavAnalysis.mat']);
    stay = behav_struct.frames_stay_cell;
    bf  = behav_struct.frames_bf_cell;
    run = behav_struct.frames_run_cell;
    dfOvF_behavStates = figure;
    for n = 1:size(dfOvF,1);
        dfOvF_stay = dfOvF(n,cell2mat(stay));
        dfOvF_bf = dfOvF(n,cell2mat(bf));
        dfOvF_run = dfOvF(n,cell2mat(run));
        ave_dfOvF_stay = mean(dfOvF_stay);
        ste_dfOvF_stay = std(dfOvF_stay)/sqrt(length(dfOvF_stay));
        ave_dfOvF_bf = mean(dfOvF_bf);
        ste_dfOvF_bf = std(dfOvF_bf)/sqrt(length(dfOvF_bf));
        ave_dfOvF_run = mean(dfOvF_run);
        ste_dfOvF_run = std(dfOvF_run)/sqrt(length(dfOvF_run)); 
        x = [1,2,3];
        scatter(x,[ave_dfOvF_stay, ave_dfOvF_bf, ave_dfOvF_run]);
        hold on;
        errorbar(x,[ave_dfOvF_stay, ave_dfOvF_bf, ave_dfOvF_run],...
            [ste_dfOvF_stay,ste_dfOvF_bf,ste_dfOvF_run]);
    end
    xlabel ('behavioral state');
    set(gca,'XTick',x,'XTicklabel',{'stationary','rigidity','run'});
    ylabel('df/f');
    title(['df/f for each beahvioral state',sessions{ii}]);
    saveas(dfOvF_behavStates, [image_dest '_dfOvF_behavStates']);
end

