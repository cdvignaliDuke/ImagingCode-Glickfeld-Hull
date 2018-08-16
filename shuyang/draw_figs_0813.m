%% SECTION - assign pathnames and datasets to be analyzed/written. 
clear;
%NEED TO UPDATE THIS SO IT ACCESSES SPREADSHEET INSTEAD OF JUST WRITING IN THE NAMES
sessions = {'180414_img1008_1'}; 
days = '1008-180414_1';
%bx_source     = ['Z:\Data\Behv_MovingDots\behavior_raw'];
%image_source_base  = ['Z:\Data\WF imaging\']; %location of permanently stored image files for retreiving meta data
image_dest_base    = ['Z:\Analysis\WF_MovingDots_Analysis\BxAndAnalysisOutputs\']; %stores the data on crash in the movingDots analysis folder
% behavior analysis results 
behav_dest = ['Z:\Analysis\WF_MovingDots_Analysis\behavioral_analysis\' days];
color_code = {'c','r','y','g'};

%% SECTION I: draw df/f &. speed
for i = 1:length(sessions)
    image_dest = [image_dest_base sessions{i} '\' sessions{i}];
    speed_struct = load([behav_dest '\' days '_behavAnalysis.mat']);
    speed = speed_struct.speed;
    cReverse_vec = speed_struct.cReverse_vec;
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
    speed_struct = load([behav_dest '\' days '_behavAnalysis.mat']);
    speed = speed_struct.speed;
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
    frames_states = load([behav_dest '\' days '_behavAnalysis.mat']);
    stay = frames_states.frames_stay_cell;
    bf  = frames_states.frames_bf_cell;
    run = frames_states.frames_run_cell;
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


%% SECTION IV: draw scatter plot for df/f of right before and right after every running window
for ii = 1: length(sessions)
    image_dest = [image_dest_base sessions{ii} '\' sessions{ii}];
    dfOvF_strct = load([image_dest, '_dfOvF_staybase.mat']);
    dfOvF = dfOvF_strct.dfOvF_staybase;
    frames_struct = load([behav_dest '\' days '_behavAnalysis.mat']);
    frames_run_cell = frames_struct.frames_run_cell;
    speed_struct = load([behav_dest '\' days '_behavAnalysis.mat']);
    speed = speed_struct.speed;
    [frames_befo_run,frames_aft_run,frames_runTrigger,frames_run_buffer_mat] = findFrames_runWindows (speed,frames_run_cell);
    
    dfOvF_run_buffer = figure;
    for k = 1: size(dfOvF,1)
        dfOvF_befoRun = dfOvF(k,cell2mat(frames_befo_run));
        dfOvF_aftRun = dfOvF(k,cell2mat(frames_aft_run));
        dfOvF_run = dfOvF(k,cell2mat(frames_run_cell));
        ave_dfOvF_befoRun = mean(dfOvF_befoRun);
        ste_dfOvF_befoRun = std(dfOvF_befoRun)/sqrt(length(dfOvF_befoRun));
        ave_dfOvF_run = mean(dfOvF_run);
        ste_dfOvF_run = std(dfOvF_run)/sqrt(length(dfOvF_run));
        ave_dfOvF_aftRun = mean(dfOvF_aftRun);
        ste_dfOvF_aftRun = std(dfOvF_aftRun)/sqrt(length(dfOvF_aftRun));
        x = [1,2,3];
        scatter(x,[ave_dfOvF_befoRun, ave_dfOvF_run, ave_dfOvF_aftRun]);
        hold on;
        errorbar(x,[ave_dfOvF_befoRun, ave_dfOvF_run, ave_dfOvF_aftRun],...
            [ste_dfOvF_befoRun,ste_dfOvF_run,ste_dfOvF_aftRun]);hold on;
    end
    %xlabel ('');
    set(gca,'XTick',x,'XTicklabel',{'right before','run','right after'});
    ylabel('df/f');
    title(['df/f right before and after running',sessions{ii}]);
    saveas(dfOvF_run_buffer, [image_dest '_dfOvF_runVsSurround']);
end
 


%% SECTION V: df/f before, after and during running, every 300ms.
for ii = 1: length(sessions)
    image_dest = [image_dest_base sessions{ii} '\' sessions{ii}];
    dfOvF_strct = load([image_dest, '_dfOvF_staybase.mat']);
    dfOvF = dfOvF_strct.dfOvF_staybase;
    behav_output = load([behav_dest '\' days '_behavAnalysis.mat']);
    %frames_run_buffer_mat = behav_output.frames_run_buffer_mat;
    speed = behav_output.speed;
    % frames run mat is a matrix for 'run buffer', which includes right
    % before run, run, and right after run.
    
    %% calculate average speed of every 300ms
    speed_run_mat = nan(size(frames_run_buffer_mat,1), size(frames_run_buffer_mat,2));
    for v = 1:size(speed_run_mat,1)
        temp = frames_run_buffer_mat(v,~isnan(frames_run_buffer_mat(v,:)));
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
    mean_run_every300ms(bin) = mean(mean_columns(length(mean_columns)-2:length(mean_columns)));
    ste_run_every300ms(bin) = std(mean_columns(length(mean_columns)-2:length(mean_columns)))/sqrt(3);
    %% calculate average fluorescence of every 300ms
    dfOvF_staybase_mat = nan(size(frames_run_buffer_mat,1), size(frames_run_buffer_mat,2));
    dfOvF_run_cell = {};
    
    dfOvF_run_300ms = figure; p = panel(); p.pack(2,1);
    for q = 1:size(dfOvF,1)
        for n = 1: size(dfOvF_staybase_mat,1)
            temp = frames_run_buffer_mat(n,~isnan(frames_run_buffer_mat(n,:)));
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
        mean_every300ms(bin) = mean(mean_columns(length(mean_columns)-2:length(mean_columns)));
        ste_every300ms(bin) = std(mean_columns(length(mean_columns)-2:length(mean_columns)))/sqrt(3);
        
        %% plot ROIs and speed
        x = (-2: 3: (bin-1)*3);
        p(1,1).select();
        hold on;
        scatter(x,mean_every300ms);hold on;
        errorbar(x,mean_every300ms,ste_every300ms);
        ylabel ('df/f');
        xlim([-2 max(x)]);
        set(gca,'XTick',x,'XTicklabel',{'before run','1','','','','5','','','','',...
          '10','','','','','15','','','','','20','','after run'});
      %'','','','25','','','', '','30',...
       % '','','','','35','','','','','40','','','','','45','','','','','50','','',...
      %   '','','55','','','','','60','','','','','65','','','','','70','','','','',...
       %  '75','','','','','80','','','','','85','','','','','90','','','','','','',...
      %  '','','','','','','','','','','','','','','','','',
          
           
           %set(gca,'XTick',x,'XTicklabel',{'-3~0','0~3','3~6','6~9','9~12','12~15',...
           % '15~18','18~21','21~24','24~27','27~30','30~33','33~36','36~39',...
           % '39~42','42~45','45~48','48~51','51~54','54~57','57~60','60~63',...
           % '63~66','66~69','69~72','72~75','75~78','78~81','81~84','84~87',...
           % '87~90','90~93','93~96','96~99','99~102','102~105','105~108',...
           % '108~111','111~114','114~117','117~120','120~123','123~126',...
           % '126~129','129~132','after run'});
        p(2,1).select();
        scatter(x,mean_run_every300ms);
        errorbar(x,mean_run_every300ms,ste_run_every300ms);
        xlabel ('frames');
        xlim([-2 max(x)]);
        ylabel ('speed');
        set(gca,'XTick',x,'XTicklabel',{'before run','1','','','','5','','','','',...
          '10','','','','','15','','','','','20','','after run'});
         % '','','','25','','','','','30',...
       % '','','','','35','','','','','40','','','','','45','','','','','50','','',...
        %'','','55','','','','','60','','','','','65','','','','','70','','','','',...
       % '75','','','','','80','','','','','85','','','','','90','',
         %
        
        
        %set(gca,'XTick',x,'XTicklabel',{'-3~0','0~3','3~6','6~9','9~12','12~15',...
           % '15~18','18~21','21~24','24~27','27~30','30~33','33~36','36~39',...
           % '39~42','42~45','45~48','48~51','51~54','54~57','57~60','60~63',...
           % '63~66','66~69','69~72','72~75','75~78','78~81','81~84','84~87',...
           % '87~90','90~93','93~96','96~99','99~102','102~105','105~108',...
           % '108~111','111~114','114~117','117~120','120~123','123~126',...
           % '126~129','129~132','after run'});
        supertitle(['df/f before, during, and after running',sessions{ii}]);
        saveas(dfOvF_run_300ms, [image_dest, '_dfOvF_runWindows']);
        
    end
    % legend(handle,'ROI1', 'ROI2');

end

