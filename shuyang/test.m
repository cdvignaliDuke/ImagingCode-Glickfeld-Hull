%% SECTION V: df/f before, after and during running, every 300ms.
for ii = 1: length(sessions)
    image_dest = [image_dest_base sessions{ii} '\' sessions{ii}];
    dfOvF_strct = load([image_dest, '_dfOvF_staybase.mat']);
    dfOvF = dfOvF_strct.dfOvF_staybase;
    behav_struct = load([behav_dest '\' days '_behavAnalysis.mat']);
    %frames_run_buffer_mat = behav_output.frames_run_buffer_mat;
    speed = behav_struct.speed;
    frames_run_mat = behav_struct.frames_run_mat;
    % frames run mat is a matrix for 'run buffer', which includes right
    % before run, run, and right after run.
    
    %% calculate average speed of every 300ms
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
    mean_run_every300ms(bin) = mean(mean_columns(length(mean_columns)-2:length(mean_columns)));
    ste_run_every300ms(bin) = std(mean_columns(length(mean_columns)-2:length(mean_columns)))/sqrt(3);
    %% calculate average fluorescence of every 300ms
    dfOvF_staybase_mat = nan(size(frames_run_mat,1), size(frames_run_mat,2));
    dfOvF_run_cell = {};
    
    dfOvF_run_300ms = figure; p = panel(); p.pack(2,1);
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