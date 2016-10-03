clear all
close all
awFSAV_eye_naive100ms
rc = behavConstsAV;
av = mouseColors_naiveEye;
mice = unique({expt.SubNum});
nMice = length(mice);
str = unique({expt.SubNum});
values = cell2mat(cellfun(@str2num,str,'UniformOutput',false));
mouse_str = ['i' strjoin(str,'_i')];
mouse_ind = find(intersect(cell2mat({av.mouse}),values));
% mouse_str = [];
% for imouse = 1:nMice
%     mouse_str = [mouse_str 'i' num2str(av(mouse_ind(imouse)).mouse) '_'];  
% end
% try
%     fnout = fullfile(rc.eyeOutputDir, ['10-Nov-2015_' mouse_str]);
%     load([fnout 'EyeSummary.mat'])
% catch
        fnout = rc.eyeOutputDir;
        load(fullfile(fnout, ['\' date '_' mouse_str 'EyeSummary.mat']));
% end

set(0,'defaultfigurepaperorientation','landscape');
set(0,'defaultfigurepapersize',[11 8.5]);
set(0,'defaultfigurepaperposition',[.25 .25 [11 8.5]-0.5]);

x = -5:0.01:5;
y = x;

%% alignment summaries Aud vis Vis
for j = [1 3];
    ialign = j;
    out_mat = [2 3];
    align = mouse(2).expt(1).align(ialign).name;

%     if ~isempty(mouse(imouse).expt(1).align(ialign).av(1).outcome(out_mat(1)).avg_rad_pre) & ~isempty(mouse(imouse).expt(1).align(ialign).av(1).outcome(out_mat(2)).avg_rad_pre)

        %transient
        transientFig = figure;
        suptitle('trans')
        subplot(1,3,1)
        for imouse = 1:nMice
            a = nan(size(mouse(imouse).expt,2),1);
            b = nan(size(mouse(imouse).expt,2),1);
            for iexp = 1:size(mouse(imouse).expt,2)
                if ~isempty(mouse(imouse).expt(iexp).align(ialign).av(1).outcome(out_mat(1)).avg_rad_pre) & ~isempty(mouse(imouse).expt(iexp).align(ialign).av(1).outcome(out_mat(2)).avg_rad_pre)
                    col = av(mouse_ind(imouse)).col_str;
                    %transient
                    a(iexp) = mouse(imouse).expt(iexp).align(ialign).av(1).outcome(2).avg_rad_trans(:,1)./mouse(imouse).expt(iexp).align(1).av(1).outcome(2).avg_rad_pre(:,1);
                    b(iexp) = mouse(imouse).expt(iexp).align(ialign).av(2).outcome(2).avg_rad_trans(:,1)./mouse(imouse).expt(iexp).align(1).av(2).outcome(2).avg_rad_pre(:,1);
                    scatter(a, b, 200,['o' col]);        
                    hold on
                    xlabel('Visual (%)')
                    ylabel('Auditory (%)')
                    if iexp == 1
                        legendInfo([imouse+(imouse-1), imouse+imouse],1) = {[num2str(av(mouse_ind(imouse)).mouse) ' Transient']};
                    end
                end
            end
            errorbarxy(nanmean(a,1), nanmean(b,1), nanstd(a,[],1)/sqrt(size(a,1)), nanstd(b,[],1)/sqrt(size(b,1)), {['o' col], col, col});
            hold on
        end
        xlim([0.9 1.1])
        ylim([0.9 1.1])
        hold on
        plot(x,y, '--k')
        hold on
        hline(1,'-k')
        hold on
        vline(1,'-k')
        axis square
        %legend(legendInfo,'Location','Southeast')
        title([align ' aligned summary- transient- normalized radius values'])
%         print([fnout '_summary_' align 'align_trans_rad_norm.pdf'], '-dpdf');
%         savefig([fnout '_summary_' align 'align_trans_rad_norm.fig']);

        %sustained
        sustainedFig = figure;
        suptitle('sust')
        subplot(1,3,1)
        for imouse = 1:nMice
            a = nan(size(mouse(imouse).expt,2),1);
            b = nan(size(mouse(imouse).expt,2),1);
            for iexp = 1:size(mouse(imouse).expt,2)
                if ~isempty(mouse(imouse).expt(iexp).align(ialign).av(1).outcome(out_mat(1)).avg_rad_pre) & ~isempty(mouse(imouse).expt(iexp).align(ialign).av(1).outcome(out_mat(2)).avg_rad_pre)
                    col = av(mouse_ind(imouse)).col_str;
                    %transient
                    a(iexp) = mouse(imouse).expt(iexp).align(ialign).av(1).outcome(2).avg_rad_sust(:,1)./mouse(imouse).expt(iexp).align(1).av(1).outcome(2).avg_rad_pre(:,1);
                    b(iexp) = mouse(imouse).expt(iexp).align(ialign).av(2).outcome(2).avg_rad_sust(:,1)./mouse(imouse).expt(iexp).align(1).av(2).outcome(2).avg_rad_pre(:,1);
                    scatter(a, b, 200,['o' col]);        
                    hold on
                    xlabel('Visual (%)')
                    ylabel('Auditory (%)')
                    if iexp == 1
                        legendInfo([imouse+(imouse-1), imouse+imouse],1) = {[num2str(av(mouse_ind(imouse)).mouse) ' sustained']};
                    end
                end
            end
            errorbarxy(nanmean(a,1), nanmean(b,1), nanstd(a,[],1)/sqrt(size(a,1)), nanstd(b,[],1)/sqrt(size(b,1)), {['o' col], col, col});
            hold on
        end
        xlim([0.9 1.1])
        ylim([0.9 1.1])
        hold on
        plot(x,y, '--k')
        hold on
        hline(1,'-k')
        hold on
        vline(1,'-k')
        axis square
        %legend(legendInfo,'Location','Southeast')
        title([align ' aligned summary- sustained- normalized radius values'])
%         print([fnout '_summary_' align 'align_sust_rad_norm.pdf'], '-dpdf');
%         savefig([fnout '_summary_' align 'align_sust_rad_norm.fig']);

        %horizontal
        %transient
        figure(transientFig);
        subplot(1,3,2)
        for imouse = 1:nMice
            a = nan(size(mouse(imouse).expt,2),1);
            b = nan(size(mouse(imouse).expt,2),1);
            for iexp = 1:size(mouse(imouse).expt,2)
                if ~isempty(mouse(imouse).expt(iexp).align(ialign).av(1).outcome(out_mat(1)).avg_rad_pre) & ~isempty(mouse(imouse).expt(iexp).align(ialign).av(1).outcome(out_mat(2)).avg_rad_pre)
                    col = av(mouse_ind(imouse)).col_str;
                    %transient
                    a(iexp) = mouse(imouse).expt(iexp).align(ialign).av(1).outcome(2).avg_hor_trans(:,1)-mouse(imouse).expt(iexp).align(1).av(1).outcome(2).avg_hor_pre(:,1);
                    b(iexp) = mouse(imouse).expt(iexp).align(ialign).av(2).outcome(2).avg_hor_trans(:,1)-mouse(imouse).expt(iexp).align(1).av(2).outcome(2).avg_hor_pre(:,1);
                    scatter(a, b, 200,['o' col]);        
                    hold on
                    xlabel('Visual (deg)')
                    ylabel('Auditory (deg)')
                    if iexp == 1
                        legendInfo([imouse+(imouse-1), imouse+imouse],1) = {[num2str(av(mouse_ind(imouse)).mouse) ' Transient']};
                    end
                end
            end
            errorbarxy(nanmean(a,1), nanmean(b,1), nanstd(a,[],1)/sqrt(size(a,1)), nanstd(b,[],1)/sqrt(size(b,1)), {['o' col], col, col});
            hold on
        end
        xlim([-5 5])
        ylim([-5 5])
        hold on
        plot(x,y, '--k')
        hold on
        hline(0,'-k')
        hold on
        vline(0,'-k')
        axis square
        %legend(legendInfo,'Location','Southeast')
        title([align ' aligned summary- transient- absolute horizontal change from baseline'])
%         print([fnout '_summary_' align 'align_trans_horiz.pdf'], '-dpdf');
%         savefig([fnout '_summary_' align 'align_trans_horiz.fig']);

        %sustained
        figure(sustainedFig);
        subplot(1,3,2)
        for imouse = 1:nMice
            a = nan(size(mouse(imouse).expt,2),1);
            b = nan(size(mouse(imouse).expt,2),1);
            for iexp = 1:size(mouse(imouse).expt,2)
                if ~isempty(mouse(imouse).expt(iexp).align(ialign).av(1).outcome(out_mat(1)).avg_rad_pre) & ~isempty(mouse(imouse).expt(iexp).align(ialign).av(1).outcome(out_mat(2)).avg_rad_pre)
                    col = av(mouse_ind(imouse)).col_str;
                    %transient
                    a(iexp) = mouse(imouse).expt(iexp).align(ialign).av(1).outcome(2).avg_hor_sust(:,1)-mouse(imouse).expt(iexp).align(1).av(1).outcome(2).avg_hor_pre(:,1);
                    b(iexp) = mouse(imouse).expt(iexp).align(ialign).av(2).outcome(2).avg_hor_sust(:,1)-mouse(imouse).expt(iexp).align(1).av(2).outcome(2).avg_hor_pre(:,1);
                    scatter(a, b, 200,['o' col]);        
                    hold on
                    xlabel('Visual (deg)')
                    ylabel('Auditory (deg)')
                    if iexp == 1
                        legendInfo([imouse+(imouse-1), imouse+imouse],1) = {[num2str(av(mouse_ind(imouse)).mouse) ' sustained']};
                    end
                end
            end
            errorbarxy(nanmean(a,1), nanmean(b,1), nanstd(a,[],1)/sqrt(size(a,1)), nanstd(b,[],1)/sqrt(size(b,1)), {['o' col], col, col});
            hold on
        end
        xlim([-5 5])
        ylim([-5 5])
        hold on
        plot(x,y, '--k')
        hold on
        hline(0,'-k')
        hold on
        vline(0,'-k')
        axis square
        %legend(legendInfo,'Location','Southeast')
        title([align ' aligned summary- sustained- absolute horizontal change from baseline'])
%         print([fnout '_summary_' align 'align_sust_horiz.pdf'], '-dpdf');
%         savefig([fnout '_summary_' align 'align_sust_horiz.fig']);

        %vertical
        %transient
        figure(transientFig);
        subplot(1,3,3);
        for imouse = 1:nMice
            a = nan(size(mouse(imouse).expt,2),1);
            b = nan(size(mouse(imouse).expt,2),1);
            for iexp = 1:size(mouse(imouse).expt,2)
                if ~isempty(mouse(imouse).expt(iexp).align(ialign).av(1).outcome(out_mat(1)).avg_rad_pre) & ~isempty(mouse(imouse).expt(iexp).align(ialign).av(1).outcome(out_mat(2)).avg_rad_pre)
                    col = av(mouse_ind(imouse)).col_str;
                    %transient
                    a(iexp) = mouse(imouse).expt(iexp).align(ialign).av(1).outcome(2).avg_ver_trans(:,1)-mouse(imouse).expt(iexp).align(1).av(1).outcome(2).avg_ver_pre(:,1);
                    b(iexp) = mouse(imouse).expt(iexp).align(ialign).av(2).outcome(2).avg_ver_trans(:,1)-mouse(imouse).expt(iexp).align(1).av(2).outcome(2).avg_ver_pre(:,1);
                    scatter(a, b, 200,['o' col]);        
                    hold on
                    xlabel('Visual (deg)')
                    ylabel('Auditory (deg)')
                    if iexp == 1
                        legendInfo([imouse+(imouse-1), imouse+imouse],1) = {[num2str(av(mouse_ind(imouse)).mouse) ' Transient']};
                    end
                end
            end
            errorbarxy(nanmean(a,1), nanmean(b,1), nanstd(a,[],1)/sqrt(size(a,1)), nanstd(b,[],1)/sqrt(size(b,1)), {['o' col], col, col});
            hold on
        end
        xlim([-5 5])
        ylim([-5 5])
        hold on
        plot(x,y, '--k')
        hold on
        hline(0,'-k')
        hold on
        vline(0,'-k')
        axis square
        %legend(legendInfo,'Location','Southeast')
        title([align ' aligned summary- transient- absolute vertical change from baseline'])
        print([fnout '_summary_' align 'align_trans.pdf'], '-dpdf');
%         savefig([fnout '_summary_' align 'align_trans_vert.fig']);

        %sustained
        figure(sustainedFig);
        subplot(1,3,3)
        for imouse = 1:nMice
            a = nan(size(mouse(imouse).expt,2),1);
            b = nan(size(mouse(imouse).expt,2),1);
            for iexp = 1:size(mouse(imouse).expt,2)
                if ~isempty(mouse(imouse).expt(iexp).align(ialign).av(1).outcome(out_mat(1)).avg_rad_pre) & ~isempty(mouse(imouse).expt(iexp).align(ialign).av(1).outcome(out_mat(2)).avg_rad_pre)
                    col = av(mouse_ind(imouse)).col_str;
                    %transient
                    a(iexp) = mouse(imouse).expt(iexp).align(ialign).av(1).outcome(2).avg_ver_sust(:,1)-mouse(imouse).expt(iexp).align(1).av(1).outcome(2).avg_ver_pre(:,1);
                    b(iexp) = mouse(imouse).expt(iexp).align(ialign).av(2).outcome(2).avg_ver_sust(:,1)-mouse(imouse).expt(iexp).align(1).av(2).outcome(2).avg_ver_pre(:,1);
                    scatter(a, b, 200,['o' col]);        
                    hold on
                    xlabel('Visual (deg)')
                    ylabel('Auditory (deg)')
                    if iexp == 1
                        legendInfo([imouse+(imouse-1), imouse+imouse],1) = {[num2str(av(mouse_ind(imouse)).mouse) ' sustained']};
                    end
                end
            end
            errorbarxy(nanmean(a,1), nanmean(b,1), nanstd(a,[],1)/sqrt(size(a,1)), nanstd(b,[],1)/sqrt(size(b,1)), {['o' col], col, col});
            hold on
        end
        xlim([-5 5])
        ylim([-5 5])
        hold on
        plot(x,y, '--k')
        hold on
        hline(0,'-k')
        hold on
        vline(0,'-k')
        axis square
        %legend(legendInfo,'Location','Southeast')
        title([align ' aligned summary- sustained- absolute vertical change from baseline'])
        print([fnout '_summary_' align 'align_sust.pdf'], '-dpdf');
%         savefig([fnout '_summary_' align 'align_sust_vert.fig']);
%     end
end

% %% hit vs miss
% ialign = 1;
% out_mat = [2 3];
% align = mouse(2).align(ialign).name;
% 
% if ~isempty(mouse(imouse).align(ialign).av(1).outcome(out_mat(1)).avg_rad_pre) & ~isempty(mouse(imouse).align(ialign).av(1).outcome(out_mat(2)).avg_rad_pre)
% 
% figure;
% for imouse = 1:nMice
%     if ~isempty(mouse(imouse).align(ialign).av(1).outcome(out_mat(1)).avg_rad_pre) & ~isempty(mouse(imouse).align(ialign).av(1).outcome(out_mat(2)).avg_rad_pre)
%     col = av(mouse_ind(imouse)).col_str;
%     %transient
%     a = mouse(imouse).align(ialign).av(1).outcome(2).avg_rad_trans(:,1)./mouse(imouse).align(1).av(1).outcome(2).avg_rad_pre(:,1);
%     b = mouse(imouse).align(ialign).av(1).outcome(3).avg_rad_trans(:,1)./mouse(imouse).align(1).av(1).outcome(3).avg_rad_pre(:,1);
%     scatter(a, b, 200,['*' col]);        
%     hold on
%     hold on
%     % sustained
%     a = mouse(imouse).align(ialign).av(1).outcome(2).avg_rad_sust(:,1)./mouse(imouse).align(1).av(1).outcome(2).avg_rad_pre(:,1);
%     b = mouse(imouse).align(ialign).av(1).outcome(3).avg_rad_sust(:,1)./mouse(imouse).align(1).av(1).outcome(3).avg_rad_pre(:,1);
%     scatter(a, b, 200, ['o' col]);
%     hold on        
% 
%     xlabel('Success (%)')
%     ylabel('Miss (%)')
%     legendInfo([imouse+(imouse-1), imouse+imouse],1) = {[num2str(av(mouse_ind(imouse)).mouse) ' Transient'], [num2str(av(mouse_ind(imouse)).mouse) ' Sustained']};
%     end
% end
%     legend(legendInfo,'Location','Southeast')
% for imouse  = 1:nMice
%         if ~isempty(mouse(imouse).align(ialign).av(1).outcome(out_mat(1)).avg_rad_pre) & ~isempty(mouse(imouse).align(ialign).av(1).outcome(out_mat(2)).avg_rad_pre)
%         col = av(mouse_ind(imouse)).col_str;
%         a = mouse(imouse).align(ialign).av(1).outcome(2).avg_rad_trans(:,1)./mouse(imouse).align(1).av(1).outcome(2).avg_rad_pre(:,1);
%         b = mouse(imouse).align(ialign).av(1).outcome(3).avg_rad_trans(:,1)./mouse(imouse).align(1).av(1).outcome(3).avg_rad_pre(:,1);
%         errorbarxy(mean(a,1), mean(b,1), std(a,[],1)/sqrt(size(a,1)), std(b,[],1)/sqrt(size(b,1)), {['*' col], col, col});
%         hold on
%         a = mouse(imouse).align(ialign).av(1).outcome(2).avg_rad_sust(:,1)./mouse(imouse).align(1).av(1).outcome(2).avg_rad_pre(:,1);
%         b = mouse(imouse).align(ialign).av(1).outcome(3).avg_rad_sust(:,1)./mouse(imouse).align(1).av(1).outcome(3).avg_rad_pre(:,1);
%         errorbarxy(mean(a,1), mean(b,1), std(a,[],1)/sqrt(size(a,1)), std(b,[],1)/sqrt(size(b,1)), {['o' col], col, col});
%         hold on
%         xlim([0.9 1.1])
%         ylim([0.9 1.1])
%         hold on
%         plot(x,y, '--k')
%         hold on
%         hline(1,'-k')
%         hold on
%         vline(1,'-k')
%         axis square
%         hold on
%         end
% end    
%     title([align ' aligned summary- normalized radius values press_SM'])
%     print([fnout '_summary_' align 'align_rad_norm_SM.pdf'], '-dpdf');
% 
% figure;
% for imouse = 1:nMice
%     if ~isempty(mouse(imouse).align(ialign).av(1).outcome(out_mat(1)).avg_rad_pre) & ~isempty(mouse(imouse).align(ialign).av(1).outcome(out_mat(2)).avg_rad_pre)
%     col = av(mouse_ind(imouse)).col_str;
%     %transient
%     a = mouse(imouse).align(ialign).av(1).outcome(2).avg_hor_trans(:,1)-mouse(imouse).align(1).av(1).outcome(2).avg_hor_pre(:,1);
%     b = mouse(imouse).align(ialign).av(1).outcome(3).avg_hor_trans(:,1)-mouse(imouse).align(1).av(1).outcome(3).avg_hor_pre(:,1);
%     scatter(a, b,200, ['*' col]);
%     hold on;
%     
%     %sustained
%     a = mouse(imouse).align(ialign).av(1).outcome(2).avg_hor_sust(:,1)-mouse(imouse).align(1).av(1).outcome(2).avg_hor_pre(:,1);
%     b = mouse(imouse).align(ialign).av(1).outcome(3).avg_hor_sust(:,1)-mouse(imouse).align(1).av(1).outcome(3).avg_hor_pre(:,1);
%     scatter(a, b, 200, ['o' col]);
%     hold on;
%     
%     xlabel('Success (mm)')
%     ylabel('Miss (mm)')
%     legendInfo([imouse+(imouse-1), imouse+imouse],1) = {[num2str(av(mouse_ind(imouse)).mouse) ' Transient'], [num2str(av(mouse_ind(imouse)).mouse) ' Sustained']};
%     end
% end
% 
%     legend(legendInfo,'Location','Southeast')
% for imouse  = 1:nMice
%     if ~isempty(mouse(imouse).align(ialign).av(1).outcome(out_mat(1)).avg_rad_pre) & ~isempty(mouse(imouse).align(ialign).av(1).outcome(out_mat(2)).avg_rad_pre)
%         col = av(mouse_ind(imouse)).col_str;
%         %transient
%         a = mouse(imouse).align(ialign).av(1).outcome(2).avg_hor_trans(:,1)-mouse(imouse).align(1).av(1).outcome(2).avg_hor_pre(:,1);
%         b = mouse(imouse).align(ialign).av(1).outcome(3).avg_hor_trans(:,1)-mouse(imouse).align(1).av(1).outcome(3).avg_hor_pre(:,1);
%         errorbarxy(mean(a,1), mean(b,1), std(a,[],1)/sqrt(size(a,1)), std(b,[],1)/sqrt(size(b,1)), {['*' col], col, col});
%         hold on
%         %sustained
%         a = mouse(imouse).align(ialign).av(1).outcome(2).avg_hor_sust(:,1)-mouse(imouse).align(1).av(1).outcome(2).avg_hor_pre(:,1);
%         b = mouse(imouse).align(ialign).av(1).outcome(3).avg_hor_sust(:,1)-mouse(imouse).align(1).av(1).outcome(3).avg_hor_pre(:,1);
%         errorbarxy(mean(a,1), mean(b,1), std(a,[],1)/sqrt(size(a,1)), std(b,[],1)/sqrt(size(b,1)), {['o' col], col, col});
%         hold on
%         xlim([-.06 .06])
%         ylim([-.06 .06])
%         hold on
%         plot(x,y, '--k')
%         hold on
%         hline(0,'-k')
%         hold on
%         vline(0,'-k')
%         hold on
%         plot(x,y, '--k')
%         hold on
%         hline(1,'-k')
%         hold on
%         vline(1,'-k')
%         axis square
%         hold on
%     end
% end    
%     title([align ' aligned summary- absolute horizontal change from baseline_target press_SM'])
%     print([fnout align 'align_horiz_subt_SM.pdf'], '-dpdf');
% 
% figure;
% for imouse = 1:nMice
%     if ~isempty(mouse(imouse).align(ialign).av(1).outcome(out_mat(1)).avg_rad_pre) & ~isempty(mouse(imouse).align(ialign).av(1).outcome(out_mat(2)).avg_rad_pre)
%     col = av(mouse_ind(imouse)).col_str;
%     %transient
%     a = mouse(imouse).align(ialign).av(1).outcome(2).avg_ver_trans(:,1)-mouse(imouse).align(1).av(1).outcome(2).avg_ver_pre(:,1);
%     b = mouse(imouse).align(ialign).av(1).outcome(3).avg_ver_trans(:,1)-mouse(imouse).align(1).av(1).outcome(3).avg_ver_pre(:,1);
%     scatter(a, b, 200,['*' col]);
%     hold on;
%     
%     %sustained
%     a = mouse(imouse).align(ialign).av(1).outcome(2).avg_ver_sust(:,1)-mouse(imouse).align(1).av(1).outcome(2).avg_ver_pre(:,1);
%     b = mouse(imouse).align(ialign).av(1).outcome(3).avg_ver_sust(:,1)-mouse(imouse).align(1).av(1).outcome(3).avg_ver_pre(:,1);
%     scatter(a, b, 200, ['o' col]);
%     hold on;
%     
%     xlabel('Success (mm)')
%     ylabel('Miss (mm)')
%     legendInfo([imouse+(imouse-1), imouse+imouse],1) = {[num2str(av(mouse_ind(imouse)).mouse) ' Transient'], [num2str(av(mouse_ind(imouse)).mouse) ' Sustained']};
%     end
% end
% 
%     legend(legendInfo,'Location','Southeast')
% for imouse  = 1:nMice
%     if ~isempty(mouse(imouse).align(ialign).av(1).outcome(out_mat(1)).avg_rad_pre) & ~isempty(mouse(imouse).align(ialign).av(1).outcome(out_mat(2)).avg_rad_pre)
%         col = av(mouse_ind(imouse)).col_str;
%         %transient
%         a = mouse(imouse).align(ialign).av(1).outcome(2).avg_ver_trans(:,1)-mouse(imouse).align(1).av(1).outcome(2).avg_ver_pre(:,1);
%         b = mouse(imouse).align(ialign).av(1).outcome(3).avg_ver_trans(:,1)-mouse(imouse).align(1).av(1).outcome(3).avg_ver_pre(:,1);
%         errorbarxy(mean(a,1), mean(b,1), std(a,[],1)/sqrt(size(a,1)), std(b,[],1)/sqrt(size(b,1)), {['*' col], col, col});
%         hold on
%         %sustained
%         a = mouse(imouse).align(ialign).av(1).outcome(2).avg_ver_sust(:,1)-mouse(imouse).align(1).av(1).outcome(2).avg_ver_pre(:,1);
%         b = mouse(imouse).align(ialign).av(1).outcome(3).avg_ver_sust(:,1)-mouse(imouse).align(1).av(1).outcome(3).avg_ver_pre(:,1);
%         errorbarxy(mean(a,1), mean(b,1), std(a,[],1)/sqrt(size(a,1)), std(b,[],1)/sqrt(size(b,1)), {['o' col], col, col});
%         hold on
%         xlim([-.06 .06])
%         ylim([-.06 .06])
%         hold on
%         plot(x,y, '--k')
%         hold on
%         hline(0,'-k')
%         hold on
%         vline(0,'-k')
%         hold on
%         plot(x,y, '--k')
%         hold on
%         hline(1,'-k')
%         hold on
%         vline(1,'-k')
%         axis square
%         hold on
%     end
% end    
%     title([align ' aligned summary- absolute vertical change from baseline_target press_SM'])
%     print([fnout align 'align_vert_subt_SM.pdf'], '-dpdf');
% end
%     
