clear all
close all

rc = behavConstsAV;
av = behavParamsAV;
mouse_str = [];
for imouse = 1:size(av,2)
    mouse_str = [mouse_str 'i' num2str(av(imouse).mouse) '_'];  
end
try
    fnout = fullfile(rc.eyeOutputDir, ['10-Nov-2015_' mouse_str]);
    load([fnout 'EyeSummary.mat'])
catch
    fnout = fullfile(rc.eyeOutputDir, [date '_' mouse_str]);
    load([fnout 'EyeSummary.mat'])
end

set(0,'defaultfigurepaperorientation','portrait');
set(0,'defaultfigurepapersize',[8.5 11]);
set(0,'defaultfigurepaperposition',[.25 .25 [8.5 11]-0.5]);

x = -2:0.01:2;
y = x;

%% alignment summaries Aud vis Vis - target
ialign = 3;
out_mat = [2 3];
align = mouse(2).align(ialign).name;

if ~isempty(mouse(imouse).align(ialign).av(1).outcome(out_mat(1)).avg_rad_pre) & ~isempty(mouse(imouse).align(ialign).av(1).outcome(out_mat(2)).avg_rad_pre)

figure;
for imouse = 1:size({av.mouse},2)
    if ~isempty(mouse(imouse).align(ialign).av(1).outcome(out_mat(1)).avg_rad_pre) & ~isempty(mouse(imouse).align(ialign).av(1).outcome(out_mat(2)).avg_rad_pre)
    col = av(imouse).col_str;
    %transient
    a = mouse(imouse).align(ialign).av(1).outcome(1).avg_rad_trans(:,1)./mouse(imouse).align(1).av(1).outcome(1).avg_rad_pre(:,1);
    b = mouse(imouse).align(ialign).av(2).outcome(1).avg_rad_trans(:,1)./mouse(imouse).align(1).av(2).outcome(1).avg_rad_pre(:,1);
    scatter(a, b, 200,['*' col]);        
    hold on
    hold on
    % sustained
    a = mouse(imouse).align(ialign).av(1).outcome(1).avg_rad_sust(:,1)./mouse(imouse).align(1).av(1).outcome(1).avg_rad_pre(:,1);
    b = mouse(imouse).align(ialign).av(2).outcome(1).avg_rad_sust(:,1)./mouse(imouse).align(1).av(2).outcome(1).avg_rad_pre(:,1);
    scatter(a, b, 200, ['o' col]);
    hold on        

    xlabel('Visual (%)')
    ylabel('Auditory (%)')
    legendInfo([imouse+(imouse-1), imouse+imouse],1) = {[num2str(av(imouse).mouse) ' Transient'], [num2str(av(imouse).mouse) ' Sustained']};
    end
end
    legend(legendInfo,'Location','Southeast')
for imouse  = 1:size({av.mouse},2)
        if ~isempty(mouse(imouse).align(ialign).av(1).outcome(out_mat(1)).avg_rad_pre) & ~isempty(mouse(imouse).align(ialign).av(1).outcome(out_mat(2)).avg_rad_pre)
        col = av(imouse).col_str;
        a = mouse(imouse).align(ialign).av(1).outcome(1).avg_rad_trans(:,1)./mouse(imouse).align(1).av(1).outcome(1).avg_rad_pre(:,1);
        b = mouse(imouse).align(ialign).av(2).outcome(1).avg_rad_trans(:,1)./mouse(imouse).align(1).av(2).outcome(1).avg_rad_pre(:,1);
        errorbarxy(mean(a,1), mean(b,1), std(a,[],1)/sqrt(size(a,1)), std(b,[],1)/sqrt(size(b,1)), {['*' col], col, col});
        hold on
        a = mouse(imouse).align(ialign).av(1).outcome(1).avg_rad_sust(:,1)./mouse(imouse).align(1).av(1).outcome(1).avg_rad_pre(:,1);
        b = mouse(imouse).align(ialign).av(2).outcome(1).avg_rad_sust(:,1)./mouse(imouse).align(1).av(2).outcome(1).avg_rad_pre(:,1);
        errorbarxy(mean(a,1), mean(b,1), std(a,[],1)/sqrt(size(a,1)), std(b,[],1)/sqrt(size(b,1)), {['o' col], col, col});
        hold on
        xlim([0.9 1.1])
        ylim([0.9 1.1])
        hold on
        plot(x,y, '--k')
        hold on
        hline(1,'-k')
        hold on
        vline(1,'-k')
        axis square
        hold on
        end
end    
    title([align ' aligned to target summary- normalized radius values'])
    print([fnout '_summary_' align 'align_rad_norm_target.pdf'], '-dpdf');

figure;
for imouse = 1:size({av.mouse},2)
    if ~isempty(mouse(imouse).align(ialign).av(1).outcome(out_mat(1)).avg_rad_pre) & ~isempty(mouse(imouse).align(ialign).av(1).outcome(out_mat(2)).avg_rad_pre)
    col = av(imouse).col_str;
    %transient
    a = mouse(imouse).align(ialign).av(1).outcome(1).avg_hor_trans(:,1)-mouse(imouse).align(1).av(1).outcome(1).avg_hor_pre(:,1);
    b = mouse(imouse).align(ialign).av(2).outcome(1).avg_hor_trans(:,1)-mouse(imouse).align(1).av(2).outcome(1).avg_hor_pre(:,1);
    scatter(a, b,200, ['*' col]);
    hold on;
    
    %sustained
    a = mouse(imouse).align(ialign).av(1).outcome(1).avg_hor_sust(:,1)-mouse(imouse).align(1).av(1).outcome(1).avg_hor_pre(:,1);
    b = mouse(imouse).align(ialign).av(2).outcome(1).avg_hor_sust(:,1)-mouse(imouse).align(1).av(2).outcome(1).avg_hor_pre(:,1);
    scatter(a, b, 200, ['o' col]);
    hold on;
    
    xlabel('Visual (mm)')
    ylabel('Auditory (mm)')
    legendInfo([imouse+(imouse-1), imouse+imouse],1) = {[num2str(av(imouse).mouse) ' Transient'], [num2str(av(imouse).mouse) ' Sustained']};
    end
end

    legend(legendInfo,'Location','Southeast')
for imouse  = 1:size({av.mouse},2)
    if ~isempty(mouse(imouse).align(ialign).av(1).outcome(out_mat(1)).avg_rad_pre) & ~isempty(mouse(imouse).align(ialign).av(1).outcome(out_mat(2)).avg_rad_pre)
        col = av(imouse).col_str;
        %transient
        a = mouse(imouse).align(ialign).av(1).outcome(1).avg_hor_trans(:,1)-mouse(imouse).align(1).av(1).outcome(1).avg_hor_pre(:,1);
        b = mouse(imouse).align(ialign).av(2).outcome(1).avg_hor_trans(:,1)-mouse(imouse).align(1).av(2).outcome(1).avg_hor_pre(:,1);
        errorbarxy(mean(a,1), mean(b,1), std(a,[],1)/sqrt(size(a,1)), std(b,[],1)/sqrt(size(b,1)), {['*' col], col, col});
        hold on
        %sustained
        a = mouse(imouse).align(ialign).av(1).outcome(1).avg_hor_sust(:,1)-mouse(imouse).align(1).av(1).outcome(1).avg_hor_pre(:,1);
        b = mouse(imouse).align(ialign).av(2).outcome(1).avg_hor_sust(:,1)-mouse(imouse).align(1).av(2).outcome(1).avg_hor_pre(:,1);
        errorbarxy(mean(a,1), mean(b,1), std(a,[],1)/sqrt(size(a,1)), std(b,[],1)/sqrt(size(b,1)), {['o' col], col, col});
        hold on
        xlim([-.06 .06])
        ylim([-.06 .06])
        hold on
        plot(x,y, '--k')
        hold on
        hline(0,'-k')
        hold on
        vline(0,'-k')
        hold on
        plot(x,y, '--k')
        hold on
        hline(1,'-k')
        hold on
        vline(1,'-k')
        axis square
        hold on
    end
end    
    title([align ' aligned summary- absolute horizontal change from baseline_target'])
    print([fnout align 'align_horiz_subt_target.pdf'], '-dpdf');

figure;
for imouse = 1:size({av.mouse},2)
    if ~isempty(mouse(imouse).align(ialign).av(1).outcome(out_mat(1)).avg_rad_pre) & ~isempty(mouse(imouse).align(ialign).av(1).outcome(out_mat(2)).avg_rad_pre)
    col = av(imouse).col_str;
    %transient
    a = mouse(imouse).align(ialign).av(1).outcome(1).avg_ver_trans(:,1)-mouse(imouse).align(1).av(1).outcome(1).avg_ver_pre(:,1);
    b = mouse(imouse).align(ialign).av(2).outcome(1).avg_ver_trans(:,1)-mouse(imouse).align(1).av(2).outcome(1).avg_ver_pre(:,1);
    scatter(a, b, 200,['*' col]);
    hold on;
    
    %sustained
    a = mouse(imouse).align(ialign).av(1).outcome(1).avg_ver_sust(:,1)-mouse(imouse).align(1).av(1).outcome(1).avg_ver_pre(:,1);
    b = mouse(imouse).align(ialign).av(2).outcome(1).avg_ver_sust(:,1)-mouse(imouse).align(1).av(2).outcome(1).avg_ver_pre(:,1);
    scatter(a, b, 200, ['o' col]);
    hold on;
    
    xlabel('Visual (mm)')
    ylabel('Auditory (mm)')
    legendInfo([imouse+(imouse-1), imouse+imouse],1) = {[num2str(av(imouse).mouse) ' Transient'], [num2str(av(imouse).mouse) ' Sustained']};
    end
end

    legend(legendInfo,'Location','Southeast')
for imouse  = 1:size({av.mouse},2)
    if ~isempty(mouse(imouse).align(ialign).av(1).outcome(out_mat(1)).avg_rad_pre) & ~isempty(mouse(imouse).align(ialign).av(1).outcome(out_mat(2)).avg_rad_pre)
        col = av(imouse).col_str;
        %transient
        a = mouse(imouse).align(ialign).av(1).outcome(1).avg_ver_trans(:,1)-mouse(imouse).align(1).av(1).outcome(1).avg_ver_pre(:,1);
        b = mouse(imouse).align(ialign).av(2).outcome(1).avg_ver_trans(:,1)-mouse(imouse).align(1).av(2).outcome(1).avg_ver_pre(:,1);
        errorbarxy(mean(a,1), mean(b,1), std(a,[],1)/sqrt(size(a,1)), std(b,[],1)/sqrt(size(b,1)), {['*' col], col, col});
        hold on
        %sustained
        a = mouse(imouse).align(ialign).av(1).outcome(1).avg_ver_sust(:,1)-mouse(imouse).align(1).av(1).outcome(1).avg_ver_pre(:,1);
        b = mouse(imouse).align(ialign).av(2).outcome(1).avg_ver_sust(:,1)-mouse(imouse).align(1).av(2).outcome(1).avg_ver_pre(:,1);
        errorbarxy(mean(a,1), mean(b,1), std(a,[],1)/sqrt(size(a,1)), std(b,[],1)/sqrt(size(b,1)), {['o' col], col, col});
        hold on
        xlim([-.06 .06])
        ylim([-.06 .06])
        hold on
        plot(x,y, '--k')
        hold on
        hline(0,'-k')
        hold on
        vline(0,'-k')
        hold on
        plot(x,y, '--k')
        hold on
        hline(1,'-k')
        hold on
        vline(1,'-k')
        axis square
        hold on
    end
end    
    title([align ' aligned summary- absolute vertical change from baseline_target'])
    print([fnout align 'align_vert_subt_target.pdf'], '-dpdf');
end

%% hit vs miss
ialign = 3;
out_mat = [2 3];
align = mouse(2).align(ialign).name;

if ~isempty(mouse(imouse).align(ialign).av(1).outcome(out_mat(1)).avg_rad_pre) & ~isempty(mouse(imouse).align(ialign).av(1).outcome(out_mat(2)).avg_rad_pre)

figure;
for imouse = 1:size({av.mouse},2)
    if ~isempty(mouse(imouse).align(ialign).av(1).outcome(out_mat(1)).avg_rad_pre) & ~isempty(mouse(imouse).align(ialign).av(1).outcome(out_mat(2)).avg_rad_pre)
    col = av(imouse).col_str;
    %transient
    a = mouse(imouse).align(ialign).av(1).outcome(2).avg_rad_trans(:,1)./mouse(imouse).align(1).av(1).outcome(2).avg_rad_pre(:,1);
    b = mouse(imouse).align(ialign).av(1).outcome(3).avg_rad_trans(:,1)./mouse(imouse).align(1).av(1).outcome(3).avg_rad_pre(:,1);
    scatter(a, b, 200,['*' col]);        
    hold on
    hold on
    % sustained
    a = mouse(imouse).align(ialign).av(1).outcome(2).avg_rad_sust(:,1)./mouse(imouse).align(1).av(1).outcome(2).avg_rad_pre(:,1);
    b = mouse(imouse).align(ialign).av(1).outcome(3).avg_rad_sust(:,1)./mouse(imouse).align(1).av(1).outcome(3).avg_rad_pre(:,1);
    scatter(a, b, 200, ['o' col]);
    hold on        

    xlabel('Success (%)')
    ylabel('Miss (%)')
    legendInfo([imouse+(imouse-1), imouse+imouse],1) = {[num2str(av(imouse).mouse) ' Transient'], [num2str(av(imouse).mouse) ' Sustained']};
    end
end
    legend(legendInfo,'Location','Southeast')
for imouse  = 1:size({av.mouse},2)
        if ~isempty(mouse(imouse).align(ialign).av(1).outcome(out_mat(1)).avg_rad_pre) & ~isempty(mouse(imouse).align(ialign).av(1).outcome(out_mat(2)).avg_rad_pre)
        col = av(imouse).col_str;
        a = mouse(imouse).align(ialign).av(1).outcome(2).avg_rad_trans(:,1)./mouse(imouse).align(1).av(1).outcome(2).avg_rad_pre(:,1);
        b = mouse(imouse).align(ialign).av(1).outcome(3).avg_rad_trans(:,1)./mouse(imouse).align(1).av(1).outcome(3).avg_rad_pre(:,1);
        errorbarxy(mean(a,1), mean(b,1), std(a,[],1)/sqrt(size(a,1)), std(b,[],1)/sqrt(size(b,1)), {['*' col], col, col});
        hold on
        a = mouse(imouse).align(ialign).av(1).outcome(2).avg_rad_sust(:,1)./mouse(imouse).align(1).av(1).outcome(2).avg_rad_pre(:,1);
        b = mouse(imouse).align(ialign).av(1).outcome(3).avg_rad_sust(:,1)./mouse(imouse).align(1).av(1).outcome(3).avg_rad_pre(:,1);
        errorbarxy(mean(a,1), mean(b,1), std(a,[],1)/sqrt(size(a,1)), std(b,[],1)/sqrt(size(b,1)), {['o' col], col, col});
        hold on
        xlim([0.9 1.1])
        ylim([0.9 1.1])
        hold on
        plot(x,y, '--k')
        hold on
        hline(1,'-k')
        hold on
        vline(1,'-k')
        axis square
        hold on
        end
end    
    title([align ' aligned to target summary- normalized radius values_SM'])
    print([fnout '_summary_' align 'align_rad_norm_target_SM.pdf'], '-dpdf');

figure;
for imouse = 1:size({av.mouse},2)
    if ~isempty(mouse(imouse).align(ialign).av(1).outcome(out_mat(1)).avg_rad_pre) & ~isempty(mouse(imouse).align(ialign).av(1).outcome(out_mat(2)).avg_rad_pre)
    col = av(imouse).col_str;
    %transient
    a = mouse(imouse).align(ialign).av(1).outcome(2).avg_hor_trans(:,1)-mouse(imouse).align(1).av(1).outcome(2).avg_hor_pre(:,1);
    b = mouse(imouse).align(ialign).av(1).outcome(3).avg_hor_trans(:,1)-mouse(imouse).align(1).av(1).outcome(3).avg_hor_pre(:,1);
    scatter(a, b,200, ['*' col]);
    hold on;
    
    %sustained
    a = mouse(imouse).align(ialign).av(1).outcome(2).avg_hor_sust(:,1)-mouse(imouse).align(1).av(1).outcome(2).avg_hor_pre(:,1);
    b = mouse(imouse).align(ialign).av(1).outcome(3).avg_hor_sust(:,1)-mouse(imouse).align(1).av(1).outcome(3).avg_hor_pre(:,1);
    scatter(a, b, 200, ['o' col]);
    hold on;
    
    xlabel('Success (mm)')
    ylabel('Miss (mm)')
    legendInfo([imouse+(imouse-1), imouse+imouse],1) = {[num2str(av(imouse).mouse) ' Transient'], [num2str(av(imouse).mouse) ' Sustained']};
    end
end

    legend(legendInfo,'Location','Southeast')
for imouse  = 1:size({av.mouse},2)
    if ~isempty(mouse(imouse).align(ialign).av(1).outcome(out_mat(1)).avg_rad_pre) & ~isempty(mouse(imouse).align(ialign).av(1).outcome(out_mat(2)).avg_rad_pre)
        col = av(imouse).col_str;
        %transient
        a = mouse(imouse).align(ialign).av(1).outcome(2).avg_hor_trans(:,1)-mouse(imouse).align(1).av(1).outcome(2).avg_hor_pre(:,1);
        b = mouse(imouse).align(ialign).av(1).outcome(3).avg_hor_trans(:,1)-mouse(imouse).align(1).av(1).outcome(3).avg_hor_pre(:,1);
        errorbarxy(mean(a,1), mean(b,1), std(a,[],1)/sqrt(size(a,1)), std(b,[],1)/sqrt(size(b,1)), {['*' col], col, col});
        hold on
        %sustained
        a = mouse(imouse).align(ialign).av(1).outcome(2).avg_hor_sust(:,1)-mouse(imouse).align(1).av(1).outcome(2).avg_hor_pre(:,1);
        b = mouse(imouse).align(ialign).av(1).outcome(3).avg_hor_sust(:,1)-mouse(imouse).align(1).av(1).outcome(3).avg_hor_pre(:,1);
        errorbarxy(mean(a,1), mean(b,1), std(a,[],1)/sqrt(size(a,1)), std(b,[],1)/sqrt(size(b,1)), {['o' col], col, col});
        hold on
        xlim([-.06 .06])
        ylim([-.06 .06])
        hold on
        plot(x,y, '--k')
        hold on
        hline(0,'-k')
        hold on
        vline(0,'-k')
        hold on
        plot(x,y, '--k')
        hold on
        hline(1,'-k')
        hold on
        vline(1,'-k')
        axis square
        hold on
    end
end    
    title([align ' aligned summary- absolute horizontal change from baseline_target_SM'])
    print([fnout align 'align_horiz_subt_target_SM.pdf'], '-dpdf');

figure;
for imouse = 1:size({av.mouse},2)
    if ~isempty(mouse(imouse).align(ialign).av(1).outcome(out_mat(1)).avg_rad_pre) & ~isempty(mouse(imouse).align(ialign).av(1).outcome(out_mat(2)).avg_rad_pre)
    col = av(imouse).col_str;
    %transient
    a = mouse(imouse).align(ialign).av(1).outcome(2).avg_ver_trans(:,1)-mouse(imouse).align(1).av(1).outcome(2).avg_ver_pre(:,1);
    b = mouse(imouse).align(ialign).av(1).outcome(3).avg_ver_trans(:,1)-mouse(imouse).align(1).av(1).outcome(3).avg_ver_pre(:,1);
    scatter(a, b, 200,['*' col]);
    hold on;
    
    %sustained
    a = mouse(imouse).align(ialign).av(1).outcome(2).avg_ver_sust(:,1)-mouse(imouse).align(1).av(1).outcome(2).avg_ver_pre(:,1);
    b = mouse(imouse).align(ialign).av(1).outcome(3).avg_ver_sust(:,1)-mouse(imouse).align(1).av(1).outcome(3).avg_ver_pre(:,1);
    scatter(a, b, 200, ['o' col]);
    hold on;
    
    xlabel('Success (mm)')
    ylabel('Miss (mm)')
    legendInfo([imouse+(imouse-1), imouse+imouse],1) = {[num2str(av(imouse).mouse) ' Transient'], [num2str(av(imouse).mouse) ' Sustained']};
    end
end

    legend(legendInfo,'Location','Southeast')
for imouse  = 1:size({av.mouse},2)
    if ~isempty(mouse(imouse).align(ialign).av(1).outcome(out_mat(1)).avg_rad_pre) & ~isempty(mouse(imouse).align(ialign).av(1).outcome(out_mat(2)).avg_rad_pre)
        col = av(imouse).col_str;
        %transient
        a = mouse(imouse).align(ialign).av(1).outcome(2).avg_ver_trans(:,1)-mouse(imouse).align(1).av(1).outcome(2).avg_ver_pre(:,1);
        b = mouse(imouse).align(ialign).av(1).outcome(3).avg_ver_trans(:,1)-mouse(imouse).align(1).av(1).outcome(3).avg_ver_pre(:,1);
        errorbarxy(mean(a,1), mean(b,1), std(a,[],1)/sqrt(size(a,1)), std(b,[],1)/sqrt(size(b,1)), {['*' col], col, col});
        hold on
        %sustained
        a = mouse(imouse).align(ialign).av(1).outcome(2).avg_ver_sust(:,1)-mouse(imouse).align(1).av(1).outcome(2).avg_ver_pre(:,1);
        b = mouse(imouse).align(ialign).av(1).outcome(3).avg_ver_sust(:,1)-mouse(imouse).align(1).av(1).outcome(3).avg_ver_pre(:,1);
        errorbarxy(mean(a,1), mean(b,1), std(a,[],1)/sqrt(size(a,1)), std(b,[],1)/sqrt(size(b,1)), {['o' col], col, col});
        hold on
        xlim([-.06 .06])
        ylim([-.06 .06])
        hold on
        plot(x,y, '--k')
        hold on
        hline(0,'-k')
        hold on
        vline(0,'-k')
        hold on
        plot(x,y, '--k')
        hold on
        hline(1,'-k')
        hold on
        vline(1,'-k')
        axis square
        hold on
    end
end    
    title([align ' aligned summary- absolute vertical change from baseline_target_SM'])
    print([fnout align 'align_vert_subt_target_SM.pdf'], '-dpdf');
end
    
