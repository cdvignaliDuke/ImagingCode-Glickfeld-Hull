rc = behavConstsAV;
av = behavParamsAV;
mouse_str = [];
for imouse = 1:size(av,2)
    mouse_str = [mouse_str 'i' num2str(av(imouse).mouse) '_'];  
end
fnout = fullfile(rc.eyeOutputDir, [date '_' mouse_str]);

load([fnout 'EyeSummary.mat'])

set(0,'defaultfigurepaperorientation','portrait');
set(0,'defaultfigurepapersize',[8.5 11]);
set(0,'defaultfigurepaperposition',[.25 .25 [8.5 11]-0.5]);
for imouse = 1:size({av.mouse},2)
    for iexp = 1:size(mouse(imouse).expt,2)
        for ialign = 1:3;
            if iexp == 1
                mouse(imouse).align(ialign).name = mouse(imouse).expt(iexp).align(ialign).name;
            end
            for iav = 1:2
                for iout = 1:5;
                    if iexp == 1
                        mouse(imouse).align(ialign).av(iav).outcome(iout).name = mouse(imouse).expt(iexp).align(ialign).av(iav).outcome(iout).name;
                    end
                    if length(mouse(imouse).expt(iexp).align(ialign).av(iav).outcome(iout).avg_rad_pre)>1
                        mouse(imouse).align(ialign).av(iav).outcome(iout).avg_rad_pre(iexp,:) = mouse(imouse).expt(iexp).align(ialign).av(iav).outcome(iout).avg_rad_pre;
                        mouse(imouse).align(ialign).av(iav).outcome(iout).avg_hor_pre(iexp,:) = mouse(imouse).expt(iexp).align(ialign).av(iav).outcome(iout).avg_hor_pre;
                        mouse(imouse).align(ialign).av(iav).outcome(iout).avg_ver_pre(iexp,:) = mouse(imouse).expt(iexp).align(ialign).av(iav).outcome(iout).avg_ver_pre;
                    end
                    if length(mouse(imouse).expt(iexp).align(ialign).av(iav).outcome(iout).avg_rad_trans)>1
                        mouse(imouse).align(ialign).av(iav).outcome(iout).avg_rad_trans(iexp,:) = mouse(imouse).expt(iexp).align(ialign).av(iav).outcome(iout).avg_rad_trans;
                        mouse(imouse).align(ialign).av(iav).outcome(iout).avg_rad_sust(iexp,:) = mouse(imouse).expt(iexp).align(ialign).av(iav).outcome(iout).avg_rad_sust;
                        mouse(imouse).align(ialign).av(iav).outcome(iout).avg_hor_trans(iexp,:) = mouse(imouse).expt(iexp).align(ialign).av(iav).outcome(iout).avg_hor_trans;
                        mouse(imouse).align(ialign).av(iav).outcome(iout).avg_hor_sust(iexp,:) = mouse(imouse).expt(iexp).align(ialign).av(iav).outcome(iout).avg_hor_sust;
                        mouse(imouse).align(ialign).av(iav).outcome(iout).avg_ver_trans(iexp,:) = mouse(imouse).expt(iexp).align(ialign).av(iav).outcome(iout).avg_ver_trans;
                        mouse(imouse).align(ialign).av(iav).outcome(iout).avg_ver_sust(iexp,:) = mouse(imouse).expt(iexp).align(ialign).av(iav).outcome(iout).avg_ver_sust;
                    end
                end
            end
        end
    end
end
save([fnout 'EyeSummary.mat'], 'mouse');
close all
x = -2:0.01:2;
y = x;

for ialign = 1:3
    align = mouse(2).align(ialign).name;
    if ialign == 1
        out_mat = [2 3];
    elseif ialign == 2
        out_mat = [2 4];
    elseif ialign == 3
        out_mat = [2 3];
    end
    out_str1 = mouse(imouse).align(ialign).av(iav).outcome(out_mat(1)).name;
    out_str2 = mouse(imouse).align(ialign).av(iav).outcome(out_mat(2)).name;
    
    figure;
    for imouse = 1:size({av.mouse},2)
        col = av(imouse).col_str;
        
        iout = 1;
        subplot(3,2,1)
        a = mouse(imouse).align(ialign).av(1).outcome(iout).avg_rad_trans(:,1)-mouse(imouse).align(1).av(1).outcome(iout).avg_rad_pre(:,1);
        b = mouse(imouse).align(ialign).av(2).outcome(iout).avg_rad_trans(:,1)-mouse(imouse).align(1).av(2).outcome(iout).avg_rad_pre(:,1);
        scatter(a, b, ['o' col]);
        title([' Transient radius change'])
        hold on;
        errorbarxy(mean(a,1), mean(b,1), std(a,[],1)/sqrt(size(a,1)), std(b,[],1)/sqrt(size(b,1)), {['o' col], col, col});
        subplot(3,2,2)
        a = mouse(imouse).align(ialign).av(1).outcome(iout).avg_rad_sust(:,1)-mouse(imouse).align(1).av(1).outcome(iout).avg_rad_pre(:,1);
        b = mouse(imouse).align(ialign).av(2).outcome(iout).avg_rad_sust(:,1)-mouse(imouse).align(1).av(2).outcome(iout).avg_rad_pre(:,1);
        scatter(a, b,  ['o' col]);
        title([' Sustained radius change'])
        hold on;
        errorbarxy(mean(a,1), mean(b,1), std(a,[],1)/sqrt(size(a,1)), std(b,[],1)/sqrt(size(b,1)), {['o' col], col, col});
        subplot(3,2,3)
        a = mouse(imouse).align(ialign).av(1).outcome(iout).avg_hor_trans(:,1)-mouse(imouse).align(1).av(1).outcome(iout).avg_hor_pre(:,1);
        b = mouse(imouse).align(ialign).av(2).outcome(iout).avg_hor_trans(:,1)-mouse(imouse).align(1).av(2).outcome(iout).avg_hor_pre(:,1);
        scatter(a, b, ['o' col]);
        title([' Transient horizontal change'])
        hold on;
        errorbarxy(mean(a,1), mean(b,1), std(a,[],1)/sqrt(size(a,1)), std(b,[],1)/sqrt(size(b,1)), {['o' col], col, col});
        subplot(3,2,4)
        a = mouse(imouse).align(ialign).av(1).outcome(iout).avg_hor_sust(:,1)-mouse(imouse).align(1).av(1).outcome(iout).avg_hor_pre(:,1);
        b = mouse(imouse).align(ialign).av(2).outcome(iout).avg_hor_sust(:,1)-mouse(imouse).align(1).av(2).outcome(iout).avg_hor_pre(:,1);
        scatter(a, b,  ['o' col]);
        title([' Sustained horizontal change'])
        hold on;
        errorbarxy(mean(a,1), mean(b,1), std(a,[],1)/sqrt(size(a,1)), std(b,[],1)/sqrt(size(b,1)), {['o' col], col, col});
        subplot(3,2,5)
        a = mouse(imouse).align(1).av(1).outcome(iout).avg_ver_pre(:,1)-mouse(imouse).align(ialign).av(1).outcome(iout).avg_ver_trans(:,1);
        b = mouse(imouse).align(1).av(2).outcome(iout).avg_ver_pre(:,1)-mouse(imouse).align(ialign).av(2).outcome(iout).avg_ver_trans(:,1);
        scatter(a, b, ['o' col]);
        title([' Transient vertical change'])
        hold on;
        errorbarxy(mean(a,1), mean(b,1), std(a,[],1)/sqrt(size(a,1)), std(b,[],1)/sqrt(size(b,1)), {['o' col], col, col});
        subplot(3,2,6)
        a = mouse(imouse).align(1).av(1).outcome(iout).avg_ver_pre(:,1)-mouse(imouse).align(ialign).av(1).outcome(iout).avg_ver_sust(:,1);
        b = mouse(imouse).align(1).av(2).outcome(iout).avg_ver_pre(:,1)-mouse(imouse).align(ialign).av(2).outcome(iout).avg_ver_sust(:,1);
        scatter(a, b,  ['o' col]);
        title(['Sustained vertical change'])
        hold on;
        errorbarxy(mean(a,1), mean(b,1), std(a,[],1)/sqrt(size(a,1)), std(b,[],1)/sqrt(size(b,1)), {['o' col], col, col});
        hold on
    for i = 1:6
        subplot(3,2,i)
        xlabel('Visual (mm)')
        ylabel('Auditory (mm)')
        xlim([-.1 .1])
        ylim([-.1 .1])
        hold on
        plot(x,y, '--k')
        hold on
        hline(0,'-k')
        hold on
        vline(0,'-k')
    end
    end
    suptitle([align ' aligned summary- absolute values'])
    print([fnout '_summary_' align 'align_AV_abs.pdf'], '-dpdf');

    figure;
    for imouse = 1:size({av.mouse},2)
        col = av(imouse).col_str;
        iout = 1;
        subplot(3,2,1)
        a = mouse(imouse).align(ialign).av(1).outcome(1).avg_rad_trans(:,1)./mouse(imouse).align(1).av(1).outcome(1).avg_rad_pre(:,1);
        b = mouse(imouse).align(ialign).av(2).outcome(1).avg_rad_trans(:,1)./mouse(imouse).align(1).av(2).outcome(1).avg_rad_pre(:,1);
        scatter(a, b, ['o' col]);
        title(['Transient all trials'])
        hold on;
        errorbarxy(mean(a,1), mean(b,1), std(a,[],1)/sqrt(size(a,1)), std(b,[],1)/sqrt(size(b,1)), {['o' col], col, col});
        subplot(3,2,2)
        a = mouse(imouse).align(ialign).av(1).outcome(1).avg_rad_sust(:,1)./mouse(imouse).align(1).av(1).outcome(1).avg_rad_pre(:,1);
        b = mouse(imouse).align(ialign).av(2).outcome(1).avg_rad_sust(:,1)./mouse(imouse).align(1).av(2).outcome(1).avg_rad_pre(:,1);
        scatter(a, b,  ['o' col]);
        title(['Sustained all trials'])
        hold on;
        errorbarxy(mean(a,1), mean(b,1), std(a,[],1)/sqrt(size(a,1)), std(b,[],1)/sqrt(size(b,1)), {['o' col], col, col});
        subplot(3,2,3)
        a = mouse(imouse).align(ialign).av(1).outcome(out_mat(:,1)).avg_rad_trans(:,1)./mouse(imouse).align(1).av(1).outcome(out_mat(:,1)).avg_rad_pre(:,1);
        b = mouse(imouse).align(ialign).av(1).outcome(out_mat(:,2)).avg_rad_trans(:,1)./mouse(imouse).align(1).av(1).outcome(out_mat(:,2)).avg_rad_pre(:,1);
        scatter(a, b, ['o' col]);
        title(['Transient visual trials'])
        hold on;
        errorbarxy(mean(a,1), mean(b,1), std(a,[],1)/sqrt(size(a,1)), std(b,[],1)/sqrt(size(b,1)), {['o' col], col, col});
        subplot(3,2,4)
        a = mouse(imouse).align(ialign).av(1).outcome(out_mat(:,1)).avg_rad_sust(:,1)./mouse(imouse).align(1).av(1).outcome(out_mat(:,1)).avg_rad_pre(:,1);
        b = mouse(imouse).align(ialign).av(1).outcome(out_mat(:,2)).avg_rad_sust(:,1)./mouse(imouse).align(1).av(1).outcome(out_mat(:,2)).avg_rad_pre(:,1);
        scatter(a, b,  ['o' col]);
        title(['ustained visual trials'])
        hold on;
        errorbarxy(mean(a,1), mean(b,1), std(a,[],1)/sqrt(size(a,1)), std(b,[],1)/sqrt(size(b,1)), {['o' col], col, col});
        subplot(3,2,5)
        a = mouse(imouse).align(ialign).av(2).outcome(out_mat(:,1)).avg_rad_trans(:,1)./mouse(imouse).align(1).av(2).outcome(out_mat(:,1)).avg_rad_pre(:,1);
        b = mouse(imouse).align(ialign).av(2).outcome(out_mat(:,2)).avg_rad_trans(:,1)./mouse(imouse).align(1).av(2).outcome(out_mat(:,2)).avg_rad_pre(:,1);
        scatter(a, b, ['o' col]);
        title(['Transient auditory trials'])
        hold on;
        errorbarxy(mean(a,1), mean(b,1), std(a,[],1)/sqrt(size(a,1)), std(b,[],1)/sqrt(size(b,1)), {['o' col], col, col});
        subplot(3,2,6)
        a = mouse(imouse).align(ialign).av(2).outcome(out_mat(:,1)).avg_rad_sust(:,1)./mouse(imouse).align(1).av(2).outcome(out_mat(:,1)).avg_rad_pre(:,1);
        b = mouse(imouse).align(ialign).av(2).outcome(out_mat(:,2)).avg_rad_sust(:,1)./mouse(imouse).align(1).av(2).outcome(out_mat(:,2)).avg_rad_pre(:,1);
        scatter(a, b,  ['o' col]);
        title(['Sustained auditory trials'])
        hold on;
        errorbarxy(mean(a,1), mean(b,1), std(a,[],1)/sqrt(size(a,1)), std(b,[],1)/sqrt(size(b,1)), {['o' col], col, col});
    
    for i = 1:6
        subplot(3,2,i)
        xlim([0.8 1.2])
        ylim([0.8 1.2])
        hold on
        plot(x,y, '--k')
        hold on
        hline(1,'-k')
        hold on
        vline(1,'-k')
        if i < 3
            xlabel('Visual (%)')
            ylabel('Auditory (%)')
        else
            xlabel([out_str1 ' (%)'])
            ylabel([out_str2 ' (%)'])
        end
    end
    end
    suptitle([align ' aligned summary- normalized radius values'])
    print([fnout '_summary_' align 'align_rad_norm.pdf'], '-dpdf');

    figure;
    for imouse = 1:size({av.mouse},2)
        col = av(imouse).col_str;
        iav = 1;
        subplot(3,2,1)
        a = mouse(imouse).align(ialign).av(iav).outcome(out_mat(1)).avg_rad_trans(:,1)-mouse(imouse).align(1).av(iav).outcome(out_mat(1)).avg_rad_pre(:,1);
        b = mouse(imouse).align(ialign).av(iav).outcome(out_mat(2)).avg_rad_trans(:,1)-mouse(imouse).align(1).av(iav).outcome(out_mat(2)).avg_rad_pre(:,1);
        scatter(a, b, ['o' col]);
        title([' Transient radius change'])
        hold on;
        errorbarxy(mean(a,1), mean(b,1), std(a,[],1)/sqrt(size(a,1)), std(b,[],1)/sqrt(size(b,1)), {['o' col], col, col});
        subplot(3,2,2)
        a = mouse(imouse).align(ialign).av(iav).outcome(out_mat(1)).avg_rad_sust(:,1)-mouse(imouse).align(1).av(iav).outcome(out_mat(1)).avg_rad_pre(:,1);
        b = mouse(imouse).align(ialign).av(iav).outcome(out_mat(2)).avg_rad_sust(:,1)-mouse(imouse).align(1).av(iav).outcome(out_mat(2)).avg_rad_pre(:,1);
        scatter(a, b,  ['o' col]);
        title(['Sustained radius change'])
        hold on;
        errorbarxy(mean(a,1), mean(b,1), std(a,[],1)/sqrt(size(a,1)), std(b,[],1)/sqrt(size(b,1)), {['o' col], col, col});
        subplot(3,2,3)
        a = mouse(imouse).align(ialign).av(iav).outcome(out_mat(1)).avg_hor_trans(:,1)-mouse(imouse).align(1).av(iav).outcome(out_mat(1)).avg_hor_pre(:,1);
        b = mouse(imouse).align(ialign).av(iav).outcome(out_mat(2)).avg_hor_trans(:,1)-mouse(imouse).align(1).av(iav).outcome(out_mat(2)).avg_hor_pre(:,1);
        scatter(a, b, ['o' col]);
        title([' Transient horizontal change'])
        hold on;
        errorbarxy(mean(a,1), mean(b,1), std(a,[],1)/sqrt(size(a,1)), std(b,[],1)/sqrt(size(b,1)), {['o' col], col, col});
        subplot(3,2,4)
        a = mouse(imouse).align(ialign).av(iav).outcome(out_mat(1)).avg_hor_sust(:,1)-mouse(imouse).align(1).av(iav).outcome(out_mat(1)).avg_hor_pre(:,1);
        b = mouse(imouse).align(ialign).av(iav).outcome(out_mat(2)).avg_hor_sust(:,1)-mouse(imouse).align(1).av(iav).outcome(out_mat(2)).avg_hor_pre(:,1);
        scatter(a, b,  ['o' col]);
        title(['Sustained horizontal change'])
        hold on;
        errorbarxy(mean(a,1), mean(b,1), std(a,[],1)/sqrt(size(a,1)), std(b,[],1)/sqrt(size(b,1)), {['o' col], col, col});
        subplot(3,2,5)
        a = mouse(imouse).align(1).av(iav).outcome(out_mat(1)).avg_ver_pre(:,1)-mouse(imouse).align(ialign).av(iav).outcome(out_mat(1)).avg_ver_trans(:,1);
        b = mouse(imouse).align(1).av(iav).outcome(out_mat(2)).avg_ver_pre(:,1)-mouse(imouse).align(ialign).av(iav).outcome(out_mat(2)).avg_ver_trans(:,1);
        scatter(a, b, ['o' col]);
        title([' Transient vertical change'])
        hold on;
        errorbarxy(mean(a,1), mean(b,1), std(a,[],1)/sqrt(size(a,1)), std(b,[],1)/sqrt(size(b,1)), {['o' col], col, col});
        subplot(3,2,6)
        a = mouse(imouse).align(1).av(iav).outcome(out_mat(1)).avg_ver_pre(:,1)-mouse(imouse).align(ialign).av(iav).outcome(out_mat(1)).avg_ver_sust(:,1);
        b = mouse(imouse).align(1).av(iav).outcome(out_mat(2)).avg_ver_pre(:,1)-mouse(imouse).align(ialign).av(iav).outcome(out_mat(2)).avg_ver_sust(:,1);
        scatter(a, b,  ['o' col]);
        title(['Sustained vertical change'])
        hold on;
        errorbarxy(mean(a,1), mean(b,1), std(a,[],1)/sqrt(size(a,1)), std(b,[],1)/sqrt(size(b,1)), {['o' col], col, col});
    
    for i = 1:6
        subplot(3,2,i)
        xlabel([out_str1 '(mm)'])
        ylabel([out_str2 '(mm)'])
        xlim([-.1 .1])
        ylim([-.1 .1])
        hold on
        plot(x,y, '--k')
        hold on
        hline(0,'-k')
        hold on
        vline(0,'-k')
    end
    end
    suptitle([align ' aligned visual summary- absolute values'])
    print([fnout '_summary_' align 'align_V_' out_str1 out_str2 '_abs.pdf'], '-dpdf');

    figure;
    for imouse = 1:size({av.mouse},2)
        col = av(imouse).col_str;
        iav = 2;
        subplot(3,2,1)
        a = mouse(imouse).align(ialign).av(iav).outcome(out_mat(1)).avg_rad_trans(:,1)-mouse(imouse).align(1).av(iav).outcome(out_mat(1)).avg_rad_pre(:,1);
        b = mouse(imouse).align(ialign).av(iav).outcome(out_mat(2)).avg_rad_trans(:,1)-mouse(imouse).align(1).av(iav).outcome(out_mat(2)).avg_rad_pre(:,1);
        scatter(a, b, ['o' col]);
        title([' Transient radius change'])
        hold on;
        errorbarxy(mean(a,1), mean(b,1), std(a,[],1)/sqrt(size(a,1)), std(b,[],1)/sqrt(size(b,1)), {['o' col], col, col});
        subplot(3,2,2)
        a = mouse(imouse).align(ialign).av(iav).outcome(out_mat(1)).avg_rad_sust(:,1)-mouse(imouse).align(1).av(iav).outcome(out_mat(1)).avg_rad_pre(:,1);
        b = mouse(imouse).align(ialign).av(iav).outcome(out_mat(2)).avg_rad_sust(:,1)-mouse(imouse).align(1).av(iav).outcome(out_mat(2)).avg_rad_pre(:,1);
        scatter(a, b,  ['o' col]);
        title(['Sustained radius change'])
        hold on;
        errorbarxy(mean(a,1), mean(b,1), std(a,[],1)/sqrt(size(a,1)), std(b,[],1)/sqrt(size(b,1)), {['o' col], col, col});
        subplot(3,2,3)
        a = mouse(imouse).align(ialign).av(iav).outcome(out_mat(1)).avg_hor_trans(:,1)-mouse(imouse).align(1).av(iav).outcome(out_mat(1)).avg_hor_pre(:,1);
        b = mouse(imouse).align(ialign).av(iav).outcome(out_mat(2)).avg_hor_trans(:,1)-mouse(imouse).align(1).av(iav).outcome(out_mat(2)).avg_hor_pre(:,1);
        scatter(a, b, ['o' col]);
        title([' Transient horizontal change'])
        hold on;
        errorbarxy(mean(a,1), mean(b,1), std(a,[],1)/sqrt(size(a,1)), std(b,[],1)/sqrt(size(b,1)), {['o' col], col, col});
        subplot(3,2,4)
        a = mouse(imouse).align(ialign).av(iav).outcome(out_mat(1)).avg_hor_sust(:,1)-mouse(imouse).align(1).av(iav).outcome(out_mat(1)).avg_hor_pre(:,1);
        b = mouse(imouse).align(ialign).av(iav).outcome(out_mat(2)).avg_hor_sust(:,1)-mouse(imouse).align(1).av(iav).outcome(out_mat(2)).avg_hor_pre(:,1);
        scatter(a, b,  ['o' col]);
        title(['Sustained horizontal change'])
        hold on;
        errorbarxy(mean(a,1), mean(b,1), std(a,[],1)/sqrt(size(a,1)), std(b,[],1)/sqrt(size(b,1)), {['o' col], col, col});
        subplot(3,2,5)
        a = mouse(imouse).align(1).av(iav).outcome(out_mat(1)).avg_ver_pre(:,1)-mouse(imouse).align(ialign).av(iav).outcome(out_mat(1)).avg_ver_trans(:,1);
        b = mouse(imouse).align(1).av(iav).outcome(out_mat(2)).avg_ver_pre(:,1)-mouse(imouse).align(ialign).av(iav).outcome(out_mat(2)).avg_ver_trans(:,1);
        scatter(a, b, ['o' col]);
        title([' Transient vertical change'])
        hold on;
        errorbarxy(mean(a,1), mean(b,1), std(a,[],1)/sqrt(size(a,1)), std(b,[],1)/sqrt(size(b,1)), {['o' col], col, col});
        subplot(3,2,6)
        a = mouse(imouse).align(1).av(iav).outcome(out_mat(1)).avg_ver_pre(:,1)-mouse(imouse).align(ialign).av(iav).outcome(out_mat(1)).avg_ver_sust(:,1);
        b = mouse(imouse).align(1).av(iav).outcome(out_mat(2)).avg_ver_pre(:,1)-mouse(imouse).align(ialign).av(iav).outcome(out_mat(2)).avg_ver_sust(:,1);
        scatter(a, b,  ['o' col]);
        title(['Sustained vertical change'])
        hold on;
        errorbarxy(mean(a,1), mean(b,1), std(a,[],1)/sqrt(size(a,1)), std(b,[],1)/sqrt(size(b,1)), {['o' col], col, col});
    
    for i = 1:6
        subplot(3,2,i)
        xlabel([out_str1 '(mm)'])
        ylabel([out_str2 '(mm)'])
        xlim([-.1 .1])
        ylim([-.1 .1])
        hold on
        plot(x,y, '--k')
        hold on
        hline(0,'-k')
        hold on
        vline(0,'-k')
    end
    end
    suptitle([align  'aligned auditory summary- absolute values'])
    print([fnout '_summary' align 'align_A_' out_str1 out_str2 '_abs.pdf'], '-dpdf');

%     figure;
%     for imouse = 1:size(av.mouse,2)
%         col = av(imouse).col_str;
%         subplot(3,3,1)
%         a = mouse(imouse).align(ialign).av(1).outcome(out_mat(1)).avg_rad_pre(:,1);
%         b = mouse(imouse).align(ialign).av(1).outcome(out_mat(2)).avg_rad_pre(:,1);
%         scatter(a, b, ['o' col]);
%         title(['Pre-' align ' radius: Visual'])
%         xlim([.4 .6])
%         ylim([0.4 0.6])
%         hold on;
%         errorbarxy(mean(a,1), mean(b,1), std(a,[],1)/sqrt(size(a,1)), std(b,[],1)/sqrt(size(b,1)), {['o' col], col, col});
%         subplot(3,3,2)
%         a = mouse(imouse).align(ialign).av(2).outcome(out_mat(1)).avg_rad_pre(:,1);
%         b = mouse(imouse).align(ialign).av(2).outcome(out_mat(2)).avg_rad_pre(:,1);
%         scatter(a, b,  ['o' col]);
%         title(['Pre-' align ' radius: Auditory'])    
%         xlim([0.4 0.6])
%         ylim([0.4 0.6])
%         hold on;
%         errorbarxy(mean(a,1), mean(b,1), std(a,[],1)/sqrt(size(a,1)), std(b,[],1)/sqrt(size(b,1)), {['o' col], col, col});
%         subplot(3,3,4)
%         a = mouse(imouse).align(ialign).av(1).outcome(out_mat(1)).avg_hor_pre(:,1);
%         b = mouse(imouse).align(ialign).av(1).outcome(out_mat(2)).avg_hor_pre(:,1);
%         scatter(a, b, ['o' col]);
%         title(['Pre-' align ' hor pos: Visual'])
%         xlim([1 2])
%         ylim([1 2])
%         hold on;
%         errorbarxy(mean(a,1), mean(b,1), std(a,[],1)/sqrt(size(a,1)), std(b,[],1)/sqrt(size(b,1)), {['o' col], col, col});
%         subplot(3,3,5)
%         a = mouse(imouse).align(ialign).av(2).outcome(out_mat(1)).avg_hor_pre(:,1);
%         b = mouse(imouse).align(ialign).av(2).outcome(out_mat(2)).avg_hor_pre(:,1);
%         scatter(a, b,  ['o' col]);
%         title(['Pre-' align ' hor pos: Auditory'])
%         xlim([1 2])
%         ylim([1 2])
%         hold on;
%         errorbarxy(mean(a,1), mean(b,1), std(a,[],1)/sqrt(size(a,1)), std(b,[],1)/sqrt(size(b,1)), {['o' col], col, col});
%         subplot(3,3,7)
%         a = mouse(imouse).align(ialign).av(1).outcome(out_mat(1)).avg_ver_pre(:,1);
%         b = mouse(imouse).align(ialign).av(1).outcome(out_mat(2)).avg_ver_pre(:,1);
%         scatter(a, b, ['o' col]);
%         title(['Pre-' align ' ver pos: Visual'])
%         xlim([1 2])
%         ylim([1 2])
%         hold on;
%         errorbarxy(mean(a,1), mean(b,1), std(a,[],1)/sqrt(size(a,1)), std(b,[],1)/sqrt(size(b,1)), {['o' col], col, col});
%         subplot(3,3,8)
%         a = mouse(imouse).align(ialign).av(2).outcome(out_mat(1)).avg_ver_pre(:,1);
%         b = mouse(imouse).align(ialign).av(2).outcome(out_mat(2)).avg_ver_pre(:,1);
%         scatter(a, b,  ['o' col]);
%         title(['Pre-' align ' ver pos: Auditory'])
%         xlim([1 2])
%         ylim([1 2])
%         hold on;
%         errorbarxy(mean(a,1), mean(b,1), std(a,[],1)/sqrt(size(a,1)), std(b,[],1)/sqrt(size(b,1)), {['o' col], col, col});
%         subplot(3,3,3)
%         a = mouse(imouse).align(ialign).av(1).outcome(1).avg_rad_pre(:,1);
%         b = mouse(imouse).align(ialign).av(2).outcome(1).avg_rad_pre(:,1);
%         scatter(a, b,  ['o' col]);
%         title(['Pre-' align ' radius'])
%         xlim([0.4 0.6])
%         ylim([0.4 0.6])
%         hold on;
%         errorbarxy(mean(a,1), mean(b,1), std(a,[],1)/sqrt(size(a,1)), std(b,[],1)/sqrt(size(b,1)), {['o' col], col, col});
%         subplot(3,3,6)
%         a = mouse(imouse).align(ialign).av(1).outcome(1).avg_hor_pre(:,1);
%         b = mouse(imouse).align(ialign).av(2).outcome(1).avg_hor_pre(:,1);
%         scatter(a, b, ['o' col]);
%         title(['Pre-' align ' hor pos'])
%         xlim([1 2])
%         ylim([1 2])
%         hold on;
%         errorbarxy(mean(a,1), mean(b,1), std(a,[],1)/sqrt(size(a,1)), std(b,[],1)/sqrt(size(b,1)), {['o' col], col, col});
%         subplot(3,3,9)
%         a = mouse(imouse).align(ialign).av(1).outcome(1).avg_ver_pre(:,1);
%         b = mouse(imouse).align(ialign).av(2).outcome(1).avg_ver_pre(:,1);
%         scatter(a, b,  ['o' col]);
%         title(['Pre-' align ' ver pos'])
%         xlim([1 2])
%         ylim([1 2])
%         hold on;
%         errorbarxy(mean(a,1), mean(b,1), std(a,[],1)/sqrt(size(a,1)), std(b,[],1)/sqrt(size(b,1)), {['o' col], col, col});
%     end
%     for i = [1 2 4 5 7 8]
%         subplot(3,3,i)
%         xlabel([out_str1 '(mm)'])
%         ylabel([out_str2 '(mm)'])
%         hold on
%         plot(x,y, '--k')
%     end
%     for i = [3 6 9]
%         subplot(3,3,i)
%         xlabel(['Visual (mm)'])
%         ylabel(['Auditory (mm)'])
%         hold on
%         plot(x,y, '--k')
%     end
%     suptitle(['Pre- ' align ' summary- absolute values'])
%     print([fnout '_summary_pre' align '_AV_' out_str1 out_str2 '_abs.pdf'], '-dpdf');
end