anal_base = '\\zoloto\bigstorlab\Lindsey\Analysis\120216';
%V1 mouse_list = {'AC45' 'M13' 'M14' 'M22' 'DR9'};
%areas = strvcat('PM', 'LM', 'AL');
%area_order = [2;3;1];
mouse_list = {'M10' 'Y28' 'M23'};
areas = strvcat('PM', 'AL', 'V1');
for iMouse  = [1 3]
    mouse = mouse_list{iMouse};
    if iMouse == 1
        date_mat = strvcat('110820', '110823', '110822');
        %expt = [2; 1; 3]; 
        expt = [2; 1; 4]; 
    end
    if iMouse == 2
        date_mat = strvcat('111001', '111002', '111004');
        expt = [3; 2; 5];
    end
    if iMouse == 3
        date_mat = strvcat('111123', '111126','111127');
        %expt = [5; 3; NaN];
        expt = [6; 4; NaN];
    end
    figure;
    for iArea = 1:3;
        area = areas(iArea,:);
        iexp = expt(iArea);

        subplot(2,3,iArea)
        n = all_fits(iArea).expt(iexp).n(1);
        for iCell = 1:n
            if all_fits(iArea).expt(iexp).bouton(iCell).goodfit ==1;
                scatter(log2(all_fits(iArea).expt(iexp).bouton(iCell).TF_fit), log2(all_fits(iArea).expt(iexp).bouton(iCell).SF_fit),3,'k');
                axis square
                hold on
            end
        end
        xlim([log2(1) log2(15)])
        title(area)

        plotfit = [];
        for iCell = 1:n
            if all_fits(iArea).expt(iexp).bouton(iCell).goodfit ==1;
                plotfit = cat(3, plotfit, all_fits(iArea).expt(iexp).bouton(iCell).plotfit);
            end
        end
        avg_tuning = mean(plotfit,3);
        subplot(2,3,iArea+3)
        imagesq(avg_tuning);
        colormap(gray);
        title(num2str(all_fits(iArea).expt(iexp).n(2)));
    end
    suptitle(mouse);
    fn_out = fullfile(anal_base, [matrix '_' num2str(P) 'P_' inj '_' mouse '_scatter_tuning_allareas_altexpt.pdf']);
    print(gcf, '-dpdf', fn_out);
end
 
% ind = [];
% n = all_fits(2).expt(1).n(1);
% for iCell = 1:n
%     if all_fits(2).expt(1).bouton(iCell).TF_fit>15
%         ind = [ind iCell];
%     end
% end
%     