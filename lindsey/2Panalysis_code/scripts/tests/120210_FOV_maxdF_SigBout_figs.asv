areas = ['PM'; 'LM'; 'AL'];
col = strvcat('c', 'k', 'r');
P = 2;
matrix = 'SF5xTF5';
inj = 'V1';
sum_base = 'G:\users\lindsey\analysisLG\experiments';
anal_base = '\\zoloto\bigstorlab\Lindsey\Analysis\120213';

fn_summary = fullfile(sum_base, [matrix '_' num2str(P) 'P_' inj '_axon_summary'], 'all_fits.mat');
load(fn_summary);


for iArea = 1:3
    nexp = all_fits(iArea).nexp;
    for iexp = 1:nexp
        mouse = all_fits(iArea).expt(iexp).mouse;
        date = all_fits(iArea).expt(iexp).date;
        userun = all_fits(iArea).expt(iexp).userun;

        base = 'G:\users\lindsey\analysisLG\active mice';
        outDir = fullfile(base, mouse,date);
        
        i = [];
        j = [];
        fn_local =  fullfile(outDir,'analysis',[date '_' mouse '_run' num2str(userun) '_local_max.mat']);
        load(fn_local);

        fn_stim = fullfile(outDir,'analysis',[date '_' mouse '_run' num2str(userun) '_allstim.mat']);
        load(fn_stim);
        FOV = mean(stim_off,3);

        fn_df =  fullfile(outDir,'analysis',[date '_' mouse '_run' num2str(userun) '_stack_dF_np.tif']);
        stack = readtiff(fn_df);
        dF_max = max(stack,[],3);

        figure;

        subplot(2,2,1)
        imagesq(dF_max);
        caxis([0 .5*MAX]);
        title('Maximum dF/F')
        colormap(gray);
        subplot(2,2,2)
        imagesq(dF_max);
        caxis([0 .5*MAX]);
        colormap(gray);
        hold on
        scatter(j,i,2.5,'b');
        title('Significantly responsive boutons')
        subplot(2,2,3)
        if all_fits(iArea).expt(iexp).n(2)>0
            n = all_fits(iArea).expt(iexp).n(1);
            for iCell = 1:n
                if all_fits(iArea).expt(iexp).bouton(iCell).goodfit ==1;
                    scatter(all_fits(iArea).expt(iexp).bouton(iCell).TF_fit, all_fits(iArea).expt(iexp).bouton(iCell).SF_fit,3,'k');
                    axis square
                    hold on
                end
            end
            ylim([0 .35])
            
            plotfit = [];
            for iCell = 1:n
                if all_fits(iArea).expt(iexp).bouton(iCell).goodfit ==1;
                    plotfit = cat(3, plotfit, all_fits(iArea).expt(iexp).bouton(iCell).plotfit);
                end
            end
            avg_tuning = mean(plotfit,3);
            subplot(2,2,4)
            imagesq(avg_tuning);
            colormap(gray);
        end
        
        suptitle([date ' ' mouse ' run' num2str(userun) ' ' areas(iArea,:)])
%         fn_out = fullfile(anal_base, [matrix '_' num2str(P) 'P_' inj '_' date '_' mouse '_run' num2str(userun) '_FOV_and_tuning.pdf']);
%         print(gcf, '-dpdf', fn_out);
    end
end