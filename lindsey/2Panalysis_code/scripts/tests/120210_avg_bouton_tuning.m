areas = ['PM'; 'LM'; 'AL'];
col = strvcat('c', 'k', 'r');
P = 2;
matrix = 'SF5xTF5';
inj = 'V1';
anal_base = '\\zoloto\bigstorlab\Lindsey\Analysis\120202';
sum_base = 'G:\users\lindsey\analysisLG\experiments';
fn_summary = fullfile(sum_base, [matrix '_' num2str(P) 'P_' inj '_axon_summary'], 'all_fits.mat');
load(fn_summary);

for iArea = 1:3
    nexp = all_fits(iArea).nexp;
    figure;
    for iexp = 1:nexp
        plotfit = [];
        if all_fits(iArea).expt(iexp).n(2)>0
            n = all_fits(iArea).expt(iexp).n(1);
            for iCell = 1:n
                if all_fits(iArea).expt(iexp).bouton(iCell).goodfit == 1
                    plotfit = cat(3,plotfit, all_fits(iArea).expt(iexp).bouton(iCell).plotfit);
                end
            end
            subplot(4,4,iexp)
            imagesq(mean(plotfit,3));
            title([all_fits(iArea).expt(iexp).mouse ' ' num2str(all_fits(iArea).expt(iexp).n(2))])
            colormap(gray)
            caxis([0 .75])
        end
    end
    suptitle([areas(iArea,:) ' average tuning'])
end

    