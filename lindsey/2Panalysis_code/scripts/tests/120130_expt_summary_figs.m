P = 2;
matrix = 'SF5xTF5';
inj = 'V1';
sum_base = 'G:\users\lindsey\analysisLG\experiments';

fn_summary = fullfile(sum_base, [matrix '_' num2str(P) 'P_' inj '_axon_summary'], 'all_fits.mat');
load(fn_summary);

anal_base = '\\Zoloto\bigstorlab\Lindsey\Analysis\120131';
areas = ['PM'; 'LM'; 'AL'];
SF_vec0 = [.32 .16 .08 .04 .02]; %flipped to have low to high SF in square  %flipud
TF_vec0 = [1 2 4 8 15];
[tftf,sfsf] = meshgrid(TF_vec0,SF_vec0); 
speed_grid = tftf./sfsf;
speed_grid_long = reshape(speed_grid', 1, 25);
for iArea = 3;
    image = areas(iArea,:);
    nexp = all_fits(iArea).nexp;
    list_fn = fullfile(sum_base, [matrix '_' num2str(P) 'P_' image], [matrix '_' num2str(P) 'P_' inj '_' image '_exp_list.mat']);
    load(list_fn);
    for iexp = 1:nexp
        mouse = char(exp_list.mouse_mat{iexp});
        date = char(exp_list.date_mat{iexp});
        userun = exp_list.runs_mat{iexp};
        
        figure;
        start =1;
        Fig = 1;
        for iCell = 1:all_fits(iArea).expt(iexp).n(1)
            if start >64
                fn_out = fullfile(anal_base, areas(iArea,:), [date '_' mouse '_run' num2str(userun) '_data_vs fit_Fig' num2str(Fig) '.pdf']);
                print(gcf, '-dpdf', fn_out);
                figure;
                start = 1;
                Fig = Fig+1;
            end
            subplot(8,8,start);
            imagesq(all_fits(iArea).expt(iexp).bouton(iCell).dFoverF);
            title(num2str(max(max(all_fits(iArea).expt(iexp).bouton(iCell).dFoverF,[],1),[],2)));
            subplot(8,8,start+1);
            imagesq(all_fits(iArea).expt(iexp).bouton(iCell).plotfit);
            start = start+2;
            if all_fits(iArea).expt(iexp).bouton(iCell).goodfit == 1;
                title('**');
            end
            colormap(gray);
        end
        fn_out = fullfile(anal_base, areas(iArea,:), [date '_' mouse '_run' num2str(userun) '_data_vs fit_Fig' num2str(Fig) '.pdf']);
        print(gcf, '-dpdf', fn_out);

        if all_fits(iArea).expt(iexp).n(2)>0
        figure;
        col = hot(200); 
        for iCell = 1:all_fits(iArea).expt(iexp).n(1)
            if all_fits(iArea).expt(iexp).bouton(iCell).goodfit == 1
                subplot(2,2,1)
                scatter(all_fits(iArea).expt(iexp).bouton(iCell).TF_fit, all_fits(iArea).expt(iexp).bouton(iCell).SF_fit, 3, col(round(all_fits(iArea).expt(iexp).bouton(iCell).dF_fit*100)));
                hold on;
                xlabel('TF')
                ylabel('SF')
                subplot(2,2,2)
                scatter(all_fits(iArea).expt(iexp).bouton(iCell).speed, all_fits(iArea).expt(iexp).bouton(iCell).dF_fit, 3, 'k');
                xlabel('speed')
                ylabel('dF/F')
                hold on;
                subplot(2,2,3)
                scatter(all_fits(iArea).expt(iexp).bouton(iCell).pos(1,2), -all_fits(iArea).expt(iexp).bouton(iCell).pos(1,1), 3, col(round(all_fits(iArea).expt(iexp).bouton(iCell).dF_fit*100)));
                hold on;
                subplot(2,2,4)
                scatter(all_fits(iArea).expt(iexp).bouton(iCell).pos(1,2), -all_fits(iArea).expt(iexp).bouton(iCell).pos(1,1), 3, col(round(all_fits(iArea).expt(iexp).bouton(iCell).speed./4)));
                hold on;
            end
        end
        subplot(2,2,3)
        title('Color = dF/F')
        subplot(2,2,4)
        title('Color = speed')
        fn_out = fullfile(anal_base, areas(iArea,:), [date '_' mouse '_bouton_summary.pdf']);
        print(gcf, '-dpdf', fn_out);
        
        
        goodfits = zeros(all_fits(iArea).expt(iexp).n(2), 25);
        start = 1;
        for iCell = 1:all_fits(iArea).expt(iexp).n(1)
            if all_fits(iArea).expt(iexp).bouton(iCell).goodfit == 1
                goodfits(start,:) = reshape(all_fits(iArea).expt(iexp).bouton(iCell).plotfit', 1, 25);
                start = start+1;
            end
        end
        goodfits_avg = reshape(mean(goodfits,1), [5 5])';
        gooddata = zeros(all_fits(iArea).expt(iexp).n(2), 25);
        start = 1;
        for iCell = 1:all_fits(iArea).expt(iexp).n(1)
            if all_fits(iArea).expt(iexp).bouton(iCell).goodfit == 1
                gooddata(start,:) = reshape(all_fits(iArea).expt(iexp).bouton(iCell).dFoverF', 1, 25);
                start = start+1;
            end
        end
        gooddata_avg = reshape(mean(gooddata,1), [5 5])';
        uspeeds = unique(speed_grid);
        speed_tuning_fit = zeros(length(uspeeds),2);
        speed_tuning_data = zeros(length(uspeeds),2);
        for ispeed = 1:length(uspeeds)
            ind = find(speed_grid_long == uspeeds(ispeed));
            ispeed_goodfits = [];
            ispeed_gooddata = [];
            for iCond = ind
                ispeed_goodfits = [ispeed_goodfits; goodfits(:,iCond)];
                ispeed_gooddata = [ispeed_gooddata; gooddata(:,iCond)];
            end
            speed_tuning_fit(ispeed,1) = mean(ispeed_goodfits,1);
            speed_tuning_fit(ispeed,2) = std(ispeed_goodfits, [],1)./size(ispeed_goodfits,1);
            speed_tuning_data(ispeed,1) = mean(ispeed_gooddata,1);
            speed_tuning_data(ispeed,2) = std(ispeed_gooddata, [],1)./size(ispeed_gooddata,1);
        end
        figure;
        subplot(2,2,1)
        imagesq(goodfits_avg);
        title('fit')
        subplot(2,2,2)
        imagesq(gooddata_avg);
        title('data')
        subplot(2,2,3)
        errorbar(uspeeds, speed_tuning_fit(:,1), speed_tuning_fit(:,2));
        title('fit')
        xlim([0 800])
        xlabel('speed')
        ylabel('dF/F')
        subplot(2,2,4)
        errorbar(uspeeds, speed_tuning_data(:,1), speed_tuning_data(:,2));
        title('data')
        xlim([0 800])
        xlabel('speed')
        ylabel('dF/F')
        colormap(gray);
        fn_out = fullfile(anal_base, areas(iArea,:), [date '_' mouse '_speed_tuning.pdf']);
        print(gcf, '-dpdf', fn_out);
        end
    end
end
    
