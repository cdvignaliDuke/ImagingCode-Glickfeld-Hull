areas = ['PM'; 'LM'; 'AL'];
for iArea = 1:3;
    P = 2;
    matrix = 'SF5xTF5';
    image = areas(iArea,:);
    inj = 'V1';
    nON=12;
    nOFF=12;
    Nshuf = 500;
    SF_vec0 = [.32 .16 .08 .04 .02]; %flipped to have low to high SF in square  %flipud
    TF_vec0 = [1 2 4 8 15];

    sum_base = 'G:\users\lindsey\analysisLG\experiments';
    base = 'G:\users\lindsey\analysisLG\active mice';
    anal_base = '\\zoloto\bigstorlab\Lindsey\Analysis\120116';

    list_fn = fullfile(sum_base, [matrix '_' num2str(P) 'P_' image], [matrix '_' num2str(P) 'P_' inj '_' image '_exp_list.mat']);
    load(list_fn);
    nexp = size(exp_list.mouse_mat,2);
    mouse_list = [];
    mice = exp_list.mouse_mat;

    for iexp = 1:nexp;
        mouse = char(exp_list.mouse_mat{iexp});
        date = char(exp_list.date_mat{iexp});
        userun = exp_list.runs_mat{iexp};
        count_prot = exp_list.prot_mat{iexp};
        run = exp_list.run_mat{iexp};
        blanks = exp_list.blanks_mat{iexp};
        dir = exp_list.dir_mat{iexp};

        base = 'G:\users\lindsey\analysisLG\active mice';    
        outDir = fullfile(base, mouse,date);

        fn_local= fullfile(outDir, 'analysis', [date '_' mouse '_run' num2str(userun) '_local_max.mat']);
        load(fn_local);
            fn_fit = fullfile(outDir, 'analysis', [date '_' mouse '_run' num2str(userun) '_Fit_struct.mat']);   
        load(fn_fit);
        fn_lbub = fullfile(outDir, 'analysis', [date '_' mouse '_run' num2str(userun) '_lbub_fits.mat']);
        load(fn_lbub);

        figure;
        start = 1;
        fig = 1;
        for iCell = 1:n_pix    
            if start>64
                fn_out = fullfile(anal_base, [date '_' mouse '_run' num2str(userun) '_Data_vs_Fit_' num2str(fig) '.pdf']);
                print(gcf, '-dpdf', fn_out);
                figure;
                fig = 1+fig;
                start = 1;
            end
            subplot(8, 8, start);
            imagesq(Fit_struct(iCell).True.s_.data);
            subplot(8, 8, start+1);
            imagesq(Fit_struct(iCell).True.s_.k2b_plot);
            if lbub_diff(iCell,4)<2 
                if lbub_diff(iCell,5)<2
                    title('**');
                end
            end
            colormap(gray);
            start= start+2;
            if iCell == n_pix
                fn_out = fullfile(anal_base, [date '_' mouse '_run' num2str(userun) '_Data_vs_Fit_' num2str(fig) '.pdf']);
                print(gcf, '-dpdf', fn_out);
            end
        end

        figure;
        start = 1;
        fig = 1;
        for iCell = 1:n_pix    
            if start>64
                fn_out = fullfile(anal_base, [date '_' mouse '_run' num2str(userun) '_GoodFits_' num2str(fig) '.pdf']);
                print(gcf, '-dpdf', fn_out);
                figure;
                fig = 1+fig;
                start = 1;
            end
            if lbub_diff(iCell,4)<2 
                if lbub_diff(iCell,5)<2
                    subplot(8, 8, start);
                    imagesq(Fit_struct(iCell).True.s_.k2b_plot);
                    start = start+1;
                end
            end
            colormap(gray);
            if iCell == n_pix
                fn_out = fullfile(anal_base, [date '_' mouse '_run' num2str(userun) '_GoodFits_' num2str(fig) '.pdf']);
                print(gcf, '-dpdf', fn_out);
            end
        end
    end
end
                    
        figure;
        for iCell = 1:n_pix
            if lbub_diff(iCell,4)<2 
                if lbub_diff(iCell,5)<2
                    loglog(2.^lbub_fits(iCell,5,4), 2.^lbub_fits(iCell,4,4),'.k')
                    hold on
                end
            end
        end
        xlim([1 15])
        ylim([0.02 0.32])
        xlabel('Temporal Frequency')
        ylabel('Spatial Frequency')
        fn_out = fullfile(anal_base, [date '_' mouse '_run' num2str(userun) '_SF_vs_TF_scatter.pdf']);
                print(gcf, '-dpdf', fn_out);



        pass = find(lbub_diff(:,4)<2 & lbub_diff(:,5)<2);

        speed = (2.^lbub_fits(:,5,4))./(2.^lbub_fits(:,4,4));
        dF = squeeze(lbub_fits(:,1,4));
        dF_pass = dF(pass);
        speed_pass = speed(pass);
        speed_vec0 = TF_vec0./SF_vec0;

        speed_avg = zeros(5,2);
        dF_avg = zeros(5,2);
        n = zeros(5,1);
        for ispeed = 1:5
            if ispeed == 1
                ind = find(speed_pass<=speed_vec0(1,ispeed)+0.00001);
            else
                ind = find(speed_pass<=speed_vec0(1,ispeed)+0.00001 & speed_pass>speed_vec0(1,ispeed-1)+0.00001);
            end
            speed_avg(ispeed,1) = mean(speed_pass(ind,:),1);
            speed_avg(ispeed,2) = std(speed_pass(ind,:),1)./sqrt(size(ind,1));
            dF_avg(ispeed,1) = mean(dF_pass(ind,:),1);
            dF_avg(ispeed,2) = std(dF_pass(ind,:),1)./sqrt(size(ind,1));
            n(ispeed,1) = length(ind);
        end
        figure;
        ploterr(speed_avg(:,1), dF_avg(:,1), speed_avg(:,2), dF_avg(:,2), 'logx')
        hold on
        for ispeed = 1:5
            text(speed_avg(ispeed,1)-.1*speed_avg(ispeed,1), 0.1, num2str(n(ispeed, 1)));
        end
        ylim([0 2]);
        xlim([0 750]);
        xlabel('Speed');
        ylabel('dF/F');

        fn_out = fullfile(anal_base, [date '_' mouse '_run' num2str(userun) '_Speed_vs_dF_errorbar.pdf']);
        print(gcf, '-dpdf', fn_out);

        figure;
        for iCell = 1:n_pix
            if lbub_diff(iCell,4)<2 
                if lbub_diff(iCell,5)<2
                    semilogx(speed(iCell,:), lbub_fits(iCell,1,4),'.k')
                    hold on
                end
            end
        end
        xlim([0 1000])
        xlabel('Speed')
        ylabel('dF/F')

        fn_out = fullfile(anal_base, [date '_' mouse '_run' num2str(userun) '_Speed_vs_dF_scatter.pdf']);
        print(gcf, '-dpdf', fn_out);
    end
end
    
