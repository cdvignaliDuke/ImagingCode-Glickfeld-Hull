anal_base = '\\zoloto\bigstorlab\Lindsey\Analysis\120130';
fn_summary = fullfile(sum_base, [matrix '_' num2str(P) 'P_' inj '_axon_summary'], 'all_fits.mat');
load(fn_summary);

areas = ['PM'; 'LM'; 'AL'];
P = 2;
inj = 'V1';
matrix = 'SF5xTF5';
Goodfits = [];
for iArea = 1:3
    TF = [];
    SF = [];
    sigma_TF = [];
    sigma_SF = [];
    dF = [];
    speed = [];
    xi = [];
    plotfit = [];
    data = [];
    nexp = all_fits(iArea).nexp;
    image = areas(iArea,:);
    sum_base = 'G:\users\lindsey\analysisLG\experiments';
    list_fn = fullfile(sum_base, [matrix '_' num2str(P) 'P_' image], [matrix '_' num2str(P) 'P_' inj '_' image '_exp_list.mat']);
    load(list_fn);
    for iexp = 1:nexp
        mouse = char(exp_list.mouse_mat{iexp});
        date = char(exp_list.date_mat{iexp});
        userun = exp_list.runs_mat{iexp};
        base = 'G:\users\lindsey\analysisLG\active mice';    
        outDir = fullfile(base, mouse,date);
        n = all_fits(iArea).expt(iexp).n(1);
        fn_out = fullfile(outDir, 'analysis', [date '_' mouse '_run' num2str(userun) '_lbub_fits.mat']);
        load(fn_out);
        for iCell = 1:n
            if lbub_diff(iCell,4)<2;
                if lbub_diff(iCell,5)<2;
                    TF = [TF; all_fits(iArea).expt(iexp).bouton(iCell).TF_fit];
                    SF = [SF; all_fits(iArea).expt(iexp).bouton(iCell).SF_fit];
                    sigma_TF = [sigma_TF; all_fits(iArea).expt(iexp).bouton(iCell).sigma_TF];
                    sigma_SF = [sigma_SF; all_fits(iArea).expt(iexp).bouton(iCell).sigma_SF];
                    dF = [dF; all_fits(iArea).expt(iexp).bouton(iCell).dF_fit];
                    speed = [speed; all_fits(iArea).expt(iexp).bouton(iCell).speed];
                    xi = [xi; all_fits(iArea).expt(iexp).bouton(iCell).xi_fit];
                    plotfit = [plotfit; reshape(all_fits(iArea).expt(iexp).bouton(iCell).plotfit', 1, 25)];
                    data = [data; reshape(all_fits(iArea).expt(iexp).bouton(iCell).dFoverF', 1, 25)];
                end
            end
        end
    end
    Goodfits(iArea).TF = TF;
    Goodfits(iArea).SF = SF;
    Goodfits(iArea).sigma_TF = sigma_TF;
    Goodfits(iArea).sigma_SF = sigma_SF;
    Goodfits(iArea).dF = dF;
    Goodfits(iArea).speed = speed;
    Goodfits(iArea).xi = xi;
    Goodfits(iArea).plotfit = plotfit;
    Goodfits(iArea).data = data;
end
fn_good = fullfile(sum_base, [matrix '_' num2str(P) 'P_' inj '_axon_summary'], 'good_fits.mat');
save(fn_good, 'Goodfits');

        
figure;
col = strvcat('c', 'k', 'r');
for iArea= 1:3;
    subplot(2,2,1);
    [H, stats, xCDF, yCDF] = cdfplot_LG(Goodfits(iArea).TF);
    semilogx(xCDF, yCDF, col(iArea,:));
    hold on
    xlim([1 15])
    xlabel('TF');
    subplot(2,2,2);
    [H, stats, xCDF, yCDF] = cdfplot_LG(Goodfits(iArea).SF);
    semilogx(xCDF, yCDF, col(iArea,:));
    hold on
    xlabel('SF');
    xlim([0.02 0.32])
    subplot(2,2,3);
    [H, stats, xCDF, yCDF] = cdfplot_LG(Goodfits(iArea).dF);
    plot(xCDF, yCDF, col(iArea,:));
    hold on
    xlabel('dF');
    xlim([.01 3])
    subplot(2,2,4);
    [H, stats, xCDF, yCDF] = cdfplot_LG(Goodfits(iArea).speed);
    semilogx(xCDF, yCDF, col(iArea,:));
    hold on
    xlim([3 750])
    xlabel('speed');
end

fn_out = fullfile(anal_base, [matrix '_' num2str(P) 'P_' inj '_allarea_log_cumdist.pdf']);
         print(gcf, '-dpdf', fn_out);
 figure;
col = strvcat('c', 'k', 'r');
for iArea= 1:3;
    subplot(2,2,1);
    [H, stats, xCDF, yCDF] = cdfplot_LG(Goodfits(iArea).sigma_TF);
    plot(xCDF, yCDF, col(iArea,:));
    hold on
    xlim([0 3])
    xlabel('sigma TF');
    subplot(2,2,2);
    [H, stats, xCDF, yCDF] = cdfplot_LG(Goodfits(iArea).sigma_SF);
    plot(xCDF, yCDF, col(iArea,:));
    hold on
    xlabel('sigma SF');
    xlim([0 3])
    subplot(2,2,3);
    [H, stats, xCDF, yCDF] = cdfplot_LG(Goodfits(iArea).xi);
    plot(xCDF, yCDF, col(iArea,:));
    hold on
    xlabel('xi');
    xlim([-2 2])
end
fn_out = fullfile(anal_base, [matrix '_' num2str(P) 'P_' inj '_allarea_log_cumdistA.pdf']);
        print(gcf, '-dpdf', fn_out);
        
        
 vec0 = zeros(2,7);
 vec0(1,:) = [.5 1 2 4 8 16 32];
 vec0(2,:) = [.01 .02 .04 .08 .16 .32 .64];    
 vec1 = interp2(vec0');
 TF_vec = zeros(8,1);
 SF_vec = zeros(8,1);
 TF_vec(1,1)= vec1(2,1);
 TF_vec(2,1)= 1.0001;
 TF_vec(3:6,1) = vec1(4:2:10,1);
 TF_vec(7,1) = 14.999;
 TF_vec(8,1) = vec1(12,1);
 SF_vec(1,1)= vec1(2,3);
 SF_vec(2,1)= .020001;
 SF_vec(3:6,1) = vec1(4:2:10,3);
 SF_vec(7,1) = .3199;
 SF_vec(8,1) = vec1(12,3);
 for iArea = 1:3
 figure; 
 subplot(2,2,1)
 scatter(Goodfits(iArea).TF,Goodfits(iArea).SF, 2);
 hold on
 text(7.5, .35, num2str(length(Goodfits(iArea).TF)));
 subplot(2,2,2)
 edges = log2(TF_vec);
 h= histc(log2(Goodfits(iArea).TF),edges);
 bar(1:8, h);
 xlabel('TF')
 subplot(2,2,3)
 edges = log2(SF_vec);
 h= histc(log2(Goodfits(iArea).SF),edges);
 bar(1:8, h);
 xlabel('SF')
 subplot(2,2,4)
 edges = log2(TF_vec./fliplr(SF_vec));
 h= histc(log2(Goodfits(iArea).speed),edges);
 bar(1:8, h);
 xlabel('Speed')
 suptitle(['All Goodfits ' areas(iArea,:)])
 fn_out = fullfile(anal_base, [matrix '_' num2str(P) 'P_' inj areas(iArea,:), '_all_goodfits_log.pdf']);
        print(gcf, '-dpdf', fn_out);
 end
 
 figure;
 start = 1;
 for iArea = 1:3
 subplot(3,3,start)
  edges = log2(TF_vec);
 h= histc(log2(Goodfits(iArea).TF),edges);
 bar(1:8, h);
  xlim([0 8])
  subplot(3,3,start+1)
  edges = log2(SF_vec);
  h= histc(log2(Goodfits(iArea).SF),edges);
 bar(1:8, h);
  xlim([0 8])
  subplot(3,3,start+2)
  edges = log2(TF_vec./flipud(SF_vec));
  h= histc(log2(Goodfits(iArea).speed),edges);
 bar(1:8, h);
 xlim([0 8])
 start= start+3;
 end
subplot(3,3,1)
title('TF');
subplot(3,3,2)
title('SF');
subplot(3,3,3)
title('Speed');

fn_out = fullfile(anal_base, [matrix '_' num2str(P) 'P_' inj '_all_areas_log_speed_SF_TF_highlowpass.pdf']);
        print(gcf, '-dpdf', fn_out);


SF_vec_ud = flipud(SF_vec);
ind_sum = zeros(5,5,3);
start =1;
figure;
for iArea = 1:3
    for iTF = 1:5
        for iSF = 1:5
            ind_both = find(Goodfits(iArea).TF>TF_vec(iTF) & Goodfits(iArea).TF<TF_vec(iTF+1) & Goodfits(iArea).SF<SF_vec_ud(iSF) & Goodfits(iArea).SF>SF_vec_ud(iSF+1));
            ind_sum(iSF,iTF,iArea) = length(ind_both);
        end
    end
    subplot(3,2,start)
    scatter(Goodfits(iArea).TF, Goodfits(iArea).SF,2);
    axis square
    ylim([0 .32])
    xlim([0 15])
    subplot(3,2,start+1)
    imagesq(ind_sum(:,:,iArea));
    title(areas(iArea,:));
    colormap(gray)
    colorbar
    start= start+2;
end
fn_out = fullfile(anal_base, [matrix '_' num2str(P) 'P_' inj '_all_areas_bouton_density.pdf']);
        print(gcf, '-dpdf', fn_out);

edges = log2(TF_vec./flipud(SF_vec));
figure;
df_avg= zeros(7,3);
df_sem= zeros(7,3);
speed_avg = zeros(7,3);
speed_sem = zeros(7,3);
for iArea = 1:3
    [n bin] = histc(log2(Goodfits(iArea).speed), edges);
    for ibin = 1:7
        ind = find(bin==ibin);
        df_avg(ibin,iArea) = mean(Goodfits(iArea).dF(ind));
        df_sem(ibin,iArea) = std(Goodfits(iArea).dF(ind),[],1)./sqrt(length(ind));
        speed_avg(ibin,iArea)=mean(Goodfits(iArea).speed(ind));
        speed_sem(ibin,iArea) = std(Goodfits(iArea).speed(ind),[],1)./sqrt(length(ind));
    end
    subplot(3,1,iArea)
    ploterr(speed_avg(:,iArea), df_avg(:,iArea), speed_sem(:,iArea), df_sem(:,iArea),'logx');
    axis square
    xlim([0 1000])
    ylim([0 .5])
end

figure;
ploterr(speed_avg, df_avg, speed_sem, df_sem,'logx');
axis square
xlim([0 1000])
ylim([0 .5])
figure;
col = strvcat('c','k','r');
for iArea = 1:3
shadedErrorBar(edges(1:7,:), df_avg(:,iArea),df_sem(:,iArea), col(iArea,:),'logx');
hold on
end
ylim([0 .6]);
xlabel('log2(Speed)');
ylabel('dF');
fn_out = fullfile(anal_base, [matrix '_' num2str(P) 'P_' inj '_allareas_speed_vs_dF.pdf']);
        print(gcf, '-dpdf', fn_out);
