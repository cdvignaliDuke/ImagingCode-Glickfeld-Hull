fig_base = '\\zmey\storlab\users\Lindsey\Projects\HVAs\Manuscript\Figures_2012';
fig = 3;

areas = strvcat('PM', 'LM', 'AL');
col = strvcat('c', 'k', 'r');
P = 2;
matrix = 'SF5xTF5';
inj = 'V1';
sum_base = 'G:\users\lindsey\analysisLG\experiments';

fn_summary = fullfile(sum_base, [matrix '_' num2str(P) 'P_' inj '_axon_summary'], 'all_fits.mat');
load(fn_summary);
fn_good = fullfile(sum_base, [matrix '_' num2str(P) 'P_' inj '_axon_summary'], 'good_fits.mat');
load(fn_good);

figure;
for iArea= 1:3;
    subplot(2,3,1);
    [H, stats, xCDF, yCDF] = cdfplot_LG(log2(Goodfits(iArea).TF));
    plot(xCDF, yCDF, col(iArea,:));
    hold on
    title('TF')
    axis square
    xlabel('log2(TF)');
    xlim([log2(1) log2(15)])
    subplot(2,3,2);
    [H, stats, xCDF, yCDF] = cdfplot_LG(log2(Goodfits(iArea).SF));
    plot(xCDF, yCDF, col(iArea,:));
    hold on
    xlabel('log2(SF)');
    axis square
    title('SF')
    xlim([log2(.02) log2(.32)])
    subplot(2,3,3);
    [H, stats, xCDF, yCDF] = cdfplot_LG(log2(Goodfits(iArea).speed));
    plot(xCDF, yCDF, col(iArea,:));
    hold on
    xlabel('log2(speed)');
    axis square
    title('Speed')
    xlim([log2(1/.32) log2(15/.02)])
end

mouse_list = {'Y13' 'X32' 'DR7' 'DR9' 'AC39' 'AC42' 'AC44' 'AC45' 'Y18' 'Y26' 'M13' 'M14' 'M22' 'M31'};
inj = 'V1';
matrix = 'SF5xTF5';
P = 2;
area_order = [1 2; 2 3; 1 3];
for iMouse = 1:length(mouse_list);
    mouse = mouse_list{iMouse};
    for ipair = 1:3
        medians = [NaN NaN];
        for iArea = 1:2
            image = areas(area_order(ipair, iArea),:);
            sum_base = 'G:\users\lindsey\analysisLG\experiments';
            list_fn = fullfile(sum_base, [matrix '_' num2str(P) 'P_' image], [matrix '_' num2str(P) 'P_' inj '_' image '_exp_list.mat']);
            load(list_fn);
            median_speed = [];
            nexp = all_fits(area_order(ipair,iArea)).nexp;
            for iexp = 1:nexp
                exp_mouse = exp_list.mouse_mat{iexp};
                if length(exp_mouse)==length(mouse)
                    if exp_mouse == mouse;
                        if all_fits(area_order(ipair,iArea)).expt(iexp).n(2)>25
                            median_speed = [median_speed all_fits(area_order(ipair,iArea)).expt(iexp).median_speed];
                        end
                    end
                end
            end
            if length(median_speed)>0
                if length(median_speed)>1
                    median_speed = mean(median_speed);
                end
                medians(:,iArea) = median_speed;
            end
        end
        subplot(2,3,ipair+3)
        semilogy(1:2, medians, '-ob')
        axis square
        hold on
        title([areas(area_order(ipair,1),:) ' vs ' areas(area_order(ipair,2),:)])
        xlim([0 3])
        ylim([1 1000])
    end    
end
for ipair = 1:3
    for iArea = 1:2
        subplot(2,3,ipair+3)
        plot(iArea, median(Goodfits(area_order(ipair,iArea)).speed), 'sk');
    end
end

fn_out = fullfile(fig_base, ['Figure' num2str(fig)], [matrix '_' num2str(P) 'P_' inj '_all_areas_cdf_medians.ps']);
        print(gcf, '-depsc2', fn_out);
        
        
 %cell body data
 fn_cells = 'G:\users\lindsey\svn\users\lindsey\analysis\Mark cell data\FOR_LG_SF_TF_Speed_AL_PM_V1cellbodies.mat';
 load(fn_cells);
 
 figure;
subplot(1,3,1);
[H, stats, xCDF, yCDF] = cdfplot_LG(log2(SF_AL));
plot(xCDF, yCDF, '-r');
hold on
[H, stats, xCDF, yCDF] = cdfplot_LG(log2(SF_PM));
plot(xCDF, yCDF, '-c');
title('SF')
axis square
xlabel('log2(SF)');
xlim([log2(.02) log2(.32)])
subplot(1,3,2);
TF_AL_adj = TF_AL;
find(TF_AL_adj<1)=1;
find(TF_AL_adj>15)=15;
TF_PM_adj = TF_PM;
find(TF_PM_adj<1)=1;
find(TF_PM_adj>15)=15;
[H, stats, xCDF, yCDF] = cdfplot_LG(log2(TF_AL_adj));
plot(xCDF, yCDF, '-r');
hold on
[H, stats, xCDF, yCDF] = cdfplot_LG(log2(TF_PM_adj));
plot(xCDF, yCDF, '-c');
title('TF')
axis square
xlabel('log2(TF)');
xlim([log2(1) log2(15)])
subplot(1,3,3);
Speed_AL_adj = Speed_AL;
find(Speed_AL_adj<1)=1;
find(Speed_AL_adj>15)=15;
Speed_PM_adj = Speed_PM;
find(Speed_PM_adj<3.12)=3.12;
find(Speed_PM_adj>750)=750;
[H, stats, xCDF, yCDF] = cdfplot_LG(log2(Speed_AL_adj));
plot(xCDF, yCDF, '-r');
hold on
[H, stats, xCDF, yCDF] = cdfplot_LG(log2(Speed_PM_adj));
plot(xCDF, yCDF, '-c');
title('Speed')
axis square
xlabel('log2(Speed)');
xlim([log2(3.12) log2(750)])

fn_out = fullfile(fig_base, ['Figure' num2str(fig)], [matrix '_' num2str(P) 'P_' inj '_cell_body_cdfs.ps']);
        print(gcf, '-depsc2', fn_out);

figure;
subplot(2,2,1)
[n_bouton_PM bin_bouton_PM] = histc(log2(Goodfits(1).speed), edges);
plot(plot_edges, n_bouton_PM(1:7,:)./sum(n_bouton_PM), col(1,:))
hold on
[n_soma_PM bin_soma_PM] = histc(log2(Speed_PM_adj), edges);
plot(plot_edges, n_soma_PM(1:7,:)./sum(n_soma_PM), 'b')
axis square
ylim([0 .6])
subplot(2,2,2)
[n_bouton_AL bin_bouton_AL] = histc(log2(Goodfits(3).speed), edges);
plot(plot_edges, n_bouton_AL(1:7,:)./sum(n_bouton_AL), col(3,:))
hold on
[n_soma_AL bin_soma_AL] = histc(log2(Speed_AL_adj), edges);
plot(plot_edges, n_soma_AL(1:7,:)./sum(n_soma_AL), 'm')
axis square
ylim([0 .6])

fn_out = fullfile(fig_base, ['Figure' num2str(fig)], [matrix '_' num2str(P) 'P_' inj '_soma_bouton_histograms.ps']);
        print(gcf, '-depsc2', fn_out);

AL_soma_median = median(Speed_AL);
PM_soma_median = median(Speed_PM);
AL_bouton_median = median(Goodfits(3).speed);
PM_bouton_median = median(Goodfits(1).speed);
LM_bouton_median = median(Goodfits(2).speed);

figure;
semilogy(1, LM_bouton_median, 'sk');
hold on
semilogy(2, AL_bouton_median, 'sk')
hold on
semilogy(3, AL_soma_median, 'sk')
hold on
semilogy(4, PM_bouton_median, 'sk')
hold on
semilogy(5, PM_soma_median, 'sk')
hold on
xlim([0 6])
ylim([1 1000])
fn_out = fullfile(fig_base, ['Figure' num2str(fig)], [matrix '_' num2str(P) 'P_' inj '_medians.ps']);
        print(gcf, '-depsc2', fn_out);