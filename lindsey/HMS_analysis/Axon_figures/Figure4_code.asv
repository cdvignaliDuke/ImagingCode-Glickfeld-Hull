fig_base = '\\zmey\storlab\users\Lindsey\Projects\HVAs\Manuscript\Figures_2012';
fig = 4;
col = strvcat('c', 'k', 'r');
P = 2;
matrix = 'SF5xTF5';
inj = 'V1';
sum_base = 'G:\users\lindsey\analysisLG\experiments';

fn_summary = fullfile(sum_base, [matrix '_' num2str(P) 'P_' inj '_axon_summary'], 'all_fits.mat');
load(fn_summary);
fn_good = fullfile(sum_base, [matrix '_' num2str(P) 'P_' inj '_axon_summary'], 'good_fits.mat');
load(fn_good);

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
edges = log2(TF_vec./flipud(SF_vec));
plot_edges = [edges(2) mean(edges(2:3)) mean(edges(3:4)) mean(edges(4:5)) mean(edges(5:6)) mean(edges(6:7)) edges(7)];
SF_plot_edges = log2([.02 .025 .045 .09 .18 .28 .32]);
TF_plot_edges = log2([1 1.25 2.25 4.5 9 13.5 15]);
%average dF by speed
df_avg= zeros(7,3,3);
df_sem= zeros(7,3,3);
figure;
for iArea = 1:3
    [n_SF bin_SF] = histc(log2(Goodfits(iArea).SF), log2(SF_vec));
    [n_TF bin_TF] = histc(log2(Goodfits(iArea).TF), log2(TF_vec));
    [n_speed bin_speed] = histc(log2(Goodfits(iArea).speed), edges);
    for ibin = 1:7
        ind_SF= find(bin_SF==ibin);
        ind_TF= find(bin_TF==ibin);
        ind_speed= find(bin_speed==ibin);
        df_avg(ibin,iArea,1) = mean(Goodfits(iArea).dF(ind_SF));
        df_sem(ibin,iArea,1) = std(Goodfits(iArea).dF(ind_SF),[],1)./sqrt(length(ind_SF));
        df_avg(ibin,iArea,2) = mean(Goodfits(iArea).dF(ind_TF));
        df_sem(ibin,iArea,2) = std(Goodfits(iArea).dF(ind_TF),[],1)./sqrt(length(ind_TF));
        df_avg(ibin,iArea,3) = mean(Goodfits(iArea).dF(ind_speed));
        df_sem(ibin,iArea,3) = std(Goodfits(iArea).dF(ind_speed),[],1)./sqrt(length(ind_speed));
    end
end
for iArea = 1:3
    subplot(2,2,1)
    errorbar(SF_plot_edges, df_avg(:,iArea,1), df_sem(:,iArea,1), col(iArea,:));
    hold on
    xlabel('log2(SF)');
    axis square
    ylabel('dF/F');
    xlim([log2(.018) log2(.35)])
    ylim([0 .7])
    subplot(2,2,2)
    errorbar(TF_plot_edges, df_avg(:,iArea,2), df_sem(:,iArea,2), col(iArea,:));
    hold on
    xlabel('log2(TF)');
    axis square
    ylabel('dF/F');
    xlim([log2(.9) log2(16)])
    ylim([0 .7])
    subplot(2,2,3)
    errorbar(plot_edges, df_avg(:,iArea,3), df_sem(:,iArea,3), col(iArea,:));
    hold on
    xlabel('log2(speed)');
    axis square
    ylabel('dF/F');
    xlim([log2(2) log2(1000)])
    ylim([0 .7])
end

fn_out = fullfile(fig_base, ['Figure' num2str(fig)], [matrix '_' num2str(P) 'P_' inj '_SF_TF_speed_vs_dF.ps']);
        print(gcf, '-depsc2', fn_out);
fn_out = fullfile(anal_base,  [matrix '_' num2str(P) 'P_' inj '_SF_TF_speed_vs_dF.pdf']);
        print(gcf, '-dpdf', fn_out);

%average speed by dF 
edges_dF = [0 .2 .4 .8 1.6 3];
plot_dF = [.1 .3 .6 1.2 2.3];
speed_avg= zeros(5,3);
speed_sem= zeros(5,3);
SF_avg= zeros(5,3);
SF_sem= zeros(5,3);
TF_avg= zeros(5,3);
TF_sem= zeros(5,3);

subplot(2,2,2)
figure;
for iArea = 1:3
    [n bin] = histc(Goodfits(iArea).dF, edges_dF);
    for ibin = 1:5
        ind = find(bin==ibin);
        speed_avg(ibin,iArea) = mean(log2(Goodfits(iArea).speed(ind)));
        speed_sem(ibin,iArea) = std(log2(Goodfits(iArea).speed(ind)),[],1)./sqrt(length(ind));
        SF_avg(ibin,iArea) = mean(log2(Goodfits(iArea).SF(ind)));
        SF_sem(ibin,iArea) = std(log2(Goodfits(iArea).SF(ind)),[],1)./sqrt(length(ind));
        TF_avg(ibin,iArea) = mean(log2(Goodfits(iArea).TF(ind)));
        TF_sem(ibin,iArea) = std(log2(Goodfits(iArea).TF(ind)),[],1)./sqrt(length(ind));
    end
    subplot(2,2,1)
    errorbar(log2(plot_dF), speed_avg(:,iArea), speed_sem(:,iArea), col(iArea,:));
    hold on
    ylim([-6.5 -.5]);
    xlim([-3.5 1.5]);
    xlabel('dF/F');
    axis square
    ylabel('log2(Speed)');
    subplot(2,2,2)
    errorbar(log2(plot_dF), SF_avg(:,iArea), SF_sem(:,iArea), col(iArea,:));
    hold on
    ylim([log2(0.02) log2(.32)]);
    xlim([-3.5 1.5]);
    xlabel('dF/F');
    axis square
    ylabel('log2(SF)');
    subplot(2,2,3)
    errorbar(log2(plot_dF), TF_avg(:,iArea), TF_sem(:,iArea), col(iArea,:));
    hold on
    ylim([log2(1) log2(15)]);
    xlim([-3.5 1.5]);
    xlabel('dF/F');
    axis square
    ylabel('log2(TF)');
end
fn_out = fullfile(anal_base,  [matrix '_' num2str(P) 'P_' inj '_SF_TF_speed_vs_binneddF.pdf']);
        print(gcf, '-dpdf', fn_out);
fn_out = fullfile(fig_base, ['Figure' num2str(fig)], [matrix '_' num2str(P) 'P_' inj '_SF_TF_speed_vs_binneddF.ps']);
        print(gcf, '-depsc2', fn_out);
% histograms by dF bin
figure
for iArea = 1:3
    [n_dF bin_dF] = histc(Goodfits(iArea).dF, edges_dF);
    for ibin = 1:length(n_dF)-1
        ind = find(bin_dF == ibin);
        [n_speed bin_speed] = histc(log2(Goodfits(iArea).speed(ind)), edges);
        [n_SF bin_SF] = histc(log2(Goodfits(iArea).SF(ind)), log2(SF_vec));
        [n_TF bin_TF] = histc(log2(Goodfits(iArea).TF(ind)), log2(TF_vec));
        subplot(3,5,ibin)
        plot(plot_edges, n_speed(1:7,:)./sum(n_speed,1), col(iArea,:));
        hold on
        axis square
        xlim([0 10]);
       ylim([0 .6]);
        xlabel('speed')
        subplot(3,5,ibin+5)
        plot(SF_plot_edges, n_SF(1:7,:)./sum(n_SF,1), col(iArea,:));
        hold on
        axis square
        xlim([log2(0.02) log2(.32)]);
        xlabel('SF')
        subplot(3,5,ibin+10)
        plot(TF_plot_edges, n_TF(1:7,:)./sum(n_TF,1), col(iArea,:));
        hold on
        axis square
        xlim([log2(1) log2(15)]);
        xlabel('TF')
    end   
end

fn_out = fullfile(fig_base, ['Figure' num2str(fig)], [matrix '_' num2str(P) 'P_' inj '_SF_TF_speed_vs_dF_histograms.ps']);
        print(gcf, '-depsc2', fn_out);
fn_out = fullfile(anal_base,  [matrix '_' num2str(P) 'P_' inj '_SF_TF_speed_vs_dF_histograms.pdf']);
        print(gcf, '-dpdf', fn_out);


figure;
for iArea = 1:3
[n_speed bin_speed] = histc(log2(Goodfits(iArea).speed), edges);
subplot(2,2,1)
plot(plot_edges, n_speed(1:7,:)./sum(n_speed,1), col(iArea,:))
hold on
end
axis square
fn_out = fullfile(fig_base, ['Figure' num2str(fig)], [matrix '_' num2str(P) 'P_' inj '_speed_histogram.ps']);
        print(gcf, '-depsc2', fn_out);

        
subplot(2,2,4)
for iArea = 1:3
    [n bin] = histc(Goodfits(iArea).dF, edges_dF);
    plot(log2(plot_dF), n(1:5,:)./(sum(n)), col(iArea,:));
    hold on
    xlabel('dF/F');
    xlim([-3.5 1.5])
    ylabel('Fraction of boutons')
    axis square
end


%classification discriminators
ErrMin_TF = zeros(1,5);
ErrMin_SF = zeros(1,5);
ErrMin_speed = zeros(1,5);
for ibin = 1:length(n_dF)-1
    ind_PM = find(bin_dF_PM == ibin);
    ind_AL = find(bin_dF_AL == ibin);

    SF_vec1 = log2(Goodfits(3).SF(ind_AL));
    SF_vec2 = log2(Goodfits(1).SF(ind_PM));
    TF_vec1 = log2(Goodfits(3).TF(ind_AL));
    TF_vec2 = log2(Goodfits(1).TF(ind_PM));
    speed_vec1 = log2(Goodfits(3).speed(ind_AL));
    speed_vec2 = log2(Goodfits(1).speed(ind_PM));

    SL = [SF_vec1; SF_vec2];
    SW = [TF_vec1; TF_vec2];
    SX = [speed_vec1; speed_vec2];

    group = [1*ones(length(SF_vec1),1); 2*ones(length(SF_vec2),1)];

    TF2 = SW;
    [i,j] = sort(TF2);
    val1 = (i(1:(end-1)) + i(2:end))/2;
    group1 = group(j);
    PercErr = zeros(size(val1));
    for count1 = 1:length(PercErr)
       val11 = val1(count1);
       ind1 = find(i<val11);
       N1_correct = length(find(i<val11 & group1==2));
       N1_false = length(find(i<val11 & group1==1));
       N2_correct = length(find(i>=val11 & group1==1));
       N2_false = length(find(i>=val11 & group1==2));
       PercErr(count1) = (N1_false + N2_false)./(length(i));
    end
    [i2,j2] = min(PercErr);
    ErrMin_TF(1,ibin) = 1-i2;
    ErrMin_TF_val = val1(j2);
    SF2 = SL;
    [i,j] = sort(SF2);
    val1 = (i(1:(end-1)) + i(2:end))/2;
    group1 = group(j);
    PercErr = zeros(size(val1));
    for count1 = 1:length(PercErr)
       val11 = val1(count1);
       ind1 = find(i<val11);
       N1_correct = length(find(i<val11 & group1==1));
       N1_false = length(find(i<val11 & group1==2));
       N2_correct = length(find(i>=val11 & group1==2));
       N2_false = length(find(i>=val11 & group1==1));
       PercErr(count1) = (N1_false + N2_false)./(length(i));
    end
    [i2,j2] = min(PercErr);
    ErrMin_SF(1,ibin) = 1-i2;
    ErrMin_SF_val = val1(j2);
    speed2 = SX;
    [i,j] = sort(speed2);
    val1 = (i(1:(end-1)) + i(2:end))/2;
    group1 = group(j);
    PercErr = zeros(size(val1));
    for count1 = 1:length(PercErr)
       val11 = val1(count1);
       ind1 = find(i<val11);
       N1_correct = length(find(i<val11 & group1==2));
       N1_false = length(find(i<val11 & group1==1));
       N2_correct = length(find(i>=val11 & group1==1));
       N2_false = length(find(i>=val11 & group1==2));
       PercErr(count1) = (N1_false + N2_false)./(length(i));
    end
    [i2,j2] = min(PercErr);
    ErrMin_speed(1,ibin) = 1-i2;
    ErrMin_speed_val = val1(j2);
end

figure;
subplot(2,2,1)
plot(log2(plot_dF), ErrMin_TF, '-k')
hold on
plot(log2(plot_dF), ErrMin_SF, '-b')
hold on
plot(log2(plot_dF), ErrMin_speed, '-r')
hold on
ylim([0.5 1])
xlim([-3.5 1.5]);

fn_out = fullfile(fig_base, ['Figure' num2str(fig)], [matrix '_' num2str(P) 'P_' inj '_SF_TF_speed_discriminators.ps']);
        print(gcf, '-depsc2', fn_out);

%median dF by speed
df_med= zeros(7,3);
df_lb= zeros(7,3);
df_ub= zeros(7,3);
ub_pct = .60;
lb_pct = .40;
for iArea = 1:3
    [n bin] = histc(log2(Goodfits(iArea).speed), edges);
    for ibin = 1:7
        ind = find(bin==ibin);
        ind_bounds = [ceil(length(ind)*lb_pct) ceil(length(ind)*ub_pct)];
        [order, index] = sort(Goodfits(iArea).dF(ind));
        lb_ind = index(ind_bounds(:,1));
        ub_ind = index(ind_bounds(:,2));
        df_med(ibin,iArea) = median(Goodfits(iArea).dF(ind));
        df_lb(ibin,iArea) = Goodfits(iArea).dF(ind(lb_ind));
        df_ub(ibin,iArea) = Goodfits(iArea).dF(ind(ub_ind));
    end
    subplot(2,2,2)
    ploterr(speed_avg(:,iArea), df_med(:,iArea),speed_sem(:,iArea), {df_lb(:,iArea), df_ub(:,iArea)}, col(iArea,:));
    hold on
    ylim([0 .6]);
    xlim([edges(2,:) edges(7,:)]);
    xlabel('log2(Speed)');
    ylabel('dF/F');
    axis square
end
title('Median') 

%90/10 dF/F
figure;
for iArea = 1:3
    [order, index] = sort(Goodfits(iArea).dF);
    percent_90 = ceil(length(index)*.9);
    [n_90 bin_90] = histc(log2(Goodfits(iArea).speed(index(1:percent_90,:))), edges);
    [n_10 bin_10] = histc(log2(Goodfits(iArea).speed(index(1+percent_90:end,:))), edges);
    for ibin = 1:7
        ind_90 = find(bin_90==ibin);
        ind_10 = find(bin_10==ibin);
        df_90_avg(ibin,iArea) = mean(Goodfits(iArea).dF(index(ind_90)));
        df_90_sem(ibin,iArea) = std(Goodfits(iArea).dF(index(ind_90)),[],1)./sqrt(length(ind_90));
        df_10_avg(ibin,iArea) = mean(Goodfits(iArea).dF(index(ind_10+percent_90)));
        df_10_sem(ibin,iArea) = std(Goodfits(iArea).dF(index(ind_10+percent_90)),[],1)./sqrt(length(ind_10));
    end
    subplot(2,2,1)
    errorbar(plot_edges, df_90_avg(:,iArea), df_90_sem(:,iArea), col(iArea,:));
    hold on
    ylim([0 .6]);
    xlim([log2(2) log2(1000)]);
    xlabel('log2(Speed)');
    axis square
    ylabel('dF/F');
    subplot(2,2,3)
    errorbar(plot_edges, df_10_avg(:,iArea), df_10_sem(:,iArea), col(iArea,:));
    hold on
    ylim([0 2]);
    xlim([log2(2) log2(1000)]);
    xlabel('log2(Speed)');
    axis square
    ylabel('dF/F');
end

%speed histograms by 90/10 dF/F
figure;
for iArea = 1:3
    [order, index] = sort(Goodfits(iArea).dF);
    percent_90 = ceil(length(index)*.9);
    [n_90_speed bin_90_speed] = histc(log2(Goodfits(iArea).speed(index(1:percent_90,:))), edges);
    [n_10_speed bin_10_speed] = histc(log2(Goodfits(iArea).speed(index(1+percent_90:end,:))), edges);
    subplot(2,2,1)
    plot(plot_edges, n_90_speed(1:7,:)./max(n_90_speed,[],1), col(iArea,:));
    hold on
    xlim([log2(2) log2(1000)]);
    ylim([0 1.2])
    axis square
    subplot(2,2,2)
    plot(plot_edges, n_10_speed(1:7,:)./max(n_10_speed,[],1), col(iArea,:));
    hold on
    xlim([log2(2) log2(1000)]);
    ylim([0 1.2])
    axis square
end
