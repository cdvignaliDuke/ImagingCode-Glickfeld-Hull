
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

edges_speed = log2(TF_vec./flipud(SF_vec));

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
end
figure;
for iArea = 1:3
    subplot(2,2,1)
    edges = log2(TF_vec);
    h= histc(log2(Goodfits(iArea).TF),edges);
    plot(edges(1:7,:), h(1:7,:)./sum(h), col(iArea,:));
    hold on
    xlabel('log2 TF')
    ylabel('fraction of boutons')
    subplot(2,2,2)
    edges = log2(SF_vec);
    h= histc(log2(Goodfits(iArea).SF),edges);
    plot(edges(1:7,:), h(1:7,:)./sum(h), col(iArea,:));
    hold on
    xlabel('log2 SF')
    ylabel('fraction of boutons')
    subplot(2,2,3)
    edges = log2(TF_vec./flipud(SF_vec));
    h= histc(log2(Goodfits(iArea).speed),edges);
    plot(edges(1:7,:), h(1:7,:)./sum(h), col(iArea,:));
    hold on
    xlabel('log2 speed')
    ylabel('fraction of boutons')
    subplot(2,2,4)
    shadedErrorBar(edges_speed(1:7,:), df_avg(:,iArea),df_sem(:,iArea), col(iArea,:));
    xlabel('log2 speed')
    ylabel('dF/F')
    ylim([0 .55])
    hold on
 end


fn_out = fullfile(anal_base, [matrix '_' num2str(P) 'P_' inj '_all_areas_log_speed_SF_TF_highlowpass.pdf']);
        print(gcf, '-dpdf', fn_out);
%difference of curves
df_diff = df_avg(:,1)-df_avg(:,3);
edges = log2(TF_vec./flipud(SF_vec));
dist_diff = histc(log2(Goodfits(1).speed),edges)./sum(histc(log2(Goodfits(1).speed),edges))- histc(log2(Goodfits(3).speed),edges)./sum(histc(log2(Goodfits(3).speed),edges));

figure;
for iArea = 1:3
    subplot(1,3,1)
    edges = log2(TF_vec./flipud(SF_vec));
    h= histc(log2(Goodfits(iArea).speed),edges);
    plot(edges(1:7,:), h(1:7,:)./sum(h), col(iArea,:));
    hold on
    xlabel('log2 speed')
    ylabel('fraction of boutons')
    axis square
    subplot(1,3,2)
    shadedErrorBar(edges_speed(1:7,:), df_avg(:,iArea),df_sem(:,iArea), col(iArea,:));
    xlabel('log2 speed')
    ylabel('dF/F')
    ylim([0 .55])
    axis square
    hold on
end
    subplot(1,3,3)
    plot(edges(1:7,:), df_diff./max(df_diff,[],1), '--k');
    hold on
    plot(edges(1:7,:), dist_diff(1:7,:)./max(dist_diff,[],1), '-k');
    axis square
    xlabel('log2 speed')
    ylabel('difference')
    ylim([-1.5 1.5])

    fn_out = fullfile(anal_base, [matrix '_' num2str(P) 'P_' inj '_AL_PM_speed_diff.pdf']);
        print(gcf, '-dpdf', fn_out);
%no high or lowpass SF/TF
for iArea = 1:3
    ind_TFbp = find(Goodfits(iArea).TF< 14.9999 & Goodfits(iArea).TF> 1.0001);
    ind_SFbp = find(Goodfits(iArea).SF< .3119 & Goodfits(iArea).SF> .0201);
    ind_bp = intersect(ind_TFbp, ind_SFbp);
    Goodfits(iArea).TFbp = Goodfits(iArea).TF(ind_bp);
    Goodfits(iArea).SFbp = Goodfits(iArea).SF(ind_bp);
    Goodfits(iArea).dFbp = Goodfits(iArea).dF(ind_bp);
    Goodfits(iArea).speedbp = Goodfits(iArea).speed(ind_bp);
    Goodfits(iArea).sigma_TFbp = Goodfits(iArea).sigma_TF(ind_bp);
    Goodfits(iArea).sigma_SFbp = Goodfits(iArea).sigma_SF(ind_bp);
end

df_BP_avg= zeros(7,3);
df_BP_sem= zeros(7,3);
speed_BP_avg = zeros(7,3);
speed_BP_sem = zeros(7,3);
figure;
edges = log2(TF_vec./flipud(SF_vec));
for iArea = 1:3
    [n bin] = histc(log2(Goodfits(iArea).speedbp), edges);
    for ibin = 1:7
        ind_ibin = find(bin==ibin);
        df_BP_avg(ibin,iArea) = mean(Goodfits(iArea).dF(ind_ibin));
        df_BP_sem(ibin,iArea) = std(Goodfits(iArea).dF(ind_ibin),[],1)./sqrt(length(ind_ibin));
        speed_BP_avg(ibin,iArea)=mean(Goodfits(iArea).speed(ind_ibin));
        speed_BP_sem(ibin,iArea) = std(Goodfits(iArea).speed(ind_ibin),[],1)./sqrt(length(ind_ibin));
    end
    subplot(2,2,1)
    edges = log2(TF_vec./flipud(SF_vec));
    h= histc(log2(Goodfits(iArea).speedbp),edges);
    plot(edges(1:7,:), h(1:7,:)./sum(h), col(iArea,:));
    hold on
    xlabel('log2 speed')
    ylabel('fraction of boutons')
    subplot(2,2,2)
    plot(edges_speed(1:7,:), df_BP_avg(:,iArea), col(iArea,:));
    xlabel('log2 speed')
    ylabel('dF/F')
    ylim([0 .55])
    hold on
end
fn_out = fullfile(anal_base, [matrix '_' num2str(P) 'P_' inj '_all_areas_speed_tuning_bandpass.pdf']);
        print(gcf, '-dpdf', fn_out);

 figure;
col = strvcat('c', 'k', 'r');
for iArea= 1:3;
    subplot(2,2,1);
    [H, stats, xCDF, yCDF] = cdfplot_LG(Goodfits(iArea).TFbp);
    semilogx(xCDF, yCDF, col(iArea,:));
    hold on
    xlim([1 15])
    xlabel('TF');
    subplot(2,2,2);
    [H, stats, xCDF, yCDF] = cdfplot_LG(Goodfits(iArea).SFbp);
    semilogx(xCDF, yCDF, col(iArea,:));
    hold on
    xlabel('SF');
    xlim([0.02 0.32])
    subplot(2,2,3);
    [H, stats, xCDF, yCDF] = cdfplot_LG(Goodfits(iArea).speedbp);
    semilogx(xCDF, yCDF, col(iArea,:));
    hold on
    xlim([3 750])
    xlabel('speed');
end

fn_out = fullfile(anal_base, [matrix '_' num2str(P) 'P_' inj '_allarea_log_cumdist_BPonly.pdf']);
         print(gcf, '-dpdf', fn_out);
 
figure;
col = strvcat('c', 'k', 'r');
for iArea= 1:3;
    subplot(2,2,1);
    [H, stats, xCDF, yCDF] = cdfplot_LG(Goodfits(iArea).sigma_TFbp);
    plot(xCDF, yCDF, col(iArea,:));
    hold on
    xlim([0 3])
    xlabel('sigma TF');
    subplot(2,2,2);
    [H, stats, xCDF, yCDF] = cdfplot_LG(Goodfits(iArea).sigma_SFbp);
    plot(xCDF, yCDF, col(iArea,:));
    hold on
    xlabel('sigma SF');
    xlim([0 3])
end
fn_out = fullfile(anal_base, [matrix '_' num2str(P) 'P_' inj '_allarea_log_cumdist_sigmas_BPonly.pdf']);
         print(gcf, '-dpdf', fn_out);
        
%group into dF

edges_dF = [0 .1  .2  .4 .8 1.6 3];
edges_speed = log2(TF_vec./flipud(SF_vec));
figure    
for iArea = 1:3
    [n bin] = histc(Goodfits(iArea).dF, edges_dF);
    for ibin = 1:length(edges_dF)-1
        ind_ibin = find(bin==ibin);
        h= histc(log2(Goodfits(iArea).speed(ind_ibin)),edges_speed);
        subplot(2,3, ibin);
        plot(edges(1:7,:), h(1:7,:)./sum(h), col(iArea,:));
        hold on
        xlabel('log2 speed')
        ylabel('fraction of boutons')
        ylim([0 .7])
        text(8, (.7-(iArea*.1)),num2str(n(ibin)))
        title(['dF/F = ' num2str(edges_dF(ibin)) ' to ' num2str(edges_dF(ibin+1))])
    end
end

fn_out = fullfile(anal_base, [matrix '_' num2str(P) 'P_' inj '_all_areas_speed_histograms_by_dF.pdf']);
        print(gcf, '-dpdf', fn_out);
        
edges_dF = [0 .1  .2  .4 .8 1.6 3];
edges_speed = log2(TF_vec./flipud(SF_vec));
figure    
for iArea = 1:3
    [n bin] = histc(Goodfits(iArea).dFbp, edges_dF);
    for ibin = 1:length(edges_dF)-1
        ind_ibin = find(bin==ibin);
        h= histc(log2(Goodfits(iArea).speedbp(ind_ibin)),edges_speed);
        subplot(2,3, ibin);
        plot(edges(1:7,:), h(1:7,:)./sum(h), col(iArea,:));
        hold on
        xlabel('log2 speed')
        ylabel('fraction of boutons')
        ylim([0 .7])
        text(8, (.7-(iArea*.1)),num2str(n(ibin)))
        title(['dF/F = ' num2str(edges_dF(ibin)) ' to ' num2str(edges_dF(ibin+1))])
    end
end

fn_out = fullfile(anal_base, [matrix '_' num2str(P) 'P_' inj '_all_areas_speed_histograms_by_dF_BPonly.pdf']);
print(gcf, '-dpdf', fn_out);


edges_dF = [.01:.01:3];
edges_speed = log2(TF_vec./flipud(SF_vec));
use_bins= [2 5];
figure    
for iArea = [3 1]
    [n bin] = histc(log2(Goodfits(iArea).speed), edges_speed);
    for ibin = 2
        ind_ibin = find(bin==ibin);
        h = histc(Goodfits(iArea).dF(ind_ibin), edges_dF);
        subplot(1,2, 1);
        plot(edges_dF, h./sum(h), col(iArea,:));
        axis square
        hold on
        xlabel('dF/F')
        ylabel('fraction of boutons')
        ylim([0 .1])
        title([num2str(2^(edges_speed(ibin))) ' deg/sec'])
    end
end

figure;
edges_dF = [0 .1  .2  .4 .8 1.6 3];
dF_inv = [7:-1:1];
edges_speed = log2(TF_vec./flipud(SF_vec));
dF_speed_density = zeros(7,7,3);
start = 1;
for iArea= 1:3
    [n bin] = histc(log2(Goodfits(iArea).speed), edges_speed);
    for ibin = 1:7
        ind_ibin = find(bin==ibin);
        [n_dF bin_dF] = histc(Goodfits(iArea).dF(ind_ibin), edges_dF);
        dF_speed_density(:,ibin,iArea) = flipud(n_dF(1:7,:));
    end
    subplot(2,3,start)
    imagesq(dF_speed_density(:,:,iArea));
    colormap(gray);
    title(areas(iArea,:));
    xlabel('speed')
    ylabel('dF/F')
    subplot(2,3,start+3)
    for ispeed = 1:7
        for idF = 1:7
            text(ispeed,dF_inv(idF), num2str(dF_speed_density(idF,ispeed,iArea)));
            hold on
        end
    end
    xlim([1 7])
    ylim([1 7])
    axis square
    axis off
    start = start+1;
end
fn_out = fullfile(anal_base, [matrix '_' num2str(P) 'P_' inj '_all_areas_speed_by_dF.pdf']);
print(gcf, '-dpdf', fn_out);

dF_speed_normspeed = zeros(7,7,3);
for iArea = 1:3
    for ispeed = 1:7
        dF_speed_normspeed(:,ispeed,iArea) = dF_speed_density(:,ispeed,iArea)./max(dF_speed_density(:,ispeed,iArea),[],1);
    end
end
dF_speed_normdf = zeros(7,7,3);
for iArea = 1:3
    for idf = 1:7
        dF_speed_normdf(idf,:,iArea) = dF_speed_density(idf,:,iArea)./max(dF_speed_density(idf,:,iArea),[],2);
    end
end

figure;
start = 1;
for iArea = 1:3
subplot(3,3,start)
    imagesq(dF_speed_density(:,:,iArea));
    colormap(gray);
    title(areas(iArea,:));
    xlabel('speed')
    ylabel('dF/F')
    subplot(3,3,start+3)
    imagesq(dF_speed_normdf(:,:,iArea));
    title('Norm dF');
    xlabel('speed')
    ylabel('dF/F')
    subplot(3,3,start+6)
    imagesq(dF_speed_normspeed(:,:,iArea));
    title('Norm speed');
    xlabel('speed')
    ylabel('dF/F')
    start= start+1;
end

fn_out = fullfile(anal_base, [matrix '_' num2str(P) 'P_' inj '_all_areas_speed_by_dF_Norm.pdf']);
print(gcf, '-dpdf', fn_out);
