%mean vs meadian dF
df_avg= zeros(7,3);
df_sem= zeros(7,3);
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

%average dF by speed

figure;
for iArea = 1:3
    [n bin] = histc(log2(Goodfits(iArea).speed), edges);
    for ibin = 1:7
        ind = find(bin==ibin);
        df_avg(ibin,iArea) = mean(Goodfits(iArea).dF(ind));
        df_sem(ibin,iArea) = std(Goodfits(iArea).dF(ind),[],1)./sqrt(length(ind));
    end
    subplot(2,2,1)
    shadedErrorBar(edges(1:7,:), df_avg(:,iArea),df_sem(:,iArea), col(iArea,:));
    hold on
    ylim([0 1]);
    xlabel('log2(Speed)');
    axis square
    ylabel('dF/F');
end
title('Mean')

for iArea = 1:3
    [n bin] = histc(log2(Goodfits(iArea).speed), edges);
    for ibin = 1:7
        ind = find(bin==ibin);
        df_avg(ibin,iArea) = median(Goodfits(iArea).dF(ind));
        df_sem(ibin,iArea) = std(Goodfits(iArea).dF(ind),[],1)./sqrt(length(ind));
    end
    subplot(2,2,2)
    shadedErrorBar(edges(1:7,:), df_avg(:,iArea),df_sem(:,iArea), col(iArea,:));
    hold on
    ylim([0 1]);
    xlabel('log2(Speed)');
    ylabel('dF/F');
    axis square
end
title('Median')

% fn_out = fullfile(anal_base, [matrix '_' num2str(P) 'P_' inj '_all_areas_speed_vs_dF_mean_median.pdf']);
%         print(gcf, '-dpdf', fn_out);
        

%average speed by dF
speed_avg= zeros(6,3);
speed_sem= zeros(6,3);
dF_avg= zeros(6,3);
dF_sem= zeros(6,3);

edges_dF = [0 .1  .2  .4 .8 1.6 3];
dF_steps = [.05 .15 .3 .6 1.2 2.3];
subplot(2,2,3)
for iArea = 1:3
    [n bin] = histc(Goodfits(iArea).dF, edges_dF);
    for ibin = 1:6
        ind = find(bin==ibin);
        speed_avg(ibin,iArea) = mean(log2(Goodfits(iArea).speed(ind)));
        speed_sem(ibin,iArea) = std(log2(Goodfits(iArea).speed(ind)),[],1)./sqrt(length(ind));
        dF_avg(ibin,iArea) = mean(Goodfits(iArea).dF(ind));
        dF_sem(ibin,iArea) = std(Goodfits(iArea).dF(ind),[],1)./sqrt(length(ind));
    end
    ploterr(dF_avg(:,iArea), speed_avg(:,iArea),dF_sem(:,iArea), speed_sem(:,iArea), col(iArea,:));
    hold on
    ylim([0 10]);
    xlabel('dF/F');
    axis square
    ylabel('log2(Speed)');
end
fn_out = fullfile(anal_base, [matrix '_' num2str(P) 'P_' inj '_all_areas_speed_vs_dF.pdf']);
        print(gcf, '-dpdf', fn_out);

figure;
for iArea = 1:3
    subplot(2,2,iArea)
    scatter(Goodfits(iArea).dF, log2(Goodfits(iArea).speed))
end

figure;
exp_col = strvcat('r', 'm', 'c', 'b', 'y', 'g', 'k','r', 'm', 'c', 'b', 'y', 'g', 'k', 'r');
for iArea = 3
    nexp = all_fits(iArea).nexp;
    subplot(2,2,iArea)
    for iexp = 1:nexp
        if iexp<6 | iexp>7
            if all_fits(iArea).expt(iexp).n(2)>0
                n = all_fits(iArea).expt(iexp).n(1);
                for iCell = 1:n
                    if all_fits(iArea).expt(iexp).bouton(iCell).goodfit == 1
                        scatter(all_fits(iArea).expt(iexp).bouton(iCell).dF_fit, log2(all_fits(iArea).expt(iexp).bouton(iCell).speed), exp_col(iexp,:));
                        hold on
                    end
                end
            end
        end
    end
end

AL_dF = [];
AL_speed = [];
for iArea = 3
    nexp = all_fits(iArea).nexp;
    subplot(2,2,iArea)
    for iexp = 1:nexp
        if iexp<6 | iexp>7
            if all_fits(iArea).expt(iexp).n(2)>0
                n = all_fits(iArea).expt(iexp).n(1);
                for iCell = 1:n
                    if all_fits(iArea).expt(iexp).bouton(iCell).goodfit == 1
                        AL_dF = [AL_dF all_fits(iArea).expt(iexp).bouton(iCell).dF_fit];
                        AL_speed = [AL_speed all_fits(iArea).expt(iexp).bouton(iCell).speed];
                    end
                end
            end
        end
    end
end


figure;

AL_df_avg = zeros(7,1);
AL_df_sem = zeros(7,1);
[n bin] = histc(log2(AL_speed), edges);
for ibin = 1:7
    ind = find(bin==ibin);
    AL_df_avg(ibin,:) = mean(AL_dF(1,ind));
    AL_df_sem(ibin,:) = std(AL_dF(1,ind),[],2)./sqrt(length(ind));
end
subplot(2,2,1)
shadedErrorBar(edges(1:7,:), AL_df_avg,AL_df_sem);
hold on
ylim([0 1]);
xlabel('log2(Speed)');
axis square
ylabel('dF/F');

subplot(2,2,4)
for iArea = 1:3
    ind = find(Goodfits(iArea).dF>1.5);
    scatter(ind, iArea*ones(1,length(ind)), col(iArea,:));
    hold on
end

%ttest
p_ttestB = zeros(6,1);
h_ttestB = zeros(6,1);
[n_PM bin_PM] = histc(Goodfits(1).dF, edges_dF);
[n_AL bin_AL] = histc(Goodfits(3).dF, edges_dF);
for ibin = 1:6
    ind_PM = find(bin_PM==ibin);
    ind_AL = find(bin_AL==ibin);
    [h_ttestB1,p_ttestB1] = ttest2(Goodfits(1).speed(ind_PM),Goodfits(3).speed(ind_AL));
    p_ttestB(ibin,:) = p_ttestB1;
    h_ttestB(ibin,:) = h_ttestB1;
end
plot(mean(dF_avg(:,[1 3]),2), p_ttestB);
hold on
ylim([0 .000001]);
xlabel('dF/F');
axis square
ylabel('p-value');



%speed vs dF space
edges_dF = [0 .1  .2  .4 .8 1.6 3];
dF_inv = [7:-1:1];
edges_speed = log2(TF_vec./flipud(SF_vec));
dF_speed_density = zeros(7,7,3);
dF_speed_normspeed = zeros(7,7,3);
dF_speed_normdf = zeros(7,7,3);

for iArea= 1:3
    [n bin] = histc(log2(Goodfits(iArea).speed), edges_speed);
    for ibin = 1:7
        ind_ibin = find(bin==ibin);
        [n_dF bin_dF] = histc(Goodfits(iArea).dF(ind_ibin), edges_dF);
        dF_speed_density(:,ibin,iArea) = flipud(n_dF(1:7,:));
    end
    for ispeed = 1:7
        dF_speed_normspeed(:,ispeed,iArea) = dF_speed_density(:,ispeed,iArea)./max(dF_speed_density(:,ispeed,iArea),[],1);
    end
    for idf = 1:7
        dF_speed_normdf(idf,:,iArea) = dF_speed_density(idf,:,iArea)./max(dF_speed_density(idf,:,iArea),[],2);
    end
end

figure
start = 1;
for iArea = 1:3
    subplot(4,3,start)
    imagesq(dF_speed_density(:,:,iArea));
    colormap(gray);
    title(areas(iArea,:));
    ylabel('dF/F')
    subplot(4,3,start+3)
    imagesq(dF_speed_normdf(:,:,iArea));
    title('Norm dF');
    ylabel('dF/F')
    subplot(4,3,start+6)
    imagesq(dF_speed_normspeed(:,:,iArea));
    title('Norm speed');
    xlabel('speed')
    ylabel('dF/F')
    subplot(4,3,start+9)
    for ispeed = 1:7
        for idF = 1:7
            text(ispeed,dF_inv(idF), num2str(dF_speed_density(idF,ispeed,iArea)), 'Fontsize', 5);
            hold on
        end
    end
    xlim([1 7])
    ylim([1 7])
    axis square
    axis off
    start= start+1;
end
%     fn_out = fullfile(anal_base, [matrix '_' num2str(P) 'P_' inj '_all_areas_speed_vs_dF_matrix.pdf']);
%         print(gcf, '-dpdf', fn_out);
        

shade = [0 0 0; 0 0 0; 0.1 0.1 0.1; 0.2 0.2 0.2; 0.4 0.4 0.4; 0.6 0.6 0.6; 0.8 0.8 0.8];
figure;
for iArea = 1:3
    subplot(2,3,iArea)
    for idF = 2:7
        plot(edges_speed(1:7,:), squeeze(dF_speed_density(idF,:,iArea))', 'Color', shade(idF,:));
        hold on
    end
    if iArea == 1;
    ylabel('Number of boutons')
    end
    title(areas(iArea,:));
    subplot(2,3,iArea+3)
    for idF = 2:7
        plot(edges_speed(1:7,:), squeeze(dF_speed_normdf(idF,:,iArea))', 'Color', shade(idF,:));
        hold on
    end
    if iArea == 1;
    ylabel('Normalized fraction of boutons')
    end
    xlabel('log2(speed)')
end
     fn_out = fullfile(anal_base, [matrix '_' num2str(P) 'P_' inj '_all_areas_speed_vs_dF.pdf']);
         print(gcf, '-dpdf', fn_out);

figure;
for iArea = 1:3
    scatter(log2(Goodfits(iArea).speed), Goodfits(iArea).dF, 3, col(iArea,:));
    hold on
end
xlabel('log2(speed)')
ylabel('dF/F')
fn_out = fullfile(anal_base, [matrix '_' num2str(P) 'P_' inj '_all_areas_speed_vs_dF_scatter.pdf']);
         print(gcf, '-dpdf', fn_out);
