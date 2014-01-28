SF_vec0 = [.32 .16 .08 .04 .02]; %flipped to have low to high SF in square  %flipud
TF_vec0 = [1 2 4 8 15];
[tftf,sfsf] = meshgrid(TF_vec0,SF_vec0); 
speed_grid = tftf./sfsf;
speed_grid_long = reshape(speed_grid', 1, 25);
SF_grid_long = reshape(sfsf', 1, 25);
TF_grid_long = reshape(tftf', 1, 25);
uspeeds = unique(speed_grid);

all_speed_fit = zeros(length(uspeeds), 2, 3);
all_speed_data = zeros(length(uspeeds), 2, 3);
all_TF_fit = zeros(length(TF_vec0), 2, 3);
all_TF_data = zeros(length(TF_vec0), 2, 3);
all_SF_fit = zeros(length(SF_vec0), 2, 3);
all_SF_data = zeros(length(SF_vec0), 2, 3);
for iArea = 1:3
    nexp = all_fits(iArea).nexp;
    goodfits = [];
    gooddata = [];
    for iexp = 1:nexp
        n = all_fits(iArea).expt(iexp).n(1);
        for iCell = 1:n
        	if all_fits(iArea).expt(iexp).bouton(iCell).goodfit == 1
                goodfits = [goodfits; reshape(all_fits(iArea).expt(iexp).bouton(iCell).plotfit', 1, 25)];
                gooddata = [gooddata; reshape(all_fits(iArea).expt(iexp).bouton(iCell).dFoverF', 1, 25)];
            end
        end
    end
    goodfits_avg = reshape(mean(goodfits,1), 5, 5)';
    gooddata_avg = reshape(mean(gooddata,1), 5, 5)';
    speed_tuning_fit = zeros(length(uspeeds),2);
    speed_tuning_data = zeros(length(uspeeds),2);
    TF_tuning_fit = zeros(length(TF_vec0),2);
    TF_tuning_data = zeros(length(TF_vec0),2);
    SF_tuning_fit = zeros(length(SF_vec0),2);
    SF_tuning_data = zeros(length(SF_vec0),2);
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
    all_speed_fit(:,:,iArea) = speed_tuning_fit;
    all_speed_data(:,:,iArea) = speed_tuning_data;
    
    for iTF = 1:length(TF_vec0)
        ind = find(TF_grid_long == TF_vec0(iTF));
        iTF_goodfits = [];
        iTF_gooddata = [];
        for iCond = ind;
            iTF_goodfits = [iTF_goodfits; goodfits(:,iCond)];
            iTF_gooddata = [iTF_gooddata; gooddata(:,iCond)];
        end
        TF_tuning_fit(iTF,1) = mean(iTF_goodfits,1);
        TF_tuning_fit(iTF,2) = std(iTF_goodfits, [],1)./size(iTF_goodfits,1);
        TF_tuning_data(iTF,1) = mean(iTF_gooddata,1);
        TF_tuning_data(iTF,2) = std(iTF_gooddata, [],1)./size(iTF_gooddata,1);
    end
    all_TF_fit(:,:,iArea) = TF_tuning_fit;
    all_TF_data(:,:,iArea) = TF_tuning_data;
    
    for iSF = 1:length(SF_vec0)
        ind = find(SF_grid_long == SF_vec0(iSF));
        iSF_goodfits = [];
        iSF_gooddata = [];
        for iCond = ind;
            iSF_goodfits = [iSF_goodfits; goodfits(:,iCond)];
            iSF_gooddata = [iSF_gooddata; gooddata(:,iCond)];
        end
        SF_tuning_fit(iSF,1) = mean(iSF_goodfits,1);
        SF_tuning_fit(iSF,2) = std(iSF_goodfits, [],1)./size(iSF_goodfits,1);
        SF_tuning_data(iSF,1) = mean(iSF_gooddata,1);
        SF_tuning_data(iSF,2) = std(iSF_gooddata, [],1)./size(iSF_gooddata,1);
    end
    all_SF_fit(:,:,iArea) = SF_tuning_fit;
    all_SF_data(:,:,iArea) = SF_tuning_data;
    
    figure;
    subplot(4,2,1)
    imagesq(goodfits_avg);
    title('fit')
    subplot(4,2,2)
    imagesq(gooddata_avg);
    title('data')
    subplot(4,2,3)
    errorbar(TF_vec0, TF_tuning_fit(:,1), TF_tuning_fit(:,2));
    axis square;    
    xlim([0 15])
    ylim([0 .2])
    xlabel('TF')
    ylabel('dF/F')
    subplot(4,2,4)
    errorbar(TF_vec0, TF_tuning_data(:,1), TF_tuning_data(:,2));
    axis square; 
    xlim([0 15])
    ylim([0 .2])
    xlabel('TF')
    ylabel('dF/F')
    subplot(4,2,5)
    errorbar(SF_vec0, SF_tuning_data(:,1), SF_tuning_data(:,2));
    axis square; 
    xlim([0 .32])
    ylim([0 .2])
    xlabel('SF')
    ylabel('dF/F')
    subplot(4,2,6)
    errorbar(SF_vec0, SF_tuning_data(:,1), SF_tuning_data(:,2));
    axis square; 
    xlim([0 .32])
    ylim([0 .2])
    xlabel('SF')
    ylabel('dF/F')
    subplot(4,2,7)
    errorbar(uspeeds, speed_tuning_fit(:,1), speed_tuning_fit(:,2));
    axis square; 
    xlim([0 800])
    ylim([0 .2])
    xlabel('speed')
    ylabel('dF/F')
    subplot(4,2,8)
    errorbar(uspeeds, speed_tuning_data(:,1), speed_tuning_data(:,2));
    axis square; 
    xlim([0 800])
    ylim([0 .2])
    xlabel('speed')
    ylabel('dF/F')
    suptitle([areas(iArea,:) ' speed tuning curves']);
    colormap(gray);
    
    fn_out = fullfile(anal_base, [matrix '_' num2str(P) 'P_' inj '_' areas(iArea,:) 'speed_tuning_curves.pdf']);
    print(gcf, '-dpdf', fn_out);
end

col = strvcat('c','k','r');
figure;
for iArea = 1:3
    subplot(3,2,1)
    errorbar(TF_vec0, squeeze(all_TF_fit(:,1,iArea)), squeeze(all_TF_fit(:,2,iArea)),col(iArea,:));
    hold on;
    axis square;    
    title('fit')
    xlim([0 15])
    ylim([0 .2])
    xlabel('TF')
    ylabel('dF/F')
    subplot(3,2,2)
    errorbar(TF_vec0, squeeze(all_TF_data(:,1,iArea)), squeeze(all_TF_data(:,2,iArea)),col(iArea,:));
    hold on;
    axis square;
    title('data')
    xlim([0 15])
    ylim([0 .2])
    xlabel('TF')
    ylabel('dF/F')
    subplot(3,2,3)
    errorbar(SF_vec0, squeeze(all_SF_fit(:,1,iArea)), squeeze(all_SF_fit(:,2,iArea)),col(iArea,:));
    hold on;
    axis square;    
    title('fit')
    xlim([0 .32])
    ylim([0 .2])
    xlabel('SF')
    ylabel('dF/F')
    subplot(3,2,4)
    errorbar(SF_vec0, squeeze(all_SF_data(:,1,iArea)), squeeze(all_SF_data(:,2,iArea)),col(iArea,:));
    hold on;
    axis square;
    title('data')
    xlim([0 .32])
    ylim([0 .2])
    xlabel('SF')
    ylabel('dF/F')
    subplot(3,2,5)
    errorbar(uspeeds, squeeze(all_speed_fit(:,1,iArea)), squeeze(all_speed_fit(:,2,iArea)),col(iArea,:));
    hold on;
    axis square;    
    title('fit')
    xlim([0 800])
    ylim([0 .2])
    xlabel('speed')
    ylabel('dF/F')
    subplot(3,2,6)
    errorbar(uspeeds, squeeze(all_speed_data(:,1,iArea)), squeeze(all_speed_data(:,2,iArea)),col(iArea,:));
    hold on;
    axis square;
    title('data')
    xlim([0 800])
    ylim([0 .2])
    xlabel('speed')
    ylabel('dF/F')
end
fn_out = fullfile(anal_base, [matrix '_' num2str(P) 'P_' inj '_' 'allareas_speed_tuning_curves.pdf']);
    print(gcf, '-dpdf', fn_out);

