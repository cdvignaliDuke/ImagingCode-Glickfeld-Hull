areas = ['PM'; 'LM'; 'AL'];
col = strvcat('c', 'k', 'r');
P = 2;
matrix = 'SF5xTF5';
inj = 'V1';
anal_base = '\\zoloto\bigstorlab\Lindsey\Analysis\120213';
sum_base = 'G:\users\lindsey\analysisLG\experiments';
fn_summary = fullfile(sum_base, [matrix '_' num2str(P) 'P_' inj '_axon_summary'], 'all_fits.mat');
load(fn_summary);
fn_good = fullfile(sum_base, [matrix '_' num2str(P) 'P_' inj '_axon_summary'], 'good_fits.mat');
load(fn_good);

SF_vec0 = [.32 .16 .08 .04 .02]; %flipped to have low to high SF in square  %flipud
TF_vec0 = [1 2 4 8 15];
[tftf,sfsf] = meshgrid(TF_vec0,SF_vec0); 
TF_vec_adj = [1 2 4 8 16];
tftf2 = meshgrid(TF_vec_adj); 
speed_grid = tftf2./sfsf;
speed_grid_long = reshape(speed_grid', 1, 25);
SF_grid_long = reshape(sfsf', 1, 25);
TF_grid_long = reshape(tftf', 1, 25);
uspeeds = unique(speed_grid);

all_speed_fit = zeros(length(uspeeds), 2, 3);
all_TF_fit = zeros(length(TF_vec0), 2, 3);
all_SF_fit = zeros(length(SF_vec0), 2, 3);
all_SFxTF_fit = zeros(5,5,3);
for iArea = 1:3
    all_SFxTF_fit(:,:,iArea) = reshape(mean(Goodfits(iArea).plotfit,1),5,5)';    
    speed_tuning_fit = zeros(length(uspeeds),2);
    TF_tuning_fit = zeros(length(TF_vec0),2);
    SF_tuning_fit = zeros(length(SF_vec0),2);
    for ispeed = 1:length(uspeeds)
        ind = find(speed_grid_long == uspeeds(ispeed));
        if length(ind)>1
            goodfits_ispeed = mean(Goodfits(iArea).plotfit(:,ind),2);
        elseif length(ind) == 1
            goodfits_ispeed = Goodfits(iArea).plotfit(:,ind);
        end
        speed_tuning_fit(ispeed,1) = mean(goodfits_ispeed,1);
        speed_tuning_fit(ispeed,2) = std(goodfits_ispeed, [],1)./sqrt(size(goodfits_ispeed,1));
    end
    all_speed_fit(:,:,iArea) = speed_tuning_fit;
    for iTF = 1:length(TF_vec0)
        ind = find(TF_grid_long == TF_vec0(iTF));
        if length(ind)>1
            goodfits_iTF = mean(Goodfits(iArea).plotfit(:,ind),2);
        elseif length(ind) == 1
            goodfits_iTF = Goodfits(iArea).plotfit(:,ind);
        end
        TF_tuning_fit(iTF,1) = mean(goodfits_iTF,1);
        TF_tuning_fit(iTF,2) = std(goodfits_iTF, [],1)./sqrt(size(goodfits_iTF,1));
    end
    all_TF_fit(:,:,iArea) = TF_tuning_fit;
    for iSF = 1:length(SF_vec0)
        ind = find(SF_grid_long == SF_vec0(iSF));
        if length(ind)>1
            goodfits_iSF = mean(Goodfits(iArea).plotfit(:,ind),2);
        elseif length(ind) == 1
            goodfits_iSF = Goodfits(iArea).plotfit(:,ind);
        end
        SF_tuning_fit(iSF,1) = mean(goodfits_iSF,1);
        SF_tuning_fit(iSF,2) = std(goodfits_iSF, [],1)./sqrt(size(goodfits_iSF,1));
    end
    all_SF_fit(:,:,iArea) = SF_tuning_fit;
end

figure;
start = 1;
for iArea = 1:3
    subplot(2,3,start)
    imagesq(all_SFxTF_fit(:,:,iArea))
    title(areas(iArea,:));
    colormap(gray)
    subplot(2,3,4)
    errorbar(log2(TF_vec0), all_TF_fit(:,1,iArea),all_SF_fit(:,2,iArea),col(iArea,:));
    hold on
    title('TF')
    ylabel('dF/F')
    xlabel('log2(TF)')
    xlim([-1 5])
    ylim([0 .2])
    subplot(2,3,5)
    errorbar(log2(SF_vec0), all_SF_fit(:,1,iArea),all_SF_fit(:,2,iArea),col(iArea,:));
    hold on
    title('SF')
    xlabel('log2(SF)')
    ylim([0 .2])
    xlim([-6 -1])
    subplot(2,3,6)
    errorbar(log2(uspeeds), all_speed_fit(:,1,iArea),all_speed_fit(:,2,iArea),col(iArea,:));
    hold on
    title('speed')
    xlabel('log2(speed)')
    ylim([0 .2])
    xlim([1 10])
    start = start+1;
end
suptitle([matrix '  ' num2str(P) 'P  ' inj 'axons-  ' areas(1,:) '(' num2str(size(Goodfits(1).plotfit,1)) ')  ' areas(2,:) '(' num2str(size(Goodfits(2).plotfit,1)) ')  ' areas(3,:) '(' num2str(size(Goodfits(3).plotfit,1)) ')'])

fn_out = fullfile(anal_base, [matrix '_' num2str(P) 'P_' inj '_all_areas_tuning_curves.pdf']);
        print(gcf, '-dpdf', fn_out);