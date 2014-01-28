areas = strvcat('PM', 'LM', 'AL', 'RL', 'AM');
col = strvcat('c', 'k', 'r', 'm', 'b');
P = 2;
matrix = 'SF5xTF5';
inj = 'V1';
sum_base = 'G:\users\lindsey\analysisLG\experiments';

fn_summary = fullfile(sum_base, [matrix '_' num2str(P) 'P_' inj '_axon_summary'], 'all_fits.mat');
load(fn_summary);
fn_good = fullfile(sum_base, [matrix '_' num2str(P) 'P_' inj '_axon_summary'], 'good_fits.mat');
load(fn_good);

figure;
leg = [];
for iArea= 1:5;
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
    leg = [leg; size(Goodfits(iArea).TF,1)];
end
legend(num2str(leg));
fn_out = fullfile('\\zmey\storlab\users\Lindsey\Projects\HVAs\2P Axon Imaging\Analysis\121004', [matrix '_' num2str(P) 'P_' inj '_TF_SF_speed_cdfs.ps']);
print(gcf, '-depsc', fn_out);


%squares and tuning curves
x = 1:5;
y = 5:-1:1;
[x_grid y_grid] = meshgrid(x,y);
x_grid_long = reshape(x_grid', [25 1]);
y_grid_long = reshape(y_grid', [25 1]);
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

all_speed_fit = zeros(length(uspeeds), 2, 5);
all_TF_fit = zeros(length(TF_vec0), 2, 5);
all_SF_fit = zeros(length(SF_vec0), 2, 5);
all_speed_fit_norm = zeros(length(uspeeds), 2, 5);
all_TF_fit_norm = zeros(length(TF_vec0), 2, 5);
all_SF_fit_norm = zeros(length(SF_vec0), 2, 5);
all_SFxTF_fit = zeros(5,5,5);
for iArea = 1:5
    all_SFxTF_fit(:,:,iArea) = reshape(mean(Goodfits(iArea).plotfit,1),5,5)';
    plotfit_norm = Goodfits(iArea).plotfit./repmat(max(Goodfits(iArea).plotfit,[],2),[1 25]);
    speed_tuning_fit = zeros(length(uspeeds),2);
    TF_tuning_fit = zeros(length(TF_vec0),2);
    SF_tuning_fit = zeros(length(SF_vec0),2);
    speed_tuning_fit_norm = zeros(length(uspeeds),2);
    TF_tuning_fit_norm = zeros(length(TF_vec0),2);
    SF_tuning_fit_norm = zeros(length(SF_vec0),2);
    for ispeed = 1:length(uspeeds)
        ind = find(speed_grid_long == uspeeds(ispeed));
        if length(ind)>1
            goodfits_ispeed = mean(Goodfits(iArea).plotfit(:,ind),2);
            goodfits_ispeed_norm = mean(plotfit_norm(:,ind),2);
        elseif length(ind) == 1
            goodfits_ispeed = Goodfits(iArea).plotfit(:,ind);
            goodfits_ispeed_norm = plotfit_norm(:,ind);
        end
        speed_tuning_fit(ispeed,1) = mean(goodfits_ispeed,1);
        speed_tuning_fit(ispeed,2) = std(goodfits_ispeed, [],1)./sqrt(size(goodfits_ispeed,1));
        speed_tuning_fit_norm(ispeed,1) = mean(goodfits_ispeed_norm,1);
        speed_tuning_fit_norm(ispeed,2) = std(goodfits_ispeed_norm, [],1)./sqrt(size(goodfits_ispeed_norm,1));
    end
    all_speed_fit(:,:,iArea) = speed_tuning_fit;
    all_speed_fit_norm(:,:,iArea) = speed_tuning_fit_norm;
    for iTF = 1:length(TF_vec0)
        ind = find(TF_grid_long == TF_vec0(iTF));
        if length(ind)>1
            goodfits_iTF = mean(Goodfits(iArea).plotfit(:,ind),2);
            goodfits_iTF_norm = mean(plotfit_norm(:,ind),2);
        elseif length(ind) == 1
            goodfits_iTF = Goodfits(iArea).plotfit(:,ind);
            goodfits_iTF_norm = plotfit_norm(:,ind);
        end
        TF_tuning_fit(iTF,1) = mean(goodfits_iTF,1);
        TF_tuning_fit(iTF,2) = std(goodfits_iTF, [],1)./sqrt(size(goodfits_iTF,1));
        TF_tuning_fit_norm(iTF,1) = mean(goodfits_iTF_norm,1);
        TF_tuning_fit_norm(iTF,2) = std(goodfits_iTF_norm, [],1)./sqrt(size(goodfits_iTF_norm,1));
    end
    all_TF_fit(:,:,iArea) = TF_tuning_fit;
    all_TF_fit_norm(:,:,iArea) = TF_tuning_fit_norm;
    for iSF = 1:length(SF_vec0)
        ind = find(SF_grid_long == SF_vec0(iSF));
        if length(ind)>1
            goodfits_iSF = mean(Goodfits(iArea).plotfit(:,ind),2);
            goodfits_iSF_norm = mean(plotfit_norm(:,ind),2);
        elseif length(ind) == 1
            goodfits_iSF = Goodfits(iArea).plotfit(:,ind);
            goodfits_iSF_norm = plotfit_norm(:,ind);
        end
        SF_tuning_fit(iSF,1) = mean(goodfits_iSF,1);
        SF_tuning_fit(iSF,2) = std(goodfits_iSF, [],1)./sqrt(size(goodfits_iSF,1));
        SF_tuning_fit_norm(iSF,1) = mean(goodfits_iSF_norm,1);
        SF_tuning_fit_norm(iSF,2) = std(goodfits_iSF_norm, [],1)./sqrt(size(goodfits_iSF_norm,1));
    end
    all_SF_fit(:,:,iArea) = SF_tuning_fit;
    all_SF_fit_norm(:,:,iArea) = SF_tuning_fit_norm;
end

all_SFxTF_long = zeros(25,5);
for iArea = 1:5
    all_SFxTF_long(:,iArea) = reshape((all_SFxTF_fit(:,:,iArea)./max(max(all_SFxTF_fit(:,:,iArea),[],1),[],2))'*50,[25 1]);        
end
area_order = [1:5];
figure;
start = 1;
for iArea = 1:5
    subplot(1,5,start)
    for iCond = 1:25
        scatter(x_grid_long(iCond,1),y_grid_long(iCond,1), all_SFxTF_long(iCond,area_order(iArea))^2,'.k')
        hold on
        axis square
        axis off
        xlim([0 6])
        ylim([0 6])
    end
    title(num2str(length(Goodfits(area_order(iArea)).dF)));
    colormap(gray)
    start = start+1;
end
fn_out = fullfile('\\zmey\storlab\users\Lindsey\Projects\HVAs\2P Axon Imaging\Analysis\121004', [matrix '_' num2str(P) 'P_' inj '_SFxTF_squares.ps']);
print(gcf, '-depsc', fn_out);


for iArea = 4:5
    image = areas(iArea,:);
    nexp = all_fits(iArea).nexp;
    all_goodfits = [];
    figure;
    for iexp = 1:nexp
        nCells = all_fits(iArea).expt(iexp).n(1);
        mouse = all_fits(iArea).expt(iexp).mouse;
        date = all_fits(iArea).expt(iexp).date;
        userun = all_fits(iArea).expt(iexp).userun;
        for iCell = 1:nCells
            if all_fits(iArea).expt(iexp).bouton(iCell).goodfit==1
                all_goodfits = cat(3,all_goodfits,all_fits(iArea).expt(iexp).bouton(iCell).plotfit);
            end
        end
        avg_goodfits = mean(all_goodfits,3);
        avg_goodfits_long = reshape((avg_goodfits./max(max(avg_goodfits,[],1),[],2))'*50,[25 1]);
        subplot(4,4,iexp)
        for iCond = 1:25
            scatter(x_grid_long(iCond,1),y_grid_long(iCond,1), avg_goodfits_long(iCond,1)^2,'.k')
            hold on
            axis square
            axis off
            xlim([0 6])
            ylim([0 6])
        end
        title([mouse ' ' date ' ' num2str(userun) ' ' num2str(all_fits(iArea).expt(iexp).n(2))]);
    end
    fn_out = fullfile('\\zmey\storlab\users\Lindsey\Projects\HVAs\2P Axon Imaging\Analysis\121004', [matrix '_' num2str(P) 'P_' inj '_' image '_allFOVs.ps']);
    print(gcf, '-depsc', fn_out);
end

mouse_list = {'M31' 'M38' 'M41' 'CM114'};
for iArea = 4:5
    nexp = all_fits(iArea).nexp;    
    figure;
    for iMouse = 1:length(mouse_list)
        all_goodfits = [];
        mouse_name = char(mouse_list(iMouse));
        for iexp = 1:nexp
            mouse = all_fits(iArea).expt(iexp).mouse;
            if length(mouse) == length(mouse_name)
                if mouse == mouse_name
                    nCells = all_fits(iArea).expt(iexp).n(1);
                    date = all_fits(iArea).expt(iexp).date;
                    userun = all_fits(iArea).expt(iexp).userun;
                    for iCell = 1:nCells
                        if all_fits(iArea).expt(iexp).bouton(iCell).goodfit==1
                            all_goodfits = cat(3,all_goodfits,all_fits(iArea).expt(iexp).bouton(iCell).plotfit);
                        end
                    end
                end
            end
        end
        avg_goodfits = mean(all_goodfits,3);
        if length(avg_goodfits)>0
            avg_goodfits_long = reshape((avg_goodfits./max(max(avg_goodfits,[],1),[],2))'*50,[25 1]);
            subplot(1,5,iMouse)
            for iCond = 1:25
                scatter(x_grid_long(iCond,1),y_grid_long(iCond,1), avg_goodfits_long(iCond,1)^2,'.k')
                hold on
                axis square
                axis off
                xlim([0 6])
                ylim([0 6])
            end
            title([mouse_name ' ' num2str(size(all_goodfits,3))]);
        end
    end
    fn_out = fullfile('\\zmey\storlab\users\Lindsey\Projects\HVAs\2P Axon Imaging\Analysis\121004', [matrix '_' num2str(P) 'P_' inj '_' image '_allmice.ps']);
    print(gcf, '-depsc', fn_out);
end

figure;
for iArea = 1:5
    subplot(2,3,1)
    errorbar(log2(TF_vec0), all_TF_fit(:,1,iArea),all_SF_fit(:,2,iArea),col(iArea,:));
    hold on
    title('TF')
    ylabel('dF/F')
    xlabel('log2(TF)')
    xlim([-1 5])
    ylim([0 .2])
    subplot(2,3,2)
    errorbar(log2(SF_vec0), all_SF_fit(:,1,iArea),all_SF_fit(:,2,iArea),col(iArea,:));
    hold on
    title('SF')
    xlabel('log2(SF)')
    ylim([0 .2])
    xlim([-6 -1])
    subplot(2,3,3)
    errorbar(log2(uspeeds), all_speed_fit(:,1,iArea),all_speed_fit(:,2,iArea),col(iArea,:));
    hold on
    title('speed')
    xlabel('log2(speed)')
    ylim([0 .2])
    xlim([1 10])
    start = start+1;
    subplot(2,3,4)
    errorbar(log2(TF_vec0), all_TF_fit_norm(:,1,iArea),all_SF_fit_norm(:,2,iArea),col(iArea,:));
    hold on
    title('TF')
    ylabel('dF/F')
    xlabel('log2(TF)')
    xlim([-1 5])
    ylim([0 .5])
    subplot(2,3,5)
    errorbar(log2(SF_vec0), all_SF_fit_norm(:,1,iArea),all_SF_fit_norm(:,2,iArea),col(iArea,:));
    hold on
    title('SF')
    xlabel('log2(SF)')
    ylim([0 .5])
    xlim([-6 -1])
    subplot(2,3,6)
    errorbar(log2(uspeeds), all_speed_fit_norm(:,1,iArea),all_speed_fit_norm(:,2,iArea),col(iArea,:));
    hold on
    title('speed')
    xlabel('log2(speed)')
    ylim([0 .5])
    xlim([1 10])
    start = start+1;
end
fn_out = fullfile('\\zmey\storlab\users\Lindsey\Projects\HVAs\2P Axon Imaging\Analysis\121004', [matrix '_' num2str(P) 'P_' inj '_TF_SF_speed_tuningcurves.ps']);
print(gcf, '-depsc', fn_out);
        
