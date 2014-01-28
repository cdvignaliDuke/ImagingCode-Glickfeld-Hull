fig_base = '\\zmey\storlab\users\Lindsey\Projects\HVAs\Manuscript\Figures_2012';
fig = 5;

P = 2;
matrix = 'SF5xTF5';
inj = 'LM';
sum_base = 'G:\users\lindsey\analysisLG\experiments';

fn_summary = fullfile(sum_base, [matrix '_' num2str(P) 'P_' inj '_axon_summary'], 'all_fits.mat');
load(fn_summary);
fn_good = fullfile(sum_base, [matrix '_' num2str(P) 'P_' inj '_axon_summary'], 'good_fits.mat');
load(fn_good);

areas = strvcat('PM', 'AL', 'V1');
col = strvcat('c', 'r', 'k');
%example mouse

figure;
mouse ='M23';
expt = [6; 4];
x = 1:5;
y = 5:-1:1;
[x_grid y_grid] = meshgrid(x,y);
x_grid_long = reshape(x_grid', [25 1]);
y_grid_long = reshape(y_grid', [25 1]);
start = 1;
for iArea = 1:2
    area = areas(iArea,:);
    iexp = expt(iArea);
    n = all_fits(iArea).expt(iexp).n(1);
    goodfit = [];
    speed = [];
    for iCell = 1:n
        if all_fits(iArea).expt(iexp).bouton(iCell).goodfit ==1;
            goodfit = [goodfit iCell];
            speed = [speed all_fits(iArea).expt(iexp).bouton(iCell).speed];
        end
    end
    [order index] = sort(speed);
    percentiles = round([1:1:6]*(length(goodfit)/7));
    for iper = 1:6
        plotfit =  reshape(all_fits(iArea).expt(iexp).bouton(goodfit(index(percentiles(iper)))).plotfit', [25 1]);
        plotfit_norm = (plotfit./max(plotfit,[],1))*20; 
        subplot(6,6,start)
        for iCond = 1:25
            scatter(x_grid_long(iCond,1),y_grid_long(iCond,1), plotfit_norm(iCond,:).^2,'.k')
            hold on
            axis square
            axis off
            xlim([0 6])
            ylim([0 6])
        end
        title(num2str(all_fits(iArea).expt(iexp).bouton(goodfit(index(percentiles(iper)))).dF_fit));
        start = start+1;
    end
end

avg_tuning_long = zeros(25,3);
for iArea = 1:2
    iexp = expt(iArea);
    plotfit = [];
    n = all_fits(iArea).expt(iexp).n(1);
    for iCell = 1:n
        if all_fits(iArea).expt(iexp).bouton(iCell).goodfit ==1;
            plotfit = cat(3, plotfit, all_fits(iArea).expt(iexp).bouton(iCell).plotfit);
        end
    end
    avg_tuning = mean(plotfit,3);
    avg_tuning_long(:,iArea) = reshape((avg_tuning./max(max(avg_tuning,[],1),[],2))'*20,[25 1]);
end

for iArea = 1:2
    iexp = expt(iArea);
    subplot(6,6,start)
    for iCond = 1:25
        scatter(x_grid_long(iCond,1),y_grid_long(iCond,1), avg_tuning_long(iCond,iArea)^2,'.k')
        hold on
        axis square
        axis off
        xlim([0 6])
        ylim([0 6])
    end
    start = start+1;
    title(num2str(all_fits(iArea).expt(iexp).n(2)));
end

%other example mice
mouse ='M10';
expt = [2; 1];

avg_tuning_long = zeros(25,3);
for iArea = 1:2
    iexp = expt(iArea);
    plotfit = [];
    n = all_fits(iArea).expt(iexp).n(1);
    for iCell = 1:n
        if all_fits(iArea).expt(iexp).bouton(iCell).goodfit ==1;
            plotfit = cat(3, plotfit, all_fits(iArea).expt(iexp).bouton(iCell).plotfit);
        end
    end
    avg_tuning = mean(plotfit,3);
    avg_tuning_long(:,iArea) = reshape((avg_tuning./max(max(avg_tuning,[],1),[],2))'*20,[25 1]);
end

for iArea = 1:2
    iexp = expt(iArea);
    subplot(6,6,start)
    for iCond = 1:25
        scatter(x_grid_long(iCond,1),y_grid_long(iCond,1), avg_tuning_long(iCond,iArea)^2,'.k')
        hold on
        axis square
        axis off
        xlim([0 6])
        ylim([0 6])
    end
    start = start+1;
    title(num2str(all_fits(iArea).expt(iexp).n(2)));
end

mouse ='Y28';
expt = [3; 2];

avg_tuning_long = zeros(25,3);
for iArea = 1:2
    iexp = expt(iArea);
    plotfit = [];
    n = all_fits(iArea).expt(iexp).n(1);
    for iCell = 1:n
        if all_fits(iArea).expt(iexp).bouton(iCell).goodfit ==1;
            plotfit = cat(3, plotfit, all_fits(iArea).expt(iexp).bouton(iCell).plotfit);
        end
    end
    avg_tuning = mean(plotfit,3);
    avg_tuning_long(:,iArea) = reshape((avg_tuning./max(max(avg_tuning,[],1),[],2))'*20,[25 1]);
end

for iArea = 1:2
    iexp = expt(iArea);
    subplot(6,6,start)
    for iCond = 1:25
        scatter(x_grid_long(iCond,1),y_grid_long(iCond,1), avg_tuning_long(iCond,iArea)^2,'.k')
        hold on
        axis square
        axis off
        xlim([0 6])
        ylim([0 6])
    end
    start = start+1;
    title(num2str(all_fits(iArea).expt(iexp).n(2)));
end
        
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

all_SFxTF_long = zeros(25,3);
for iArea = 1:3
    all_SFxTF_long(:,iArea) = reshape((all_SFxTF_fit(:,:,iArea)./max(max(all_SFxTF_fit(:,:,iArea),[],1),[],2))'*20,[25 1]);        
end

for iArea = 1:2
    subplot(6,6,start)
    for iCond = 1:25
        scatter(x_grid_long(iCond,1),y_grid_long(iCond,1), all_SFxTF_long(iCond,iArea)^2,'.k')
        hold on
        axis square
        axis off
        xlim([0 6])
        ylim([0 6])
    end
    start = start+1;
    title(num2str(length(Goodfits(iArea).dF)));
end
fn_out = fullfile(fig_base, ['Figure' num2str(fig)], [matrix '_' num2str(P) 'P_' inj '_all_areas_example_tuning.ps']);
        print(gcf, '-depsc2', fn_out);

figure;
for iArea = 1:3
    title(num2str(length(Goodfits(iArea).dF)));
    colormap(gray)
    subplot(1,3,1)
    errorbar(log2(TF_vec0), all_TF_fit(:,1,iArea),all_SF_fit(:,2,iArea),col(iArea,:));
    hold on
    title('TF')
    ylabel('dF/F')
    xlabel('log2(TF)')
    xlim([-1 5])
    ylim([0 .25])
    subplot(1,3,2)
    errorbar(log2(SF_vec0), all_SF_fit(:,1,iArea),all_SF_fit(:,2,iArea),col(iArea,:));
    hold on
    title('SF')
    xlabel('log2(SF)')
    ylim([0 .25])
    xlim([-6 -1])
    subplot(1,3,3)
    errorbar(log2(uspeeds), all_speed_fit(:,1,iArea),all_speed_fit(:,2,iArea),col(iArea,:));
    hold on
    title('speed')
    xlabel('log2(speed)')
    ylim([0 .25])
    xlim([1 10])
    start = start+1;
end
suptitle([matrix '  ' num2str(P) 'P  ' inj 'axons-  ' areas(1,:) '(' num2str(size(Goodfits(1).plotfit,1)) ')  ' areas(2,:) '(' num2str(size(Goodfits(2).plotfit,1)) ')  ' areas(3,:) '(' num2str(size(Goodfits(3).plotfit,1)) ')'])

fn_out = fullfile(fig_base, ['Figure' num2str(fig)], [matrix '_' num2str(P) 'P_' inj '_all_areas_tuning_curves.ps']);
        print(gcf, '-depsc', fn_out);
        
%scatters
for iArea = 1:2;
    area = areas(iArea,:);
    iexp = expt(iArea);

    subplot(2,3,iArea)
    n = all_fits(iArea).expt(iexp).n(1);
    for iCell = 1:n
        if all_fits(iArea).expt(iexp).bouton(iCell).goodfit ==1;
            scatter(log2(all_fits(iArea).expt(iexp).bouton(iCell).TF_fit), log2(all_fits(iArea).expt(iexp).bouton(iCell).SF_fit),3,'k');
            axis square
            box on
            hold on
            xlim([log2(1) log2(15)])
            ylim([log2(0.02) log2(0.32)])
        end
    end
    title(area)
end
        
