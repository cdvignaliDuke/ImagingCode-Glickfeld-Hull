fig_base = '\\zmey\storlab\users\Lindsey\Projects\HVAs\Manuscript\Resubmission_Figures';
fig = 7;

P = 2;
matrix = 'SF5xTF5';
inj = 'LM';
sum_base = 'G:\users\lindsey\analysisLG\experiments';

fn_summary = fullfile(sum_base, [matrix '_' num2str(P) 'P_' inj '_axon_summary'], 'all_fits.mat');
load(fn_summary);
fn_good = fullfile(sum_base, [matrix '_' num2str(P) 'P_' inj '_axon_summary'], 'good_fits.mat');
load(fn_good);

%individual mouse tuning
mouse_list = {'M23' 'M10' 'Y28'};
areas = strvcat('PM', 'AL');
fn_summary = fullfile(sum_base, [matrix '_' num2str(P) 'P_' inj '_axon_summary'], 'all_fits.mat');
load(fn_summary);
figure;
for iMouse = 1:3
    mouse = char(mouse_list(iMouse));
    for iArea = 1:2
        image = areas(iArea,:);
        sum_base = 'G:\users\lindsey\analysisLG\experiments';
        list_fn = fullfile(sum_base, [matrix '_' num2str(P) 'P_' image], [matrix '_' num2str(P) 'P_' inj '_' image '_exp_list.mat']);
        load(list_fn);
        plotfit = [];
        nexp = all_fits(iArea).nexp;
        for iexp = 1:nexp
            exp_mouse = all_fits(iArea).expt(iexp).mouse;
            if length(exp_mouse)==length(mouse)
                if exp_mouse == mouse;
                    userun = all_fits(iArea).expt(iexp).userun;
                    date = all_fits(iArea).expt(iexp).date;
                    base = 'G:\users\lindsey\analysisLG\active mice';    
                    outDir = fullfile(base, mouse,date);  
                    fn_out = fullfile(outDir, 'analysis', [date '_' mouse '_run' num2str(userun) '_lbub_fits.mat']);
                    load(fn_out)
                    for iCell = 1:all_fits(iArea).expt(iexp).n(1)
                        if find(goodfit_ind == iCell)
                            plotfit_ind = reshape(all_fits(iArea).expt(iexp).bouton(iCell).plotfit',[25 1]);
                            plotfit = [plotfit plotfit_ind];
                        end
                    end
                end
            end
        end
        plotfit_avg = mean(plotfit,2);
        plotfit_norm = (plotfit_avg./max(plotfit_avg,[],1))*55; 
        subplot(3,3,iArea+(iMouse-1)*3)
        for iCond = 1:25
            scatter(x_grid_long(iCond,1),y_grid_long(iCond,1), plotfit_norm(iCond,:).^2,'.k')
            hold on
            axis square
            axis off
            xlim([0 6])
            ylim([0 6])
        end
        title(num2str(size(plotfit,2)));
    end
end
fn_out = fullfile(fig_base, ['Figure' num2str(fig)], [matrix '_' num2str(P) 'P_' inj '_tuning_example_mice.ps']);        
print(gcf, '-depsc2', fn_out);

figure;
inj_col = strvcat('g', 'y','c', 'r');
inj_list = strvcat('V1', 'LM');
start = 1;
for i_inj = 1:2
    inj = inj_list(i_inj,:);
    fn_good = fullfile(sum_base, [matrix '_' num2str(P) 'P_' inj '_axon_summary'], 'good_fits.mat');
    load(fn_good);
    if inj == 'V1'
        area_order = [1,3,2];
    elseif inj =='LM'
        area_order = [1,2,3];
    end

    for iArea= 1:2;
        subplot(2,3,1);
        [H, stats, xCDF, yCDF] = cdfplot_LG(log2(Goodfits(area_order(iArea)).TF));
        plot(xCDF, yCDF, inj_col(start,:));
        hold on
        title('TF')
        axis square
        xlabel('log2(TF)');
        xlim([log2(1) log2(15)])
        subplot(2,3,2);
        [H, stats, xCDF, yCDF] = cdfplot_LG(log2(Goodfits(area_order(iArea)).SF));
        plot(xCDF, yCDF, inj_col(start,:));
        hold on
        xlabel('log2(SF)');
        axis square
        title('SF')
        xlim([log2(.02) log2(.32)])
        subplot(2,3,3);
        [H, stats, xCDF, yCDF] = cdfplot_LG(log2(Goodfits(area_order(iArea)).speed));
        plot(xCDF, yCDF, inj_col(start,:));
        hold on
        xlabel('log2(speed)');
        axis square
        title('Speed')
        xlim([log2(1/.32) log2(15/.02)])
        start = start+1;
    end
end
fn_out = fullfile(fig_base, ['Figure' num2str(fig)], [matrix '_' num2str(P) 'P_' inj '_all_areas_cdfs.ps']);
        print(gcf, '-depsc', fn_out);

figure;        
mouse_list = {'M10' 'Y28' 'M23' 'M15' 'Y25' 'M21'}; 
inj = 'LM';
matrix = 'SF5xTF5';
P = 2;
areas = strvcat('PM', 'AL', 'V1');
for iMouse = 1:length(mouse_list);
    mouse = mouse_list{iMouse};
    medians = [NaN NaN];
    for iArea = 1:2
        image = areas(iArea,:);
        sum_base = 'G:\users\lindsey\analysisLG\experiments';
        list_fn = fullfile(sum_base, [matrix '_' num2str(P) 'P_' image], [matrix '_' num2str(P) 'P_' inj '_' image '_exp_list.mat']);
        load(list_fn);
        median_speed = [];
        nexp = all_fits(iArea).nexp;
        for iexp = 1:nexp
            exp_mouse = exp_list.mouse_mat{iexp};
            if length(exp_mouse)==length(mouse)
                if exp_mouse == mouse;
                    if all_fits(iArea).expt(iexp).n(2)>25
                        median_speed = [median_speed all_fits(iArea).expt(iexp).median_speed];
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
    subplot(1,3,1)
    semilogy(1:2, medians, '-ob')
    axis square
    hold on
    title([areas(1,:) ' vs ' areas(2,:)])
    xlim([0 3])
    ylim([1 1000])
end
semilogy(1, median(Goodfits(1).speed), 'sk');
hold on
semilogy(2, median(Goodfits(2).speed), 'sk');
fn_out = fullfile(fig_base, ['Figure' num2str(fig)], [matrix '_' num2str(P) 'P_' inj '_all_mice_medians.ps']);
        print(gcf, '-depsc', fn_out);

%separate AL/PM with V1/LM injections
inj_col = strvcat('k', 'g');
inj_list = strvcat('V1', 'LM');
for i_inj = 1:2
    inj = inj_list(i_inj,:);
    fn_good = fullfile(sum_base, [matrix '_' num2str(P) 'P_' inj '_axon_summary'], 'good_fits.mat');
    load(fn_good);
    for iArea = 1:2
        if inj == 'V1'
            area_order = [1,3,2];
        elseif inj =='LM'
            area_order = [1,2,3];
        end
        subplot(2,3,4+iArea);
        [H, stats, xCDF, yCDF] = cdfplot_LG(log2(Goodfits(area_order(iArea)).speed));
        plot(xCDF, yCDF, inj_col(i_inj,:));
        hold on
        xlabel('log2(speed)');
        axis square
        title('Speed')
        xlim([log2(1/.32) log2(15/.02)])
    end
end

%average individual mice
mouse_list = {'M10' 'Y28' 'M23' 'M15' 'Y25' 'M21'}; 
inj = 'LM';
areas = strvcat('PM', 'AL');
matrix = 'SF5xTF5';
col = strvcat('c', 'r');
P = 2;
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
avg_tuning = [];
figure;
for iArea = 1:2
    image = areas(iArea,:);
    sum_base = 'G:\users\lindsey\analysisLG\experiments';
    list_fn = fullfile(sum_base, [matrix '_' num2str(P) 'P_' image], [matrix '_' num2str(P) 'P_' inj '_' image '_exp_list.mat']);
    load(list_fn);
    plotfits = [];
    plotfits_long = [];
    for iMouse = 1:length(mouse_list);
        mouse = mouse_list{iMouse};
        plotfit = [];
        nexp = all_fits(iArea).nexp;
        for iexp = 1:nexp
            exp_mouse = all_fits(iArea).expt(iexp).mouse;
            if length(exp_mouse)==length(mouse)
                if exp_mouse == mouse;
                    userun = all_fits(iArea).expt(iexp).userun;
                    date = all_fits(iArea).expt(iexp).date;
                    base = 'G:\users\lindsey\analysisLG\active mice';    
                    outDir = fullfile(base, mouse,date);  
                    fn_out = fullfile(outDir, 'analysis', [date '_' mouse '_run' num2str(userun) '_lbub_fits.mat']);
                    load(fn_out)
                    for iCell = 1:all_fits(iArea).expt(iexp).n(1)
                        if find(goodfit_ind == iCell)
                        plotfit = cat(3,plotfit, all_fits(iArea).expt(iexp).bouton(iCell).plotfit);
                        end
                    end
                end
            end
        end
        if size(plotfit,3)>25
            plotfit_avg = mean(plotfit,3);
            plotfits = cat(3, plotfits, plotfit_avg./max(max(plotfit_avg,[],2),[],1));
            plotfit_long_norm = reshape(plotfit_avg', [25 1])./max(max(plotfit_avg,[],2),[],1);
            plotfits_long = [plotfits_long plotfit_long_norm]; 
        end
        avg_tuning(iArea).mouse(iMouse).name = mouse;
        avg_tuning(iArea).mouse(iMouse).number = size(plotfit,3);
        avg_tuning(iArea).mouse(iMouse).tuning = plotfit_avg;
        avg_tuning(iArea).mouse(iMouse).tuning_norm = plotfit_avg./max(max(plotfit_avg,[],2),[],1);
    end
    subplot(2,3,iArea)
    plotfit_all =  reshape(mean(plotfits,3)', [25 1]);
    plotfit_norm = (plotfit_all./max(plotfit_all,[],1))*75; 
    for iCond = 1:25
        scatter(x_grid_long(iCond,1),y_grid_long(iCond,1), plotfit_norm(iCond,:).^2,'.k')
        hold on
        axis square
        axis off
        xlim([0 6])
        ylim([0 6])
    end
    title(num2str(size(plotfits,3)));
    
    speed_tuning_fit = zeros(length(uspeeds),2);
    TF_tuning_fit = zeros(length(TF_vec0),2);
    SF_tuning_fit = zeros(length(SF_vec0),2);
    for iTF = 1:length(TF_vec0)
        ind = find(TF_grid_long == TF_vec0(iTF));
        if length(ind)>1
            goodfits_iTF = mean(plotfits_long(ind,:),1);
        elseif length(ind) == 1
            goodfits_iTF = plotfits_long(ind,:);
        end
        TF_tuning_fit(iTF,1) = mean(goodfits_iTF,2);
        TF_tuning_fit(iTF,2) = std(goodfits_iTF, [],2)./sqrt(size(goodfits_iTF,2));
    end
    for iSF= 1:length(SF_vec0)
        ind = find(SF_grid_long == SF_vec0(iSF));
        if length(ind)>1
            goodfits_iSF = mean(plotfits_long(ind,:),1);
        elseif length(ind) == 1
            goodfits_iSF = plotfits_long(ind,:);
        end
        SF_tuning_fit(iSF,1) = mean(goodfits_iSF,2);
        SF_tuning_fit(iSF,2) = std(goodfits_iSF, [],2)./sqrt(size(goodfits_iSF,2));
    end
    for ispeed = 1:length(uspeeds)
        ind = find(speed_grid_long == uspeeds(ispeed));
        if length(ind)>1
            goodfits_ispeed = mean(plotfits_long(ind,:),1);
        elseif length(ind) == 1
            goodfits_ispeed = plotfits_long(ind,:);
        end
        speed_tuning_fit(ispeed,1) = mean(goodfits_ispeed,2);
        speed_tuning_fit(ispeed,2) = std(goodfits_ispeed, [],2)./sqrt(size(goodfits_ispeed,2));
    end
    subplot(2,3,4)
    errorbar(log2(TF_vec0), TF_tuning_fit(:,1),TF_tuning_fit(:,2),col(iArea,:));
    hold on
    title('TF')
    ylabel('dF/F')
    xlabel('log2(TF)')
    xlim([-1 5])
    ylim([0 1])
    subplot(2,3,5)
    errorbar(log2(SF_vec0),SF_tuning_fit(:,1),SF_tuning_fit(:,2),col(iArea,:));
    hold on
    title('SF')
    xlabel('log2(SF)')
    ylim([0 1])
    xlim([-6 -1])
    subplot(2,3,6)
    errorbar(log2(uspeeds), speed_tuning_fit(:,1),speed_tuning_fit(:,2),col(iArea,:));
    hold on
    title('speed')
    xlabel('log2(speed)')
    ylim([0 1])
    xlim([1 10])
end

fn_summary = fullfile(sum_base, [matrix '_' num2str(P) 'P_' inj '_axon_summary'], 'avg_tuning.mat');
save(fn_summary, 'avg_tuning');

fn_out = fullfile(fig_base, ['Figure' num2str(fig)], [matrix '_' num2str(P) 'P_' inj '_tuning_mouse_avg.ps']);        
print(gcf, '-depsc2', fn_out); 
        
%stats
r_same = zeros(length(mouse_list), length(mouse_list), 3);
for iArea = 1:2
    for iMouse1 = 1:length(mouse_list)
        if avg_tuning(iArea).mouse(iMouse1).number>25
            for iMouse2 = 1:length(mouse_list)
                if avg_tuning(iArea).mouse(iMouse2).number>25
                    if iMouse2 ~= iMouse1
                        r_same(iMouse1,iMouse2,iArea) = corr(reshape(avg_tuning(iArea).mouse(iMouse1).tuning_norm',[25,1]), reshape(avg_tuning(iArea).mouse(iMouse2).tuning_norm',[25,1]));
                    else
                        r_same(iMouse1,iMouse2,iArea) = NaN;
                    end
                end
            end
        else
            r_same(iMouse1,:,iArea) = NaN;
            r_same(:,iMouse1,iArea) = NaN;
        end
    end
    r_same_long(:,iArea) = triu2vec(r_same(:,:,iArea));
end
r_same_avg = nanmean(r_same_long,1);
r_same_sem = nanstd(r_same_long,[],1)./sqrt(sum(isnan(r_same_long),1));

r_same_all = reshape(r_same_long, [size(r_same_long,1)*size(r_same_long,2) 1]);
r_same_all_avg = nanmean(r_same_all,1);
r_same_all_sem = nanstd(r_same_all,[],1)./sqrt(sum(isnan(r_same_all),1));

r_pair = zeros(length(mouse_list), 1);
area_order = [1 2];
for ipair = 1
    for iMouse = 1:length(mouse_list)
        if avg_tuning(area_order(ipair, 1)).mouse(iMouse).number>25
            if avg_tuning(area_order(ipair, 2)).mouse(iMouse).number>25
                r_pair(iMouse,ipair) = corr(reshape(avg_tuning(area_order(ipair, 1)).mouse(iMouse).tuning_norm',[25,1]), reshape(avg_tuning(area_order(ipair, 2)).mouse(iMouse).tuning_norm',[25,1]));
            else
                r_pair(iMouse,ipair) = NaN;
            end
        else
            r_pair(iMouse,ipair) = NaN;
        end
    end
end
r_pair_avg = nanmean(r_pair,1);
r_pair_sem = nanstd(r_pair,[],1)./sqrt(sum(isnan(r_pair),1));

[h_PM_AL p_PM_AL] = ttest2(r_same_long(:,1), r_pair(:,1));
[h_AL_PM p_AL_PM] = ttest2(r_same_long(:,2), r_pair(:,1));

r_diff = zeros(length(mouse_list), length(mouse_list), 3);
area_order = [1 2];
for ipair = 1
    for iMouse1 = 1:length(mouse_list)
        if avg_tuning(area_order(ipair, 1)).mouse(iMouse1).number>25
            for iMouse2 = 1:length(mouse_list)
                if avg_tuning(area_order(ipair, 2)).mouse(iMouse2).number>25
                    r_diff(iMouse1,iMouse2,ipair) = corr(reshape(avg_tuning(area_order(ipair, 1)).mouse(iMouse1).tuning_norm',[25,1]), reshape(avg_tuning(area_order(ipair, 2)).mouse(iMouse2).tuning_norm',[25,1]));
                else
                    r_diff(iMouse1,iMouse2,ipair) = NaN;
                end
            end
        else
            r_diff(iMouse1,:,ipair) = NaN;
        end
    end
    r_diff_long(:,ipair) = triu2vec(r_diff(:,:,ipair));
end
r_diff_avg = nanmean(r_diff_long,1);
r_diff_sem = nanstd(r_diff_long,[],1)./sqrt(sum(isnan(r_diff_long),1));

r_diff_all = reshape(r_diff_long, [size(r_diff_long,1)*size(r_diff_long,2) 1]);
r_diff_all_avg = nanmean(r_diff_all,1);
r_diff_all_sem = nanstd(r_diff_all,[],1)./sqrt(sum(isnan(r_diff_all),1));

[h_PM_AL_all p_PM_AL_all] = ttest2(r_same_long(:,1), r_diff_long(:,1));
[h_AL_PM_all p_AL_PM_all] = ttest2(r_same_long(:,2), r_diff_long(:,1));

[h_all p_all] = ttest2(r_same_all, r_diff_all);

fn_summary = fullfile(sum_base, [matrix '_' num2str(P) 'P_' inj '_axon_summary'], 'avg_tuning_correlations.mat');
save(fn_summary, 'r_same_long', 'r_same_all', 'r_same_all_avg', 'r_same_all_sem', 'r_diff_long', 'r_diff_all','r_diff_all_avg', 'r_diff_all_sem', 'p_all');