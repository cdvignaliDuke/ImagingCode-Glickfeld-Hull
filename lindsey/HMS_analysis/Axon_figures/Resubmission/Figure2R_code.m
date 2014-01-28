fig_base = '\\zmey\storlab\users\Lindsey\Projects\HVAs\Manuscript\Resubmission_Figures';
fig = '2';

areas = strvcat('PM', 'LM', 'AL');
col = strvcat('c', 'k', 'r');
P = 2;
matrix = 'SF5xTF5';
inj = 'V1';
sum_base = 'G:\users\lindsey\analysisLG\experiments';

%individual mouse tuning
mouse_list = {'AC45' 'M14' 'Y26'};
fn_summary = fullfile(sum_base, [matrix '_' num2str(P) 'P_' inj '_axon_summary'], 'all_fits.mat');
load(fn_summary);
figure;
for iMouse = 1:3
    mouse = char(mouse_list(iMouse));
    for iArea = 1:3
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


fn_summary = fullfile(sum_base, [matrix '_' num2str(P) 'P_' inj '_axon_summary'], 'all_fits.mat');
load(fn_summary);
fn_good = fullfile(sum_base, [matrix '_' num2str(P) 'P_' inj '_axon_summary'], 'good_fits.mat');
load(fn_good);
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
mouse_list = {'Y13' 'X32' 'DR7' 'DR9' 'AC39' 'AC42' 'AC44' 'AC45' 'Y18' 'Y26' 'M13' 'M14' 'M22' 'M31'};
figure;
for iArea = 1:3
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

for iArea = 1:3
    avg_tuning(iArea).all_tuning_long = [];
    for iMouse = 1:length(mouse_list);
        if avg_tuning(iArea).mouse(iMouse).number>25;
        tuning_long = reshape(avg_tuning(iArea).mouse(iMouse).tuning_norm', [25 1]);
        avg_tuning(iArea).all_tuning_long = [avg_tuning(iArea).all_tuning_long tuning_long];
        end
    end
end
for iArea = 1:3
    tuning_fit_all(iArea).TF_tuning = [];
    tuning_fit_all(iArea).SF_tuning = [];
    tuning_fit_all(iArea).speed_tuning = [];
    for iTF = 1:length(TF_vec0)
        ind = find(TF_grid_long == TF_vec0(iTF));
        if length(ind)>1
            goodfits_iTF = mean(avg_tuning(iArea).all_tuning_long(ind,:),1);
        elseif length(ind) == 1
            goodfits_iTF = avg_tuning(iArea).all_tuning_long(ind,:);
        end
        tuning_fit_all(iArea).TF_tuning = [tuning_fit_all(iArea).TF_tuning; goodfits_iTF];
    end
    for iSF = 1:length(SF_vec0)
        ind = find(SF_grid_long == SF_vec0(iSF));
        if length(ind)>1
            goodfits_iSF = mean(avg_tuning(iArea).all_tuning_long(ind,:),1);
        elseif length(ind) == 1
            goodfits_iSF = avg_tuning(iArea).all_tuning_long(ind,:);
        end
        tuning_fit_all(iArea).SF_tuning = [tuning_fit_all(iArea).SF_tuning; goodfits_iSF];
    end
    for ispeed = 1:length(uspeeds)
        ind = find(speed_grid_long == uspeeds(ispeed));
        if length(ind)>1
            goodfits_ispeed = mean(avg_tuning(iArea).all_tuning_long(ind,:),1);
        elseif length(ind) == 1
            goodfits_ispeed = avg_tuning(iArea).all_tuning_long(ind,:);
        end
        tuning_fit_all(iArea).speed_tuning = [tuning_fit_all(iArea).speed_tuning; goodfits_ispeed];
    end
end
for iArea = 1:3
    tuning_fit_all(iArea).TF_tuning_high = mean(tuning_fit_all(iArea).TF_tuning(4:5,:),1);
    tuning_fit_all(iArea).TF_tuning_low= mean(tuning_fit_all(iArea).TF_tuning(1:2,:),1);
    tuning_fit_all(iArea).TF_tuning_ratio= tuning_fit_all(iArea).TF_tuning_high./tuning_fit_all(iArea).TF_tuning_low;
    tuning_fit_all(iArea).TF_tuning_metric= (tuning_fit_all(iArea).TF_tuning_high-tuning_fit_all(iArea).TF_tuning_low)./(tuning_fit_all(iArea).TF_tuning_high+tuning_fit_all(iArea).TF_tuning_low);
    tuning_fit_all(iArea).SF_tuning_high = mean(tuning_fit_all(iArea).SF_tuning(4:5,:),1);
    tuning_fit_all(iArea).SF_tuning_low= mean(tuning_fit_all(iArea).SF_tuning(1:2,:),1);
    tuning_fit_all(iArea).SF_tuning_ratio= tuning_fit_all(iArea).SF_tuning_high./tuning_fit_all(iArea).SF_tuning_low;
    tuning_fit_all(iArea).SF_tuning_metric= (tuning_fit_all(iArea).SF_tuning_high-tuning_fit_all(iArea).SF_tuning_low)./(tuning_fit_all(iArea).SF_tuning_high+tuning_fit_all(iArea).SF_tuning_low);
    tuning_fit_all(iArea).speed_tuning_high = mean(tuning_fit_all(iArea).speed_tuning(7:9,:),1);
    tuning_fit_all(iArea).speed_tuning_low= mean(tuning_fit_all(iArea).speed_tuning(1:3,:),1);
    tuning_fit_all(iArea).speed_tuning_ratio= tuning_fit_all(iArea).speed_tuning_high./tuning_fit_all(iArea).speed_tuning_low;
    tuning_fit_all(iArea).speed_tuning_metric= (tuning_fit_all(iArea).speed_tuning_high-tuning_fit_all(iArea).speed_tuning_low)./(tuning_fit_all(iArea).speed_tuning_high+tuning_fit_all(iArea).speed_tuning_low);
end

TF_tuning_ratio_avg = zeros(2,3);
SF_tuning_ratio_avg = zeros(2,3);
speed_tuning_ratio_avg = zeros(2,3);

for iArea = 1:3
    TF_tuning_ratio_avg(1,iArea) = mean(tuning_fit_all(iArea).TF_tuning_ratio,2);
    TF_tuning_ratio_avg(2,iArea) = std(tuning_fit_all(iArea).TF_tuning_ratio,[],2)./sqrt(size(tuning_fit_all(iArea).TF_tuning_ratio,2));
    SF_tuning_ratio_avg(1,iArea) = mean(tuning_fit_all(iArea).SF_tuning_ratio,2);
    SF_tuning_ratio_avg(2,iArea) = std(tuning_fit_all(iArea).SF_tuning_ratio,[],2)./sqrt(size(tuning_fit_all(iArea).SF_tuning_ratio,2));
    speed_tuning_ratio_avg(1,iArea) = mean(tuning_fit_all(iArea).speed_tuning_ratio,2);
    speed_tuning_ratio_avg(2,iArea) = std(tuning_fit_all(iArea).speed_tuning_ratio,[],2)./sqrt(size(tuning_fit_all(iArea).speed_tuning_ratio,2));
end
area_order = [1 2; 2 3; 1 3];
h_TF_SF_speed_ratio = zeros(3,3);
p_TF_SF_speed_ratio = zeros(3,3);
h_TF_SF_speed_metric = zeros(3,3);
p_TF_SF_speed_metric = zeros(3,3);
for ipair = 1:3
    [h_TF_SF_speed_ratio(1,ipair) p_TF_SF_speed_ratio(1,ipair)] = ttest2(tuning_fit_all(area_order(ipair, 1)).TF_tuning_ratio, tuning_fit_all(area_order(ipair, 2)).TF_tuning_ratio);
    [h_TF_SF_speed_ratio(2,ipair) p_TF_SF_speed_ratio(2,ipair)] = ttest2(tuning_fit_all(area_order(ipair, 1)).SF_tuning_ratio, tuning_fit_all(area_order(ipair, 2)).SF_tuning_ratio);
    [h_TF_SF_speed_ratio(3,ipair) p_TF_SF_speed_ratio(3,ipair)] = ttest2(tuning_fit_all(area_order(ipair, 1)).speed_tuning_ratio, tuning_fit_all(area_order(ipair, 2)).speed_tuning_ratio);
    [h_TF_SF_speed_metric(1,ipair) p_TF_SF_speed_metric(1,ipair)] = ttest2(tuning_fit_all(area_order(ipair, 1)).TF_tuning_metric, tuning_fit_all(area_order(ipair, 2)).TF_tuning_metric);
    [h_TF_SF_speed_metric(2,ipair) p_TF_SF_speed_metric(2,ipair)] = ttest2(tuning_fit_all(area_order(ipair, 1)).SF_tuning_metric, tuning_fit_all(area_order(ipair, 2)).SF_tuning_metric);
    [h_TF_SF_speed_metric(3,ipair) p_TF_SF_speed_metric(3,ipair)] = ttest2(tuning_fit_all(area_order(ipair, 1)).speed_tuning_metric, tuning_fit_all(area_order(ipair, 2)).speed_tuning_metric);
end

%tuning correlations
r_same = zeros(length(mouse_list), length(mouse_list), 3);
for iArea = 1:3
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
r_same_sem = nanstd(r_same_long,[],1)./sqrt(sum(~isnan(r_same_long),1));

r_pair = zeros(length(mouse_list), 3);
area_order = [1 2; 2 3; 1 3];
for ipair = 1:3
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
r_pair_sem = nanstd(r_pair,[],1)./sqrt(sum(~isnan(r_pair),1));

[h_PM_LM p_PM_LM] = ttest2(r_same_long(:,1), r_pair(:,1));
[h_PM_AL p_PM_AL] = ttest2(r_same_long(:,1), r_pair(:,3));
[h_AL_LM p_AL_LM] = ttest2(r_same_long(:,3), r_pair(:,2));
[h_AL_PM p_AL_PM] = ttest2(r_same_long(:,3), r_pair(:,3));
[h_LM_PM p_LM_PM] = ttest2(r_same_long(:,2), r_pair(:,1));
[h_LM_AL p_LM_AL] = ttest2(r_same_long(:,2), r_pair(:,2));

r_diff = zeros(length(mouse_list), length(mouse_list), 3);
area_order = [1 2; 2 3; 1 3];
for ipair = 1:3
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
r_diff_sem = nanstd(r_diff_long,[],1)./sqrt(sum(~isnan(r_diff_long),1));

[h_PM_LM_all p_PM_LM_all] = ttest2(r_same_long(:,1), r_diff_long(:,1));
[h_PM_AL_all p_PM_AL_all] = ttest2(r_same_long(:,1), r_diff_long(:,3));
[h_AL_LM_all p_AL_LM_all] = ttest2(r_same_long(:,3), r_diff_long(:,2));
[h_AL_PM_all p_AL_PM_all] = ttest2(r_same_long(:,3), r_diff_long(:,3));
[h_LM_PM_all p_LM_PM_all] = ttest2(r_same_long(:,2), r_diff_long(:,1));
[h_LM_AL_all p_LM_AL_all] = ttest2(r_same_long(:,2), r_diff_long(:,2));

fn_summary = fullfile(sum_base, [matrix '_' num2str(P) 'P_' inj '_axon_summary'], 'avg_tuning_correlations.mat');
save(fn_summary, 'r_same_long', 'r_same_avg', 'r_same_sem', 'r_diff_long', 'r_diff_avg', 'r_diff_sem', 'p_PM_LM_all', 'p_PM_AL_all', 'p_LM_PM_all', 'p_LM_AL_all','p_AL_LM_all', 'p_AL_PM_all');
