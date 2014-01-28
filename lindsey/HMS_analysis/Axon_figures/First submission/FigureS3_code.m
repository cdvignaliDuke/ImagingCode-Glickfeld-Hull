fig_base = '\\zmey\storlab\users\Lindsey\Projects\HVAs\Manuscript\Figures_2012';
fig = 'S3';
areas = strvcat('PM', 'LM', 'AL');
area_order = [2;3;1];
col = strvcat('c', 'k', 'r');
P = 2;
matrix = 'SF5xTF5';
inj = 'V1';
sum_base = 'G:\users\lindsey\analysisLG\experiments';

fn_summary = fullfile(sum_base, [matrix '_' num2str(P) 'P_' inj '_axon_summary'], 'all_fits.mat');
load(fn_summary);
fn_good = fullfile(sum_base, [matrix '_' num2str(P) 'P_' inj '_axon_summary'], 'good_fits.mat');
load(fn_good);


%example mouse
mouse ='AC45';
date_list = strvcat('110822', '110820', '110823');
userun = 1:4;
expt = [8; 4; 6];
dirs = 2;
base = 'G:\users\lindsey\analysisLG\active mice';

x = 1:5;
y = 5:-1:1;
[x_grid y_grid] = meshgrid(x,y);
x_grid_long = reshape(x_grid', [25 1]);
y_grid_long = reshape(y_grid', [25 1]);

figure;
begin = 1;
for iArea = 1:3
    date = date_list(iArea,:);
    outDir = fullfile(base, mouse,date);
    fn_df =  fullfile(outDir,'analysis',[date '_' mouse '_run' num2str(userun) '_stack_dF_all.tif']);
    stack = readtiff(fn_df);
    fn_mask =  fullfile(outDir,'analysis',[date '_' mouse '_run' num2str(userun) '_ttest_mask.mat']);
    load(fn_mask);
    local_max = [];
    fn_local =  fullfile(outDir,'analysis',[date '_' mouse '_run' num2str(userun) '_local_max.mat']);
    load(fn_local);
    siz = size(stack);
    dF_nodir = zeros(siz(1), siz(2), 25);
    start = 1;
    for iCond = 1:25
        dF_nodir(:,:,iCond) = mean(stack(:,:,start:start+1),3);
        start = start+2;
    end
    dF_avg = squeeze(mean(mean(dF_nodir,1),2));
    dF_norm= (dF_avg./max(dF_avg,[],1)*20);
    dF_nodir_long = reshape(dF_nodir, [siz(1)*siz(2) 25]);
    dF_max = max(dF_nodir,[],3);
    mask = find(ttest_mask == 1);
    dF_mask = mean(dF_nodir_long(mask,:),1)';
    dF_mask_norm = (dF_mask./max(dF_mask,[],1)*20);
    subplot(6,6,begin+3)
    imagesq(dF_max);
    colormap(gray)
    caxis([0 .5])
    subplot(6,6,begin)
    for iCond = 1:25
        scatter(x_grid_long(iCond,1),y_grid_long(iCond,1), dF_norm(iCond,:).^2,'.k')
        hold on
        axis square
        axis off
        xlim([0 6])
        ylim([0 6])
    end
    title(num2str(max(dF_avg,[],1)));
    subplot(6,6,begin+9)
    imagesq(ttest_mask);
    colormap(gray)
    subplot(6,6,begin+6)
    for iCond = 1:25
        scatter(x_grid_long(iCond,1),y_grid_long(iCond,1), dF_mask_norm(iCond,:).^2,'.k')
        hold on
        axis square
        axis off
        xlim([0 6])
        ylim([0 6])
    end
    title(num2str(max(dF_mask,[],1)));
    iexp = expt(iArea);
    plotfit = [];
    all_data = [];
    good_data = [];
    n = all_fits(iArea).expt(iexp).n(1);
    for iCell = 1:n
        all_data = cat(3, all_data, reshape(all_fits(iArea).expt(iexp).bouton(iCell).dFoverF, [5 5]));
        if all_fits(iArea).expt(iexp).bouton(iCell).goodfit ==1;
            good_data = cat(3, good_data, reshape(all_fits(iArea).expt(iexp).bouton(iCell).dFoverF, [5 5]));
            plotfit = cat(3, plotfit, all_fits(iArea).expt(iexp).bouton(iCell).plotfit);
        end
    end
    avg_data_all = mean(all_data,3);
    avg_data_long = reshape((avg_data_all./max(max(avg_data_all,[],1),[],2))'*20,[25 1]);
    avg_data_good = mean(all_data,3);
    avg_data_good_long = reshape((avg_data_good./max(max(avg_data_good,[],1),[],2))'*20,[25 1]);
    avg_tuning = mean(plotfit,3);
    avg_tuning_long = reshape((avg_tuning./max(max(avg_tuning,[],1),[],2))'*20,[25 1]);
    subplot(6,6,begin+12)
    for iCond = 1:25
        scatter(x_grid_long(iCond,1),y_grid_long(iCond,1), avg_data_long(iCond,:).^2,'.k')
        hold on
        axis square
        axis off
        xlim([0 6])
        ylim([0 6])
    end
    title(num2str(max(max(avg_data_all,[],1),[],2)));
    subplot(6,6,begin+18)
    for iCond = 1:25
        scatter(x_grid_long(iCond,1),y_grid_long(iCond,1), avg_data_good_long(iCond,:).^2,'.k')
        hold on
        axis square
        axis off
        xlim([0 6])
        ylim([0 6])
    end
    title(num2str(max(max(avg_data_good,[],1),[],2)));
    bouton_pix_all = zeros(size(local_max));
    bouton_pix_good = zeros(size(local_max));
    n = all_fits(iArea).expt(iexp).n(1);
    for iCell = 1:n
        pos = all_fits(iArea).expt(iexp).bouton(iCell).pos;
        suby = pos(1)-1:pos(1)+1;
        subx = pos(2)-1:pos(2)+1;
        bouton_pix_all(suby,subx) = 1;
        if all_fits(iArea).expt(iexp).bouton(iCell).goodfit == 1
            bouton_pix_good(suby,subx) = 1;
        end
    end
    subplot(6,6,begin+15)
    imagesq(bouton_pix_all);
    colormap(gray)
    subplot(6,6,begin+24)
    for iCond = 1:25
        scatter(x_grid_long(iCond,1),y_grid_long(iCond,1), avg_tuning_long(iCond,:).^2,'.k')
        hold on
        axis square
        axis off
        xlim([0 6])
        ylim([0 6])
    end
    title(num2str(max(max(avg_tuning,[],1),[],2)));
    subplot(6,6,begin+27)
    imagesq(bouton_pix_good);
    begin = begin+1;
end
fn_out = fullfile(fig_base, ['Figure' fig], [matrix '_' num2str(P) 'P_' inj '_all_areas_' mouse '_tuning_FOV_pix_boutons.ps']);
        print(gcf, '-depsc2', fn_out);
        
%for all FOVs
figure;    
begin =1;
for iArea = 1:3
    image = areas(iArea,:);
    sum_base = 'G:\users\lindsey\analysisLG\experiments';
    list_fn = fullfile(sum_base, [matrix '_' num2str(P) 'P_' image], [matrix '_' num2str(P) 'P_' inj '_' image '_exp_list.mat']);
    load(list_fn);
    nexp = size(exp_list.mouse_mat,2);
    stack_dF_avg_all = zeros(25,nexp);
    stack_dF_mask_avg_all = zeros(25,nexp);
    plotfit = [];
    good_data = [];
    all_data = [];
    for iexp = 1:nexp;   
        mouse = char(exp_list.mouse_mat{iexp});
        date = char(exp_list.date_mat{iexp});
        userun = exp_list.runs_mat{iexp};
        count_protocol = exp_list.prot_mat{iexp};
        run = exp_list.run_mat{iexp};
        blanks = exp_list.blanks_mat{iexp};
        dirs = exp_list.dir_mat{iexp};

        base = 'G:\users\lindsey\analysisLG\active mice';    
        outDir = fullfile(base, mouse,date);
        
        fn_out = fullfile(outDir, 'analysis', [date '_' mouse '_run' num2str(userun) '_stack_dF_all.tif']);
        stack_dF = readtiff(fn_out);
        fn_ttest= fullfile(outDir,'analysis',[date '_' mouse '_run' num2str(userun) '_ttest_mask.mat']);
        load(fn_ttest);
        
        siz = size(stack_dF);
        stack_dF_long = reshape(stack_dF, [siz(1)*siz(2) siz(3)]);
        stack_dF_avg = mean(stack_dF_long,1)';
        ind = find(ttest_mask==1);
        stack_dF_mask_avg = mean(stack_dF_long(ind,:),1)';
        
        if dirs ==1
            nCond = 25;
        elseif dirs ==2
            nCond = 50;
        end

        stack_dF_avg_temp = zeros(25,1);
        stack_dF_mask_avg_temp = zeros(25,1);
        if nCond == 50;
            start = 1;
            for iCond = 1:25
                stack_dF_avg_temp(iCond,:) = mean(stack_dF_avg(start:start+1,:),1);
                stack_dF_mask_avg_temp(iCond,:) = mean(stack_dF_avg(start:start+1,:),1);
                start = start+2;
            end
            stack_dF_avg = stack_dF_avg_temp;
            stack_dF_mask_avg = stack_dF_mask_avg_temp;
        end        
        stack_dF_avg_all(:,iexp) =stack_dF_avg(1:25,:);
        stack_dF_mask_avg_all(:,iexp) =stack_dF_mask_avg(1:25,:);
        
        n = all_fits(iArea).expt(iexp).n(1);
        for iCell = 1:n
            all_data = cat(3, all_data, reshape(all_fits(iArea).expt(iexp).bouton(iCell).dFoverF, [5 5]));
            if all_fits(iArea).expt(iexp).bouton(iCell).goodfit == 1
                plotfit = cat(3,plotfit,all_fits(iArea).expt(iexp).bouton(iCell).plotfit);
               good_data = cat(3, good_data, reshape(all_fits(iArea).expt(iexp).bouton(iCell).dFoverF, [5 5]));
            end
        end        
    end
    stack_dF = nanmean(stack_dF_avg_all,2)';
    stack_dF_norm = (stack_dF./max(stack_dF,[],2))*20;
    subplot(6,6,begin)
    for iCond = 1:25
        scatter(x_grid_long(iCond,1),y_grid_long(iCond,1), stack_dF_norm(:,iCond).^2,'.k')
        hold on
        axis square
        axis off
        xlim([0 6])
        ylim([0 6])
    end
    title(num2str(max(stack_dF,[],2)));
    stack_mask = nanmean(stack_dF_mask_avg_all,2)';
    stack_mask_norm = (stack_mask./max(stack_mask,[],2))*20;
    subplot(6,6,begin+6)
    for iCond = 1:25
        scatter(x_grid_long(iCond,1),y_grid_long(iCond,1), stack_mask_norm(:,iCond).^2,'.k')
        hold on
        axis square
        axis off
        xlim([0 6])
        ylim([0 6])
    end
    title(num2str(max(stack_mask,[],2)));
    all_data_avg = mean(all_data,3);
    all_data_long = reshape(all_data_avg',[25 1]);
    all_data_norm = (all_data_long./max(all_data_long,[],1))*20;
    subplot(6,6,begin+12)
    for iCond = 1:25
        scatter(x_grid_long(iCond,1),y_grid_long(iCond,1), all_data_norm(iCond,:).^2,'.k')
        hold on
        axis square
        axis off
        xlim([0 6])
        ylim([0 6])
    end
    title(num2str(max(all_data_long,[],1)));
    good_data_avg = mean(good_data,3);
    good_data_long = reshape(good_data_avg',[25 1]);
    good_data_norm = (good_data_long./max(good_data_long,[],1))*20;
    subplot(6,6,begin+18)
    for iCond = 1:25
        scatter(x_grid_long(iCond,1),y_grid_long(iCond,1), good_data_norm(iCond,:).^2,'.k')
        hold on
        axis square
        axis off
        xlim([0 6])
        ylim([0 6])
    end
    title(num2str(max(good_data_long,[],1)));
    plotfit_avg = mean(plotfit,3);
    plotfit_long = reshape(plotfit_avg',[25 1]);
    plotfit_norm = (plotfit_long./max(plotfit_long,[],1))*20;
    subplot(6,6,begin+24)
    for iCond = 1:25
        scatter(x_grid_long(iCond,1),y_grid_long(iCond,1), plotfit_norm(iCond,:).^2,'.k')
        hold on
        axis square
        axis off
        xlim([0 6])
        ylim([0 6])
    end
    title(num2str(max(plotfit_long,[],1)));
    begin = begin+1;
end
       
fn_out = fullfile(fig_base, ['Figure' fig], [matrix '_' num2str(P) 'P_' inj '_all_areas_all_FOV_tuning.ps']);
        print(gcf, '-depsc2', fn_out);

%speed tuning binned by speed
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
speed_tuning_all = zeros(6,length(uspeeds),2,3);
speed_tuning_mask = zeros(6,length(uspeeds),2,3);
speed_tuning_all_data = zeros(6,length(uspeeds),2,3);
speed_tuning_good_data = zeros(6,length(uspeeds),2,3);
speed_tuning_fit = zeros(6,length(uspeeds),2,3);
edges_dF= [0    0.2000    0.4000    0.8000    1.6000    5.0000];
n_all_dF = zeros(6,3);
n_mask_dF = zeros(6,3);
n_all_data_dF = zeros(6,3);
n_good_data_dF = zeros(6,3);
n_fit_dF = zeros(6,3);
h_speed_all = zeros(1,3);
p_speed_all = zeros(1,3);
h_speed_mask = zeros(1,3);
p_speed_mask = zeros(1,3);
h_speed_alldata = zeros(1,3);
p_speed_alldata = zeros(1,3);
h_speed_gooddata = zeros(1,3);
p_speed_gooddata = zeros(1,3);
for iArea = [1 3]
    image = areas(iArea,:);
    sum_base = 'G:\users\lindsey\analysisLG\experiments';
    list_fn = fullfile(sum_base, [matrix '_' num2str(P) 'P_' image], [matrix '_' num2str(P) 'P_' inj '_' image '_exp_list.mat']);
    load(list_fn);
    nexp = size(exp_list.mouse_mat,2);
    stack_all = [];
    stack_mask_all = [];
    all_data = [];
    good_data = [];
    for iexp= 1:nexp 
        mouse = char(exp_list.mouse_mat{iexp});
        date = char(exp_list.date_mat{iexp});
        userun = exp_list.runs_mat{iexp};
        count_protocol = exp_list.prot_mat{iexp};
        run = exp_list.run_mat{iexp};
        blanks = exp_list.blanks_mat{iexp};
        dirs = exp_list.dir_mat{iexp};

        base = 'G:\users\lindsey\analysisLG\active mice';    
        outDir = fullfile(base, mouse,date);
        
        fn_out = fullfile(outDir, 'analysis', [date '_' mouse '_run' num2str(userun) '_stack_dF_all.tif']);
        stack_dF = readtiff(fn_out);
        fn_ttest= fullfile(outDir,'analysis',[date '_' mouse '_run' num2str(userun) '_ttest_mask.mat']);
        load(fn_ttest);
        
        siz = size(stack_dF);
        ind = find(ttest_mask==1);
        if dirs ==1
            nCond = 25;
        elseif dirs ==2
            nCond = 50;
        end
        stack_dF_temp = zeros(siz(1),siz(2),25);
        if nCond == 50;
            start = 1;
            for iCond = 1:25
                stack_dF_temp(:,:,iCond) = mean(stack_dF(:,:,start:start+1),3);
                start = start+2;
            end
            stack_dF = stack_dF_temp;
        end        
        stack_long = reshape(stack_dF(:,:,1:25),[siz(1)*siz(2), 25]);
        stack_all = [stack_all; stack_long];
        stack_mask = stack_long(ind, :);
        stack_mask_all = [stack_mask_all; stack_mask];
        
        n = all_fits(iArea).expt(iexp).n(1);
        for iCell = 1:n
            all_data = [all_data; reshape(reshape(all_fits(iArea).expt(iexp).bouton(iCell).dFoverF, [5 5])',[1 25])];
            if all_fits(iArea).expt(iexp).bouton(iCell).goodfit == 1
               good_data = [good_data; reshape(reshape(all_fits(iArea).expt(iexp).bouton(iCell).dFoverF, [5 5])',[1 25])];
            end
        end    
    end
                    
    %division into df bins
    max_dF_all = max(stack_all,[],2);
    [n_all_dF(:,iArea) bin_all_dF] = histc(max_dF_all,edges_dF);
    max_dF_mask = max(stack_mask_all,[],2);
    [n_mask_dF(:,iArea) bin_mask_dF] = histc(max_dF_mask,edges_dF);
    max_dF_all_data = max(all_data,[],2);
    [n_all_data_dF(:,iArea) bin_all_data_dF] = histc(max_dF_all_data,edges_dF);
    max_dF_good_data = max(good_data,[],2);
    [n_good_data_dF(:,iArea) bin_good_data_dF] = histc(max_dF_good_data,edges_dF);
    [n_fit_dF(:,iArea) bin_fit_dF] = histc(Goodfits(iArea).dF, edges_dF);
    ind_all = find(bin_all_dF>0);
    %ttests
    [h_speed_all(:,iArea), p_speed_all(:,iArea)] = ttest(mean(stack_all(ind_all,[1 2 3 6 7 11]),2), mean(stack_all(ind_all,[15 19 20 23 24 25]),2));
    [h_speed_mask(:,iArea), p_speed_mask(:,iArea)] = ttest(mean(stack_mask_all(:,[1 2 3 6 7 11]),2), mean(stack_mask_all(:,[15 19 20 23 24 25]),2));
    [h_speed_alldata(:,iArea), p_speed_alldata(:,iArea)] = ttest(mean(all_data(:,[1 2 3 6 7 11]),2), mean(all_data(:,[15 19 20 23 24 25]),2));
    [h_speed_gooddata(:,iArea), p_speed_gooddata(:,iArea)] = ttest(mean(good_data(:,[1 2 3 6 7 11]),2), mean(good_data(:,[15 19 20 23 24 25]),2));

    for ibin = 1:length(n_dF)
        if ibin<6
            ind_all_dF = find(bin_all_dF == ibin);
            ind_mask_dF = find(bin_mask_dF == ibin);
            ind_all_data_dF = find(bin_all_data_dF == ibin);
            ind_good_data_dF = find(bin_good_data_dF == ibin);
            ind_fit_dF = find(bin_fit_dF == ibin);
        end
        for ispeed = 1:length(uspeeds)
            ind_speed = find(speed_grid_long == uspeeds(ispeed));
            if ibin<6
                if length(ind_speed)>1
                    all_ispeed = mean(stack_all(ind_all_dF,ind_speed),2);
                    mask_ispeed = mean(stack_mask_all(ind_mask_dF,ind_speed),2);
                    alldata_ispeed = mean(all_data(ind_all_data_dF,ind_speed),2);
                    gooddata_ispeed = mean(good_data(ind_good_data_dF,ind_speed),2);
                    goodfits_ispeed = mean(Goodfits(iArea).plotfit(ind_fit_dF,ind_speed),2);
                elseif length(ind_speed) == 1
                    all_ispeed = stack_all(ind_all_dF,ind_speed);
                    mask_ispeed = stack_mask_all(ind_mask_dF,ind_speed);
                    alldata_ispeed = all_data(ind_all_data_dF,ind_speed);
                    gooddata_ispeed = good_data(ind_good_data_dF,ind_speed);
                    goodfits_ispeed = Goodfits(iArea).plotfit(ind_fit_dF,ind_speed);
                end
            elseif ibin == 6
                if length(ind_speed)>1
                    all_ispeed = nanmean(stack_all(ind_all,ind_speed),2);
                    mask_ispeed = mean(stack_mask_all(:,ind_speed),2);
                    alldata_ispeed = mean(all_data(:,ind_speed),2);
                    gooddata_ispeed = mean(good_data(:,ind_speed),2);
                    goodfits_ispeed = mean(Goodfits(iArea).plotfit(:,ind_speed),2);
                elseif length(ind_speed) == 1
                    all_ispeed = stack_all(ind_all,ind_speed);
                    mask_ispeed = stack_mask_all(:,ind_speed);
                    alldata_ispeed = all_data(:,ind_speed);
                    gooddata_ispeed = good_data(:,ind_speed);
                    goodfits_ispeed = Goodfits(iArea).plotfit(:,ind_speed);
                end
            end
            speed_tuning_all(ibin,ispeed,1,iArea) = nanmean(all_ispeed,1);
            speed_tuning_all(ibin,ispeed,2,iArea) = nanstd(all_ispeed, [],1)./sqrt(size(all_ispeed,1));
            speed_tuning_mask(ibin,ispeed,1,iArea) = mean(mask_ispeed,1);
            speed_tuning_mask(ibin,ispeed,2,iArea) = std(mask_ispeed, [],1)./sqrt(size(mask_ispeed,1));
            speed_tuning_all_data(ibin,ispeed,1,iArea) = mean(alldata_ispeed,1);
            speed_tuning_all_data(ibin,ispeed,2,iArea) = std(alldata_ispeed, [],1)./sqrt(size(alldata_ispeed,1));
            speed_tuning_good_data(ibin,ispeed,1,iArea) = mean(gooddata_ispeed,1);
            speed_tuning_good_data(ibin,ispeed,2,iArea) = std(gooddata_ispeed, [],1)./sqrt(size(gooddata_ispeed,1));
            speed_tuning_fit(ibin,ispeed,1,iArea) = mean(goodfits_ispeed,1);
            speed_tuning_fit(ibin,ispeed,2,iArea) = std(goodfits_ispeed, [],1)./sqrt(size(goodfits_ispeed,1));
        end
    end
end

figure;
for iArea  = [1 3]
    begin = 0;
    for ibin = 1:5;
        subplot(5,5,ibin+begin)
        errorbar(log2(uspeeds), speed_tuning_all(ibin,:,1,iArea),speed_tuning_all(ibin,:,2,iArea),col(iArea,:));
        hold on
        xlim([1 10])
        axis square
        subplot(5,5,ibin+1+begin)
        errorbar(log2(uspeeds), speed_tuning_mask(ibin,:,1,iArea),speed_tuning_mask(ibin,:,2,iArea),col(iArea,:));
        hold on
        xlim([1 10])
        axis square
        subplot(5,5,ibin+2+begin)
        errorbar(log2(uspeeds), speed_tuning_all_data(ibin,:,1,iArea),speed_tuning_all_data(ibin,:,2,iArea),col(iArea,:));
        hold on
        xlim([1 10])
        axis square
        subplot(5,5,ibin+3+begin)
        errorbar(log2(uspeeds), speed_tuning_good_data(ibin,:,1,iArea),speed_tuning_good_data(ibin,:,2,iArea),col(iArea,:));
        hold on
        xlim([1 10])
        axis square
        subplot(5,5,ibin+4+begin)
        errorbar(log2(uspeeds), speed_tuning_fit(ibin,:,1,iArea),speed_tuning_fit(ibin,:,2,iArea),col(iArea,:));
        hold on
        xlim([1 10])
        axis square
        begin = begin+4;
    end
end
for iplot = 16:20
subplot(5,5,iplot)
ylim([0 .75])
end
fn_out = fullfile(fig_base, ['Figure' num2str(fig)], [matrix '_' num2str(P) 'P_' inj '_speed_tuning_by_dF.ps']);
        print(gcf, '-depsc2', fn_out);

figure;
for iArea  = [1 3]
    for ibin = 6;
        subplot(1,5,1)
        errorbar(log2(uspeeds), speed_tuning_all(ibin,:,1,iArea),speed_tuning_all(ibin,:,2,iArea),col(iArea,:));
        hold on
        xlim([1 10])
        axis square
        subplot(1,5,2)
        errorbar(log2(uspeeds), speed_tuning_mask(ibin,:,1,iArea),speed_tuning_mask(ibin,:,2,iArea),col(iArea,:));
        hold on
        xlim([1 10])
        axis square
        subplot(1,5,3)
        errorbar(log2(uspeeds), speed_tuning_all_data(ibin,:,1,iArea),speed_tuning_all_data(ibin,:,2,iArea),col(iArea,:));
        hold on
        xlim([1 10])
        axis square
        subplot(1,5,4)
        errorbar(log2(uspeeds), speed_tuning_good_data(ibin,:,1,iArea),speed_tuning_good_data(ibin,:,2,iArea),col(iArea,:));
        hold on
        xlim([1 10])
        axis square
        subplot(1,5,5)
        errorbar(log2(uspeeds), speed_tuning_fit(ibin,:,1,iArea),speed_tuning_fit(ibin,:,2,iArea),col(iArea,:));
        hold on
        xlim([1 10])
        axis square
    end
end
fn_out = fullfile(fig_base, ['Figure' num2str(fig)], [matrix '_' num2str(P) 'P_' inj '_speed_tuning_all_dF.ps']);
        print(gcf, '-depsc2', fn_out);

for iplot = 1
subplot(1,5,iplot)
ylim([0 .1])
end