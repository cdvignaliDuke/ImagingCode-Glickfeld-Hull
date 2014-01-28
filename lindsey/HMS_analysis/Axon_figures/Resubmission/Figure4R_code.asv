fig_base = '\\zmey\storlab\users\Lindsey\Projects\HVAs\Manuscript\Resubmission_Figures';
fig = '4';
areas = strvcat('PM', 'LM', 'AL');
area_order = [2;3;1];
col = strvcat('c', 'k', 'r');
P = 2;
matrix = 'SF5xTF5';
inj = 'V1L5';
sum_base = 'G:\users\lindsey\analysisLG\experiments';

fn_summary = fullfile(sum_base, [matrix '_' num2str(P) 'P_' inj '_axon_summary'], 'all_fits.mat');
load(fn_summary);
fn_good = fullfile(sum_base, [matrix '_' num2str(P) 'P_' inj '_axon_summary'], 'good_fits.mat');
load(fn_good);

x = 1:5;
y = 5:-1:1;
[x_grid y_grid] = meshgrid(x,y);
x_grid_long = reshape(x_grid', [25 1]);
y_grid_long = reshape(y_grid', [25 1]);

%individual mouse tuning
mouse_list = {'CM94' 'VC1'};
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

%average tuning of FOV
figure;
for iArea = 1:3
    nexp = all_fits(iArea).nexp;
    image = areas(iArea,:);
    sum_base = 'G:\users\lindsey\analysisLG\experiments';
    list_fn = fullfile(sum_base, [matrix '_' num2str(P) 'P_' image], [matrix '_' num2str(P) 'P_' inj '_' image '_exp_list.mat']);
    load(list_fn);
    x = 1:5;
    y = 5:-1:1;
    [x_grid y_grid] = meshgrid(x,y);
    x_grid_long = reshape(x_grid', [25 1]);
    y_grid_long = reshape(y_grid', [25 1]);
    plotfit_all = mean(Goodfits(iArea).plotfit,1)';
    plotfit_norm_all = (plotfit_all./max(plotfit_all,[],1))*55; 
    subplot(3,3,iArea)
    for iCond = 1:25
        scatter(x_grid_long(iCond,1),y_grid_long(iCond,1), plotfit_norm_all(iCond,:).^2,'.k')
        hold on
        axis square
        axis off
        xlim([0 6])
        ylim([0 6])
    end
    title([areas(iArea,:) ' ' num2str(size(Goodfits(iArea).plotfit,1))]);
end
fn_out = fullfile(fig_base, ['Figure' num2str(fig)], [matrix '_' num2str(P) 'P_' inj '_all_expts_avg_tuning.ps']);
    print(gcf, '-depsc2', fn_out);

%Tuning curves
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
    all_SFxTF_long(:,iArea) = reshape((all_SFxTF_fit(:,:,iArea)./max(max(all_SFxTF_fit(:,:,iArea),[],1),[],2))'*75,[25 1]);        
end
area_order = [2;3;1];
figure;
start = 1;
for iArea = 1:3
    subplot(2,3,start)
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
    subplot(2,3,4)
    errorbar(log2(TF_vec0), all_TF_fit(:,1,iArea),all_SF_fit(:,2,iArea),col(iArea,:));
    hold on
    title('TF')
    ylabel('dF/F')
    xlabel('log2(TF)')
    xlim([-1 5])
    ylim([0 .4])
    subplot(2,3,5)
    errorbar(log2(SF_vec0), all_SF_fit(:,1,iArea),all_SF_fit(:,2,iArea),col(iArea,:));
    hold on
    title('SF')
    xlabel('log2(SF)')
    ylim([0 .4])
    xlim([-6 -1])
    subplot(2,3,6)
    errorbar(log2(uspeeds), all_speed_fit(:,1,iArea),all_speed_fit(:,2,iArea),col(iArea,:));
    hold on
    title('speed')
    xlabel('log2(speed)')
    ylim([0 .5])
    xlim([1 10])
    start = start+1;
end
suptitle([matrix '  ' num2str(P) 'P  ' inj 'axons-  ' areas(1,:) '(' num2str(size(Goodfits(1).plotfit,1)) ')  ' areas(2,:) '(' num2str(size(Goodfits(2).plotfit,1)) ')  ' areas(3,:) '(' num2str(size(Goodfits(3).plotfit,1)) ')'])

fn_out = fullfile(fig_base, ['Figure' num2str(fig)], [matrix '_' num2str(P) 'P_' inj '_all_areas_tuning_curves.ps']);
        print(gcf, '-depsc', fn_out);

%cumulative distributions 
figure;
for iArea= 1:3;
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
end

fn_out = fullfile(fig_base, ['Figure' num2str(fig)], [matrix '_' num2str(P) 'P_' inj '_all_areas_cdfs.ps']);
        print(gcf, '-depsc', fn_out);
        
%add no shuffle data from AL
areas = ['AL'];
inj = 'V1L5';
P = 2;
for iArea = 1
    matrix = 'SF5xTF5';
    image = areas(iArea,:);
    sum_base = 'G:\users\lindsey\analysisLG\experiments';
    list_fn = fullfile(sum_base, [matrix '_' num2str(P) 'P_' image], [matrix '_' num2str(P) 'P_' inj '_' image '_exp_list.mat']);
    load(list_fn);
    nexp = size(exp_list.mouse_mat,2);  
    x = [];
    for iexp = 1:nexp
        mouse = char(exp_list.mouse_mat{iexp});
        date = char(exp_list.date_mat{iexp});
        userun = exp_list.runs_mat{iexp};
        count_protocol = exp_list.prot_mat{iexp};
        run = exp_list.run_mat{iexp};
        blanks = exp_list.blanks_mat{iexp};
        dirs = exp_list.dir_mat{iexp};

        if dirs ==1
            nCond = 25;
        elseif dirs ==2
            nCond = 50;
        end

        base = 'G:\users\lindsey\analysisLG\active mice';    
        outDir = fullfile(base, mouse,date);
        
        fn_out = fullfile(outDir, 'analysis', [date '_' mouse '_run' num2str(userun) '_Fit_struct_ttest_pt01_noshuf.mat']);   
        load(fn_out)
        fn_out = fullfile(outDir, 'analysis', [date '_' mouse '_run' num2str(userun) '_local_max_ttest_pt01.mat']);   
        load(fn_out)
        
        for iCell = 1:n_pix
            x = [x; Fit_struct(iCell).True.s_.x];
        end
    end
end

subplot(2,3,1)
[H, stats, xCDF, yCDF] = cdfplot_LG(TF');
plot(xCDF, yCDF, 'g');
hold on
subplot(2,3,2)
[H, stats, xCDF, yCDF] = cdfplot_LG(SF');
plot(xCDF, yCDF, 'g');
subplot(2,3,3)
[H, stats, xCDF, yCDF] = cdfplot_LG(log2(((2.^TF')./(2.^SF'))));
plot(xCDF, yCDF, 'g');

fn_out = fullfile(fig_base, ['Figure' num2str(fig)], [matrix '_' num2str(P) 'P_' inj '_all_areas_cdfs_lowerthresh_AL.ps']);
        print(gcf, '-depsc', fn_out);
