fig_base = '\\zmey\storlab\users\Lindsey\Projects\HVAs\Manuscript\Resubmission_Figures';
fig = 5;

areas = strvcat('PM', 'LM', 'AL');
col = strvcat('c', 'k', 'r');
P = 2;
matrix = 'SF5xTF5';
inj = 'V1L5';
sum_base = 'G:\users\lindsey\analysisLG\experiments';

fn_summary = fullfile(sum_base, [matrix '_' num2str(P) 'P_' inj '_axon_summary'], 'all_fits.mat');
load(fn_summary);
fn_good = fullfile(sum_base, [matrix '_' num2str(P) 'P_' inj '_axon_summary'], 'good_fits.mat');
load(fn_good);

figure;
inj_list = {'V1' 'V1L5'};
matrix = 'SF5xTF5';
P = 2;

area_order = [1 2; 2 3; 1 3];
for iMouse = 1:length(mouse_list);
    mouse = mouse_list{iMouse};
    for ipair = 1:3
        medians = [NaN NaN];
        for iArea = 1:2
            for iinj = 1:2
                inj = char(inj_list(iinj));
                if iinj == 1
                    mouse_list = {'Y13' 'X32' 'DR7' 'DR9' 'AC39' 'AC42' 'AC44' 'AC45' 'Y18' 'Y26' 'M13' 'M14' 'M22' 'M31'};
                elseif iinj == 2
                    mouse_list = {'VC1' 'CM94'};
                end
                fn_summary = fullfile(sum_base, [matrix '_' num2str(P) 'P_' inj '_axon_summary'], 'all_fits.mat');
                load(fn_summary);
                fn_good = fullfile(sum_base, [matrix '_' num2str(P) 'P_' inj '_axon_summary'], 'good_fits.mat');
                load(fn_good);
                image = areas(area_order(ipair, iArea),:);
                sum_base = 'G:\users\lindsey\analysisLG\experiments';
                list_fn = fullfile(sum_base, [matrix '_' num2str(P) 'P_' image], [matrix '_' num2str(P) 'P_' inj '_' image '_exp_list.mat']);
                load(list_fn);
                median_speed = [];
                nexp = all_fits(area_order(ipair,iArea)).nexp;
                for iexp = 1:nexp
                    exp_mouse = exp_list.mouse_mat{iexp};
                    if length(exp_mouse)==length(mouse)
                        if exp_mouse == mouse;
                            if all_fits(area_order(ipair,iArea)).expt(iexp).n(2)>25
                                median_speed = [median_speed all_fits(area_order(ipair,iArea)).expt(iexp).median_speed];
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
                subplot(2,3,ipair+3)
                if iinj == 1
                    semilogy(1:2, medians, '-ob')
                elseif iinj == 2
                    semilogy(1:2, medians, '-^b')
                end
                axis square
                hold on
                title([areas(area_order(ipair,1),:) ' vs ' areas(area_order(ipair,2),:)])
                xlim([0 3])
                ylim([1 1000])
            end
        end    
    end
end

fn_out = fullfile(fig_base, ['Figure' num2str(fig)], [matrix '_' num2str(P) 'P_' inj '_all_mice_medians.ps']);
        print(gcf, '-depsc2', fn_out);

figure;        
for iArea = 1:3
    inj = 'V1';
    fn_good = fullfile(sum_base, [matrix '_' num2str(P) 'P_' inj '_axon_summary'], 'good_fits.mat');
    load(fn_good);
    Goodfits_WT = Goodfits(iArea).speed;
    inj = 'V1L5';
    fn_good = fullfile(sum_base, [matrix '_' num2str(P) 'P_' inj '_axon_summary'], 'good_fits.mat');
    load(fn_good);
    Goodfits_Rbp = Goodfits(iArea).speed;
    subplot(2,3,iArea)
    semilogy(iArea, median([Goodfits_WT; Goodfits_Rbp]), 'ok');
    ylim([1 1000])
end

fn_out = fullfile(fig_base, ['Figure' num2str(fig)], [matrix '_' num2str(P) 'P_' inj '_all_areas_medians.ps']);
        print(gcf, '-depsc2', fn_out);
        
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
        
%cumulative distributions with WT + RBP4cre
inj_list = {'V1' 'V1L5'};
P = 2;
matrix = 'SF5xTF5';
figure;
title_build = [];
for iArea = 1:3
    TF_all = [];
    SF_all = [];
    Speed_all = [];
    for iinj = 1:2
        inj = char(inj_list(iinj));
        fn_good = fullfile(sum_base, [matrix '_' num2str(P) 'P_' inj '_axon_summary'], 'good_fits.mat');
        load(fn_good);
        TF_all = [TF_all; Goodfits(iArea).TF];
        SF_all = [SF_all; Goodfits(iArea).SF];
        Speed_all = [Speed_all; Goodfits(iArea).speed];
    end
    subplot(1,3,1);
    [H, stats, xCDF, yCDF] = cdfplot_LG(log2(TF_all));
    plot(xCDF, yCDF, col(iArea,:));
    hold on
    title('TF')
    axis square
    xlabel('log2(TF)');
    xlim([log2(1) log2(15)])
    subplot(1,3,2);
    [H, stats, xCDF, yCDF] = cdfplot_LG(log2(SF_all));
    plot(xCDF, yCDF, col(iArea,:));
    hold on
    xlabel('log2(SF)');
    axis square
    title('SF')
    xlim([log2(.02) log2(.32)])
    subplot(1,3,3);
    [H, stats, xCDF, yCDF] = cdfplot_LG(log2(Speed_all));
    plot(xCDF, yCDF, col(iArea,:));
    hold on
    xlabel('log2(speed)');
    axis square
    title('Speed')
    xlim([log2(1/.32) log2(15/.02)])   
    title_build = [title_build ' ' areas(iArea,:) ' ' num2str(size(TF_all,1))];
end
suptitle(title_build);
fn_out = fullfile(fig_base, ['Figure' num2str(fig)], [matrix '_' num2str(P) 'P_all_areas_cdf_V1&V1L5.ps']);
        print(gcf, '-depsc2', fn_out);
    