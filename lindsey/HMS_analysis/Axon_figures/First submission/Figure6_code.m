fig_base = '\\zmey\storlab\users\Lindsey\Projects\HVAs\Manuscript\Figures_2012';
fig = 6;

P = 2;
matrix = 'SF5xTF5';
inj = 'LM';
sum_base = 'G:\users\lindsey\analysisLG\experiments';

fn_summary = fullfile(sum_base, [matrix '_' num2str(P) 'P_' inj '_axon_summary'], 'all_fits.mat');
load(fn_summary);
fn_good = fullfile(sum_base, [matrix '_' num2str(P) 'P_' inj '_axon_summary'], 'good_fits.mat');
load(fn_good);

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
    subplot(2,3,4)
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
fn_out = fullfile(fig_base, ['Figure' num2str(fig)], [matrix '_' num2str(P) 'P_' inj '_all_areas_cdf_medians.ps']);
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


        
 