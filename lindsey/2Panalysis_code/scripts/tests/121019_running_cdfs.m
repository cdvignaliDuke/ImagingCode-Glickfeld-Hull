clear all
areas = ['PM'; 'LM'; 'AL'; 'RL'; 'AM'];
inj = 'V1';
P = 2;
all_goodfits = [];
for iArea = 1:5
    matrix = 'SF5xTF5';
    image = areas(iArea,:);
    sum_base = 'G:\users\lindsey\analysisLG\experiments';
    list_fn = fullfile(sum_base, [matrix '_' num2str(P) 'P_' image], [matrix '_' num2str(P) 'P_' inj '_' image '_exp_list.mat']);
    load(list_fn);
    nexp = size(exp_list.mouse_mat,2);
    all_goodfits(iArea).fits = [];
%     all_goodfits(iArea).plotfit_run = [];
%     all_goodfits(iArea).plotfit_norun = [];
%     all_goodfits(iArea).data_run = [];
%     all_goodfits(iArea).data_norun = [];
%     run_mat(iArea).ratio = [];
%     run_mat(iArea).good = [];
    run_ratio(iArea).mat = [];
    for iexp = 1:nexp
        mouse = char(exp_list.mouse_mat{iexp});
        date = char(exp_list.date_mat{iexp});
        userun = exp_list.runs_mat{iexp};
        count_protocol = exp_list.prot_mat{iexp};
        run = exp_list.run_mat{iexp};
        blanks = exp_list.blanks_mat{iexp};
        dirs = exp_list.dir_mat{iexp};
        
        if run == 1
            if dirs ==1
                nCond = 25;
            elseif dirs ==2
                nCond = 50;
            end

            base = 'G:\users\lindsey\analysisLG\active mice';    
            outDir = fullfile(base, mouse,date);
            
            fn_out = fullfile(outDir, 'analysis', [date '_' mouse '_run' num2str(userun) '_wheel.mat']);   
            load(fn_out);
            
            x = length(find(Wheel_sorted ==0));
            y = length(find(Wheel_sorted >0 & Wheel_sorted<=0.1));
            z = length(find(Wheel_sorted >0.1));
            
            run_ratio(iArea).mat = [run_ratio(iArea).mat; [x y z]];
        end
    end
end
            
            if x >49
                if y >49
                    run_mat(iArea).good = [run_mat(iArea).good iexp];
                    
                    fn_out = fullfile(outDir, 'analysis', [date '_' mouse '_run' num2str(userun) '_lbub_fits_run_norun.mat']);   
                    load(fn_out);
                    fn_out = fullfile(outDir, 'analysis', [date '_' mouse '_run' num2str(userun) '_Fit_struct_run.mat']);   
                    load(fn_out);

                    goodfits = intersect(goodfit_ind.run, goodfit_ind.norun);
                    all_goodfits(iArea).fits = cat(1,all_goodfits(iArea).fits, lbub_fits(goodfits,:,:,:));
                    for iCell = goodfits
                        all_goodfits(iArea).plotfit_run = cat(3,all_goodfits(iArea).plotfit_run, Fit_struct_run(iCell).True_run.s_.k2b_plot);
                        all_goodfits(iArea).plotfit_norun = cat(3,all_goodfits(iArea).plotfit_norun, Fit_struct_run(iCell).True_norun.s_.k2b_plot);
                        all_goodfits(iArea).data_run = cat(3,all_goodfits(iArea).data_run, Fit_struct_run(iCell).True_run.s_.data);
                        all_goodfits(iArea).data_norun = cat(3,all_goodfits(iArea).data_norun, Fit_struct_run(iCell).True_norun.s_.data);
                    end
                end
            end
        end
    end
end

figure;
for iArea = 1:5
    subplot(3,5,iArea)
    imagesq(mean(all_goodfits(iArea).plotfit_run,3));
    title([areas(iArea,:) ' ' num2str(size(all_goodfits(iArea).plotfit_run,3))])
    caxis([0 .4])
    text(2,6.5,num2str(max(max(mean(all_goodfits(iArea).plotfit_run,3),[],1),[],2)))
    subplot(3,5,iArea+5)
    imagesq(mean(all_goodfits(iArea).plotfit_norun,3));
    colormap(jet)
    caxis([0 .4])
    text(2,6.5,num2str(max(max(mean(all_goodfits(iArea).plotfit_norun,3),[],1),[],2)))
end

for iArea= 1:5;
    figure;
    subplot(2,2,1);
    [H, stats, xCDF, yCDF] = cdfplot_LG(squeeze(all_goodfits(iArea).fits(:,1,4,3)));
    plot(xCDF, yCDF, 'k');
    hold on
    [H, stats, xCDF, yCDF] = cdfplot_LG(squeeze(all_goodfits(iArea).fits(:,1,4,2)));
    plot(xCDF, yCDF, 'r');
    title('Amplitude')
    axis square
    xlabel('dF/F');
    xlim([0 3])
    subplot(2,2,2);
    [H, stats, xCDF, yCDF] = cdfplot_LG(squeeze(all_goodfits(iArea).fits(:,5,4,3)));
    plot(xCDF, yCDF, 'k');
    hold on
    [H, stats, xCDF, yCDF] = cdfplot_LG(squeeze(all_goodfits(iArea).fits(:,5,4,2)));
    plot(xCDF, yCDF, 'r');
    title('TF')
    axis square
    xlabel('log2(TF)');
    xlim([log2(1) log2(15)])
    subplot(2,2,3);
    [H, stats, xCDF, yCDF] = cdfplot_LG(squeeze(all_goodfits(iArea).fits(:,4,4,3)));
    plot(xCDF, yCDF, 'k');
    hold on
    [H, stats, xCDF, yCDF] = cdfplot_LG(squeeze(all_goodfits(iArea).fits(:,4,4,2)));
    plot(xCDF, yCDF, 'r');
    xlabel('log2(SF)');
    axis square
    title('SF')
    xlim([log2(.02) log2(.32)])
    subplot(2,2,4);
    [H, stats, xCDF, yCDF] = cdfplot_LG(log2(squeeze(2.^(all_goodfits(iArea).fits(:,5,4,3))./2.^(all_goodfits(iArea).fits(:,4,4,3)))));
    plot(xCDF, yCDF, 'k');
    hold on
    [H, stats, xCDF, yCDF] = cdfplot_LG(log2(squeeze(2.^(all_goodfits(iArea).fits(:,5,4,2))./2.^(all_goodfits(iArea).fits(:,4,4,2)))));    
    plot(xCDF, yCDF, 'r');
    xlabel('log2(speed)');
    axis square
    title('Speed')
    xlim([log2(1/.32) log2(15/.02)])
    suptitle(areas(iArea,:))
end

for iArea= 1:5;
    figure;
    subplot(2,2,1);
    scatter(all_goodfits(iArea).fits(:,1,4,2),all_goodfits(iArea).fits(:,1,4,3),'.k');
    hold on
    plot(median(all_goodfits(iArea).fits(:,1,4,2),1), median(all_goodfits(iArea).fits(:,1,4,3),1), 'ob');
    hold on
    x = 1:.1:3;
    y = x;
    plot(x,y);
    xlim([0 3])
    ylim([0 3])
    title('Amplitude')
    axis square
    xlabel('run dF/F');
    ylabel('norun dF/F');
    subplot(2,2,2);
    scatter(all_goodfits(iArea).fits(:,5,4,2),all_goodfits(iArea).fits(:,5,4,3),'.k');
    hold on
    plot(median(all_goodfits(iArea).fits(:,5,4,2),1), median(all_goodfits(iArea).fits(:,5,4,3),1), 'ob');
    hold on
    x = log2(1):.1:log2(15);
    y = x;
    plot(x,y);
    title('TF')
    axis square
    xlabel('run Hz');
    ylabel('norun Hz');
    xlim([log2(1) log2(15)])
    ylim([log2(1) log2(15)])
    subplot(2,2,3);
    scatter(all_goodfits(iArea).fits(:,4,4,2),all_goodfits(iArea).fits(:,4,4,3),'.k');
    hold on
    plot(median(all_goodfits(iArea).fits(:,4,4,2),1), median(all_goodfits(iArea).fits(:,4,4,3),1), 'ob');
    hold on
    x = log2(.02):.1:log2(.32);
    y = x;
    plot(x,y);
    xlabel('run cpd');
    xlabel('norun cpd');
    axis square
    title('SF')
    xlim([log2(.02) log2(.32)])
    ylim([log2(.02) log2(.32)])
    subplot(2,2,4);
    scatter(log2(2.^all_goodfits(iArea).fits(:,5,4,2)./2.^all_goodfits(iArea).fits(:,4,4,2)),log2(2.^all_goodfits(iArea).fits(:,5,4,3)./2.^all_goodfits(iArea).fits(:,4,4,3)),'.k');
    hold on
    plot(median(log2(2.^all_goodfits(iArea).fits(:,5,4,2)./2.^all_goodfits(iArea).fits(:,4,4,2))), median(log2(2.^all_goodfits(iArea).fits(:,5,4,3)./2.^all_goodfits(iArea).fits(:,4,4,3))), 'ob');
    hold on
    x = log2(1/.32):.1:log2(15/.02);
    y = x;
    plot(x,y);
    xlabel('run speed');
    xlabel('norun speed');
    axis square
    title('Speed')
    xlim([log2(1/.32) log2(15/.02)])
    ylim([log2(1/.32) log2(15/.02)])
    suptitle(areas(iArea,:))
end

amp_change = [];
SF_change = [];
TF_change = [];
speed_change = [];
for iArea  = 1:5
    amp_change(iArea).all = squeeze((all_goodfits(iArea).fits(:,1,4,2)-all_goodfits(iArea).fits(:,1,4,3))./(all_goodfits(iArea).fits(:,1,4,3)+all_goodfits(iArea).fits(:,1,4,2)));
    amp_change(iArea).avg = mean((all_goodfits(iArea).fits(:,1,4,2)-all_goodfits(iArea).fits(:,1,4,3))./(all_goodfits(iArea).fits(:,1,4,3)+all_goodfits(iArea).fits(:,1,4,2)),1);
    amp_change(iArea).sem = std((all_goodfits(iArea).fits(:,1,4,2)-all_goodfits(iArea).fits(:,1,4,3))./(all_goodfits(iArea).fits(:,1,4,3)+all_goodfits(iArea).fits(:,1,4,2)),[],1)./sqrt(size(all_goodfits(iArea).fits,1));
    TF_change(iArea).all = squeeze((2.^all_goodfits(iArea).fits(:,5,4,2)-2.^all_goodfits(iArea).fits(:,5,4,3))./(2.^all_goodfits(iArea).fits(:,5,4,3)+2.^all_goodfits(iArea).fits(:,5,4,2)));
    TF_change(iArea).avg = mean((2.^all_goodfits(iArea).fits(:,5,4,2)-2.^all_goodfits(iArea).fits(:,5,4,3))./(2.^all_goodfits(iArea).fits(:,5,4,3)+2.^all_goodfits(iArea).fits(:,5,4,2)),1);
    TF_change(iArea).sem = std((2.^all_goodfits(iArea).fits(:,5,4,2)-2.^all_goodfits(iArea).fits(:,5,4,3))./(2.^all_goodfits(iArea).fits(:,5,4,3)+2.^all_goodfits(iArea).fits(:,5,4,2)),[],1)./sqrt(size(all_goodfits(iArea).fits,1));
    SF_change(iArea).all = squeeze((2.^all_goodfits(iArea).fits(:,4,4,2)-2.^all_goodfits(iArea).fits(:,4,4,3))./(2.^all_goodfits(iArea).fits(:,4,4,3)+2.^all_goodfits(iArea).fits(:,4,4,2)));
    SF_change(iArea).avg = mean((2.^all_goodfits(iArea).fits(:,4,4,2)-2.^all_goodfits(iArea).fits(:,4,4,3))./(2.^all_goodfits(iArea).fits(:,4,4,3)+2.^all_goodfits(iArea).fits(:,4,4,2)),1);
    SF_change(iArea).sem = std((2.^all_goodfits(iArea).fits(:,4,4,2)-2.^all_goodfits(iArea).fits(:,4,4,3))./(2.^all_goodfits(iArea).fits(:,4,4,3)+2.^all_goodfits(iArea).fits(:,4,4,2)),[],1)./sqrt(size(all_goodfits(iArea).fits,1));
    speed_change(iArea).all = squeeze(((2.^all_goodfits(iArea).fits(:,5,4,2)./2.^all_goodfits(iArea).fits(:,4,4,2))-(2.^all_goodfits(iArea).fits(:,5,4,3)./2.^all_goodfits(iArea).fits(:,4,4,3)))./(((2.^all_goodfits(iArea).fits(:,5,4,2)./2.^all_goodfits(iArea).fits(:,4,4,2))+(2.^all_goodfits(iArea).fits(:,5,4,3)./2.^all_goodfits(iArea).fits(:,4,4,3)))));
    speed_change(iArea).avg = mean(((2.^all_goodfits(iArea).fits(:,5,4,2)./2.^all_goodfits(iArea).fits(:,4,4,2))-(2.^all_goodfits(iArea).fits(:,5,4,3)./2.^all_goodfits(iArea).fits(:,4,4,3)))./(((2.^all_goodfits(iArea).fits(:,5,4,2)./2.^all_goodfits(iArea).fits(:,4,4,2))+(2.^all_goodfits(iArea).fits(:,5,4,3)./2.^all_goodfits(iArea).fits(:,4,4,3)))),1);
    speed_change(iArea).sem = std(((2.^all_goodfits(iArea).fits(:,5,4,2)./2.^all_goodfits(iArea).fits(:,4,4,2))-(2.^all_goodfits(iArea).fits(:,5,4,3)./2.^all_goodfits(iArea).fits(:,4,4,3)))./(((2.^all_goodfits(iArea).fits(:,5,4,2)./2.^all_goodfits(iArea).fits(:,4,4,2))+(2.^all_goodfits(iArea).fits(:,5,4,3)./2.^all_goodfits(iArea).fits(:,4,4,3)))),1)./sqrt(size(all_goodfits(iArea).fits,1));
end
figure;
for iArea = 1:5
    subplot(2,2,1)
    bar(iArea, amp_change(iArea).avg)
    hold on
    errorbar(iArea, amp_change(iArea).avg,amp_change(iArea).sem,'d')
    title('Amplitude')
    subplot(2,2,2)
    bar(iArea, TF_change(iArea).avg)
    hold on
    errorbar(iArea, TF_change(iArea).avg,TF_change(iArea).sem,'d')
    title('TF')
    subplot(2,2,3)
    bar(iArea, SF_change(iArea).avg)
    hold on
    errorbar(iArea, SF_change(iArea).avg,SF_change(iArea).sem,'d')
    title('SF')
    subplot(2,2,4)
    bar(iArea, speed_change(iArea).avg)
    hold on
    errorbar(iArea, speed_change(iArea).avg,speed_change(iArea).sem,'d')
    title('Speed')
end

for iArea = 1:5
[h_amp(1,iArea) p_amp(1,iArea)] = ttest(amp_change(iArea).all);
[h_TF(1,iArea) p_TF(1,iArea)] = ttest(TF_change(iArea).all);
[h_SF(1,iArea) p_SF(1,iArea)] = ttest(SF_change(iArea).all);
[h_speed(1,iArea) p_speed(1,iArea)] = ttest(speed_change(iArea).all);
end