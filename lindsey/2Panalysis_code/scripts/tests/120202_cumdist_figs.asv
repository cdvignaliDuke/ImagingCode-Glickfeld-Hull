figure;
for iArea= 1:3;
    subplot(2,2,1);
    [H, stats, xCDF, yCDF] = cdfplot_LG(log2(Goodfits(iArea).TF));
    plot(xCDF, yCDF, col(iArea,:));
    hold on
    title('TF')
    xlabel('log2(TF)');
    xlim([log2(1) log2(15)])
    subplot(2,2,2);
    [H, stats, xCDF, yCDF] = cdfplot_LG(log2(Goodfits(iArea).SF));
    plot(xCDF, yCDF, col(iArea,:));
    hold on
    xlabel('log2(SF)');
    title('SF')
    xlim([log2(.02) log2(.32)])
    subplot(2,2,3);
    [H, stats, xCDF, yCDF] = cdfplot_LG(Goodfits(iArea).dF);
    plot(xCDF, yCDF, col(iArea,:));
    hold on
    xlabel('dF/F');
    title('dF/F');
    subplot(2,2,4);
    [H, stats, xCDF, yCDF] = cdfplot_LG(log2(Goodfits(iArea).speed));
    plot(xCDF, yCDF, col(iArea,:));
    hold on
    xlabel('log2(speed)');
    title('Speed')
    xlim([log2(1/.32) log2(15/.02)])
end
suptitle([matrix '  ' num2str(P) 'P  ' inj 'axons-  ' areas(1,:) '(' num2str(size(Goodfits(1).plotfit,1)) ')  ' areas(2,:) '(' num2str(size(Goodfits(2).plotfit,1)) ')  ' areas(3,:) '(' num2str(size(Goodfits(3).plotfit,1)) ')'])
fn_out = fullfile(anal_base, [matrix '_' num2str(P) 'P_' inj '_allarea_log_TF_SF_dF_speed_cumdist.pdf']);
         print(gcf, '-dpdf', fn_out);
figure;
for iArea= 1:3;
    subplot(2,2,1);
    [H, stats, xCDF, yCDF] = cdfplot_LG(Goodfits(iArea).sigma_TF);
    plot(xCDF, yCDF, col(iArea,:));
    hold on
    xlim([0 3])
    xlabel('sigma TF');
    title(' ')
    subplot(2,2,2);
    [H, stats, xCDF, yCDF] = cdfplot_LG(Goodfits(iArea).sigma_SF);
    plot(xCDF, yCDF, col(iArea,:));
    hold on
    xlabel('sigma SF');
    title(' ')
    xlim([0 3])
    subplot(2,2,3);
    [H, stats, xCDF, yCDF] = cdfplot_LG(Goodfits(iArea).xi);
    plot(xCDF, yCDF, col(iArea,:));
    hold on
    xlabel('xi');
    title(' ')
    xlim([-2 2])
end
suptitle([matrix '  ' num2str(P) 'P  ' inj 'axons-  ' areas(1,:) '(' num2str(size(Goodfits(1).plotfit,1)) ')  ' areas(2,:) '(' num2str(size(Goodfits(2).plotfit,1)) ')  ' areas(3,:) '(' num2str(size(Goodfits(3).plotfit,1)) ')'])
fn_out = fullfile(anal_base, [matrix '_' num2str(P) 'P_' inj '_allarea_sigmas_xi_cumdist.pdf']);
        print(gcf, '-dpdf', fn_out);
 
        
 %distributions
 figure; 
 for iArea = 1:3
    subplot(2,2,1);
    [counts, bins] = hist(log2(Goodfits(iArea).TF), 100);
    plot(bins, counts./size(Goodfits(iArea).TF,1), col(iArea,:))
    hold on
    title('TF')
    xlabel('log2(TF)');
    xlim([log2(1) log2(15)])
    subplot(2,2,2);
    [counts, bins] = hist(log2(Goodfits(iArea).SF), 100);
    plot(bins, counts./size(Goodfits(iArea).TF,1), col(iArea,:))
    hold on
    xlabel('log2(SF)');
    title('SF')
    xlim([log2(.02) log2(.32)])
    subplot(2,2,3);
    [counts, bins] = hist(Goodfits(iArea).dF, 100);
    plot(bins, counts./size(Goodfits(iArea).TF,1), col(iArea,:))
    hold on
    xlabel('dF/F');
    title('dF/F');
    subplot(2,2,4);
    [counts, bins] = hist(log2(Goodfits(iArea).speed), 100);
    plot(bins, counts./size(Goodfits(iArea).TF,1), col(iArea,:))
    hold on
    xlabel('log2(speed)');
    title('Speed')
    xlim([log2(1/.32) log2(15/.02)])
 end
 %histograms       
 vec0 = zeros(2,7);
 vec0(1,:) = [.5 1 2 4 8 16 32];
 vec0(2,:) = [.01 .02 .04 .08 .16 .32 .64];    
 vec1 = interp2(vec0');
 TF_vec = zeros(8,1);
 SF_vec = zeros(8,1);
 TF_vec(1,1)= vec1(2,1);
 TF_vec(2,1)= 1.0001;
 TF_vec(3:6,1) = vec1(4:2:10,1);
 TF_vec(7,1) = 14.999;
 TF_vec(8,1) = vec1(12,1);
 SF_vec(1,1)= vec1(2,3);
 SF_vec(2,1)= .020001;
 SF_vec(3:6,1) = vec1(4:2:10,3);
 SF_vec(7,1) = .3199;
 SF_vec(8,1) = vec1(12,3);
 for iArea = 1:3
 figure; 
 subplot(2,2,1)
 scatter(Goodfits(iArea).TF,Goodfits(iArea).SF, 2);
 hold on
 text(7.5, .35, num2str(length(Goodfits(iArea).TF)));
 subplot(2,2,2)
 edges = log2(TF_vec);
 h= histc(log2(Goodfits(iArea).TF),edges);
 bar(1:8, h);
 xlabel('TF')
 subplot(2,2,3)
 edges = log2(SF_vec);
 h= histc(log2(Goodfits(iArea).SF),edges);
 bar(1:8, h);
 xlabel('SF')
 subplot(2,2,4)
 edges = log2(TF_vec./fliplr(SF_vec));
 h= histc(log2(Goodfits(iArea).speed),edges);
 bar(1:8, h);
 xlabel('Speed')
 suptitle(['All Goodfits ' areas(iArea,:)])
 fn_out = fullfile(anal_base, [matrix '_' num2str(P) 'P_' inj areas(iArea,:), '_all_goodfits_log.pdf']);
        print(gcf, '-dpdf', fn_out);
 end
 
 figure;
 start = 1;
 for iArea = 1:3
 subplot(3,3,start)
  edges = log2(TF_vec);
 h= histc(log2(Goodfits(iArea).TF),edges);
 bar(1:8, h);
  xlim([0 8])
  subplot(3,3,start+1)
  edges = log2(SF_vec);
  h= histc(log2(Goodfits(iArea).SF),edges);
 bar(1:8, h);
  xlim([0 8])
  subplot(3,3,start+2)
  edges = log2(TF_vec./flipud(SF_vec));
  h= histc(log2(Goodfits(iArea).speed),edges);
 bar(1:8, h);
 xlim([0 8])
 start= start+3;
 end
subplot(3,3,1)
title('TF');
subplot(3,3,2)
title('SF');
subplot(3,3,3)
title('Speed');

fn_out = fullfile(anal_base, [matrix '_' num2str(P) 'P_' inj '_all_areas_log_speed_SF_TF_highlowpass.pdf']);
        print(gcf, '-dpdf', fn_out);


