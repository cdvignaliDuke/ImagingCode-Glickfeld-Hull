%all OI for 2Hz and 8Hz
figure;
subplot(2,2,1)
col = strvcat('c', 'k', 'r');
x = [];
for iArea = 1:3
    [H, stats, xCDF, yCDF] = cdfplot_LG(all_fits_dir(iArea).x_2Hz_bestSF(:,3));
    plot(xCDF, yCDF, col(iArea,:));
    hold on
    x= [x; size(all_fits_dir(iArea).x_2Hz_bestSF(:,3),1)];
end
legend(num2str(x));
title('2Hz Best SF')
subplot(2,2,2)
col = strvcat('c', 'k', 'r');
x = [];
for iArea = 1:3
    [H, stats, xCDF, yCDF] = cdfplot_LG(all_fits_dir(iArea).x_8Hz_bestSF(:,3));
    plot(xCDF, yCDF, col(iArea,:));
    hold on
    x= [x; size(all_fits_dir(iArea).x_8Hz_bestSF(:,3),1)];
end
title('8Hz Best SF')
legend(num2str(x));
suptitle('Ori Index')
fn_out = fullfile('\\zmey\storlab\users\Lindsey\Projects\HVAs\2P Axon Imaging', ['DIR_' num2str(P) 'P_' inj], ['DIR_' num2str(P) 'P_' inj '_OI_2Hz_8Hz_bestSF.ps']);
print(gcf, '-depsc', fn_out);

%all OI for 0.04cpd and 0.16cpd
figure;
subplot(2,2,1)
col = strvcat('c', 'k', 'r');
x = [];
for iArea = 1:3
    [H, stats, xCDF, yCDF] = cdfplot_LG(all_fits_dir(iArea).x_pt04cpd_bestTF(:,3));
    plot(xCDF, yCDF, col(iArea,:));
    hold on
    x= [x; size(all_fits_dir(iArea).x_pt04cpd_bestTF(:,3),1)];
end
legend(num2str(x));
title('0.04cpd Best TF')
subplot(2,2,2)
col = strvcat('c', 'k', 'r');
x = [];
for iArea = 1:3
    [H, stats, xCDF, yCDF] = cdfplot_LG(all_fits_dir(iArea).x_pt16cpd_bestTF(:,3));
    plot(xCDF, yCDF, col(iArea,:));
    hold on
    x= [x; size(all_fits_dir(iArea).x_pt16cpd_bestTF(:,3),1)];
end
title('0.16cpd Best TF')
legend(num2str(x));
suptitle('Ori Index')
fn_out = fullfile('\\zmey\storlab\users\Lindsey\Projects\HVAs\2P Axon Imaging', ['DIR_' num2str(P) 'P_' inj], ['DIR_' num2str(P) 'P_' inj '_OI_pt04cpd_pt16cpd_bestTF.ps']);
print(gcf, '-depsc', fn_out);

%all DI for 2Hz and 8Hz
figure;
subplot(2,2,1)
col = strvcat('c', 'k', 'r');
x = [];
for iArea = 1:3
    [H, stats, xCDF, yCDF] = cdfplot_LG(all_fits_dir(iArea).x_orituned_2Hz_bestSF(:,8));
    plot(xCDF, yCDF, col(iArea,:));
    hold on
    x= [x; size(all_fits_dir(iArea).x_orituned_2Hz_bestSF(:,8),1)];
end
legend(num2str(x));
title('2Hz Best SF')
subplot(2,2,2)
col = strvcat('c', 'k', 'r');
x = [];
for iArea = 1:3
    [H, stats, xCDF, yCDF] = cdfplot_LG(all_fits_dir(iArea).x_orituned_8Hz_bestSF(:,8));
    plot(xCDF, yCDF, col(iArea,:));
    hold on
    x= [x; size(all_fits_dir(iArea).x_orituned_8Hz_bestSF(:,8),1)];
end
title('8Hz Best SF')
legend(num2str(x));
suptitle('Dir Index')
fn_out = fullfile('\\zmey\storlab\users\Lindsey\Projects\HVAs\2P Axon Imaging', ['DIR_' num2str(P) 'P_' inj], ['DIR_' num2str(P) 'P_' inj '_DI_2Hz_8Hz_bestSF.ps']);
print(gcf, '-depsc', fn_out);

%all DI for 0.04cpd and 0.16cpd
figure;
subplot(2,2,1)
col = strvcat('c', 'k', 'r');
x = [];
for iArea = 1:3
    [H, stats, xCDF, yCDF] = cdfplot_LG(all_fits_dir(iArea).x_orituned_pt04cpd_bestTF(:,8));
    plot(xCDF, yCDF, col(iArea,:));
    hold on
    x= [x; size(all_fits_dir(iArea).x_orituned_pt04cpd_bestTF(:,8),1)];
end
legend(num2str(x));
title('0.04cpd Best TF')
subplot(2,2,2)
col = strvcat('c', 'k', 'r');
x = [];
for iArea = 1:3
    [H, stats, xCDF, yCDF] = cdfplot_LG(all_fits_dir(iArea).x_orituned_pt16cpd_bestTF(:,8));
    plot(xCDF, yCDF, col(iArea,:));
    hold on
    x= [x; size(all_fits_dir(iArea).x_orituned_pt04cpd_bestTF(:,8),1)];
end
title('0.16cpd Best TF')
legend(num2str(x));
suptitle('Dir Index')
fn_out = fullfile('\\zmey\storlab\users\Lindsey\Projects\HVAs\2P Axon Imaging', ['DIR_' num2str(P) 'P_' inj], ['DIR_' num2str(P) 'P_' inj '_DI_pt04cpd_pt16cpd_bestTF.ps']);
print(gcf, '-depsc', fn_out);

%OI for 2Hz and 8Hz BestSFTF
figure;
subplot(2,2,1)
col = strvcat('c', 'k', 'r');
x = [];
for iArea = 1:3
    [H, stats, xCDF, yCDF] = cdfplot_LG([all_fits_dir(iArea).x_bestSFTF_2Hz_pt04cpd(:,3); all_fits_dir(iArea).x_bestSFTF_2Hz_pt16cpd(:,3)]);
    plot(xCDF, yCDF, col(iArea,:));
    hold on
    x= [x; size([all_fits_dir(iArea).x_bestSFTF_2Hz_pt04cpd(:,3); all_fits_dir(iArea).x_bestSFTF_2Hz_pt16cpd(:,3)],1)];
end
legend(num2str(x));
title('2Hz Best SFTF')
subplot(2,2,2)
col = strvcat('c', 'k', 'r');
x = [];
for iArea = 1:3
    [H, stats, xCDF, yCDF] = cdfplot_LG([all_fits_dir(iArea).x_bestSFTF_8Hz_pt04cpd(:,3); all_fits_dir(iArea).x_bestSFTF_8Hz_pt16cpd(:,3)]);
    plot(xCDF, yCDF, col(iArea,:));
    hold on
    x= [x; size([all_fits_dir(iArea).x_bestSFTF_8Hz_pt04cpd(:,3); all_fits_dir(iArea).x_bestSFTF_8Hz_pt16cpd(:,3)],1)];
end
title('8Hz Best SFTF')
legend(num2str(x));
suptitle('Ori Index')
fn_out = fullfile('\\zmey\storlab\users\Lindsey\Projects\HVAs\2P Axon Imaging', ['DIR_' num2str(P) 'P_' inj], ['DIR_' num2str(P) 'P_' inj '_OI_2Hz_8Hz_bestSFTF.ps']);
print(gcf, '-depsc', fn_out);

%OI for 0.04 and 0.16 cpd BestSFTF
figure;
subplot(2,2,1)
col = strvcat('c', 'k', 'r');
x = [];
for iArea = 1:3
    [H, stats, xCDF, yCDF] = cdfplot_LG([all_fits_dir(iArea).x_bestSFTF_2Hz_pt04cpd(:,3); all_fits_dir(iArea).x_bestSFTF_8Hz_pt04cpd(:,3)]);
    plot(xCDF, yCDF, col(iArea,:));
    hold on
    x= [x; size([all_fits_dir(iArea).x_bestSFTF_2Hz_pt04cpd(:,3); all_fits_dir(iArea).x_bestSFTF_8Hz_pt04cpd(:,3)],1)];
end
legend(num2str(x));
title('0.04cpd Best SFTF')
subplot(2,2,2)
col = strvcat('c', 'k', 'r');
x = [];
for iArea = 1:3
    [H, stats, xCDF, yCDF] = cdfplot_LG([all_fits_dir(iArea).x_bestSFTF_2Hz_pt16cpd(:,3); all_fits_dir(iArea).x_bestSFTF_8Hz_pt16cpd(:,3)]);
    plot(xCDF, yCDF, col(iArea,:));
    hold on
    x= [x; size([all_fits_dir(iArea).x_bestSFTF_2Hz_pt16cpd(:,3); all_fits_dir(iArea).x_bestSFTF_8Hz_pt16cpd(:,3)],1)];
end
title('0.16cpd Best SFTF')
legend(num2str(x));
suptitle('Ori Index')
fn_out = fullfile('\\zmey\storlab\users\Lindsey\Projects\HVAs\2P Axon Imaging', ['DIR_' num2str(P) 'P_' inj], ['DIR_' num2str(P) 'P_' inj '_OI_pt04cpd_pt16cpd_bestSFTF.ps']);
print(gcf, '-depsc', fn_out);

%DI for 2Hz and 8Hz BestSFTF
figure;
subplot(2,2,1)
col = strvcat('c', 'k', 'r');
x = [];
for iArea = 1:3
    [H, stats, xCDF, yCDF] = cdfplot_LG([all_fits_dir(iArea).x_bestSFTF_2Hz_pt04cpd(:,8); all_fits_dir(iArea).x_bestSFTF_2Hz_pt16cpd(:,8)]);
    plot(xCDF, yCDF, col(iArea,:));
    hold on
    x= [x; size([all_fits_dir(iArea).x_bestSFTF_2Hz_pt04cpd(:,8); all_fits_dir(iArea).x_bestSFTF_2Hz_pt16cpd(:,8)],1)];
end
legend(num2str(x));
title('2Hz Best SFTF')
subplot(2,2,2)
col = strvcat('c', 'k', 'r');
x = [];
for iArea = 1:3
    [H, stats, xCDF, yCDF] = cdfplot_LG([all_fits_dir(iArea).x_bestSFTF_8Hz_pt04cpd(:,8); all_fits_dir(iArea).x_bestSFTF_8Hz_pt16cpd(:,8)]);
    plot(xCDF, yCDF, col(iArea,:));
    hold on
    x= [x; size([all_fits_dir(iArea).x_bestSFTF_8Hz_pt04cpd(:,8); all_fits_dir(iArea).x_bestSFTF_8Hz_pt16cpd(:,8)],1)];
end
title('8Hz Best SFTF')
legend(num2str(x));
suptitle('Dir Index')
fn_out = fullfile('\\zmey\storlab\users\Lindsey\Projects\HVAs\2P Axon Imaging', ['DIR_' num2str(P) 'P_' inj], ['DIR_' num2str(P) 'P_' inj '_DI_2Hz_8Hz_bestSFTF.ps']);
print(gcf, '-depsc', fn_out);

%DI for 0.04 and 0.16 cpd BestSFTF
figure;
subplot(2,2,1)
col = strvcat('c', 'k', 'r');
x = [];
for iArea = 1:3
    [H, stats, xCDF, yCDF] = cdfplot_LG([all_fits_dir(iArea).x_bestSFTF_2Hz_pt04cpd(:,8); all_fits_dir(iArea).x_bestSFTF_8Hz_pt04cpd(:,8)]);
    plot(xCDF, yCDF, col(iArea,:));
    hold on
    x= [x; size([all_fits_dir(iArea).x_bestSFTF_2Hz_pt04cpd(:,8); all_fits_dir(iArea).x_bestSFTF_8Hz_pt04cpd(:,8)],1)];
end
legend(num2str(x));
title('0.04cpd Best SFTF')
subplot(2,2,2)
col = strvcat('c', 'k', 'r');
x = [];
for iArea = 1:3
    [H, stats, xCDF, yCDF] = cdfplot_LG([all_fits_dir(iArea).x_bestSFTF_2Hz_pt16cpd(:,8); all_fits_dir(iArea).x_bestSFTF_8Hz_pt16cpd(:,8)]);
    plot(xCDF, yCDF, col(iArea,:));
    hold on
    x= [x; size([all_fits_dir(iArea).x_bestSFTF_2Hz_pt16cpd(:,8); all_fits_dir(iArea).x_bestSFTF_8Hz_pt16cpd(:,8)],1)];
end
title('0.16cpd Best SFTF')
legend(num2str(x));
suptitle('Dir Index')
fn_out = fullfile('\\zmey\storlab\users\Lindsey\Projects\HVAs\2P Axon Imaging', ['DIR_' num2str(P) 'P_' inj], ['DIR_' num2str(P) 'P_' inj '_DI_pt04cpd_pt16cpd_bestSFTF.ps']);
print(gcf, '-depsc', fn_out);

%Mark's analysis
figure;
subplot(2,2,1)
col = strvcat('c', 'k', 'r');
x = [];
for iArea = 1:2
    [H, stats, xCDF, yCDF] = cdfplot_LG(all_fits_dir(iArea).x_2Hz_bestSF(:,3));
    plot(xCDF, yCDF, col(iArea,:));
    hold on
    x= [x; size(all_fits_dir(iArea).x_2Hz_bestSF(:,3),1)];
end
for iArea = 3
    [H, stats, xCDF, yCDF] = cdfplot_LG(all_fits_dir(iArea).x_8Hz_bestSF(:,3));
    plot(xCDF, yCDF, col(iArea,:));
    hold on
    x= [x; size(all_fits_dir(iArea).x_8Hz_bestSF(:,3),1)];
end
legend(num2str(x));
title('Ori Index')
subplot(2,2,2)
col = strvcat('c', 'k', 'r');
x = [];
for iArea = 1:2
    [H, stats, xCDF, yCDF] = cdfplot_LG(all_fits_dir(iArea).x_2Hz_bestSF(:,8));
    plot(xCDF, yCDF, col(iArea,:));
    hold on
    x= [x; size(all_fits_dir(iArea).x_2Hz_bestSF(:,3),1)];
end
for iArea = 3
    [H, stats, xCDF, yCDF] = cdfplot_LG(all_fits_dir(iArea).x_8Hz_bestSF(:,8));
    plot(xCDF, yCDF, col(iArea,:));
    hold on
    x= [x; size(all_fits_dir(iArea).x_8Hz_bestSF(:,3),1)];
end
title('Dir Index')
legend(num2str(x));
suptitle('2Hz in LM/PM 8Hz in AL')
fn_out = fullfile('\\zmey\storlab\users\Lindsey\Projects\HVAs\2P Axon Imaging', ['DIR_' num2str(P) 'P_' inj], ['DIR_' num2str(P) 'P_' inj '_OI_DI_2Hz_LMPM_8Hz_AL.ps']);
print(gcf, '-depsc', fn_out);