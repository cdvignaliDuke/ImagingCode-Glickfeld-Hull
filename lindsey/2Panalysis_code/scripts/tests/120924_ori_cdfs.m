
%all OI and DI for bestSFTF
figure;
subplot(2,2,1)
col = strvcat('c', 'k', 'r');
x = [];
for iArea = 1:3
    [H, stats, xCDF, yCDF] = cdfplot_LG(all_fits_dir(iArea).x_bestSFTF(:,3));
    plot(xCDF, yCDF, col(iArea,:));
    hold on
    x= [x; size(all_fits_dir(iArea).x_bestSFTF(:,3),1)];
end
title('Ori Index')
legend(num2str(x));
subplot(2,2,2)
col = strvcat('c', 'k', 'r');
x = [];
for iArea = 1:3
    [H, stats, xCDF, yCDF] = cdfplot_LG(all_fits_dir(iArea).x_bestSFTF_orituned(:,8));
    plot(xCDF, yCDF, col(iArea,:));
    hold on
    x= [x; size(all_fits_dir(iArea).x_bestSFTF_orituned(:,3),1)];
end
title('Dir Index')
legend(num2str(x));
suptitle('BestSFTF')

fn_out = fullfile('\\zmey\storlab\users\Lindsey\Projects\HVAs\2P Axon Imaging', ['DIR_' num2str(P) 'P_' inj], ['DIR_' num2str(P) 'P_' inj '_OI_DI_bestSFTF.ps']);
print(gcf, '-depsc', fn_out);


%all OI for each SFTF
figure;
subplot(2,2,1)
col = strvcat('c', 'k', 'r');
x = [];
for iArea = 1:3
    [H, stats, xCDF, yCDF] = cdfplot_LG(all_fits_dir(iArea).x_2Hz_pt04cpd(:,3));
    plot(xCDF, yCDF, col(iArea,:));
    hold on
    x= [x; size(all_fits_dir(iArea).x_2Hz_pt04cpd(:,3),1)];
end
legend(num2str(x));
title('2Hz 0.04 cpd')
subplot(2,2,2)
col = strvcat('c', 'k', 'r');
x = [];
for iArea = 1:3
    [H, stats, xCDF, yCDF] = cdfplot_LG(all_fits_dir(iArea).x_2Hz_pt16cpd(:,3));
    plot(xCDF, yCDF, col(iArea,:));
    hold on
    x= [x; size(all_fits_dir(iArea).x_2Hz_pt16cpd(:,3),1)];
end
title('2Hz 0.16 cpd')
legend(num2str(x));
subplot(2,2,3)
col = strvcat('c', 'k', 'r');
x = [];
for iArea = 1:3
    [H, stats, xCDF, yCDF] = cdfplot_LG(all_fits_dir(iArea).x_8Hz_pt04cpd(:,3));
    plot(xCDF, yCDF, col(iArea,:));
    hold on
    x= [x; size(all_fits_dir(iArea).x_8Hz_pt04cpd(:,3),1)];
end
title('8Hz 0.04 cpd')
legend(num2str(x));
subplot(2,2,4)
x = [];
col = strvcat('c', 'k', 'r');
for iArea = 1:3
    [H, stats, xCDF, yCDF] = cdfplot_LG(all_fits_dir(iArea).x_8Hz_pt16cpd(:,3));
    plot(xCDF, yCDF, col(iArea,:));
    hold on
    x= [x; size(all_fits_dir(iArea).x_8Hz_pt16cpd(:,3),1)];
end
title('8Hz 0.16 cpd')
legend(num2str(x));
suptitle('Ori Index')
fn_out = fullfile('\\zmey\storlab\users\Lindsey\Projects\HVAs\2P Axon Imaging', ['DIR_' num2str(P) 'P_' inj], ['DIR_' num2str(P) 'P_' inj '_OI_iSFTF.ps']);
print(gcf, '-depsc', fn_out);

%all DI for each SFTF
figure;
subplot(2,2,1)
col = strvcat('c', 'k', 'r');
x = [];
for iArea = 1:3
    [H, stats, xCDF, yCDF] = cdfplot_LG(all_fits_dir(iArea).x_orituned_2Hz_pt04cpd(:,8));
    plot(xCDF, yCDF, col(iArea,:));
    hold on
    x= [x; size(all_fits_dir(iArea).x_orituned_2Hz_pt04cpd(:,3),1)];
end
title('2Hz 0.04 cpd')
legend(num2str(x));
subplot(2,2,2)
col = strvcat('c', 'k', 'r');
x = [];
for iArea = 1:3
    [H, stats, xCDF, yCDF] = cdfplot_LG(all_fits_dir(iArea).x_orituned_2Hz_pt16cpd(:,8));
    plot(xCDF, yCDF, col(iArea,:));
    hold on
    x= [x; size(all_fits_dir(iArea).x_orituned_2Hz_pt16cpd(:,3),1)];
end
title('2Hz 0.16 cpd')
legend(num2str(x));
subplot(2,2,3)
col = strvcat('c', 'k', 'r');
x = [];
for iArea = 1:3
    [H, stats, xCDF, yCDF] = cdfplot_LG(all_fits_dir(iArea).x_orituned_8Hz_pt04cpd(:,8));
    plot(xCDF, yCDF, col(iArea,:));
    hold on
    x= [x; size(all_fits_dir(iArea).x_orituned_8Hz_pt04cpd(:,3),1)];
end
title('8Hz 0.04 cpd')
legend(num2str(x));
subplot(2,2,4)
col = strvcat('c', 'k', 'r');
x = [];
for iArea = 1:3
    [H, stats, xCDF, yCDF] = cdfplot_LG(all_fits_dir(iArea).x_orituned_8Hz_pt16cpd(:,8));
    plot(xCDF, yCDF, col(iArea,:));
    hold on
    x= [x; size(all_fits_dir(iArea).x_orituned_8Hz_pt16cpd(:,3),1)];
end
title('8Hz 0.16 cpd')
legend(num2str(x));
suptitle('Dir Index')
fn_out = fullfile('\\zmey\storlab\users\Lindsey\Projects\HVAs\2P Axon Imaging', ['DIR_' num2str(P) 'P_' inj], ['DIR_' num2str(P) 'P_' inj '_DI_iSFTF.ps']);
print(gcf, '-depsc', fn_out);

%all OI for bestSFTF for each SFTF
figure;
subplot(2,2,1)
col = strvcat('c', 'k', 'r');
x = [];
for iArea = 1:3
    [H, stats, xCDF, yCDF] = cdfplot_LG(all_fits_dir(iArea).x_bestSFTF_2Hz_pt04cpd(:,3));
    plot(xCDF, yCDF, col(iArea,:));
    hold on
    x= [x; size(all_fits_dir(iArea).x_bestSFTF_2Hz_pt04cpd(:,3),1)];
end
legend(num2str(x));
title('2Hz 0.04 cpd')
subplot(2,2,2)
col = strvcat('c', 'k', 'r');
x = [];
for iArea = 1:3
    [H, stats, xCDF, yCDF] = cdfplot_LG(all_fits_dir(iArea).x_bestSFTF_2Hz_pt16cpd(:,3));
    plot(xCDF, yCDF, col(iArea,:));
    hold on
    x= [x; size(all_fits_dir(iArea).x_bestSFTF_2Hz_pt16cpd(:,3),1)];
end
title('2Hz 0.16 cpd')
legend(num2str(x));
subplot(2,2,3)
col = strvcat('c', 'k', 'r');
x = [];
for iArea = 1:3
    [H, stats, xCDF, yCDF] = cdfplot_LG(all_fits_dir(iArea).x_bestSFTF_8Hz_pt04cpd(:,3));
    plot(xCDF, yCDF, col(iArea,:));
    hold on
    x= [x; size(all_fits_dir(iArea).x_bestSFTF_8Hz_pt04cpd(:,3),1)];
end
title('8Hz 0.04 cpd')
legend(num2str(x));
subplot(2,2,4)
x = [];
col = strvcat('c', 'k', 'r');
for iArea = 1:3
    [H, stats, xCDF, yCDF] = cdfplot_LG(all_fits_dir(iArea).x_bestSFTF_8Hz_pt16cpd(:,3));
    plot(xCDF, yCDF, col(iArea,:));
    hold on
    x= [x; size(all_fits_dir(iArea).x_bestSFTF_8Hz_pt16cpd(:,3),1)];
end
title('8Hz 0.16 cpd')
legend(num2str(x));
suptitle('Ori Index BestSFTF')
fn_out = fullfile('\\zmey\storlab\users\Lindsey\Projects\HVAs\2P Axon Imaging', ['DIR_' num2str(P) 'P_' inj], ['DIR_' num2str(P) 'P_' inj '_OI_bestSFTF_iSFTF.ps']);
print(gcf, '-depsc', fn_out);

%all DI for bestSFTF for each SFTF
figure;
subplot(2,2,1)
col = strvcat('c', 'k', 'r');
x = [];
for iArea = 1:3
    [H, stats, xCDF, yCDF] = cdfplot_LG(all_fits_dir(iArea).x_orituned_bestSFTF_2Hz_pt04cpd(:,8));
    plot(xCDF, yCDF, col(iArea,:));
    hold on
    x= [x; size(all_fits_dir(iArea).x_orituned_bestSFTF_2Hz_pt04cpd(:,3),1)];
end
title('2Hz 0.04 cpd')
legend(num2str(x));
subplot(2,2,2)
col = strvcat('c', 'k', 'r');
x = [];
for iArea = 1:3
    [H, stats, xCDF, yCDF] = cdfplot_LG(all_fits_dir(iArea).x_orituned_bestSFTF_2Hz_pt16cpd(:,8));
    plot(xCDF, yCDF, col(iArea,:));
    hold on
    x= [x; size(all_fits_dir(iArea).x_orituned_bestSFTF_2Hz_pt16cpd(:,3),1)];
end
title('2Hz 0.16 cpd')
legend(num2str(x));
subplot(2,2,3)
col = strvcat('c', 'k', 'r');
x = [];
for iArea = 1:3
    [H, stats, xCDF, yCDF] = cdfplot_LG(all_fits_dir(iArea).x_orituned_bestSFTF_8Hz_pt04cpd(:,8));
    plot(xCDF, yCDF, col(iArea,:));
    hold on
    x= [x; size(all_fits_dir(iArea).x_orituned_bestSFTF_8Hz_pt04cpd(:,3),1)];
end
title('8Hz 0.04 cpd')
legend(num2str(x));
subplot(2,2,4)
col = strvcat('c', 'k', 'r');
x = [];
for iArea = 1:3
    [H, stats, xCDF, yCDF] = cdfplot_LG(all_fits_dir(iArea).x_orituned_bestSFTF_8Hz_pt16cpd(:,8));
    plot(xCDF, yCDF, col(iArea,:));
    hold on
    x= [x; size(all_fits_dir(iArea).x_orituned_bestSFTF_8Hz_pt16cpd(:,3),1)];
end
title('8Hz 0.16 cpd')
legend(num2str(x));
suptitle('Dir Index BestSFTF')
fn_out = fullfile('\\zmey\storlab\users\Lindsey\Projects\HVAs\2P Axon Imaging', ['DIR_' num2str(P) 'P_' inj], ['DIR_' num2str(P) 'P_' inj '_DI_bestSFTF_iSFTF.ps']);
print(gcf, '-depsc', fn_out);

figure;
for iArea = 1:3
    subplot(2,3,iArea)
    [H, stats, xCDF, yCDF] = cdfplot_LG(all_fits_dir(iArea).x_bestSFTF_2Hz_pt16cpd(:,3));
    plot(xCDF, yCDF, 'b');
    hold on
    [H, stats, xCDF, yCDF] = cdfplot_LG(all_fits_dir(iArea).x_bestSFTF_2Hz_pt04cpd(:,3));
    plot(xCDF, yCDF, 'c');
    hold on
    [H, stats, xCDF, yCDF] = cdfplot_LG(all_fits_dir(iArea).x_bestSFTF_8Hz_pt16cpd(:,3));
    plot(xCDF, yCDF, 'g');
    hold on
    [H, stats, xCDF, yCDF] = cdfplot_LG(all_fits_dir(iArea).x_bestSFTF_8Hz_pt04cpd(:,3));
    plot(xCDF, yCDF, 'r');
    hold on
    title(areas(iArea,:))
    xlabel('OSI')
end
fn_out = fullfile('\\zmey\storlab\users\Lindsey\Projects\HVAs\2P Axon Imaging', ['DIR_' num2str(P) 'P_' inj], ['DIR_' num2str(P) 'P_' inj '_OI_iArea_bestSFTF_iSFTF.ps']);
print(gcf, '-depsc', fn_out);
for iArea = 1:3
    subplot(2,3,iArea+3)
    [H, stats, xCDF, yCDF] = cdfplot_LG(all_fits_dir(iArea).x_orituned_bestSFTF_2Hz_pt16cpd(:,8));
    plot(xCDF, yCDF, 'b');
    hold on
    [H, stats, xCDF, yCDF] = cdfplot_LG(all_fits_dir(iArea).x_orituned_bestSFTF_2Hz_pt04cpd(:,8));
    plot(xCDF, yCDF, 'c');
    hold on
    [H, stats, xCDF, yCDF] = cdfplot_LG(all_fits_dir(iArea).x_orituned_bestSFTF_8Hz_pt16cpd(:,8));
    plot(xCDF, yCDF, 'g');
    hold on
    [H, stats, xCDF, yCDF] = cdfplot_LG(all_fits_dir(iArea).x_orituned_bestSFTF_8Hz_pt04cpd(:,8));
    plot(xCDF, yCDF, 'r');
    hold on
    title(areas(iArea,:))
    xlabel('DSI')
end

figure;
for iArea = 1:3
    subplot(2,3,iArea)
    [H, stats, xCDF, yCDF] = cdfplot_LG(all_fits_dir(iArea).x_2Hz_pt16cpd(:,3));
    plot(xCDF, yCDF, 'b');
    hold on
    [H, stats, xCDF, yCDF] = cdfplot_LG(all_fits_dir(iArea).x_2Hz_pt04cpd(:,3));
    plot(xCDF, yCDF, 'c');
    hold on
    [H, stats, xCDF, yCDF] = cdfplot_LG(all_fits_dir(iArea).x_8Hz_pt16cpd(:,3));
    plot(xCDF, yCDF, 'g');
    hold on
    [H, stats, xCDF, yCDF] = cdfplot_LG(all_fits_dir(iArea).x_8Hz_pt04cpd(:,3));
    plot(xCDF, yCDF, 'r');
    hold on
    title(areas(iArea,:))
    xlabel('OSI')
end
for iArea = 1:3
    subplot(2,3,iArea+3)
    [H, stats, xCDF, yCDF] = cdfplot_LG(all_fits_dir(iArea).x_orituned_2Hz_pt16cpd(:,8));
    plot(xCDF, yCDF, 'b');
    hold on
    [H, stats, xCDF, yCDF] = cdfplot_LG(all_fits_dir(iArea).x_orituned_2Hz_pt04cpd(:,8));
    plot(xCDF, yCDF, 'c');
    hold on
    [H, stats, xCDF, yCDF] = cdfplot_LG(all_fits_dir(iArea).x_orituned_8Hz_pt16cpd(:,8));
    plot(xCDF, yCDF, 'g');
    hold on
    [H, stats, xCDF, yCDF] = cdfplot_LG(all_fits_dir(iArea).x_orituned_8Hz_pt04cpd(:,8));
    plot(xCDF, yCDF, 'r');
    hold on
    title(areas(iArea,:))
    xlabel('DSI')
end
fn_out = fullfile('\\zmey\storlab\users\Lindsey\Projects\HVAs\2P Axon Imaging', ['DIR_' num2str(P) 'P_' inj], ['DIR_' num2str(P) 'P_' inj '_DI_iArea_bestSFTF_iSFTF.ps']);
print(gcf, '-depsc', fn_out);