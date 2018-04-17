 tt = (1:size(data_tc_spont(:,sampleN(ic)),1)) * 33;
    hold on;
 subplot(2,3,3);   plot(tt, data_tc_spont(:,sampleN(ic)) + 10, 'k');
%     plot( data_tc_spont(:,sampleN(ic)) + 2, 'k');
    hold on;
subplot(2,3,3);    plot(tt(1:end-1), events_diff_x1(:,sampleN(ic)) - 3, 'b');
%     plot( events_diff_x1(:,sampleN(ic)), 'b');
    scale = max(events_diff_x1(:,sampleN(ic)));
subplot(2,3,3);    plot( events_ind{sampleN(ic)}*33, scale*0.5, 'k.');


hold on; subplot(1,2,1);line([700 776], [250 250], 'color', 'white', 'LineWidth', 2)
subplot(1,2,1);text('position', [230 81], 'fontsize', 6, 'string', '100 µm', 'color', 'white')

figure; ax1=subplot(2,2,1);imagesq(img_mean);colormap(ax1,gray)
hold on; subplot(2,2,1);line([700 776], [250 250], 'color', 'white', 'LineWidth', 2)
subplot(2,2,1);text('position', [230 81], 'fontsize', 6, 'string', '100 µm', 'color', 'white')

m1 = maskTmp;
m1(m1>100) = 0;
mask_rgb(:,:,1) = m1*0.01;
m1(m1==100) = 0.6863;
mm = reshape(mask_final, 264, 796);
mm(mm == 15) = 0;
mm(mm>0)= 1;
m2 = m1+mm*0.8745;
mask_rgb(:,:,2) = m2;
m4 = mm;
m4(m4<1) = 0.5608;
mask_rgb(:,:,3) = m4;
ax2=subplot(2,2,2);imagesq(mask_rgb);


m1 = maskTmp;
m1(m1==100 | m1 == 102) = 0;
mask_rgb(:,:,1) = m1*0.01;
m1(m1==101) = 0.6863;
mm = reshape(mask_final, 264, 796);
mm(mm == 43) = 0;
mm(mm>0)= 1;
m2 = m1+mm*0.8745;
mask_rgb(:,:,2) = m2;
m4 = mm;
m4(m4<1) = 0.5608;
mask_rgb(:,:,3) = m4;
ax3=subplot(2,2,3);imagesq(mask_rgb);

m1 = maskTmp;
m1(m1==101 | m1 == 100) = 0;
mask_rgb(:,:,1) = m1*0.01;
m1(m1==102) = 0.6863;
mm = reshape(mask_final, 264, 796);
mm(mm == 55) = 0;
mm(mm>0)= 1;
m2 = m1+mm*0.8745;
mask_rgb(:,:,2) = m2;
m4 = mm;
m4(m4<1) = 0.5608;
mask_rgb(:,:,3) = m4;
ax4=subplot(2,2,4);imagesq(mask_rgb);