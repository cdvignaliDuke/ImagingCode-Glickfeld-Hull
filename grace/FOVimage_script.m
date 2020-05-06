CD = 'Z:\All_Staff\home\grace\2P_Imaging\190801_i1306\001';
cd(CD);
load('001_000_001.mat')
nframes = info.config.frames;
data = sbxread('001_000_001', 0, nframes);
if size(data,1) == 1
    data_g = squeeze(data);
else
    data_r = squeeze(data(2,:,:,:));
    data_g = squeeze(data(1,:,:,:));
end
figure; subplot(2,1,1); imagesc(mean(data_g,3))
subplot(2,1,2); imagesc(mean(data_r,3));

data_avg = mean(data_r(:,:,101:200),3);
figure; imagesc(data_avg)
[out data_reg] = stackRegister(data_r,data_avg);
data_r_reg_avg = mean(data_reg,3);
figure; imagesc(data_r_reg_avg);
axis off
truesize
title('190801 i1306 1040 nm - Red Channel')
colormap gray
print('Z:\All_Staff\home\grace\Analysis\2P\190801_i1306\190801_i1306_run-001_FOVred.pdf','-dpdf','-bestfit')
save('Z:\All_Staff\home\grace\Analysis\2P\190801_i1306\190801_i1306_run-001_avgFOVred.mat','data_r_reg_avg')
writetiff(data_r_reg_avg, 'Z:\All_Staff\home\grace\Analysis\2P\190801_i1306\190801_i1306_run-001_avgFOVred.tiff')

[outs data_g_reg] = stackRegister_MA(data_g, [],[], out);
data_g_reg_avg = mean(data_g_reg,3);
figure; imagesc(data_g_reg_avg);
axis off
truesize
title('190801 i1306 1040 nm - Green Channel')
colormap gray
print('Z:\All_Staff\home\grace\Analysis\2P\190801_i1306\190801_i1306_run-001_FOVgreen.pdf','-dpdf','-bestfit')
save('Z:\All_Staff\home\grace\Analysis\2P\190801_i1306\190801_i1306_run-001_avgFOVgreen.mat','data_g_reg_avg')
writetiff(data_g_reg_avg, 'Z:\All_Staff\home\grace\Analysis\2P\190801_i1306\190801_i1306_run-001_avgFOVgreen.tiff')

sz = size(data_reg);
rgb = nan(sz(1),sz(2),3);
rgb(:,:,1) = data_r_reg_avg./max(max(data_r_reg_avg,[],1),[],2);
rgb(:,:,2) = data_g_reg_avg./(max(max(data_r_reg_avg,[],1),[],2)./4);
figure; imshow(double(rgb))
axis off
truesize
title('190801 i1306 1040 nm')
print('Z:\All_Staff\home\grace\Analysis\2P\190801_i1306\190801_i1306_run-001_FOV_RGB.pdf','-dpdf','-bestfit')
