cd('Z:\Data\2P_imaging\170501_img92\img92');
fname = 'img92_000_000';
for ii = 1%:20
    data = squeeze(sbxread(fname,ii*1,900));
    data_avg = mean(data,3);
    data_max= max(data,[],3);
    figure; imagesc(data_max); colormap gray;
    title(['img92 170501 max proj frames ', num2str(ii*5000) ':' num2str(ii*5000+100)]);
end