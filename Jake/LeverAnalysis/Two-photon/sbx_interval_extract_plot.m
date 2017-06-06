cd('Z:\Data\2P_imaging\170529_img94\img94');
fname = 'img94_000_000';
for ii = 1%:20
    data = squeeze(sbxread(fname,10000,5000));
    data_avg = mean(data,3);
    data_max= max(data,[],3);
    figure; imagesc(data_max); colormap gray;
    title(['30HZ img94 170529 max proj frames 10000-15000']);
end