

clear;
cd('Z:\Data\2P_imaging\170426_img90\img90');
%cd('Z:\Data\2P_imaging\img93\img93');
fname = 'img90_000_000';
for ii = 1%:20
    data = squeeze(sbxread(fname,1,1000));
    data_avg = mean(data,3);
    data_max= max(data,[],3);
    figure; imagesc(data_max); colormap gray;
    title(['30HZ img90 170426 recon max proj frames 1:1000']);
end


%% write a tiff movie to analyze diff. in sessions

cd('Z:\Data\2P_imaging\170518_img93\img93');
fname = 'img93_000_000';
for ii = 1%:20
    data = squeeze(sbxread(fname,0,5000));
    data_avg = mean(data,3);
    data_max= max(data,[],3);
    figure; imagesc(data_max); colormap gray;
    title(['30HZ img93 170518 recon max proj frames 1:5000']);
end
writetiff(data, 'Z:\Data\2P_imaging\170518_img93\img93_tiff_1_5000');