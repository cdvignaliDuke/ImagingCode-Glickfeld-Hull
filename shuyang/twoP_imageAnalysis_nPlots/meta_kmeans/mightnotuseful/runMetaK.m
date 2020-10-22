function runMetaK(downsampled_movie, outdir, figName, fileName)
z = size(downsampled_movie,3);
% downsampled_movie = downsampled_movie(:,:,1:2:10000);
[x,y,z2] = size(downsampled_movie);
if z > 10000       % reduce frames if too large
    downsampled_movie = downsampled_movie(:, :, 1:round(z*0.7));
    %     downsampled_movie = downsampled_movie(:, :, 1:4236);
    %     z = 4236;
     z2 = round(z * 0.7);
end

movie2d = double(reshape(downsampled_movie, x*y, z2));

% normalize df/f
meanF = mean(movie2d, 2);

df_f = bsxfun(@minus, movie2d, meanF);

df_f = bsxfun(@rdivide, df_f, meanF);

max_dff = max(df_f, [], 2);

df_f_norm = bsxfun(@rdivide, df_f, max_dff);

% clear downsampled_movie
 df_f_frame = reshape(df_f_norm, x, y, z2);
%  sc = [df_f_frame(10,25,:);df_f_frame(28,45,:);df_f_frame(37,40,:);df_f_frame(40,65,:);df_f_frame(48,77,:);df_f_frame(51,81,:)]; % img46
% sc = [df_f_frame(18,15,:);df_f_frame(30,15,:);df_f_frame(9,39,:);df_f_frame(32,47,:);df_f_frame(33,34,:);df_f_frame(33,34,:)];
%  sc = [df_f_frame(27,19,:);df_f_frame(33,23,:);df_f_frame(75,46,:);df_f_frame(75,46,:);df_f_frame(84,52,:)]; %img 30 151009
%  sc = [df_f_frame(54,40,:);df_f_frame(61,54,:);df_f_frame(68,60,:);df_f_frame(65,73,:);df_f_frame(6,25,:)]; %img 38
% sc = [df_f_frame(53,42,:);df_f_frame(53,60,:);df_f_frame(64,64,:);df_f_frame(67,43,:);df_f_frame(6,22,:)]; %img 38
%  sc = squeeze(sc);
 sc = 0;
% run meta k means
[allclusters, centroidcorr, dendmem, dunnsinitial]=meta_k_means(df_f_norm, 'correlation', sc);

% figure; fig = imagesc(coocMatrix);
% colorbar; title('Co-occurence Matrix'); colormap jet;
% saveas(fig, [outdir,'co-occurence',date,'.fig']);

% assign color to different cluster
coocMatrix = zeros(x*y, x*y);
indx = zeros(x*y, 1);
cc = reshape(clusters, x*y, 1);
for i = 1:max(max(cc))
%     indx((allclusters{i, 1})) = i;
    temp = cc==i;
    for f = 1:length(temp)
        coocMatrix(temp(f), temp) = i;
        coocMatrix(temp, temp(f)) = i;
    end
end


clusters = reshape(indx, x, y);
figure; fig = imagesc(clusters); 
title(['#frame = ', num2str(z)]); 
saveas(fig, [outdir, figName]);

%clear downsampled_movie
save([outdir, fileName], 'clusters', 'centroidcorr', 'dendmem', 'dunnsinitial', 'clusters', 'coocMatrix');

end