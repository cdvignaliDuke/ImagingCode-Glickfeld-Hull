data_gcamp = double(squeeze(data(1,:,:,:)));
data_tdtomato = double(squeeze(data(1,:,:,:)));

dF_gcamp = bsxfun(@minus,data_gcamp,mean(data_gcamp,3));
dFoverF_gcamp = bsxfun(@rdivide,dF_gcamp,mean(data_gcamp,3));

img_tdtomato = mean(data_tdtomato,3);

images = figure;
colormap gray
subplot(2,1,1)
imagesc(img_tdtomato)
title('Red PMT, 920nm')
subplot(2,1,2)
imagesc(mean(data_gcamp,3))
title('Green PMT, 920nm')

figure;
imagesc(max(dFoverF_gcamp,[],3))
colormap gray