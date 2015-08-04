%% scale bar
siz1 = size(data_reg,1);
siz2 = size(data_reg,2);
sizImg = zeros(siz1,siz2);

umperpixel = 550/size(data_reg,2);
pixels50um = ceil(50/umperpixel);
pixels25um = ceil(25/umperpixel);


sizImg(siz1-50:siz1-45,siz2-(pixels25um+30):siz2-30) = 1;
figure;imagesq(sizImg);colormap(gray)

figure; imagesq(maxDFoverF); colormap(gray)

maxDFoverFImg = maxDFoverF;

maxDFoverFImg = maxDFoverF(10:siz1-10,200:siz2-100);

figure; imagesq(maxDFoverFImg); colormap(gray)

siz3 = size(maxDFoverFImg,1);
siz4 = size(maxDFoverFImg,2);
sizImg2 = zeros(siz3,siz4);
sizImg2(siz3-50:siz3-47,siz4-(pixels25um+30):siz4-30) = 1;
figure;imagesq(sizImg2);colormap(gray)

scalebarnumber = max(max(maxDFoverFImg));

maxDFoverFImg(siz3-50:siz3-47,siz4-(pixels25um+30):siz4-30) = scalebarnumber;
figure; imagesq(maxDFoverFImg); colormap(gray)

writetiff(maxDFoverFImg,'maxDFoverFscbr')

%% cells mask img
cell1 = 111;
cell1mask = zeros(siz1,siz2);
cell1mask(mask_cell == cell1) = 1;
figure;imagesq(cell1mask);colormap(gray)
cell1mask = cell1mask(10:siz1-10,200:siz2-100);
figure;imagesq(cell1mask);colormap(gray)

writetiff(cell1mask, ['cell' num2str(cell1) 'mask'])