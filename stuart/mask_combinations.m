normImage = mat2gray(data_dfof_max_lind);
normImage1 = mat2gray(data_dfof_max{4});
figure; imagesq(normImage)
figure; imagesq(normImage1)
figure; imagesq(normImage+normImage1)
figure; imagesq(normImage+normImage1*2)%this one looks strongest with most cells

data_dfof_max{5}=(normImage+normImage1*2);


