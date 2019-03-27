%b = bwboundaries(mask_cell(:,:,1));

img = imread('Z:\Analysis\2P_MovingDots_Analysis\imaging_analysis\190131_img1018\AVG_190131_img1018_001_jpeg_1_1000.jpg');
imshow('Z:\Analysis\2P_MovingDots_Analysis\imaging_analysis\190131_img1018\AVG_190131_img1018_001_jpeg_1_1000.jpg'); hold on;
%figure;
for i  = 1:size(mask3D,3)
    bound = cell2mat(bwboundaries(mask3D(:,:,i)));
    randcolor = rand(1,4);
    plot(bound(:,2),bound(:,1),'.','color',randcolor); hold on;
end