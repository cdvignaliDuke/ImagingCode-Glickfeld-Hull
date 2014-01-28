function [neurons,glia,neuropil]=scriptFindCells(green_img,red_img,params);
%SCRIPTFINDCELLS Runs IMFINDCELLS to segment neurons and glia
%  [neurons,glia,neuropil]=scriptFindCells(green_img,red_img);
%

if nargin < 3
    [params]  = imFindCellsParamsSet   % displays default parameters
end
    
close all;

%% contrast adaptation
green_adapt = imScale(stackLocalContrastAdaptation(green_img,params.sigma,params.offset));
red_adapt = imScale(stackLocalContrastAdaptation(red_img,params.sigma,params.offset));

figure;imagesc(green_adapt);axis square;colormap gray;
title('Green adapt');

figure;imagesc(red_adapt);axis square;colormap gray;
title('Red adapt');

% assumes images have same scale
neurons_adapt= max(green_adapt - red_adapt,0);
glia_adapt = red_adapt;

figure;
imagesc(neurons_adapt);axis square;colorbar;colormap gray;
title('Neurons - Contrast adapted')
figure;
imagesc(glia_adapt);axis square;colorbar;colormap gray;
title('Glia - Contrast adapted')

%% image segmentation
[neurons, bwneurons] = imFindCells (neurons_adapt, params);
[glia, bwglia] = imFindCells (glia_adapt, params);

green_log = imLog(green_adapt);
red_log =  imLog(red_adapt);
neurons_log = imLog(neurons_adapt);
glia_log = imLog(glia_adapt);

figure;imshow(imScale(neurons_log,'uint16'));title('Neurons - Input image');
figure;imshow(imShade(neurons_log,bwneurons));title('Neurons - Segmentation');
figure;imshow(imScale(glia_log,'uint16'));title('Glia - Input image');
figure;imshow(imShade(glia_log,bwglia));title('Glia - Segmentation');

overlay = cat(3,bwglia,bwneurons,zeros(size(bwneurons)));
figure;imshow(overlay);title('Neurons (green) and glia (red) segmentation');
figure;imshow(green_log);
figure;imshow(imShade(green_log,bwglia,bwneurons))

rgb = cat(3,red_log,green_log,zeros(size(red_log)));
figure;imshow(rgb);

%% create neuropil mask
se = strel('disk',3);
neuropil = ~imdilate(bwneurons|bwglia,se);
figure
imagesc(neuropil);axis square;
title('Neuropil mask');

return;