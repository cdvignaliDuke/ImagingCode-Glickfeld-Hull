imRGB = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\grace\Analysis\2P\200114_i1314\200114_i1314_runs-003\untitled.tif';
uiopen(imRGB);
imGRAY = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\grace\Analysis\2P\200114_i1314\200114_i1314_runs-003\grey.tif';
uiopen(imGRAY);
imRGB2 = cat(3,imGRAY,imGRAY,imGRAY);
reg_shifts = load('\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\grace\Analysis\2P\200114_i1314\200114_i1314_runs-003\200114_i1314_runs-003_reg_shifts.mat');
reg2 = reg_shifts.data_reg_avg;

%%Red Wash
figure;
imagesc(reg2);
% colormap(map)
% colormap(flipud(jet(256)));
clim([800 900]);
truesize;
axis off


%%Tone Scale
figure;
imHSV = rgb2hsv(imRGB);
imHSV(:,:,1) = 0./360;
imHSV(:,:,2) = 1 - imHSV(:,:,2);
result = hsv2rgb(imHSV);
imshow(result)

sz_target = size(reg2);
rgb_reg = zeros(sz_target(1), sz_target(2), 3);
rgb_reg(:,:,1) = reg2;
rgb_reg(:,:,2) = reg2;
rgb_reg(:,:,3) = reg2;
inHSV = rgb2hsv(rgb_reg(:,:,1));