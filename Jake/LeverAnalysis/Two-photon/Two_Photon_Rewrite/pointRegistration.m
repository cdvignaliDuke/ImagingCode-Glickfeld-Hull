%% control point method
%alternate method where you pick three points on the two images you think
%are the same and the algorithm alings them. May be altering the image too
%much. 
% data_avg = mean(img2,3);
% target_avg = mean(img2,3);
% AVG = double(data_avg);
% AVG(find(AVG>1e5)) = 0;
% AVG = (AVG./max(max(abs(AVG)))); 
% target = double(target_avg);
% target(find(target>1e5)) = 0;
% target = (target./max(max(abs(target)))); 
% 
% AVG = mask2;
% target = mask1;
% [input_points, base_points] = cpselect(AVG,target,'Wait', true) 
% sz_target  = size(target);
% mytform    = maketform('affine',input_points(1:3,:), base_points(1:3,:));
% registered2 = imtransform(mean(img_reg_1000_2,3),mytform,'XData',[1 sz_target(2)],'YData',[1 sz_target(1)]); 

%% Thresholding Method
% session 1, do not need many frames
% homomorphic filtering remove multiplicative noise
% [data, pre] = homomorphicFilter(img(:,:,1:10:end/2));
data = NaN(size(img_reg,1), size(img_reg,2), floor(size(img_reg,3)/5));
for ff = 1:size(data,3)
    if ismember(ff/100, [0:100:5000]/100);
        ff
    end
    data(:,:,ff) = mean(img_reg(:,:,(ff-1)*5+1:ff*5),3);
    %data = cat(3, data, mean(img_reg(:,:,(ff-1)*5+1:ff*5),3));
end

data2 = [];
for ff = 1:floor(size(img_reg2,3)/5)
    ff
    data2 = cat(3, data2, mean(img_reg2(:,:,(ff-1)*5+1:ff*5),3));
end

[data_filtered, pre] = homomorphicFilter(img_reg(:,:,1:4:end)); % helps to remove background noise
% select rois based on intensity
roi = detectSingleFrameRois(data);
% combine/remove overlapped cells 
mask_final = processMaskpfgc(roi);
% combine highly correlated cells
threshold = 0.8;  %signal correlation threshold
[ ~, mask3D, ~] = finalMask(img_reg(:,:,1:10000), mask_final, threshold, out_dir);

% session 2

[data2, pre] = homomorphicFilter(img_reg2(:,:,1:4:end));
roi_2 = detectSingleFrameRois(data2);
mask_final_2 = processMaskpfgc(roi_2);
threshold = 0.8;
[ ~, mask3D, ~] = finalMask(img_reg2(:,:,1:10000), mask_final_2, threshold, out_dir);

% some ways to show cells like contour
figure;
imshow(mat2gray(mean(overlay,3)));truesize
for i = 1:size(mask3D,3) 
    hold on
    contour(mask3D(:,:,i),'g')
    
end
for i = 1:size(mask3D_2,3)
    hold on
    mask_shift = imtranslate(mask3D_2(:,:,i), [xshift,yshift]);
    contour(mask_shift,'r')
end

% the amount of xy shift based on visual inspection, can adjust after 
% viewing the merged mask in ImageJ
xshift = 12; %img92 ((87 - 107) + (99 - 120))/2 or -24; img90 -6;
yshift = 10;%img92 ((184 - 192) + (60 - 66))/2 ; img90 12;
img_reg_1000_1 = img_reg(:,:,1:1000);
img_reg_1000_2 = img_reg2(:,:,1:1000);

% using Matlab imtranslate first to move image and then stackRegister to
% register and upsampling
registered = imtranslate(img_reg2(:,:,1), [xshift,yshift]);
[~,registered] = stackRegister(img_reg_1000_2, registered, 10);

% just some ways to show before and after
% session 1 sample image
[data, ~] = homomorphicFilter(img_reg_1000_1(:,:,2:end));
mean_data1 = mean(data,3);
% figure;imagesc(mean_data1);truesize;
% figure;imshow(mat2gray(mean_data1)); truesize
% session 2 sample image
[data2_ori, ~] = homomorphicFilter(img_reg_1000_2(:,:,2:end));
mean_data2_ori = mean(data2_ori,3);
figure;imagesc(mean_data2_ori);truesize
% cross session registration example of session 2
[data2, ~] = homomorphicFilter(registered(:,:,2:end));
mean_data2 = mean(data2,3);
figure;subplot(1,2,1);imagesc(mean_data1); axis off
subplot(1,2,2);imagesc(img_comb_ref); axis off; truesize

figure;subplot(1,2,1);imshow(histeq(mat2gray(mean_data1))); axis off
subplot(1,2,2);imshow(mat2gray(img_comb_ref)); axis off;truesize

[D, img_comb_ref] = imregdemons(mean_data2, mean_data1, [32 16 8 4],'AccumulatedFieldSmoothing',2.5,'PyramidLevels',4);

overlay = imfuse(mean_data1, img_comb_ref,'falsecolor','Scaling','joint','ColorChannels',[1 2 0]);figure;imagesc(overlay);truesize

mask1 = mask_raw; mask1(mask1>0) = 1; 
mask2 = warp_test; mask2(mask2>0) = 2;
mask3 = mask1 + mask2;
mask1 = mask_raw; mask2 = warp_test;
mask1(mask3==3) = 200;
mask2(mask3==3) = 0;
mask3 = mask1 + mask2;
overlay = imfuse(mask1, mask2,'falsecolor','Scaling','independent','ColorChannels',[2 1 0]);figure;imagesc(overlay);truesize

overlay = imfuse(mask1,mask2,'blend','Scaling','joint');figure;imagesc(overlay);truesize
% boundaries = rangefilt(mask2,ones(3)) > 0;
% figure;imoverlay(mask1,boundaries,[0 0 0]);