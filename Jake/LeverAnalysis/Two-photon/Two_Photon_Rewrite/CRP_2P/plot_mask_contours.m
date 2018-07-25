function plot_mask_contours(mask3D, avg_reg_img)
%function to draw contours of individual masks over a single FoV.
%either draw contours over a black background or an average image
if isempty(avg_reg_img)
    canvas_img = zeros(size(mask3D,1), size(mask3D,2));
    color_spec = [1,0.4,0.7;   1,1,0;   0,1,0;  1,0,1;   1,0.5,0;   0.2,0.6,1,;   1,0,0;  1,1,1];
else
    canvas_img = avg_reg_img; 
    color_spec = [1,0.4,0.7;   1,1,0;   0,1,0;  1,0,1;   1,0.5,0;   0.2,0.6,1,;   1,0,0;]; %dont use white if plotting on grayscale image
end

spec_ind = repmat([1:size(color_spec,1)], 1, 100);
figure; 
%plot the background
imagesc(canvas_img); colormap gray;
hold on; 

%plot each contour
for mask_num = 1:size(mask3D,3)
    contour(mask3D(:,:,mask_num), 'Color', [color_spec(spec_ind(mask_num),:)]);
end

end