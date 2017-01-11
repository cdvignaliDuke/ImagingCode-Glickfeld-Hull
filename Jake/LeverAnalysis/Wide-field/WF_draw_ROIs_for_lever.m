function WF_draw_ROIs_for_lever_full_size(days, image_source, image_dest);
%Use this function to manually draw ROIs. 

% ---------------
image_source  = [image_source '\' days '_MMStack.ome.tif'];
%----loads the a sample of the imaging data
img = [];
for ii = 215:265;
    img_temp = double(imread(image_source,ii));
    img = cat(3,img,img_temp);
end

%alternative method for selecting image sample
% for i = 1:20
%     for ii = 51:80;
%         img_temp = double(imread(image_source,ii));
%         img = cat(3,img,img_temp);
%     end
%     if min(mean(mean(img,2),1))>0
%         break
%     end
% end

% calculate average and other states
avg_img = squeeze(mean(img,3));
std_img = std(img,[], 2);
all_sd = std(img(:));
clear img;

%plot the avg_img to draw the ROIs
f = figure; im = imagesc(avg_img); axis ij; colormap jet; colormap jet; shading flat;
num_cluster = input('Enter no of clusters:');

for i=1:num_cluster
    roi(i) =  impoly;
    disp(['Finished ROI ' num2str(i)]);
end

%save size/location information for each ROI
roi_mask = [];
roi_position = {};
for i=1:length(roi)
    c_mask = roi(i).createMask(im);
    roi_mask(i,:) = c_mask(:);
    roi_position{i} = roi(i).getPosition;
    roi_position{i} = [roi_position{i}; roi_position{i}(1,:)];
end
cluster.roi_mask =roi_mask;
cluster.roi_mask_combined = sum(roi_mask,1);
cluster.roi_position = roi_position;
cluster.num_cluster = num_cluster;
close gcf;
save([image_dest '_cluster'], 'cluster');
return
