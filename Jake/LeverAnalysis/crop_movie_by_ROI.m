function [img, sz]  = crop_movie_by_ROI(image_source, meta_data2, ROI_x, ROI_y, frame_info)
%loads the tiff files, crops them based on the ROI and then concatenates them into img. 

%allocate memory and set variables
img = [];
sz = [];
first_frame = 1; 
last_frame = length(frame_info.times);

%download raw movie
img = [];
for ii = 1:length(meta_data2)
    img_subset = readtiff(meta_data2{ii}(1).Filename);
    img = cat(3,img,img_subset);
    disp(['file #', num2str(ii), ' download complete']);
end

% crop the movie according to ROI coords and calculate frame size
img = img(ROI_x, ROI_y, :);
sz = size(img);
return

