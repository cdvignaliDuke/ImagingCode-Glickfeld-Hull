function cropped_movie = crop_image_according_to_blackout_mask(full_movie_reg, full_frame_mask, dir_for_cropped_image_coords );

%check for existing x/y min/max coords in PCA_output_dir and
%kmeans_output_dir. If neither exist then determine new coords and save
%them.

%check for existing file
if exist(dir_for_cropped_image_coords) == 2;
    load(dir_for_cropped_image_coords);
    coords_already_exist = 1;
else
    coords_already_exist = 0;
    %find coordinates to crop image based on the min and max x/y values of non-zero values
    x_max = max(find(sum(full_frame_mask,1)));
    x_min = min(find(sum(full_frame_mask,1)));
    y_max = max(find(sum(full_frame_mask,2)));
    y_min = min(find(sum(full_frame_mask,2)));
end
cropped_movie = full_movie_reg([y_min:y_max], [x_min:x_max], :);
if coords_already_exist == 0;
    save(dir_for_cropped_image_coords, 'x_min', 'x_max', 'y_min', 'y_max');
end
end 