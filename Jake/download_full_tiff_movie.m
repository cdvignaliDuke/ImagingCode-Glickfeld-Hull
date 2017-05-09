function full_movie = download_full_tiff_movie(meta_data2)
%uses the information contained in the meta_data of the tiff file to
%download all the frames of a tiff movie. If the movie is large it will be
%divided up amongst separate files. This function reads those files in
%order and concatenates the frames into one matrix. 
full_movie = []; 
    for frame_file_num = 1:length(meta_data2)-1 %cycle through all of the full sized frame containing files. Each one should be 4000+frames
        img_subset = readtiff(meta_data2{frame_file_num}(1).Filename);
        full_movie = cat(3, full_movie, img_subset);
        disp(['file number ', num2str(frame_file_num), ' download complete']);
    end
end
