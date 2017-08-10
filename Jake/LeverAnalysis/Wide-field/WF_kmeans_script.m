%%WF_kmeans_script
clear;
days = {'151021_img29', '151022_img29', '151009_img30', '151011_img30', '151211_img32', '151212_img32', '160129_img35', '160131_img35', '160129_img36', ...
    '160131_img36', '160314_img38', '160315_img38', '160319_img41', '160320_img41', '160606_img46', '160722_img53', '160904_img55'}; %'150718_img27', '150719_img27', '150716_img28', '150717_img28', 
bx_source     = ['Z:\Data\WidefieldImaging\GCaMP\behavior\'];
image_source_base  = ['Z:\Data\WidefieldImaging\GCaMP\']; %location of permanently stored image files for retreiving meta data
image_dest_base    = ['Z:\Analysis\WF Lever Analysis\BxAndAnalysisOutputs\']; %stores the data on crash in the lever analysis folder
bx_outputs_dir =     ['Z:\Analysis\WF Lever Analysis\BxAndAnalysisOutputs\BxOutputs\'];
PCA_output_dir_base =    ['Z:\Analysis\WF Lever Analysis\PCA_output_dir\'];
kmeans_inputs_dir =      ['Z:\Analysis\WF Lever Analysis\Meta-kmeans_inputs\'];
kmeans_output_dir_base = ['Z:\Analysis\WF Lever Analysis\kmeans_output_dir\'];
old_cd = cd; %save old cd so I can restore it later
downsample_factor = 0.2;


for ii = [3,4,5,6,11,12]  
    if ~exist([kmeans_inputs_dir, days{ii}, '\', days{ii}, '_stable_movie_avg.mat'])
        stable_frame_int = input(['please open first movie file for ' days{ii} ' in imagej and select a frame interval in the form of a vector with no/few motion artifacts']);
        meta_data_dir = [image_dest_base days{ii} '\' days{ii} '_meta_data'];
        load(meta_data_dir, 'meta_data2');
        stable_movie = [];
        for iii = stable_frame_int
            img_subset = imread(meta_data2{1}(1).Filename,iii);
            stable_movie = cat(3, stable_movie, img_subset);
        end
        stable_movie_avg = mean(stable_movie, 3);
        save([kmeans_inputs_dir, days{ii}, '\', days{ii}, '_stable_movie_avg'], 'stable_movie_avg');
    end
end
for ii = [3,4,5,6,11,12]  %original three datasets = [15, 16, 17]     Additionaly four datasets which could be analyzed = [9,10,13,14]
    PCA_output_dir = [PCA_output_dir_base days{ii} '\'];
    kmeans_output_dir = [kmeans_output_dir_base days{ii} '\'];
    meta_data_dir = [image_dest_base days{ii} '\' days{ii} '_meta_data'];
    load(meta_data_dir);
    
    %check for existing downsampled, cropped, mask-applied movie
    %     if exist([kmeans_output_dir, days{ii}, '_downsampled_cropped_motionCorr.mat'], 'file') == 0;
    
    %download the entire movie
    full_movie = download_full_tiff_movie(meta_data2);
    
    %motion correction of the entire movie
    stable_movie_avg = load([kmeans_inputs_dir, days{ii}, '\', days{ii}, '_stable_movie_avg']);
    stable_movie_avg = stable_movie_avg.stable_movie_avg;
    [out full_movie_reg] = stackRegister(full_movie, stable_movie_avg);  %does stackRegister alter the size of each frame?
    clear full_movie; %the motion correction does not alter the size of the movie
    
    %Check for existing black out mask of full frame. If none then manually generate one in imageJ
    if exist([kmeans_inputs_dir, days{ii}, '_movies_full_frame_mask.tif'], 'file') == 2;
        full_frame_mask = readtiff([kmeans_inputs_dir, days{ii}, '_movies_full_frame_mask.tif']);
    else
        disp(['analysis paused. Please check to make sure a mask exist for ', days{ii}, ' if not then manually generate one in imageJ now']);
        pause
    end
    
    %apply full frame mask
    full_movie_size = size(full_movie_reg);
    full_frame_mask = reshape(full_frame_mask, 1, full_movie_size(1)*full_movie_size(2));
    full_movie_reg = reshape(full_movie_reg, full_movie_size(1)*full_movie_size(2), full_movie_size(3) );
    zeros_in_mask = find(full_frame_mask == 0);
    full_movie_reg(zeros_in_mask, :) = 0;
    full_movie_reg = reshape(full_movie_reg, full_movie_size(1), full_movie_size(2), full_movie_size(3));
    full_frame_mask = reshape(full_frame_mask, full_movie_size(1), full_movie_size(2));
    
    %crop movie
    PCA_dir_for_cropped_coords = [PCA_output_dir days{ii} '_cropped_image_coords.mat'];
    cropped_movie = crop_image_according_to_blackout_mask(full_movie_reg, full_frame_mask, PCA_dir_for_cropped_coords); %automatically saves
    cropped_movie = double(cropped_movie);
    
    %down sample movie to average activity over 5x5pixel blocks
    tic
    downsampled_movie = double(imresize(cropped_movie, downsample_factor));
    toc
    
    %save the downsampled, cropped, mask applied movie to kmeans_output_dir
    save([kmeans_output_dir, days{ii}, '_downsampled_cropped_motionCorr_blackout'], 'downsampled_movie');
    
    
    %     elseif exist([kmeans_output_dir, days{ii}, '_downsampled_cropped_motionCorr.mat'], 'file') == 2;
    %             load([kmeans_output_dir, days{ii}, '_downsampled_cropped_motionCorr.mat']);
    %     end
end


    %%this is the code for actually performing the kmeans analysis=======================================================
    %needs to be put inside the forloop
    %check size of the movie
    %     movie_size = size(downsampled_movie);
    %     if length(movie_size)==3;
    %         downsampled_movie = reshape(downsampled_movie, movie_size(1)*movie_size(2), movie_size(3));
    %     end
    %
    %     %Run kmeans on df/f data
    %     %downsampled_movie 2D, consists of zeros in masked areas and integers
    %     %from 1 to 7 in F areas.
    %     downsampled_movie_mean = mean(downsampled_movie,2);  %dim1=frames dim2=pixels
    %     downsampled_movie_mean = repmat(downsampled_movie_mean, 1,size(downsampled_movie,2));
    %     dfof_downsampled = (downsampled_movie-downsampled_movie_mean)./downsampled_movie_mean;
    %     avg_f_for_sub = reshape(downsampled_movie_mean(:,1), movie_size(1), movie_size(2));
    %     downsampled_movie_max = max(dfof_downsampled,[],2);
    %     downsampled_movie_max = repmat(downsampled_movie_max, 1, size(downsampled_movie,2));
    %     downsampled_movie_min = min(dfof_downsampled,[],2);
    %     downsampled_movie_min = repmat(downsampled_movie_min, 1, size(downsampled_movie,2));
    %     downsampled_norm_max = dfof_downsampled./downsampled_movie_max;
    %     downsampled_norm_min_max = (dfof_downsampled-downsampled_movie_min)./(downsampled_movie_max-downsampled_movie_min);
    %
    %     %run a pair wise correlation to determine number of clusters
    %     pixel_pairwise_correlation = corr(downsampled_movie');
    %     save([kmeans_output_dir, days{ii}, '_pixel_pairwise_correlation'],'pixel_pairwise_correlation');
    %
    %     %Remove pairwise correlations of the blacked out pixels with each other
    %     %and F pixels. Look for rows of all NaNs, index them and remove them from both axes.
    %     index_of_nan_rows = [];
    %     for iv = 1:length(pixel_pairwise_correlation);
    %         if isempty(find(~isnan(pixel_pairwise_correlation(iv,:))));
    %             index_of_nan_rows = [index_of_nan_rows, iv];
    %         end
    %     end
    %     pixel_pairwise_correlation(:,index_of_nan_rows) = [];
    %     pixel_pairwise_correlation(index_of_nan_rows, :) = [];
    %     figure;imagesc(pixel_pairwise_correlation);
    %
    %     %run kmeans in order to cluster the downsample pixels into groups
    %     %time on y-axis and pixels on x-axis? then reshape to view the frame
    %     for num_clusters = 8:10;
    %         kmeans_inx = kmeans(dfof_downsampled,num_clusters);
    %         kmeans_inx = reshape(kmeans_inx, movie_size(1), movie_size(2));
    %         figure; subplot(1,2,1); imagesc(kmeans_inx); colorbar;
    %         title([days{ii} ' dfof kmeans clustering # clusters = ' num2str(num_clusters)]);
    %         num_pixels = [];
    %         for iv = 1:num_clusters
    %             num_pixels = [num_pixels; length(find(kmeans_inx==iv))];
    %         end
    %         subplot(1,2,2); bar([1:num_clusters], num_pixels');
    %     end
    %     for num_clusters = 8:10;
    %         kmeans_inx = kmeans(downsampled_norm_max,num_clusters);
    %         kmeans_inx = reshape(kmeans_inx, movie_size(1), movie_size(2));
    %         figure; subplot(1,2,1); imagesc(kmeans_inx); colorbar;
    %         title([days{ii} ' normalized to max: kmeans clustering # clusters = ' num2str(num_clusters)]);
    %         num_pixels = [];
    %         for iv = 1:num_clusters
    %             num_pixels = [num_pixels; length(find(kmeans_inx==iv))];
    %         end
    %         subplot(1,2,2); bar([1:num_clusters], num_pixels');
    %     end
    %     for num_clusters = 8:10;
    %         kmeans_inx = kmeans(downsampled_norm_min_max,num_clusters);
    %         kmeans_inx = reshape(kmeans_inx, movie_size(1), movie_size(2));
    %         figure; subplot(1,2,1); imagesc(kmeans_inx); colorbar;
    %         title([days{ii} ' normalized to max and min sub: kmeans clustering # clusters = ' num2str(num_clusters)]);
    %         num_pixels = [];
    %         for iv = 1:num_clusters
    %             num_pixels = [num_pixels; length(find(kmeans_inx==iv))];
    %         end
    %         subplot(1,2,2); bar([1:num_clusters], num_pixels');
    %     end

