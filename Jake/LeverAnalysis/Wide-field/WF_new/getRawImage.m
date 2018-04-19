%%WF_kmeans_script
clear;
days = {'151021_img29', '151022_img29', '151009_img30', '151011_img30', '151211_img32', '151212_img32', '160129_img35', '160131_img35', ...
    '160129_img36', '160131_img36', '160314_img38', '160315_img38', '160319_img41', '160320_img41', '160606_img46', '160722_img53', '160904_img55'}; %'150718_img27', '150719_img27', '150716_img28', '150717_img28',
bx_source     = ['Z:\home\jake\Data\WidefieldImaging\GCaMP\behavior\'];
image_source_base  = ['Z:\home\jake\Data\WidefieldImaging\GCaMP\']; %location of permanently stored image files for retreiving meta data
image_dest_base    = ['Z:\home\jake\Analysis\WF Lever Analysis\BxAndAnalysisOutputs\']; %stores the data on crash in the lever analysis folder
bx_outputs_dir =     ['Z:\home\jake\Analysis\WF Lever Analysis\BxAndAnalysisOutputs\BxOutputs\'];
PCA_output_dir_base =    ['Z:\home\jake\Analysis\WF Lever Analysis\PCA_output_dir\'];
kmeans_inputs_dir =      ['Z:\home\jake\Analysis\WF Lever Analysis\kmeans_output_dir\'];
kmeans_output_dir_base = ['Z:\home\jake\Analysis\WF Lever Analysis\Meta-kmeans_output_dir\'];
downsample_factor = 0.2;

for ii = [4,5,6,11,12,15,16,17]  %original three datasets = [15, 16, 17]     Additionaly four datasets which could be analyzed = [9,10,13,14]
    usFacs =100;
    PCA_output_dir = [PCA_output_dir_base days{ii} '\'];
    kmeans_output_dir = [kmeans_output_dir_base days{ii} '\'];
    kmeans_inputs_dir = [kmeans_inputs_dir days{ii} '\'];
    prepMovie_dir1 = [kmeans_inputs_dir, days{ii}, '_downsampled_cropped_motionCorr.mat'];
    prepMovie_dir2 = [kmeans_output_dir, days{ii}, '_downsampled_cropped_motionCorr_blackout.mat'];
    outdir = [kmeans_output_dir];
    
    NUM_FRAMES = 5;
    [~,~, last] = size(downsampled_movie);
    avg = mean((downsampled_movie(:,:,1:NUM_FRAMES:last)),3);
    fig = figure;
    imagesc(avg);axis ij; colormap jet;
    saveas(fig, 'meanIMG.fig');
    print(['../Final_mask/', 'meanIMG.pdf'], '-dpdf');
    
    meta_data_dir = [image_dest_base days{ii} '\' days{ii} '_meta_data'];
    load(meta_data_dir);
    
    %download the entire movie
    full_movie = download_full_tiff_movie(meta_data2);
    
    fig = figure;
    imagesc(full_movie(:,:,1));colormap jet;
    saveas(fig, [outdir, 'Final_mask\', 'fullIMG.fig']);
    print([outdir, 'Final_mask\', 'fullIMG.pdf'], '-dpdf');
    % load behavior data
    %     bxName = days{ii};
    %     subName = bxName(end-1:end); subDay = bxName(1:6);
    %     bxDir = dir([bx_source '*' subName '*' subDay '*.mat']);
    %     load([bx_source bxDir.name]);
    %
    %     subID = days{ii};
    %     WF_behav;
    
    
    
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

