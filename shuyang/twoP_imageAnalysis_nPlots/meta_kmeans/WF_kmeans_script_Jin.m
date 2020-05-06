%%WF_kmeans_script
clear;
days = {'151021_img29', '151022_img29', '151009_img30', '151011_img30', '151211_img32', '151212_img32', '160129_img35', '160131_img35', ...
    '160129_img36', '160131_img36', '160314_img38', '160315_img38', '160319_img41', '160320_img41', '160606_img46', '160722_img53', '160904_img55'}; %'150718_img27', '150719_img27', '150716_img28', '150717_img28',
bx_source     = 'Y:\home\jake\Data\WidefieldImaging\GCaMP\behavior\';
image_source_base  = 'Y:\home\jake\Data\WidefieldImaging\GCaMP\'; %location of permanently stored image files for retreiving meta data
image_dest_base    = 'Y:\home\jake\Analysis\WF Lever Analysis\BxAndAnalysisOutputs\'; %stores the data on crash in the lever analysis folder
bx_outputs_dir =     'Y:\home\jake\Analysis\WF Lever Analysis\BxAndAnalysisOutputs\BxOutputs\';
PCA_output_dir_base =    'Y:\home\jake\Analysis\WF Lever Analysis\PCA_output_dir\';
kmeans_inputs_dir =      'Y:\home\jake\Analysis\WF Lever Analysis\kmeans_output_dir\';
kmeans_output_dir_base = 'Y:\home\jake\Analysis\WF Lever Analysis\Meta-kmeans_output_dir\';
downsample_factor = 0.2;

for ii = [17] %[3,4,5,6,11,12]  %original three datasets = [15, 16, 17]     Additionaly four datasets which could be analyzed = [9,10,13,14]
    usFacs =100;
    PCA_output_dir = [PCA_output_dir_base days{ii} '\'];
    kmeans_output_dir = [kmeans_output_dir_base days{ii} '\'];
    kmeans_inputs_dir = [kmeans_inputs_dir days{ii} '\'];
    prepMovie_dir1 = [kmeans_inputs_dir, days{ii}, '_downsampled_cropped_motionCorr.mat'];
    prepMovie_dir2 = [kmeans_output_dir, days{ii}, '_downsampled_cropped_motionCorr_blackout.mat'];
    outdir = [kmeans_output_dir];
    if ~exist(outdir)
        mkdir(outdir);
    end
    if exist(prepMovie_dir1, 'file') == 2
        load(prepMovie_dir1);
    elseif exist(prepMovie_dir2, 'file') == 2
        load(prepMovie_dir2);
    else
        
        meta_data_dir = [image_dest_base days{ii} '\' days{ii} '_meta_data'];
        load(meta_data_dir);
        
        %download the entire movie
        full_movie = download_full_tiff_movie(meta_data2);
        %     %crop movie
            [ROI_x, ROI_y] = get_2P_ROI(full_movie);
            full_movie = full_movie(ROI_x, ROI_y, :);
        
        %motion correction of the entire movie
        if exist([kmeans_output_dir, days{ii}, 'reg_out.mat'],'file') == 2
            load([kmeans_output_dir, days{ii}, 'reg_out.mat']);
            [~,full_movie_reg]=stackRegister_MA(full_movie,full_movie(:,:,1),usFacs,out);
        else
            rf = 30*100;
            ref30 = full_movie(:,:,10:100:rf);
            sf = 100*100+10;
            samp100 = full_movie(:,:,11:100:sf);
            dshift = [];
            for r = 1:size(ref30,3)
                
                [reg_out,aa] = stackRegister(samp100, ref30(:,:,r));
                dshift = [dshift;mean(((reg_out(:,3).^2)+(reg_out(:,4).^2)))];
                
            end
            
            min_f = find(dshift == min(dshift));
            img_ref = ref30(:,:,min_f);
            [out, full_movie_reg] = stackRegister(full_movie, img_ref);
        end
        clear full_movie; %the motion correction does not alter the size of the movie
        
%         full_movie_reg = MaskBackground(full_movie);
        full_movie_reg = MaskBackground(full_movie_reg); %make movie into binary masks?
        full_movie_reg = cropImage(full_movie_reg); %only save the part that has fluorescence
        
        %down sample movie to average activity over 5x5pixel blocks
        
        %     downsampled_movie = double(imresize(cropped_movie, downsample_factor));
        downsampled_movie = double(imresize(full_movie_reg, downsample_factor)); %imresize is a Matlab function that makes the image smaller.
        
        %clear full_movie_reg
        
        %save the downsampled, cropped, mask applied movie to kmeans_output_dir
        if ~exist(kmeans_output_dir)
            mkdir(kmeans_output_dir)
        end
        save([kmeans_output_dir, days{ii}, '_downsampled_cropped_motionCorr_blackout'], 'downsampled_movie');
        save([kmeans_output_dir, days{ii}, 'reg_out.mat'], 'out', 'ROI_x', 'ROI_y');
    end
    % load behavior data
    bxName = days{ii};
    subName = bxName(end-1:end); subDay = bxName(1:6);
    bxDir = dir([bx_source '*' subName '*' subDay '*.mat']);
    load([bx_source bxDir.name]);
    
    subID = days{ii};
    load([outdir, 'Correlation.mat']);
    WF_behav;
    
    save([outdir, 'Correlation.mat'], 'outSCorr', 'outFCorr', 'outSCorr_avg', 'outFCorr_avg', 'inSCorr', 'inFCorr', ...
        'inSCorr_avg', 'inFCorr_avg', 'sFrame', 'fFrame', 'cSizeNorm', 'outSCorrN', 'outSCorrNN', 'outFCorrN', 'outFCorrNN');
    
    %metaK means
    figName = 'cluster_wholeMovie.fig';
    fileName = 'clustering_wholeMovie.mat';
    runMetaK(downsampled_movie, outdir, figName, fileName)
    
%          elseif exist([kmeans_output_dir, days{ii}, '_downsampled_cropped_motionCorr.mat'], 'file') == 2;
%                  load([kmeans_output_dir, days{ii}, '_downsampled_cropped_motionCorr.mat']);
%          end
    %this is the code for actually performing the kmeans analysis=======================================================
    %needs to be put inside the forloop
    %check size of the movie
    movie_size = size(downsampled_movie); % downsampled_movie: length*width*frame
    if length(movie_size)==3 % if movie is a 3 dimension thing, reshape it into a 2D thing: pixel*frame
        downsampled_movie = reshape(downsampled_movie, movie_size(1)*movie_size(2), movie_size(3));
    end

    %Run kmeans on df/f data
    %downsampled_movie 2D, consists of zeros in masked areas and integers
    %from 1 to 7 in F areas.
    %calculate the df/F using average F across time as F, normalize df/f in different ways 
    downsampled_movie_mean = mean(downsampled_movie,2);  %dim1=pixels dim2=frames
    downsampled_movie_mean = repmat(downsampled_movie_mean, 1,size(downsampled_movie,2));
    dfof_downsampled = (downsampled_movie-downsampled_movie_mean)./downsampled_movie_mean;
    avg_f_for_sub = reshape(downsampled_movie_mean(:,1), movie_size(1), movie_size(2));
    downsampled_movie_max = max(dfof_downsampled,[],2); % get the biggest value for each pixel across time
    downsampled_movie_max = repmat(downsampled_movie_max, 1, size(downsampled_movie,2)); %repeat for n=nframe times
    downsampled_movie_min = min(dfof_downsampled,[],2);
    downsampled_movie_min = repmat(downsampled_movie_min, 1, size(downsampled_movie,2));
    downsampled_norm_max = dfof_downsampled./downsampled_movie_max;
    downsampled_norm_min_max = (dfof_downsampled-downsampled_movie_min)./(downsampled_movie_max-downsampled_movie_min);

    %run a pair wise correlation to determine number of clusters
    pixel_pairwise_correlation = corr(downsampled_movie');% run coorelation on frame*pixel. coor is doing correlation across column (the 2nd dimension), so to decide spacial synchrony, your second dimension should be pixel/dendrite, not time
    %output is n total pixel*n total pixel
    save([kmeans_output_dir, days{ii}, '_pixel_pairwise_correlation'],'pixel_pairwise_correlation');

    %Remove pairwise correlations of the blacked out pixels with each other
    %and F pixels. Look for rows of all NaNs, index them and remove them from both axes.
    index_of_nan_rows = [];
    for iv = 1:length(pixel_pairwise_correlation)
        if isempty(find(~isnan(pixel_pairwise_correlation(iv,:))))
            index_of_nan_rows = [index_of_nan_rows, iv];
        end
    end
    pixel_pairwise_correlation(:,index_of_nan_rows) = [];
    pixel_pairwise_correlation(index_of_nan_rows, :) = [];
    figure;imagesc(pixel_pairwise_correlation);

    %run kmeans in order to cluster the downsample pixels into groups
    %time on y-axis and pixels on x-axis? then reshape to view the frame
    for num_clusters = 8:10
        kmeans_inx = kmeans(dfof_downsampled,num_clusters); % df/f: pixel*frames k-means returns a 1*pixel vector indicating which cluster that pixel belongs to
        kmeans_inx = reshape(kmeans_inx, movie_size(1), movie_size(2)); %reshape the index into length*width (2D)
        figure; subplot(1,2,1); imagesc(kmeans_inx); colorbar;
        title([days{ii} ' dfof kmeans clustering # clusters = ' num2str(num_clusters)]);
        num_pixels = [];
        for iv = 1:num_clusters %how big (many pixels) is each cluster?
            num_pixels = [num_pixels; length(find(kmeans_inx==iv))];
        end
        subplot(1,2,2); bar([1:num_clusters], num_pixels');
    end
    for num_clusters = 8:10
        kmeans_inx = kmeans(downsampled_norm_max,num_clusters);
        kmeans_inx = reshape(kmeans_inx, movie_size(1), movie_size(2));
        figure; subplot(1,2,1); imagesc(kmeans_inx); colorbar;
        title([days{ii} ' normalized to max: kmeans clustering # clusters = ' num2str(num_clusters)]);
        num_pixels = [];
        for iv = 1:num_clusters
            num_pixels = [num_pixels; length(find(kmeans_inx==iv))];
        end
        subplot(1,2,2); bar([1:num_clusters], num_pixels');
    end
    for num_clusters = 8:10
        kmeans_inx = kmeans(downsampled_norm_min_max,num_clusters);
        kmeans_inx = reshape(kmeans_inx, movie_size(1), movie_size(2));
        figure; subplot(1,2,1); imagesc(kmeans_inx); colorbar;
        title([days{ii} ' normalized to max and min sub: kmeans clustering # clusters = ' num2str(num_clusters)]);
        num_pixels = [];
        for iv = 1:num_clusters
            num_pixels = [num_pixels; length(find(kmeans_inx==iv))];
        end
        subplot(1,2,2); bar([1:num_clusters], num_pixels');
    end


end


