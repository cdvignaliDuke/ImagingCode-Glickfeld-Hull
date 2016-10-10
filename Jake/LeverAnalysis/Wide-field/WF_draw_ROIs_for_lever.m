function WF_draw_ROIs_for_lever(days, image_dest);
%Use this function to manually draw ROIs after running ShrinkMovie. 

BIN_SIZE = 1;  
% ---------------
    image_dest  = [image_dest '\' days '_ROIshrink.tif'];   %modified by JH
    info = imfinfo(image_dest);
    if(~exist('first_frame', 'var')) 
        f_frame =1;
    else
        f_frame = first_frame;
    end
    if(~exist('last_frame', 'var'))
        l_frame =length(info);
    else
        l_frame = last_frame;
    end
    %----loads the imaging data
    [img, sz]  = get_movie_by_ROI(image_dest, info,[], [], BIN_SIZE, f_frame, l_frame);
    % remove average
    avg_img = mean(img,2);
    std_img = std(img,[], 2);
    all_sd = std(img(:));
    
    %figure; imagesc(reshape(std_img', sz)); axis ij; colormap jet; colormap jet; shading flat;
    %title('std_img just for comparison');
    f = figure; im = imagesc(reshape(avg_img', sz)); axis ij; colormap jet; colormap jet; shading flat;
    num_cluster = input('Enter no of clusters:');
    
    for i=1:num_cluster
        roi(i) =  impoly;
        disp(['Finished ROI ' num2str(i)]);
    end
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
    save([image_dest(1:end-4) 'cluster'], 'cluster');
