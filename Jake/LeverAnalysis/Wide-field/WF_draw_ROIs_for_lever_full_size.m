function WF_draw_ROIs_for_lever_full_size(days, image_source, image_dest);
%Use this function to manually draw ROIs after running ShrinkMovie. 

BIN_SIZE = 1;  
% ---------------
    image_source  = [image_source '\' days '_MMStack.ome.tif'];   %modified by JH
    %----loads the imaging data
    %need to import the meta_data and get the number of frames? 
    img = [];
    for ii = 51:250;
        img_temp = double(imread(image_source,ii));
        img = cat(3,img,img_temp);
    end
    % remove average
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
