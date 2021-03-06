%This script is based off of cluster_vermis_movies.m  

BIN_SIZE = 1;  
DATA_DIR =  'S:\VermisImaging\';
FRAME_TIME_DIR = 'C:\Users\jake\TempData\';
BEHAVE_DIR = 'C:\Users\jake\TempData\behavior\';
% ---------------
days = {'160613_img40_optic_2'};
%days = {'150716_img27', '150717_img27', '150718_img27', '150719_img27'};  
%days = {'150719_img28', '150717_img28', '150716_img28'};
%shrink_movie_b(days{1},20)

for kk=1:length(days)
    ROI_name  =  days{kk};
    %temporarily deleted "days{kk} '\'" after DATA_DIR 
    image_dest  = [DATA_DIR days{kk} '\160613_img40sub_optic_2_MMStack.ome.tif'];   %modified by JH
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
    [img, sz]  = get_movie_by_ROI(image_dest, info,[], [], BIN_SIZE, f_frame, l_frame);
    % remove avergae
    avg_img = mean(img,2);
    std_img = std(img,[], 2);
    all_sd = std(img(:));
    
    f = figure; im = imagesc(reshape(std_img', sz)); axis ij; colormap jet; colormap jet; shading flat;
    num_cluster = input('Enter no of clusters:')
    
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
    cluster.allmask = reshape(sum(roi_mask,1),[1002,1004]); %Converts row vector of mask to matrix. HARDCODED for 1002 x 1004 pixel frames.
    for ii=1:num_cluster
        cluster.(sprintf('mask%d',ii))=reshape(roi_mask(ii,:),[1002,1004]); %Creates a separate matrix for each mask and stores into corresponding variables.
    end
    close gcf;
    save([image_dest(1:end-8) 'cluster'], 'cluster');
end
