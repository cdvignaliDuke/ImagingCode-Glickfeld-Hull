%Use this function to draw ROIs after running ShrinkMovie. 
clear

BIN_SIZE = 1;  
DATA_DIR =  'C:\Users\jake\TempData\';
FRAME_TIME_DIR = 'C:\Users\jake\TempData\';
BEHAVE_DIR = 'Z:\Data\WidefieldImaging\GCaMP\behavior';
% ---------------
%days = {'160131_img36','160131_img35','160129_img36','160129_img35'};
%days = {'150716_img27', '150718_img27', '150719_img27', '150716_img28', '150717_img28', '151021_img29', '151022_img29', '151009_img30', '151011_img30', '151211_img32', '151212_img32'};
%days = {'150716_img27', '150717_img27', '150718_img27', '150719_img27'};  
%days = {'150719_img28', '150717_img28', '150716_img28'};
%days = {'160208_img36', '160207_img35', '160207_img36', '160205_img35'};
days = {'160315_img38'};

for kk=1:length(days)
    ROI_name  =  days{kk};
    image_dest  = [DATA_DIR days{kk} '\' ROI_name '_ROIshrink.tif'];   %modified by JH
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
end
