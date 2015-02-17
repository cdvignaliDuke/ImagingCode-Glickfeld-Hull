

function [img, sz]  = get_movie_by_ROI(image_dest, info,ROI_x, ROI_y, BIN_SIZE, first_frame, last_frame)
img = [];
sz = [];
counter  = 0;
% --- main loop, read data from .tiff file
for k = first_frame:last_frame
    counter = counter+1;
    c_img =  double(imread(image_dest, k, 'info', info));
    c_img = c_img(ROI_x, ROI_y);
    
    
    % ----bin image 
    % --- cut edge 
    n_sz = size(c_img) - mod( size(c_img), BIN_SIZE);
    c_img = c_img(1:n_sz(1), 1:n_sz(2));
    % ----- bin data 
    b_sz =  n_sz/BIN_SIZE;
    bin_img=reshape(mean(mean(reshape(c_img,[BIN_SIZE b_sz(1) BIN_SIZE b_sz(2) ]),1),3),[b_sz(1) b_sz(2)]);
    img(:,counter) = bin_img(:); % save everything as a vectors, in the ende reshape to matrix 
    
    assert(median(img(:,counter) ) > 1); % just to be sure that camera was on
    if(counter ==1)
        sz = size(bin_img);     
    else
        assert(isequal(sz, size(bin_img))); % all images must be in the same size
    end
    
    
end
