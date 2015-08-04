
function [img, sz]  = get_movie_by_ROI(image_dest, info,ROI_x, ROI_y, BIN_SIZE, first_frame, last_frame)
img = [];
sz = [];
if(~exist('first_frame', 'var') || isempty(first_frame))
    first_frame =1;
end
if(~exist('last_frame', 'var') || isempty(last_frame))
    last_frame = length(info);
end
if(~exist('BIN_SIZE', 'var') || isempty(BIN_SIZE))
    BIN_SIZE =1;
end


c_img =  double(imread(image_dest, first_frame, 'info', info));
if(exist('ROI_x' , 'var') && exist('ROI_y' , 'var') &&   ~isempty(ROI_x) && ~isempty(ROI_y))
    c_img = c_img(ROI_x, ROI_y);
else
    ROI_x = [];
    ROI_y = [];
end
n_sz = size(c_img) - mod( size(c_img), BIN_SIZE);
sz =  n_sz/BIN_SIZE;
vec_sz = sz(1)*sz(2);
x_inx = 1:vec_sz;
n_frame= last_frame - first_frame+1;
img = zeros(vec_sz, n_frame);

% ------ main loop, read data from .tiff file
for k = 1:n_frame
    
    bin_img = get_frame(image_dest, k+first_frame-1, info, ROI_x, ROI_y, BIN_SIZE);
    img(x_inx,k) = bin_img(:); % save everything as a vectors, in the ende reshape to matrix   
    
end

