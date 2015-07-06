clear;
%  go over all the movies in a session
% 1. select an ROI
% 2. stack all the movies to one tiff file

DATA_DIR =  'C:\Users\jake\TempData\';
%DATA_DIR = '/Volumes/Promise\ RAID\ 1/mati/imaging_data/150206_img16/'
%DATA_DIR =  '/Users/mati/Documents/matlab/mouse_data/';
day = '150309_img16';
%day = '150209_img18';
session = '';
image_dest  = [DATA_DIR day '\'];
BIN_SIZE =1;

old_cd = cd(image_dest);   
%------ assume the file has a a unique string [_MMStack_[???].ome] and ?? is the number of the file -1
img_str_indicator = 'MMStack';
all_files = dir(['*' img_str_indicator '*'] );
% ---- get order of imaging
file_order = [];
for i=1:length(all_files)
    beg_inx = findstr(all_files(i).name, img_str_indicator) + length(img_str_indicator);
    end_no_inx = strfind(all_files(i).name , '.ome');
    if(end_no_inx == beg_inx)
        file_order(i) = 1;
    else
        file_order(i) = str2num(all_files(i).name(beg_inx+1:end_no_inx-1))+1;
    end
end
all_files(file_order) = all_files;

for i=1:length(all_files)
    info{i} = imfinfo(all_files(i).name);    
    if(i==1)
        all_info= info{i};
    else
        all_info = cat(1, all_info,info{i});
    end
end
frame_times = get_frame_time_by_movie_info(all_info);

dest =  [image_dest day '_ROI'];
save([dest '_frame_times'],  'frame_times');

% use first file to calculate ROI
cd('C:\Users\jake\Documents\MATLAB\MatiCode\Lever');
[ROI_x, ROI_y] = get_movie_ROI(info{1}, 1, 1000); % get the ROI -  must be a rectangle   
img = {};
sz = {};
for i=1:length(all_files)
    [img, sz] = get_movie_by_ROI(info{i}(1).Filename, info{i}, ROI_x, ROI_y, BIN_SIZE, 1, length(info{i}));
    
    % ----- write tif file with ROI only
    for j=1:size(img,2)
        c_img = img(:,j);
        s_img = uint8(reshape(c_img,sz));
        mode= 'append';
        if(i==1 && j==1)
            mode = 'overwrite';
        end
        imwrite(s_img, [dest '.tif'], 'WriteMode', mode);
        if(mod(j,200) ==1)
            disp(['movie #' num2str(i) ' frame # ' num2str(j)]);
        end
    end
end
cd(old_cd);