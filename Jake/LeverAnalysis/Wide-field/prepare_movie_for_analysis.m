clear;
tic
%  go over all the movies in a session
% 1. select an ROI
% 2. stack all the movies to one tiff file
DATA_DIR =  'Z:\Data\WidefieldImaging\GCaMP\';
day =  '160725_img53' %'160208_img35_2' % '160209_img36_NoLever_1'};   %'150717_img28' '150719_img28'
motion_correction =1;           %1 to register image   0 to take ROI without motion correction
stable_int = [1:50];         %use imageJ to find a series of frames (~100) in which little movement occurs. Use this during motion correction
session = '';
image_dest  = ['C:\Users\jake\TempData\' day '\'];
image_source  = ['Z:\Data\WidefieldImaging\GCaMP\' day '\'];
BIN_SIZE =1;
%MUST ALSO SET ncores for parpool

old_cd = cd(image_source);   
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
%------------------------------------------------------
%------------------------------------------------------
%%

% use first file to calculate ROI
cd('C:\Users\jake\Documents\Repositories\ImagingCode-Glickfeld-Hull\Jake\LeverAnalysis');
disp('please select desired ROI and then press "s"');
[ROI_x, ROI_y] = get_movie_ROI(info{1}, 1, 1000); % get the ROI -  must be a rectangle   
disp('THANKS');
img = {};
sz = {};

%% MOTION CORRECTION AND ROI EXTRACTION
if motion_correction == 1;
    % register image   motion correction
    data_stable = [];
    data = readtiff([DATA_DIR day '\'], [], num2str(all_files(1).name(1:end-4)));
    data_stable = data(:,:,stable_int);
    data_reg = mean(data_stable,3);
    figure; imagesq(data_reg); colormap(gray)
    disp('saving data_reg');
    dest =  [image_dest day '_reg'];
    save(dest, 'data_reg');
    
    data_reg_ROI = [];
    for i=1:length(all_files);
        if i > 1;
            disp(['reading tiff file' i]);
            data = readtiff([DATA_DIR day '\'], [], num2str(all_files(i).name(1:end-4)));
        end
        disp(['registering image file ' i]);
        [out data_reg_sec] = stackRegister(data, data_reg);     %could make stack regist parallel to decrease time to process. 
        data_reg_sec = data_reg_sec(ROI_x,ROI_y,:);
        if i == 1;
            data_reg_ROI = data_reg_sec;
        elseif i > 1;
            data_reg_ROI = cat(3, data_reg_ROI, data_reg_sec);
        end
    end
    disp('motion corrected ROI calculated. Now writing tif');
    %code for writing tiff after motion correction
    nframes = size(data_reg_ROI,3); %16920(9cores)    or  16936(8cores)
    ncores = 10;
    parpool = ncores;
    %set the size of the batch that each stream will read
    nbatch = nframes./ncores;
    %create a parallel loop to load the data
    parfor i=1:ncores                              %will have to write something into HAD_cmp in order to stitch the batches back together. 
        writetiff(data_reg_ROI(:,:,1+((i-1)*nbatch):i*nbatch),[image_dest day '_ROI_' num2str(i) '.tif']); 
    end
    %writetiff(data_reg_ROI,[dest '_ROI.tif']);
    
    
elseif motion_correction == 0;
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
end
cd(old_cd);
toc