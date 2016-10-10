function [all_files, all_info, info] = obtain_tif_meta_data(image_source)
%This function looks at a specific folder and identifies all the tiff files
%within it. It ASSUMES they were all collected in sequence during a single
%imaging session. It will output the time of each frame collected
%determined from the tiff metadata and the info file for each set of
%frames. 

cd(image_source);  %point MATLAB at the folder containing the image files. 
%------ assume the file has a a unique string [_MMStack_[???].ome] and ??? is the number of the file -1
img_str_indicator = 'MMStack';
all_files = dir(['*' img_str_indicator '*'] );
% ---- get order of imaging files  (often micromanager stores ~4k frames in a single file then creates a new one for sessions with more than ~4k frames)
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
        all_info = info{i};
    else
        all_info = cat(1,all_info, info{i});
        %all_info(i) = info{i};
    end
end
end