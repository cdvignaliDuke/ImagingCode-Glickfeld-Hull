function shrink_movie(day, image_source, image_dest, meta_data2, BIN_SIZE)
%Bins the data so you can run ClusterVermis WF_lever_Bx_DFoF. Reduces the
%resolution. Saves the shrunken movie. Uses a default bin size of 20 if
%BIN_SIZE is not specified.

% cd(image_source);
% mati_code_cd = 'C:\Users\jake\Documents\Repositories\ImagingCode-Glickfeld-Hull\Jake\LeverAnalysis';
%------ set path. Check to see if BIN_SIZE or meta_data were specified.
dest =  [image_dest '\'  day '_ROI'];
if(~exist('BIN_SIZE', 'var'))
    BIN_SIZE =20;
end
if(~exist('meta_data2', 'var'))
    meta_data2 = load([image_dest '\' day '_ROI' '_meta_data']);  %if meta_data does not exist then look to the analysis folder and load it
end
ROI_x =[]; %just serve as place holders for get_movie_by_ROI inputs
ROI_y = [];
shrink_img = [];
sz = {};
%------ Load and shrink images
%cd(mati_code_cd);
for i=1:length(meta_data2)
    %meta_data does not currently have info about # of
    %frames and it needs that for get_movie to work==========================================================================================================
    [img, sz] = get_movie_by_ROI(meta_data2{i}(1).Filename, meta_data2{i}, ROI_x, ROI_y, BIN_SIZE, 1, length(meta_data2{i})); % this is the function which actually does the shrinking
    %=====================================================================================
    shrink_img = [shrink_img, img];
end
%----- write tif file with ROI only
for i=1:size(shrink_img,2)
    c_img = shrink_img(:,i);
    s_img = uint8(reshape(c_img,sz));
    mode= 'append';
    if(i==1)
        mode = 'overwrite';
    end
    imwrite(s_img, [dest 'shrink.tif'], 'WriteMode', mode);
end
end