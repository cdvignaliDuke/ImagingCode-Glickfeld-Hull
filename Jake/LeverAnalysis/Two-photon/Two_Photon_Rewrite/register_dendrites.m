%% script for registering dendrites across sessions/days

%% Register masks across days 
file_info;
img94 = [];
img94.day1 = [];
img94.postlearning = [];

%load data
day = '170524_000_img94';
mask_dir = 'Z:\Analysis\Cue_reward_pairing_analysis\2P\';
load([mask_dir, day, '\Results']);
mask_final_binary = mask_final;
mask_final_binary = reshape(mask_final_binary,size(mask_raw));
mask_final_binary(mask_final_binary>0) = 1;
img94.day1 = mask_final_binary;

day = '170529_000_img94';
load([mask_dir, day, '\Results']);
mask_final_binary = mask_final;
mask_final_binary = reshape(mask_final_binary,size(mask_raw));
mask_final_binary(mask_final_binary>0) = 1;
img94.postlearning = mask_final_binary;

mask_tiff = cat(3, img94.day1, img94.postlearning);
writetiff(mask_tiff, [mask_dir, day, '\mask_tif']);

%% Calculate centroid for each dendrite

%allocate memory and define variables
num_cells = size(mask3D,3);
centroid_coords = nan(2,num_cells);

%extract centroid for each dendrite and store
for ii = 1:num_cells
    temp_coords = regionprops(mask3D(:,:,ii), 'Centroid'); %obtain centroid for each neuron
    centroid_coords(:,ii) = temp_coords.Centroid; %top row =x   bottom row=y
end
centroid_coords = round(centroid_coords);

%plot centroids onto mask to verify they are a good fit
mask_final = reshape(mask_final, size(mask_raw));

mask_flat_test = mask_flat;
for ii = 1:num_cells
    mask_flat_test([centroid_coords(2,ii)-1:centroid_coords(2,ii)+1],[centroid_coords(1,ii)-1:centroid_coords(1,ii)+1], :) = 1;
end
figure; imagesc(mask_flat_test);
  

