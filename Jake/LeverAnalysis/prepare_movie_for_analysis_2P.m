clear;
%  go over all the movies in a session
% 1. select an ROI
% 2. stack all the movies to one tiff file

base_dir =  '\\CRASH.dhe.duke.edu\data\home\jake\';
SubNum = '924';
date = '150703';
run = '_000_000';
time = '1821';
mouse = 'img24';
image_dest  = fullfile(base_dir, [date '_' mouse], mouse);
BIN_SIZE =1;

out_base = 'Z:\home\lindsey\Analysis\2P\Jake';
run_name = [date '_' mouse '_run' run(length(run)-2:end)];
out_path = fullfile(out_base,run_name);
% load MWorks file
behav_dir = [base_dir '\Data\2P_imaging\behavior'];
cd(behav_dir);
mworks = ['data-' 'i' SubNum '-' date '-' time '.mat']; 
load (mworks);

% find sbx info files
data_dir = [base_dir '\Data\2P_imaging\' date '_' mouse '\' mouse];
cd(data_dir);
fName = [mouse run];
imgMatFile = [fName '.mat'];
load(imgMatFile);

[frame_times frame_count input] = get_frame_time_by_counters(input, info);

dest =  fullfile(out_path,run_name);
save([dest '_frame_times.mat'],  'frame_times', 'input');

%load and register dataset
nframes = info.config.frames;
tic
data = sbxread(fName,0,nframes);
toc
data = squeeze(data);

%remove negative data (by addition)
data_sub = data-min(min(min(data,[],1),[],2),[],3);
clear data

% register
data_avg = mean(data_sub(:,:,90:190),3);
figure; imagesq(data_avg); colormap(gray)

[out data_reg] = stackRegister(data_sub, data_avg);
clear data_sub

dest =  fullfile(out_path,run_name);
save([dest '_data_reg.mat'],  'data_reg');

% use first file to calculate ROI
[ROI_x, ROI_y] = get_2P_ROI(data_reg); % get the ROI -  must be a rectangle   
data_reg = data_reg(ROI_x,ROI_y,:);
writetiff(data_reg,[dest '_ROI.tif']);
