
clear all; clear global; clear cd;
cd('Z:\Data\2P_imaging\170911_img036\img036');
%cd('Z:\Data\2P_imaging\img040\170912_img039\img039');
%cd('Z:\Data\2P_imaging\WindowOutcomes\img044\170921_img044\img044');
fname = 'img036_000_002';
data = squeeze(sbxread(fname,0,1000));
if length(size(data)) ==4
    data1 = squeeze(data(1,:,:,:));
    data2 = squeeze(data(2,:,:,:));
    data_avg = mean(data1,3);
else
    data_avg = mean(data,3);
end
data_max= max(data,[],3);
figure; imagesc(data_avg); colormap gray; truesize;
title(['30HZ img036 170911 recon avg frames 1:1000 001']);


%% write a tiff movie to analyze diff. in sessions

cd('Z:\Data\2P_imaging\170912_img039\img039');
%cd('Z:\Data\2P_imaging\img039\170910_img039\img039');
fname = 'img039_000_000';
for ii = 1%:20
    data = squeeze(sbxread(fname,0,1000));
    data_avg = mean(data,3);
    data_max= max(data,[],3);
    data2 = data(:,:,[1:3:end]);
    figure; imagesc(data_max); colormap gray;
    title(['30HZ img039 170912 recon max proj frames 1:1000']);
end
writetiff(data, 'Z:\Data\2P_imaging\170912_img039\img039_tiff_1_1000');


%% motion registration 

%determine frame nums to use for selecting stable frame
sz = size(data); 
laser_on_ind = [1:sz(3)];
frame_nums_for_ref30 = laser_on_ind(randi([1,size(laser_on_ind,2)],1,30));
frame_nums_for_samp100 = laser_on_ind(round(linspace(1,size(laser_on_ind,2))));

% select 30 random frames from throughout the movie
ref30 = data(:,:,frame_nums_for_ref30);

%motion register each of the 30 random frames to 100 frames from the movie. Find the one with the lowerst dshift
samp100 = data(:,:,frame_nums_for_samp100);
dshift = [];
for r = 1:size(ref30,3)
    [reg_out,aa] = stackRegister(samp100, ref30(:,:,r));
    dshift = [dshift;mean(((reg_out(:,3).^2)+(reg_out(:,4).^2)))];
end

%pick the frame which had the lowest dshift and motion register full movie to that frame
min_f = find(dshift == min(dshift));
if length(min_f) > 1
    min_f = min_f(1,1);
end
img_ref = ref30(:,:,min_f);
[reg_out, img_reg] = stackRegister(data, img_ref);
%img_reg = img_reg(:,:,[1:3:end])
writetiff(img_reg, 'Z:\Data\2P_imaging\160416_img90\img90_tiff_motion_reg');
    
    
    
    
    