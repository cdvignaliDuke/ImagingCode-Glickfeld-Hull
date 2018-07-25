
clear all; clear global; clear cd;
subjNum = '59';
session_date = '161103';
file_num = '_000_000';
frame_rate = '30Hz ';
frames_plotted = '1:1500';
cd(['Z:\Data\2P_imaging\', session_date, '_img', subjNum '\img', subjNum]);
%cd('Z:\Data\2P_imaging\171121_img050\img050');
%cd('Z:\Data\2P_imaging\img040\170912_img039\img039');
%cd('Z:\Data\2P_imaging\180507_img085\img085');
%cd('Z:\Data\2P_imaging\WindowOutcomes\180509_img085\img085');
fname = ['img', subjNum, file_num];
data = squeeze(sbxread(fname,0,1500));
load([fname, '.mat']);
if length(size(data)) ==4
    data = squeeze(data(1,:,:,:));
    data2 = squeeze(data(2,:,:,:));
    data_avg = mean(data,3);
else
    data_avg = mean(data,3);
end
data_max= max(data,[],3);
figure; imagesc(data_avg); colormap gray; %truesize;
title([frame_rate, subjNum, ' ', session_date, ' recon avg frames ', frames_plotted, ' ', file_num]);
writetiff(data, ['Z:\Data\2P_imaging\', session_date, '_img', subjNum, '\img', subjNum, '_tiff_ 1_1500']);

%% write a tiff movie to analyze diff. in sessions
% 
% subjNum = '070';
% session_date = '180104';
% cd(['Z:\Data\2P_imaging\', session_date, '_img', subjNum '\img', subjNum]);
% %cd(['Z:\Data\2P_imaging\WindowOutcomes\', session_date, '_img', subjNum '\img', subjNum]);
% fname = ['img', subjNum, '_000_000'];
% for ii = 1%:20
%     data = squeeze(sbxread(fname,0,800));
%     data_avg = mean(data,3);
%     data_max= max(data,[],3);
%     data2 = data(:,:,[1:3:end]);
%     figure; imagesc(data_max); colormap gray;
%     title(['30HZ ', subjNum, ' ', session_date, '  recon max proj frames  1:800']);
%     truesize;
% end
% writetiff(data, ['Z:\Data\2P_imaging\', session_date, '_img', subjNum, '\img', subjNum, '_tiff_ 1_801']);
% %writetiff(data, ['Z:\Data\2P_imaging\WindowOutcomes\', session_date, '_img', subjNum, '\img', subjNum, '_tiff_ 1_800']);

%% motion registration 
% 
% %determine frame nums to use for selecting stable frame
% sz = size(data); 
% laser_on_ind = [1:sz(3)];
% frame_nums_for_ref30 = laser_on_ind(randi([1,size(laser_on_ind,2)],1,30));
% frame_nums_for_samp100 = laser_on_ind(round(linspace(1,size(laser_on_ind,2))));
% 
% % select 30 random frames from throughout the movie
% ref30 = data(:,:,frame_nums_for_ref30);
% 
% %motion register each of the 30 random frames to 100 frames from the movie. Find the one with the lowerst dshift
% samp100 = data(:,:,frame_nums_for_samp100);
% dshift = [];
% for r = 1:size(ref30,3)
%     [reg_out,aa] = stackRegister(samp100, ref30(:,:,r));
%     dshift = [dshift;mean(((reg_out(:,3).^2)+(reg_out(:,4).^2)))];
% end
% 
% %pick the frame which had the lowest dshift and motion register full movie to that frame
% min_f = find(dshift == min(dshift));
% if length(min_f) > 1
%     min_f = min_f(1,1);
% end
% img_ref = ref30(:,:,min_f);
% [reg_out, img_reg] = stackRegister(data, img_ref);
% %img_reg = img_reg(:,:,[1:3:end])
% writetiff(img_reg, 'Z:\Data\2P_imaging\160416_img90\img90_tiff_motion_reg');



