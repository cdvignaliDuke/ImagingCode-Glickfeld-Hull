
clear all; clear global; clear cd;
subjNum = '1504';
session_date = '191016';
file_num = '_000_000';
frame_rate = '15Hz ';
frame_start = 0; % 1386:1639
num_frames = 700;
%cd('Y:\home\elizabeth\2Pdata\190622_19_1uMDART\x19');
cd(['Y:\home\jake\Data\2P_imaging\', session_date, '_img', subjNum '\img', subjNum]);
%cd(['Y:\home\jake\Data\2P_imaging\WindowOutcomes\', session_date, '_img', subjNum '\img', subjNum]);
%cd('Y:\home\jake\Data\2P_imaging\WindowOutcomes\EF_X_G7\EF_X');
fname = ['img', subjNum, file_num];
%fname = ['EF_X_000_000'];
data = squeeze(sbxread(fname,frame_start,num_frames, '.sbx'));
%data = squeeze(sbxread('img036_000_000',0,15, '.sbx'));
load([fname, '.mat']);
if length(size(data)) ==4
    data2 = squeeze(data(2,:,:,:));
    data = squeeze(data(1,:,:,:));
    data_avg = mean(data,3);
    data_avg2 = mean(data2,3);
else
    data_avg = mean(data,3);
end
data_max= max(data,[],3);
% figure; imagesc(data_avg); colormap gray; %truesize;
% title([frame_rate, subjNum, ' ', session_date, ' recon avg frames ', num2str(frame_start), ':', num2str(num_frames), ' ', file_num]);
writetiff(data, ['Y:\home\jake\Data\2P_imaging\', session_date, '_img', subjNum, '\img', subjNum, '_tiff_', num2str(frame_start), '_', num2str(num_frames),]);
%writetiff(data, ['Y:\home\jake\Data\2P_imaging\WindowOutcomes\', session_date, '_img', subjNum, '\img', subjNum, '_tiff_', num2str(frame_start), '_', num2str(num_frames),]);
%writetiff(data, ['Y:\home\jake\Data\2P_imaging\WindowOutcomes\img098\190612\', subjNum, '_tiff_', num2str(frame_start), '_', num2str(num_frames)]);
%writetiff(data, ['Y:\home\jake\Data\2P_imaging\WindowOutcomes\EF_96_G7\EF_X_tiff_', num2str(frame_start), '_', num2str(num_frames)]);


%% write a tiff movie to analyze diff. in sessions
% 
% subjNum = '070';
% session_date = '180104';
% cd(['Y:\home\jake\Data\2P_imaging\', session_date, '_img', subjNum '\img', subjNum]);
% %cd(['Y:\home\jake\Data\2P_imaging\WindowOutcomes\', session_date, '_img', subjNum '\img', subjNum]);
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
% writetiff(data, ['Y:\home\jake\Data\2P_imaging\', session_date, '_img', subjNum, '\img', subjNum, '_tiff_ 1_801']);
% %writetiff(data, ['Y:\home\jake\Data\2P_imaging\WindowOutcomes\', session_date, '_img', subjNum, '\img', subjNum, '_tiff_ 1_800']);

%% motion registration 

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
% img_ref = mean(data(:,:,[500:617]),3);
% min_f = find(dshift == min(dshift));
% if length(min_f) > 1
%     min_f = min_f(1,1);
% end
% img_ref = ref30(:,:,min_f);
% 
% img_ref = mean(data(:,:,175:210),3);
% [reg_out, img_reg] = stackRegister(data, img_ref);
% %img_reg = img_reg(:,:,[1:3:end])
% %writetiff(img_reg, 'Y:\home\jake\Data\2P_imaging\180322_img077\img077_tiff_motion_reg_4886');
% writetiff(img_reg, ['Y:\home\jake\Data\2P_imaging\WindowOutcomes\EF_96_G7\EF_X_tiff_motion_reg']);
% 
% %average every third frame to denoise
% img_reg_proj = stackGroupProject(img_reg, 3);
% writetiff(img_reg_proj, ['Y:\home\jake\Data\2P_imaging\WindowOutcomes\EF_96_G7\EF_X_tiff_motion_reg_3rd_avg']);
% 
% 
% % apply a gaussian filter to denoise the data
% img_reg_filter = stackFilter(img_reg); % filter IC with gaussian filter
% writetiff(img_reg_proj, ['Y:\home\jake\Data\2P_imaging\WindowOutcomes\EF_X_G7\EF_X_tiff_motion_reg_gauss_filt']);
% 
% %calculate df/f 
% avg_f = mean(img_reg(:,:,[2:100]),3);
% diff_f = bsxfun(@minus, double(img_reg), avg_f);   
% dfof = bsxfun(@rdivide, diff_f, avg_f); 
% writetiff(dfof, ['Y:\home\jake\Data\2P_imaging\WindowOutcomes\EF_X_G7\EF_X_tiff_motion_reg_dfof']);
% 
% 
% 
% 



