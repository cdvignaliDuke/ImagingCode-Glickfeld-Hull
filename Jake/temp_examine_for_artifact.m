% 
% clear all; clear global; clear cd;
% subjNum = '666';
% session_date = '190612';
% file_num = '_000_000';
% frame_rate = '30Hz ';
% frame_start = 0;
% num_frames = 10000;
% cd(['Y:\home\jake\Data\2P_imaging\Troubleshooting\', session_date,'\', subjNum '\']);
% fname = [subjNum, file_num];
% data = squeeze(sbxread(fname,frame_start,num_frames, '.sbx'));
% load([fname, '.mat']);
% if length(size(data)) ==4
%     data = squeeze(data(1,:,:,:));
%     data2 = squeeze(data(2,:,:,:));
%     data_avg = mean(data,3);
% else
%     data_avg = mean(data,3);
% end
% data_max= max(data,[],3);
% title([frame_rate, subjNum, ' ', session_date, ' recon avg frames ', num2str(frame_start), ':', num2str(num_frames), ' ', file_num]);
% %writetiff(data, ['Y:\home\jake\Data\2P_imaging\WindowOutcomes\', session_date, '_img', subjNum, '\img', subjNum, '_tiff_', num2str(frame_start), '_', num2str(num_frames),]);


%% check pcokels modulation

clear all; clear global; clear cd;
subjNum = '1032';
session_date = '190615';
file_num = '_000_000';
frame_rate = '30Hz ';
num_frames = 3000;
%cd(['Y:\home\jake\Data\2P_imaging\190615_img1030\img1030\']);
cd(['Y:\home\jake\Data\2P_imaging\190615_img1032\img1032\']);

for frame_start = [0, 30000, 60000, 90000]
    fname = ['img', subjNum, file_num];
    data = squeeze(sbxread(fname,frame_start,num_frames, '.sbx'));
    load([fname, '.mat']);
    data_avg = squeeze(mean(mean(data,2),1))';
    figure; plot(data_avg); title([session_date, ' ', subjNum, ' ', num2str(frame_start)]);
end





