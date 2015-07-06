edit Load_SBXdataset_fast.m

down = 10;

data_down = stackGroupProject(data,down);
clear data

%remove negative data (by addition)
data_sub = data_down-min(min(min(data_down,[],1),[],2),[],3);
clear data_down

% register
data_avg = mean(data_sub(:,:,100:110),3);
figure; imagesq(data_avg); colormap(gray)

[out data_reg] = stackRegister(data_sub, data_avg);
clear data_sub

% %roi found by pixel correlation then selection with imCellEditInteractive
% b = 5;
% siz = size(data_reg);
% corr_map = zeros(siz(1),siz(2));
% for ix = b:siz(2)-b
%     for iy = b:siz(1)-b
%         TC = data_reg(iy,ix,:);
%         surround = (data_reg(iy-1,ix-1,:)+data_reg(iy-1,ix,:)+data_reg(iy-1,ix+1,:)+data_reg(iy,ix-1,:)+data_reg(iy,ix+1,:)+data_reg(iy+1,ix-1,:)+data_reg(iy+1,ix,:)+data_reg(iy+1,ix+1,:))/8;
%         R = corrcoef(TC,surround);
%         corr_map(iy,ix) = R(1,2);
%     end
% end
% 
% figure; imagesq(corr_map); colormap(gray)

maxF = max(data_reg,[],3);
figure; imagesq(maxF); colormap(gray)

bwout = imCellEditInteractive(maxF); %using 0.9 and 3
mask_cell = bwlabel(bwout);

% save directory
fileSave = fullfile('Z:\analysis\',mouse,'two-photon imaging',date,ImgFolder);
cd(fileSave);

save('mask.mat','mask_cell');

