edit Load_SBXdataset_fast.m

data_sub = data-min(min(min(data,[],1),[],2),[],3);
data = data_sub;
clear data_sub

data_avg = mean(data(:,:,2000:2100),3);
figure; imagesq(data_avg); colormap(gray)

[out data_reg] = stackRegister(data, data_avg);
clear data

% get cells
masksFolder = '005'; %changes according to dataset
fileDirMasks = fullfile('Z:\analysis\',mouse,'two-photon imaging', date, masksFolder);

cd(fileDirMasks);
load('mask&TCDir.mat');

% get timecourses
dataTC = stackGetTimeCourses(data_reg, mask_cell);

%save in
try
    fileSave = fullfile('Z:\analysis\',mouse,'two-photon imaging', date, ImgFolder);
    cd(fileSave);
catch
    fileSave = fullfile('Z:\analysis\',mouse,'two-photon imaging', date);
    cd(fileSave);
    mkdir(ImgFolder);
    fileSave = fullfile('Z:\analysis\',mouse,'two-photon imaging', date, ImgFolder);
    cd(fileSave);
end
dataTimecourse.dataTC = dataTC;
dataTimecourse.mouse = mouse;
dataTimecourse.date = date;
dataTimecourse.dataset = ImgFolder;
save('dataTC.mat','dataTimecourse')